from itertools import product
from typing import List, Tuple, Union
import numpy as np
from numpy.linalg import svd
import sys

from gwonsquantum.chemistry.molecule import electronic_structure
from gwonsquantum.polynomial.transform.bravyi_kitaev import beta_matrix
from gwonsquantum.chemistry.vqe.coupled_cluster_vqe.qcc import QCC
from gwonsquantum.chemistry.vqe.coupled_cluster_vqe.uccsd import UCCSD
from gwonsquantum.chemistry.molecule import ElectronicStructure

DEBUG = False

def parse_input_file(inp_file):
    def input_matrix(line, key, result, cnt, triangle=False):
        if result[key] is None:
            result[key] = np.zeros((result['nbasis'], result['nbasis']), dtype=float)
        for i, v in enumerate(line.split()):
            if DEBUG:
                print(i, cnt, v)
            result[key] [cnt][i] = float(v)
            if triangle:
                result[key] [i][cnt] = float(v)

    with open(inp_file, 'r') as of:
        linecnt = 0
        mode = 0
        result = {
            "n_electrons": 0,
            "nbasis": 0,
            "oei": None,
            "s": None,
            "c": None,
            "lst_tei": None,
            "tei": None,
            "constant": None,
            "mo_energies": None,
        }
        while True:
            line = of.readline()
            if not line:
                break
            if len(line)<=1:
                continue
            if line.startswith('NELEC'):
                result['n_electrons'] = int(line.split()[1])
            elif line.startswith("ENUC"):
                result['constant'] = float(line.split()[1])
            elif line.startswith("NBASIS"):
                result['nbasis'] = int(line.split()[1])
            elif line.startswith('OEI'):
                mode = 1
                linecnt = 0
            elif line.startswith("ENERGIES"):
                mode = 5
                linecnt = 0
            elif line.startswith('OVERLAP'):
                mode = 2
                linecnt = 0
            elif line.startswith('MOCOEFF'):
                mode = 3
                linecnt = 0
            elif line.startswith("ERI"):
                mode = 4
                linecnt = 0
                result['eri'] = list()
            elif mode == 1:
                input_matrix(line, 'oei', result, linecnt, triangle=True)
                linecnt += 1
            elif mode == 2:
                input_matrix(line, 's', result, linecnt, triangle=True)
                linecnt += 1
            elif mode == 3:
                input_matrix(line, 'c', result, linecnt)
                linecnt += 1
            elif mode == 4:
                spt = line.split()
                assert len(spt) == 5
                if result['lst_tei'] is None :
                    result['lst_tei'] = list()
                i,j,k,l = int(spt[0]), int(spt[1]), int(spt[2]), int(spt[3])
                result['lst_tei'].append(((i, j, k, l), float(spt[4])))
            elif mode == 5:
                spt = line.split()
                if result['mo_energies'] is None:
                    result['mo_energies'] = list()
                result['mo_energies'].append(float(spt[0]))
    if DEBUG:
        print(result['mo_energies'])
        print(result['oei'])
        print(result['lst_tei'][0])
        print(result['constant'])
        print(result['n_electrons'])
    return result

def lst_eri_to_mat(nbasis, lst_eri):
    def _check_or_store(mat_eri, _i, _j, _k, _l, val):
        if abs(mat_eri[_i][_j][_k][_l]) > 1e-7:
            assert abs(mat_eri[_i][_j][_k][_l] - val) < 1e-7
        else:
            mat_eri[_i][_j][_k][_l] = val

    mat_eri = np.zeros((nbasis, nbasis, nbasis, nbasis))
    # mat_eri[i][j][l][k] = <ij|kl> = (ik|jl)

    # (ij|kl) = <ik|jl>   = mat_eri[i][k][l][j]
    # (ji|kl) = <jk|il>   = mat_eri[j][k][l][i]
    # (ij|lk) = <il|jk>   = mat_eri[i][l][k][j]
    # (ji|lk) = <jl|ik>   = mat_eri[j][l][k][i]

    # (kl|ij) = <ki|lj>   = mat_eri[k][i][j][l]
    # (lk|ij) = <li|kj>   = mat_eri[l][i][j][k]
    # (kl|ji) = <kj|li>   = mat_eri[k][j][i][l]
    # (lk|ji) = <lj|ki>   = mat_eri[l][j][i][k]

    for (i,j,k,l), v in lst_eri:
        _check_or_store(mat_eri, i, k, l, j, v)
        _check_or_store(mat_eri, j, k, l, i, v)
        _check_or_store(mat_eri, i, l, k, j, v)
        _check_or_store(mat_eri, j, l, k, i, v)
        
        _check_or_store(mat_eri, k, i, j, l, v)
        _check_or_store(mat_eri, l, i, j, k, v)
        _check_or_store(mat_eri, k, j, i, l, v)
        _check_or_store(mat_eri, l, j, i, k, v)

    return mat_eri

def call_vqe(mo_contents):
    es: ElectronicStructure ={
        "n_electrons": mo_contents['n_electrons'],
        "oei": mo_contents['oei'],
        "tei": mo_contents['tei'],
        "constant": mo_contents['constant']
    }
    n_basis = mo_contents['nbasis']
    n_spin_basis = 2 * n_basis
    n_electron = mo_contents['n_electrons']
    homo_idx = 2
    lumo_idx = 2
    frozen = [x for x in range(max(n_electron//2 - homo_idx, 0))]
    active = [x for x in range(max(n_electron//2 - homo_idx, 0), min(n_electron//2 + lumo_idx, n_basis))]
    n_frozen = len(frozen)
    n_active = len(active)
    n_virtual = n_basis - n_frozen - n_active

    vqe_obj = QCC(electronic_structure=es,
    #vqe_obj = UCCSD(electronic_structure=es,
                  mapping="JW",
                  frozen_indices=frozen,
                  active_indices=active,
                  hamiltonian_optimization=False,
                  num_entanglers=5
                  )
    opt_ret = vqe_obj.run_optimization(
        backend="statevector_simulator",
        qiskit_kwargs={"optimization_level": 3},
        trot_implementation="star",
        trot_opt_level=2
    )
    energy = opt_ret['min_val']
    xopt = opt_ret['xopt']
    #xopt = [0.0 for _ in range(5)]
    _energy, sv = vqe_obj.single_experiment(
        xopt,
        backend="statevector_simulator",
        qiskit_kwargs={"optimization_level": 3},
        trot_implementation="star",
        trot_opt_level=2,
        return_sv = True
    )
    amp_idx, amps = get_amp_from_sv(sv, "JW")
    for i, a in enumerate(amp_idx):
        amp_idx[i] = "1"*2*n_frozen + a + "0"*2*n_virtual
    return amp_idx, amps, energy 

#TODO: Consider electron number breaking.

def get_amp_from_sv(
    final_state:Union[np.ndarray, List[complex]],
    mapping="JW")\
    -> Tuple[List[str], float]:

    num_qubits = int(np.log2(len(final_state)))
    truncation = 10
    eps = 1e-6
    amps = list()

    for i, c in enumerate(final_state):
        if abs(c) < 1e-7 : continue
        vec = [int(i) for i in bin(i)[2:].rjust(num_qubits, '0')][::-1]
        if mapping == "BK":
            vec = (np.matmul(beta_matrix(len(vec), inv=True), vec) % 2)
        elif mapping == "Bosonic_Direct":
            vec = [vec[i//2] for i in range(len(vec) * 2)]
        elif mapping != "JW":
            raise ValueError
        vec_str = "".join([str(i) for i in vec])  #LSB first
        c = abs(c)**2
        if c > eps: amps.append([vec_str, c])

    amps = sorted(amps, key=lambda x: x[1], reverse=True)[:truncation]

    # Renormalization
    a_sum = sum([a[1] for a in amps])
    for i in range(len(amps)):
        amps[i][1] = amps[i][1] / a_sum

    return ([x[0] for x in amps], [x[1] for x in amps])

#### Need to set 

# hamiltonian_optmization = False
# backend = statevector_simulator

def output_file(out_path, ai, a, e):
    with open(out_path, "w") as of:
        of.write(f"{e}\n")
        of.write(f"{len(ai)}\n")
        for i, c in zip(ai, a):
            of.write(f"{i}\t{c}\n")


if __name__ == "__main__":
    ifpath = sys.argv[1]
    ofpath = sys.argv[2]
    contents = parse_input_file(ifpath)
    #mo_contents = ao_to_orth_mo(ao_contents)
    contents['tei'] = lst_eri_to_mat(contents['nbasis'], contents['lst_tei'])
    amp_idx, amps, energy = call_vqe(contents)
    output_file(ofpath, amp_idx, amps, energy)
