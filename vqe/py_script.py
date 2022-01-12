from itertools import product
from typing import List, Tuple, Union, final
import numpy as np
from numpy.linalg import svd
import sys, os

from gwonsquantum.chemistry.molecule import electronic_structure
from gwonsquantum.polynomial.transform.bravyi_kitaev import beta_matrix, bravyi_kitaev
from gwonsquantum.polynomial.transform.jordan_wigner import jordan_wigner
from gwonsquantum.chemistry.vqe.coupled_cluster_vqe.qcc import QCC
from gwonsquantum.chemistry.vqe.coupled_cluster_vqe.uccsd import UCCSD
from gwonsquantum.chemistry.molecule import ElectronicStructure
from gwonsquantum.chemistry.molecule.electronic_structure import get_hamiltonian_from_ei
from gwonsquantum.qiskitextension.pauli_measurement import measurement_statevector

DEBUG = False
mapping = "JW"

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
            "n_monomers": 0,
            "n_electrons": 0,
            "nbasis": 0,
            "oei": None,
            "s": None,
            "c": None,
            "lst_tei": None,
            "tei": None,
            "constant": None,
            "mo_energies": None,
            "env_oei" : None,
            "homo" : None,
            "lumo" : None,
            "ent" : None,
        }
        while True:
            line = of.readline()
            if not line:
                break
            if len(line)<=1:
                continue
            if line.startswith("NMONOMER"):
                result['n_monomers'] = int(line.split()[1])
            elif line.startswith('NELEC'):
                result['n_electrons'] = int(line.split()[1])
            elif line.startswith("ENUC"):
                result['constant'] = float(line.split()[1])
            elif line.startswith("NBASIS"):
                result['nbasis'] = int(line.split()[1])
            elif line.startswith("HOMO"):
                result["homo"] = int(line.split()[1])
            elif line.startswith("LUMO"):
                result["lumo"] = int(line.split()[1])
            elif line.startswith("ENT"):
                result["ent"] = int(line.split()[1])
            elif line.startswith('OEI'):
                mode = 1
                linecnt = 0
            elif line.startswith("ENV_OEI"):
                mode = 6
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
            elif mode == 6:
                input_matrix(line, 'env_oei', result, linecnt, triangle=True)
                linecnt += 1
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

def get_spin_spaces(frozen, active, n_basis):
    frozen_spin = list()
    active_spin = list()
    virtual_spin = list()
    for i in frozen:
        frozen_spin += [2*i, 2*i+1]
    for i in active:
        active_spin += [2*i, 2*i+1]
    for i in range(n_basis):
        if i not in frozen and i not in active:
            virtual_spin += [2*i, 2*i+1]
    frozen_spin = sorted(frozen_spin)
    active_spin = sorted(active_spin)
    virtual_spin = sorted(virtual_spin)
    return frozen_spin, active_spin, virtual_spin

def get_spaces(n_electron, homo_cnt, lumo_cnt, n_basis):
    frz_idx = max(n_electron//2 - homo_cnt, 0) if homo_cnt > 0 else 0
    vir_idx = min(n_electron//2 + lumo_cnt, n_basis) if lumo_cnt > 0 else n_basis
    # All lists are asserted to be sorted
    frozen = [x for x in range(frz_idx)]
    active = [x for x in range(frz_idx, vir_idx)]
    virtual = [x for x in range(vir_idx, n_basis)]
    return frozen, active, virtual

def call_vqe(mo_contents):
    # Extract input
    es: ElectronicStructure ={
        "n_electrons": mo_contents['n_electrons'],
        "oei": mo_contents['oei'],
        "tei": mo_contents['tei'],
        "constant": mo_contents['constant']
    }
    n_basis = mo_contents['nbasis']
    n_spin_basis = 2 * n_basis
    n_electron = mo_contents['n_electrons']
    nmonomers = mo_contents['n_monomers']
    env_oei = mo_contents['env_oei']
    n_entang = mo_contents['ent']
    homo_cnt = mo_contents['homo']
    lumo_cnt = mo_contents['lumo']

    # Extract spaces
    frozen, active, virtual = get_spaces(n_electron, homo_cnt, lumo_cnt, n_basis)
    frozen_spin, active_spin, virtual_spin = get_spin_spaces(frozen, active, n_basis)
    n_frozen = len(frozen)
    n_active = len(active)
    n_virtual = len(virtual)
    assert n_frozen+n_active+n_virtual == n_basis

    # Define and run VQE
    vqe_obj = QCC(electronic_structure=es,
    #vqe_obj = UCCSD(electronic_structure=es,
                  mapping=mapping,
                  frozen_indices=frozen,
                  active_indices=active,
                  hamiltonian_optimization=False,
                  num_entanglers=n_entang
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

    # Calculate environmental potential
    es_env: ElectronicStructure = {
        "n_electrons": mo_contents['n_electrons'],
        "oei": env_oei,
        "tei": np.zeros((n_basis, n_basis, n_basis, n_basis), dtype=float),
        "constant": 0.0
    }
    env_fham, n_orb = get_hamiltonian_from_ei(es_env, frozen, active)
    assert n_orb == 2*n_active
    if mapping == "JW":
        env_pham = jordan_wigner(env_fham, 2*n_active)
    elif mapping == "BK":
        env_pham = bravyi_kitaev(env_fham, 2*n_active)
    else:
        raise NotImplementedError
    dv = measurement_statevector(sv, env_pham.co_flipping_set_group) + env_pham.constant()
    assert abs(dv.imag) < 1e-7
    dv = float(dv.real)

    # Get data for postprocessing to density
    amp_idx, amps = get_amp_from_sv(sv, n_basis, frozen_spin, virtual_spin, mapping) # Diagonal terms
    corr_mat_t = get_corr_from_sv(sv, n_basis, frozen, active, virtual, mapping)

    return amp_idx, amps, energy, dv, corr_mat_t

#TODO: Consider electron number breaking.

def get_amp_from_sv(
    final_state:Union[np.ndarray, List[complex]],
    n_basis:int,
    frozen_spin:List[int],
    virtual_spin:List[int],
    mapping="JW")\
    -> Tuple[List[str], float]:

    num_qubits = int(np.log2(len(final_state)))
    truncation = 100
    eps = 1e-7
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

    amp_idx, amp_val = [x[0] for x in amps], [x[1] for x in amps]

    # Add non-active orbitals
    for i, a in enumerate(amp_idx):
        lst_idx = list()
        k=0
        for j in range(n_basis * 2):
            if j in frozen_spin:
                lst_idx.append("1")
            elif j in virtual_spin:
                lst_idx.append("0")
            else:
                lst_idx.append(a[k])
                k += 1
        amp_idx[i] = "".join(lst_idx)
        # amp_idx[i] = "1"*2*n_frozen + a + "0"*2*n_virtual

    return amp_idx, amp_val

def get_corr_from_sv(
     final_state:Union[np.ndarray, List[complex]],
     n_basis:int,
     frozen:List[int],
     active:List[int],
     virtual:List[int],
     mapping="JW")->np.ndarray:

    def _make_mask(_m, _n, _N):
        _lower_mask = int(2 ** _m - 1)
        _middle_mask = int(2 ** (_n - _m - 1) - 1) << _m
        _high_mask = 2 ** (_N - 2) - 1 - _middle_mask - _lower_mask
        return _high_mask, _middle_mask, _lower_mask


    def _recover_idx(_m, _n, _N, _i):
        _hm, _mm, _lm = _make_mask(_m, _n, _N)
        _i10 = (_i & _lm) ^ (1 << _m) ^ ((_i & _mm) << 1) ^ ((_i & _hm) << 2)
        _i01 = _i10 ^ (1 << _m) ^ (1 << _n)
        return _i10, _i01

    if mapping != "JW":
        raise NotImplementedError
    N = int(np.log2(len(final_state)))
    corr_mat = np.zeros((N//2, N//2))
    for m, n in product(range(N//2), repeat=2):
        if n<=m:
            continue
        tmp_val = 0
        for i in range(len(final_state)//4):
            i_10, j_01 = _recover_idx(2*m, 2*n, N, i) # I(2m)=1, I(2n)=0, J(2m)=0, J(2n)=1
            tmp_val += 2*(final_state[i_10].conjugate() * final_state[j_01]).real
            i_10, j_01 = _recover_idx(2*m+1, 2*n+1, N, i) # I(2m+1)=1, I(2n+1)=0, J(2m1+)=0, J(2n+1)=1
            tmp_val += 2*(final_state[i_10].conjugate() * final_state[j_01]).real
        corr_mat[m][n] = tmp_val

    corr_mat_full_t = np.zeros((n_basis, n_basis))
    for m, n in product(range(n_basis), repeat=2):
        if m in frozen or m in virtual:
            continue
        if n in frozen or n in virtual:
            continue
        corr_mat_full_t[n][m] = corr_mat[active.index(m)][active.index(n)]
    return corr_mat_full_t

#### Need to set 

# hamiltonian_optmization = False
# backend = statevector_simulator

def output_file(out_path, ai, a, e, dv, corr_mat:np.ndarray):
    try:
        tmp_f = open(out_path, "w")
        tmp_f.close()
    except FileNotFoundError:
        dir_path = '/'.join(out_path.split('/')[:-1])
        os.mkdir(dir_path)
    with open(out_path, "w") as of:
        of.write(f"{e}\n")
        of.write(f"{dv}\n")
        of.write(f"{len(ai)}\n")
        for i, c in zip(ai, a):
            of.write(f"{i}\t{c}\n")
        size_corr = corr_mat.shape[0]
        for n in range(0, size_corr):
            for m in range(0, size_corr):
                if m>n:
                    break
                of.write(f"{corr_mat[m][n]}\t")
            of.write('\n')


if __name__ == "__main__":
    ifpath = sys.argv[1]
    ofpath = sys.argv[2]
    # print(f"pyscript : VQE for {ifpath} spawned.")
    assert len(sys.argv) == 3
    contents = parse_input_file(ifpath)
    #mo_contents = ao_to_orth_mo(ao_contents)
    contents['tei'] = lst_eri_to_mat(contents['nbasis'], contents['lst_tei'])
    amp_idx, amps, energy, dv, corr_mat = call_vqe(contents)
    output_file(ofpath, amp_idx, amps, energy, dv, corr_mat)
    # print(f"pyscript : output {ofpath} Generated.")
