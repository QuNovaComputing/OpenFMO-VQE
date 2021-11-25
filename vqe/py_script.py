from itertools import product
from typing import List, Tuple, Union
import numpy as np
from numpy.linalg import svd
import sys

from gwonsquantum.polynomial.transform.bravyi_kitaev import beta_matrix
from gwonsquantum.chemistry.vqe import QCC

def parse_input_file(inp_file):
    def input_matrix(line, key, result, cnt):
        _nbasis = len(line.split())
        if result['nbasis'] == 0 :
            result['nbasis'] = _nbasis
        else :
            assert _nbasis == result['nbasis']
        if result[key] is None:
            result[key] = np.zeros((_nbasis, _nbasis), dtype=float)

        for i, v in enumerate(line.split()):
            result[key] [cnt][i] = float(v)


    with open(inp_file, 'r'):
        linecnt = 0
        mode = 0
        result = {
            "nelec": 0,
            "nbasis": 0,
            "oei": None,
            "s": None,
            "c": None,
            "mo_eri": None,
            "enuc": None,
        }
        while True:
            line = inp_file.readline()
            if not line:
                break
            if line.startswith('NELEC'):
                result['nelec'] = int(line.split()[1])
            elif line.startswith("ENUC"):
                result['enuc'] = float(line.split()[1])
            elif line.startswith('OEI'):
                mode = 1
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
                input_matrix(line, 'oei', result, linecnt)
            elif mode == 2:
                input_matrix(line, 's', result, linecnt)
            elif mode == 3:
                input_matrix(line, 'c', result, linecnt)
            elif mode == 4:
                spt = line.split()
                assert len(spt) == 5
                i,j,k,l = int(spt[0]), int(spt[1]), int(spt[2]), int(spt[3])
                result['eri'].append(((i, j, k, l), float(spt[4])))
    return result

def lst_eri_to_mat(nbasis, lst_eri):
    mat_eri = [[[[0.0 for _ in range(nbasis)]
                      for _ in range(nbasis)]
                      for _ in range(nbasis)]
                      for _ in range(nbasis)]
    for (i,j,k,l), v in lst_eri:
        mat_eri[i][j][k][l] = v
        mat_eri[i][j][l][k] = v
        mat_eri[j][i][k][l] = v
        mat_eri[j][i][l][k] = v
        
        mat_eri[k][l][i][j] = v
        mat_eri[k][l][j][i] = v
        mat_eri[l][k][i][j] = v
        mat_eri[l][k][j][i] = v

    return mat_eri
'''
def generate_mo_eri(nbasis, coeff, aoeri):
    mat_eri = lst_eri_to_mat(nbasis, aoeri)
    mo_eri = list()
    for i in range(nbasis):
        for j in range(i+1):
            for k in range(i+1):
                for l in range(k+1):
                    for a,b,c,d in product(range(nbasis), repeat=4):
                        c = coeff[i][a]*coeff[j][b]*coeff[k][c]*coeff[l][d]
                        c*= mat_eri[a][b][c][d]
                        mo_eri.append(((i,j,k,l), c))
    return mo_eri
'''

# symmetric orthogonalization
def trans_sym(S):
    U, D, U_dag = svd(S)
    D_inv_sqrt = np.sqrt(np.reciprocal(D))
    D_inv_sqrt_mat = np.diag(D_inv_sqrt)
    X = np.dot(np.dot(U, D_inv_sqrt_mat), U_dag)
    return X

# canonical orthogonalization
def trans_can(S):
    U, D, U_dag = svd(S)
    D_inv_sqrt = np.sqrt(np.reciprocal(D))
    D_inv_sqrt_mat = np.diag(D_inv_sqrt)
    X = np.dot(U, D_inv_sqrt_mat)
    return X

def ao_to_orth_mo(ao_contents):
    nelec = ao_contents['nelec']
    nbasis = ao_contents['nbasis']
    ao_oei = ao_contents['oei']
    ao_ovlp = ao_contents['s']
    coeff = ao_contents['c']
    mo_eri = ao_contents['mo_eri']
    enuc = ao_contents['enuc']

    mo_oei = np.zeros((nbasis, nbasis), dtype=float)
    mo_overlap = np.zeros((nbasis, nbasis), dtype=float)
    # Assumming MO_i = sum C_ii' AO_i'
    for i, j in product(range(nbasis), repeat=2):
        mo_oei[i][j] = sum([ao_oei[_i][_j]*coeff[i][_i]*coeff[j][_j]
        for _i, _j in product(range(nbasis), repeat=2)])
    for i, j in product(range(nbasis), repeat=2):
        mo_overlap[i][j] = sum([ao_ovlp[_i][_j]*coeff[i][_i]*coeff[j][_j]
        for _i, _j in product(range(nbasis), repeat=2)])
    # mo_eri = generate_mo_eri(nbasis, coeff, ao_eri)
    
    # TODO:Need orthogonalization


    ret = {
        "n_electron": nelec,
        "nbasis": nbasis,
        "oei": mo_oei,
        "tei": mo_eri,
        "enuc": enuc
    }
    return ret

def call_vqe(mo_contents):
    vqe_obj = QCC()
    # TODO:Finish this

def get_amp_from_sv(
    final_state:Union[np.ndarray, List[complex]],
    mapping="JW")\
    -> Tuple[List[str], float]:

    num_qubits = int(np.log2(len(final_state)))
    truncation = 10
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
        vec_str = "".join([str(i) for i in vec[::-1]])  #MSB first
        amps.append((vec_str, abs(c)**2))

    amps = sorted(amps, key=lambda x: x[1], reverse=True)[:truncation]

    # Renormalization
    a_sum = sum([a[1] for a in amps])
    for i in range(len(amps)):
        amps[i][1] = amps[i][1] / a_sum

    return [x[0] for x in amps], [x[1] for x in amps]

#### Need to set 

# hamiltonian_optmization = False
# backend = statevector_simulator


if __name__ == "__main__":
    ifpath = sys.argv[1]
    ofpath = sys.argv[2]
    ao_contents = parse_input_file(ifpath)
    mo_contents = ao_to_orth_mo(ao_contents)
    #TODO: Finish here
