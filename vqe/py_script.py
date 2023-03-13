from itertools import product
from typing import List, Tuple, Union, final
import numpy as np
from numpy.linalg import svd
import sys, os

import json

from qunova_vqe import call_vqe, establish_pulsar_connection

establish_pulsar_connection("username", "api_key")

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
            "moons": None,
            "ansatz": None,
            "threshold": None
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
            elif line.startswith("MOONS"):
                mode = 7
                linecnt = 0
            elif line.startswith("ANSATZ"):
                result["ansatz"] = line.split()[1]
            elif line.startswith("THRESHOLD"):
                result["threshold"] = line.split()[1]
            elif line.startswith("DENSITY"):
                mode = 8
                linecnt = 0
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
            elif mode == 7:
                spt = line.split()
                if result['moons'] is None:
                    result['moons'] = list()
                result['moons'].append(float(spt[0]))
            else:
                pass

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
