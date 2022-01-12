from itertools import product
from math import sqrt
import numpy as np
from typing import List, Optional
from pyscf import gto

def euclidian_dist(x1, y1, z1, x2, y2, z2):
    return sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


def max_dist(geometry):
    dst = list()
    for i, j in product(range(len(geometry)), repeat=2):
        if i == j:
            pass
        new_tuple = *geometry[i][1], *geometry[j][1]
        dst.append(euclidian_dist(*new_tuple))
    return max(dst)


def geometry_to_str(geometry):
    lst_str = list()
    for atom in geometry:
        at = atom[0]
        x, y, z = atom[1]
        x, y, z = "{:.8f}".format(x), "{:.8f}".format(y), "{:.8f}".format(z)
        lst_str.append(f" {at}   {at}   {x}   {y}   {z}")
    return '\n'.join(lst_str)


def gen_inputs(nfrag: int,
               indat1: List[int],
               geometry,
               desc: str,
               rhf: bool=False,
               icharge: Optional[List[int]] = None,
               mult1: Optional[List[int]] = None,
               respap: Optional[float] = None,
               resppc: Optional[float] = None,
               resdim: Optional[float] = None,
               mono_homo: Optional[int] = None,
               mono_lumo: Optional[int] = None,
               dim_homo: Optional[int] = None,
               dim_lumo: Optional[int] = None,
               mono_ent: Optional[int] = None,
               dim_ent: Optional[int] = None,
               fmobnd: Optional[str] = None
               ):
    if icharge is None:
        icharge = [0] * nfrag
    if mult1 is None:
        mult1 = [1] * nfrag
    max_dst = max_dist(geometry)
    if respap is None:
        respap = round(max_dst) + 1.5
    if resppc is None:
        resppc = round(max_dst) + 1.5
    if resdim is None:
        resdim = round(max_dst) + 1.5
    if mono_homo is None:
        mono_homo = -1
    if mono_lumo is None:
        mono_lumo = -1
    if dim_homo is None:
        dim_homo = -1
    if dim_lumo is None:
        dim_lumo = -1
    if mono_ent is None:
        mono_ent = -1
    if dim_ent is None:
        dim_ent = -1
    icharge_str = ",".join([str(x) for x in icharge])
    mult1_str = ",".join([str(x) for x in mult1])
    indat1_str = ",".join([str(x) for x in indat1])
    geometry_str = geometry_to_str(geometry)
    if fmobnd is None:
        fmobnd = ""
    form = '''
 $gddi niogroup=1 nioprocs=1 $end
 $fmo nfrag={0} ICHARG(1)={1}
 indat(1)={2} mult(1)={3}
 respap={4} resppc={5} resdim={6} $end
 $basis gbasis=sto ngauss=3 $end
 $fmoxyz
{7}
 $end
 $fmoprp
modorb=3
 $end
 $fmobnd
{8}
 $end
 $FMOHYB
sto-3g 5 5
1 0  -0.117784    0.542251    0.000000    0.000000    0.850774
0 1  -0.117787    0.542269    0.802107    0.000000   -0.283586
0 1  -0.117787    0.542269   -0.401054   -0.694646   -0.283586
0 1  -0.117787    0.542269   -0.401054    0.694646   -0.283586
0 1   1.003621   -0.015003    0.000000    0.000000    0.000000
 $END
'''.format(nfrag, icharge_str, indat1_str, mult1_str, respap, resppc, resdim, geometry_str, fmobnd)
    tail_vqe = '''
$contrl method=vqe vqescr=./vqe/py_script.py desc={0} $end
$vqeprp monhomo={1} monlumo={2} dimhomo={3} dimlumo={4} monent={5} diment={6} $end
    '''.format(desc, mono_homo, mono_lumo, dim_homo, dim_lumo, mono_ent, dim_ent)
    tail_rhf = '''
$contrl method=rhf $end
    '''
    if rhf:
        return form +'\n'+tail_rhf
    else:
        return form +'\n'+tail_vqe


def build_alkane(n):
    ang_CCC = 109.5
    ang_CCH = 109.5
    ang_dihedral = 120
    len_CC = 1.533350
    len_CH = 1.093511
    coord_list = list()
    at_list = list()
    for i in range(1, n + 1):
        if i == 1:
            coord_list.append("C")
            at_list.append("C")
        elif i == 2:
            coord_list.append(f"C {i - 1} {len_CC}")
            at_list.append("C")
        elif i == 3:
            coord_list.append(f"C {i - 1} {len_CC} {i - 2} {ang_CCC}")
            at_list.append("C")
        else:
            coord_list.append(f"C {i - 1} {len_CC} {i - 2} {ang_CCC} {i - 3} {180}")
            at_list.append("C")

    for i in range(1, n + 1):
        '''
        if i+2 <= n:
            coord = f"H {i} {len_CH} {i+1} {ang_CCH} {i+2} {ang_dihedral}"
        else:
            coord = f"H {i} {len_CH} {i-1} {ang_CCH} {i-2} {ang_dihedral}"
        n_append = 2 if i not in [1, n] else 3
        coord_list += [coord] * n_append
        at_list += ["H"] * n_append
        '''
        if i == 1:
            coord_list.append(f"H {i} {len_CH} {i + 1} {ang_CCH} {i + 2} {120}")
            at_list.append("H")
            coord_list.append(f"H {i} {len_CH} {i + 1} {ang_CCH} {i + 2} {0}")
            at_list.append("H")
            coord_list.append(f"H {i} {len_CH} {i + 1} {ang_CCH} {i + 2} {-120}")
            at_list.append("H")
        elif i == n:
            coord_list.append(f"H {i} {len_CH} {i - 1} {ang_CCH} {i - 2} {120}")
            at_list.append("H")
            coord_list.append(f"H {i} {len_CH} {i - 1} {ang_CCH} {i - 2} {0}")
            at_list.append("H")
            coord_list.append(f"H {i} {len_CH} {i - 1} {ang_CCH} {i - 2} {-120}")
            at_list.append("H")
        else:
            if i + 2 <= n:
                coord_list.append(f"H {i} {len_CH} {i + 1} {ang_CCH} {i + 2} {120}")
                at_list.append("H")
                coord_list.append(f"H {i} {len_CH} {i + 1} {ang_CCH} {i + 2} {-120}")
                at_list.append("H")
            else:
                coord_list.append(f"H {i} {len_CH} {i - 1} {ang_CCH} {i - 2} {120}")
                at_list.append("H")
                coord_list.append(f"H {i} {len_CH} {i - 1} {ang_CCH} {i - 2} {-120}")
                at_list.append("H")

    mol = gto.Mole()
    mol.atom = "\n".join(coord_list)
    mol.build()

    return [[at, tuple(g)] for at, g in zip(at_list, mol.atom_coords("ang"))]


if __name__ == "__main__":
    mon_lumo = 2
    mon_homo = 2
    mon_ent = 5
    dim_lumo = 4
    dim_homo = 4
    dim_ent = 10
    
    for num_c in range(3, 10):
        desc = f"C{num_c}H{num_c*2+2}_{mon_homo}-{mon_lumo}-{mon_ent}_{dim_homo}-{dim_lumo}-{dim_ent}"
        geo = build_alkane(num_c)
        fragmentation = [i for i in range(1, num_c+1)]
        for i in range(1, num_c+1):
            fragmentation += [i]*(3 if i in [1, num_c] else 2)
        fmobnd = '\n'.join([f"-{i} {i+1} sto-3g" for i in range(1, num_c)])
        input_data = gen_inputs(num_c, fragmentation, geo, desc,
                                mono_homo=mon_homo, mono_lumo=mon_lumo, mono_ent=mon_ent,
                                dim_homo=dim_homo, dim_lumo=dim_lumo, dim_ent=dim_ent,
                                fmobnd=fmobnd, 
                                rhf=False)
        with open(f"./{desc}.inp", 'w') as of:
            of.write(input_data)
    
    