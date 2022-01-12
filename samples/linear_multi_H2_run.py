from itertools import product
from math import sqrt
import numpy as np
from typing import List, Optional

def linear_geometry(atoms: List[str],
                    distances):
    """Return geometry of atoms for linear geometry.

    Args:
        atoms: list of atoms with length of n. example) ['H', 'Be', 'H']
        distances: list of interatomic distances with length of n-1. example) [1.5, 1.5] or 1.5

    Returns:
        geometry: geometry fed to the ``Molecule`` object initialization.

    """

    if len(distances) != len(atoms) - 1:
        raise ValueError
    position = [sum(distances[:i]) for i in range(len(atoms))]
    return [[a, (0, 0, p)] for a, p in zip(atoms, position)]


def biatomic_molecule_multimer_geometry(atom, n, d1, d2):
    atom = [atom] * (2 * n)
    dst = [d1 if i % 2 == 0 else d2 for i in range(2 * n - 1)]
    geo = linear_geometry(atom, dst)
    return geo


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
 
 $end
 $FMOHYB
sto-3g 5 5
1 0  -0.117784    0.542251    0.000000    0.000000    0.850774
0 1  -0.117787    0.542269    0.802107    0.000000   -0.283586
0 1  -0.117787    0.542269   -0.401054   -0.694646   -0.283586
0 1  -0.117787    0.542269   -0.401054    0.694646   -0.283586
0 1   1.003621   -0.015003    0.000000    0.000000    0.000000
 $END
'''.format(nfrag, icharge_str, indat1_str, mult1_str, respap, resppc, resdim, geometry_str)
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


if __name__ == "__main__":
    nfrag = 6
    d1 = 0.7237340367
    # d2 = 2.95
    for d2 in [0.72374]+list(np.arange(0.50, 4.51, 0.25)):
        d_str = str(d2)[:5].replace('.','_')
        desc = f"linear_H2_6__{d_str}"
        geo = biatomic_molecule_multimer_geometry("H", nfrag, d1, d2)
        fragmentation = [i // 2 + 1 for i in range(2 * nfrag)]
        input_data = gen_inputs(nfrag, fragmentation, geo, desc, rhf=True)
        with open(f"./linear_H2_6/linear_H2_6_{d_str}_rhf.inp", 'w') as of:
            of.write(input_data)
    
    
