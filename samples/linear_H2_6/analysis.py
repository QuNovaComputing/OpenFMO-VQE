import sys, os


if __name__ == '__main__':
    log_dir = sys.argv[1]
    prefix = "linear_H2_6_"

    fmo_rhf = list()
    fmo_vqe = list()

    for in_file in sorted(os.listdir(log_dir)):
        if in_file.split('.')[-1] != 'log':
            continue
        if not in_file.startswith(prefix):
            continue
        rhf = 'rhf' in in_file
        param = in_file[len(prefix):-1*len('.log')].replace('_', '.').replace('.rhf', '').ljust(5, '0')
        
        full_infile = os.path.join(log_dir, in_file)
        converged = False
        value = None
        with open(full_infile, 'r') as of:
            for l in of.readlines():
                if '==== SCC converged ====' in l:
                    converged = True
                if 'total energy(FMO) =' in l:
                    value = float(l.split()[-1])
        if rhf:
            fmo_rhf.append((param, value, converged))
        else:
            fmo_vqe.append((param, value, converged))
    
    fmo_rhf = sorted(fmo_rhf, key=lambda x : x[0])
    fmo_vqe = sorted(fmo_vqe, key=lambda x : x[0])

    print("RHF")
    for k in fmo_rhf:
        print(str(k)+',')
    
    print("VQE")
    for k in fmo_vqe:
        print(str(k)+',')
