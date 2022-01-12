import sys
import os
import subprocess

if __name__ =="__main__":
    inp_dir = sys.argv[1]
    # nfrag = 6
    # np = 10
    buf = 2048
    sequential = False

    job_list = list()
    of_list = list()
    
    for inp_file in sorted(os.listdir(inp_dir)):
        # if 'rhf' not in inp_file:
        #     # Check -q option when you use rhf
        #     continue
        if inp_file.split('.')[-1] != 'inp':
            continue
        full_inp_file = os.path.join(inp_dir, inp_file)
        with open(full_inp_file, 'r') as of:
            for l in of.readlines():
                if 'nfrag' in l:
                    ng = int(l.split("nfrag")[-1].replace("=", "").split(' ')[0])
                    break
            np = ng + 3

        cmd = ["./ofmo-master", "-ng", str(ng), "-np", str(np), "-q", "-B", str(buf)]
        full_oup_file = '.'.join(full_inp_file.split('.')[:-1]+["log"])
        cmd_tail = [full_inp_file]
        print(full_inp_file+" > "+full_oup_file)
        of = open(full_oup_file, "w")
        of_list.append(of)
        if sequential:
            subprocess.run(cmd + cmd_tail, stdout=of)
            of.close()
        else:
            p = subprocess.Popen(cmd + cmd_tail, stdout=of)
            job_list.append(p)
    while(not sequential):
        if len(job_list) == 0:
            _ = [of.close() for of in of_list]
            break
        _job_list = list(job_list)
        _of_list = list(of_list)
        for f, p in zip(_of_list, _job_list):
            if p.poll() is not None:
                # Terminated
                f.close()
                print(f"{f.name} terminated.")
                of_list.remove(f)
                job_list.remove(p)
