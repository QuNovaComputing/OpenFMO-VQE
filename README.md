GPU-accelerated OpenFMO
=======================

**GPU-accelerated OpenFMO version 1.0 relaesed. (2018-04-01)**

**OpenFMO official site:** http://www.openfmo.org/ ;  

What is OpenFMO?
----------------

**Open-architecture Implementation of 
Fragment Molecular Orbital Method
for Peta-scale Computing (OpenFMO)** 
is an open-architecture program
targeting the effective FMO calculations 
on massive parallel computers,
now including post-Peta Scale systems.


Fragment molecular orbital (FMO) 
method is a representative method to solve 
the electronic structures of large bio-molecules
including protein, DNA, sugar chain, and so on.


OpenFMO was written from scratch
by Inadomi and co-workers 
in C-language (ca. 54,000 lines)
with OpneMP and MPI hybrid programming
model.
MPI dynamic process management
is used for OpenFMO program
on the basis of the master-worker
execution model involving
master process, worker groups, and
data server.


In addition to FMO calculations,
the users can do conventional restricted Hartree-Fock (RHF) 
calculations using the "skeleton-RHF" code of OpenFMO,
which is also OpenMP and MPI hybrid program.


**In this released version OpenFMO 1.0,
GPU-accelerated FMO-RHF and RHF calculations can be performed
on GPU cluster.** 


Capabilities and Limitations
----------------------------

* Single-point Ground-state Energy Calculation

* RHF with *skeleton-RHF* Code

* FMO2-RHF (FMO2: FMO Method with Two-body Correction)

* Minimum and Double-zeta Gaussian Basis Functions; Up to Third-row Atoms (namely, H - Ar) 

  * STO-3G, 6-31G, 6-31G(d), 6-31G(d,p)

* MPI + OpenMP Parallelization for RHF and FMO2-RHF

* GPU-accelerated RHF and FMO-RHF with Fermi or Kepler microarchitecture
  supporting double-precision floating-point operations

Citing OpenFMO
--------------
 
1. **OpenFMO**

   Please cite the following references 
   for any publication including the scientific works
   obtained from OpenFMO application.

   * Takami, T.; Maki, J.; Ooba, J.; Inadomi, Y.; Honda, H.; Kobayashi, T.; 
     Nogita, R.; Aoyagi, M.; 
     Open-Architecture Implementation of Fragment Molecular Orbital Method 
     for Peta-Scale Computing.
     *CoRR*, 2007, URL: http://arxiv.org/abs/cs/0701075

   * Inadomi, Y.;, Maki, J.; Honda, H.; Takami, T.; Kobayashi T.;, 
     Aoyagi, M.; Minami, K., Performance Tuning of Parallel Fragment 
     Molecular Orbital Program (OpenFmo) for Effective Execution 
     on k-Computer. *J. Comput. Chem., Jpn.*, **12** (2):145–155, 2013. 
     doi:10.2477/jccj.2012-0031.

2. **OpenFMO with Falanx Middleware**

   In addition, please cite the following reference for 
   any publication including the scientific works obtained from 
   OpenFMO application with Falanx middleware.

   * Takefusa, A; Ikegami, T; Nakada, H.; Takano, R.; Tozawa, T.; and 
     Tanaka, Y. Scalable and Highly Available Fault Resilient 
     Programming Middleware for Exascale Computing. 
     The International Conference for High Performance Computing, 
     Networking, Storage and Analysis: SC14. 
     URL: http://sc14.supercomputing.org/sites/all/themes/sc14/files/archive/tech_poster/poster_files/post217s2-file2.pdf.

3. **GPU-accelerated RHF**

   Please cite the following reference 
   when publishing the results obtained from
   the GPU-accelerated *skeleton RHF* code in OpenFMO.

   * Umeda, H.; Hanawa, T.; Shoji, M.; Boku, T.;  Shigeta, Y., 
     Performance Benchmark of FMO Calculation with 
     GPU-accelerated Fock Matrix Preparation Routine. 
     *J. Comput. Chem., Jpn.*, **13** (6):323–324, 
     2015. doi:10.2477/jccj.2014-0053.

   * Umeda, H.; Hanawa, T.; Shoji, M.; Boku, T.; Shigeta, Y., 
     Large-Scale MO Calculation with GPU-accelerated FMO Program. 
     The International Conference for High Performance Computing, 
     Networking, Storage and Analysis: SC15. 
     URL: http://sc15.supercomputing.org/sites/all/themes/SC15images/tech_poster/tech_poster_pages/post198.html.


4. **GPU-accelerated FMO-RHF**

   The following reference should
   be cited when publishing the results utilizing 
   the GPU-accelerated FMO-RHF code in OpenFMO.

   * Umeda, H.; Hanawa, T.; Shoji, M.; Boku, T.; Shigeta, Y., 
     GPU-accelerated FMO Calculation with OpenFMO: Four-center 
     Inter-fragment Coulomb Interaction. 
     *J. Comput. Chem., Jpn.*, **14** (3):69–70, 2015. 
     doi:10.2477/jccj.2015-0041.

   * Umeda, H.; Hanawa, T.; Shoji, M.; Boku, T.; Shigeta, Y., 
     Large-Scale MO Calculation with GPU-accelerated FMO Program. 
     The International Conference for High Performance Computing, 
     Networking, Storage and Analysis: SC15. 
     URL: http://sc15.supercomputing.org/sites/all/themes/SC15images/tech_poster/tech_poster_pages/post198.html.


Licenses
--------

1. OpenFMO:

   **Copyright (c) 2007 Yuichi Inadomi, Toshiya Takami, Hiroaki Honda, 
   Jun Maki, and Mutsumi Aoyagi**

   Released under the MIT license

   http://opensource.org/licenses/mit-license.php

   The copy of the license is also included in the distribution 
   as "LICENSE_OpenFMO".

2. GPGPU parts:

   **Copyright (c) 2013 Hiroaki Umeda, Yuichi Inadomi, Toshihiro Hanawa, Mitsuo Shoji, Taisuke Boku, and Yasuteru Shigeta**

   Released under the MIT license

   http://opensource.org/licenses/mit-license.php

   The copy of the license is also included in the distribution 
   as "LICENSE_GPGPU_parts".

3. Falanx Middleware parts:

   The Falanx middleware is copyrighted by National Institute Advanced 
   Industrial Science and Technology (AIST), and is licensed under the 
   Apache License, Version 2.0. You may obtain a copy of the License at
   
   http://www.apache.org/licenses/LICENSE-2.0

   which is also included in the distribution as "LICENSE_Falanx".

4. Released Version 1.0:

   **Copyright (c) 2018 Hiroaki Umeda, Yuichi Inadomi, Hirotaka Kitoh-Nishioka, and Yasuteru Shigeta**

   Released under the MIT license

   http://opensource.org/licenses/mit-license.php

   The copy of the license is also included in the distribution 
   as "LICENSE_Released_Version_1.0".

Contact Person
--------------

Hirotaka KITOH-NISHIOKA Dr.

| JST PRESTO Researcher  
| Japan Science and Technology Agency (JST)  
| CCS Collaborative Fellow  
| Center for Computational Sciences (CCS),  
| University of Tsukuba  
| 1-1-1 Tennodai, Tsukuba, Ibaraki, JAPAN  
| E-mail: hkito{at}ccs.tsukuba.ac.jp  

At Biological Function and Information Group
of Professor Shigeta in Center for Computational Sciences, 
University of Tsukuba;
https://www.ccs.tsukuba.ac.jp/eng/research-divisions/division-of-life-sciences/biological-function-and-information-group


Acknowledgments
---------------
JST-CREST: "Development of System Software Technologies for
post-Peta Scale High Performance Computing"

Prerequisites
-------------

* LINUX/UNIX Cluster Machines
* GNU C Compiler
* Intel C Compiler
* MPI Libraries (Default: Intel MPI Library) supproting MPI_Comm_spawn functions.
* Intel MKL(Math Kernel Library)

In addition, **GPU-accelerated OpenFMO** requires:

* NVIDIA Graphics card (Fermi or Kepler microarchitecture)
  supporting double precision floating point operations
* NVIDIA drivers for GPU

How to Get
----------

OpenFMO program is available through
the repositories hosted on
github ( https://github.com/OpenFMO/OpenFMO ).

To check out the latest OpenFMO sources:

   `$ git clone https://github.com/OpenFMO/OpenFMO.git OpenFMO`

How to Compile
--------------
After checking out the release archive of OpenFMO,
you should move to its top directory:


   `$ cd OpenFMO`
   
Makefile located in the top directory
is used to build the OpenFMO executables.
By typing make command with the "help" target option,
the usage of the make command is printed:

   `$ make help`

    usage: make <target>
    original    build ofmo-master, ofmo-worker, ofmo-mserv.
    falanx      build ofmo-falanx.
    rhf         build ofmo-rhf.
    clean       remove build files.
    help        print this message.


In line with the printed explanation,
the following command yields
the "skeleton-RHF" executable, *ofmo-rhf* ,
in the top directory:

   `$ make rhf`



The following command yields
the three executables, *ofmo-master* , *ofmo-master* ,
and *ofmo-mserv* , in the top directory,
which run the FMO calculations
with the master-worker execute model.

   `$ make original`


If it is difficult to run with MPI_Comm_Spawn
for your system, you can use Falanx programming middleware 
( https://sites.google.com/site/spfalanx/ )


To build **GPU-accelerated OpenFMO** executables,
one should modify the following parts of Makefile;

   code-block:: Makefile

     xcCUDAA = KEPLER
     xcCUDAA = FERMI
     xcCUDA = 0


Default value of *xcCUDA* variable
is set to zero, which turns off nvcc compilation.
To build the codes with nvcc for the Fermi microarchitecture,
Makefile should be modified as follows;
  
   code-block:: Makefile

     #xcCUDAA = KEPLER
     xcCUDAA = FERMI
     #xcCUDA = 0


Similarly, for the Kepler microarchitecture,
Makefile should be modified as follows;

   code-block:: Makefile

     xcCUDAA = KEPLER
     #xcCUDAA = FERMI
     #xcCUDA = 0


For getting optimal performance on your system,
you may change dim2e[][] array in cuda/cuda-integ.cu.



Executing OpenFMO
-----------------

1. Setup


   * Set **OMP_NUM_THREADS** environment variable:

        code-block:: tcsh

           # csh, tcsh:
           setenv OMP_NUM_THREADS (Number_of_Threads)

        code-block:: bash

           # sh, bash:
           export OMP_NUM_THREADS=(Number_of_Threads)

   * Set **LIBRARY_PATH** environment variable:

        code-block:: tcsh

           # csh, tcsh: add to shell
           setenv LIBRARY_PATH $LD_LIBRARY_PATH

        code-block:: bash

           # sh, bash: add to shell
           export LIBRARY_PATH=$LD_LIBRARY_PATH

   * Set **OFMOPATH** environment variable that points 
     the directory storing the OpenFMO executables 
     ( *ofmo-master* , *ofmo-worker* , and *ofmo-mserv* ) 
     or "skeleton-RHF" executable ( *skel-rhf* ), 
     which is usually the directory where you compile OpenFMO 
     programs; you have to tell OpenFMO the path
     to this directory using the execution option of *ofmo-master* 
     with *-bindir* .

        code-block:: tcsh

           # csh, tcsh: 
           setenv OFMOPATH /OpenFMO/executables/install/directory

        code-block:: bash

           # sh, bash:
           export OFMOPATH=/OpenFMO/executables/install/directory

   * Set **SCRDIR** environment variable that points the directory 
     storing the temporary "scratch" files 
     for OpenFMO executables; you have to tell OpenFMO the path
     to this directory using the execution option of *ofmo-master* 
     with *-scrdir* . 

        code-block:: tcsh

           # csh, tcsh: 
           setenv SCRDIR /Pass/To/OpenFMO/Scratch/Files/Directory

        code-block:: bash

           # sh, bash:
           export SCRDIR=/Pass/To/OpenFMO/Scratch/Files/Directory

   * Prepare the input files.



2. Command Line Options


   By running the "skeleton-RHF" program, *skel-rhf* , 
   with a help command-line argument, *-h* ,
   its usage is printed:


     `$ ${OFMOPATH}/skel-rhf -h`

       Usage: skel-rhf [-snvh][-B buffer] input [density]
         -B buf: # buffer size (MB, default: 0)
         -s: sync
         -n: dryrun
         -v: verbose
         -h: show this help
       Options for GPGPU:
         -d ndev: # devices (default:0)

   Similarly, by running the OpenFMO program, *ofmo-master* ,
   with a help command-line argument, *-h* ,
   you can see some of its command-line arguments:


     `$ ${OFMOPATH}/ofmo-master -h`

       Usage: ofmo-master [options] [input [InitDens]]
         -ng #: # groups
         -np #: # total MPI procs
          -B #: buffer size / proc (MB, default: 512)
            -v: verbose
            -h: show this help
       Options for GPGPU:
          -d #: # devices (default:0)

   **Note that OpenFMO should be invoked with 
   -ng and -np command-line arguments.**

   The details of the command-line arguments to *ofmo-master* 
   are listed in Table 1 at OpenFMO official site ( http://www.openfmo.org/ ). 


3. Multi-thread Execution of "skeleton-RHF"


   1. First set **OMP_NUM_THREADS** environment variable (see :ref:`Setup`).

   2. Next set **OFMOPATH** environment variable (see :ref:`Setup`). 
      Then, execute the "skeleton-RHF" program within a single cluster node:


      `$ ${OFMOPATH}/skel-rhf Input_File_Name > Log_File_Name`

   3. Execute the GPGPU-accelerated "skeleton-RHF" program within a single 
      cluster node:

      `$ ${OFMOPATH}/skel-rhf -d 1 Input_File_Name > Log_File_Name`


4. Hybrid Execution of "skeleton-RHF"


   1. First set **OMP_NUM_THREADS** , **OFMOPATH** ,
      and  **SCRDIR** environment variables.

   2. Execute the "skeleton-RHF" program with *N* MPI processes:

      `$ mpiexec.hydra -np N ${OFMOPATH}/skel-rhf Input_File_Name > Log_File_Name`

   3. To perform GPGPU-accelerated RHF/RKS calculations with *N* MPI processes:
   
      `$ mpiexec.hydra -np N ${OFMOPATH}/skel-rhf -d 1 Input_File_Name > Log_File_Name`


5. Execution of OpenFMO


   Here, we demonstrate an example of
   the way OpenFMO is executed by using the 
   cluster including the GPU nodes; 
   one node is comprised of Intel Xeon E-5-2680 
   (2.6 GHz 8 cores) 2 CPUs
   and NVIDIA Tesla M2090 (Fermi) 4 units.


   If the GPU-accelerated FMO-RHF calculation
   is performed using 8 nodes (2 x 8 x 8 = 128 cores and 4 x 8 = 32 GPU units ) 
   with 1 data server of 1 rank and 15 worker groups of 2 ranks,
   you should run OpenFMO as follows:

   1. First, set up  **OFMOPATH** and  **SCRDIR** environment variables.

   2. **OMP_NUM_THREADS** should be set to 4 (128 cores / 32 MPI processes).

   3. Then, execute *ofmo-master* with the proper 
      command-line arguments:

      `$ mpiexec.hydra -np 1 ${OFMOPATH}/ofmo-master -np 32 -ng 15 -d 1 -bindir ${OFMOPATH} -scrdir ${SCRDIR} Input_File_Name > Log_File_Name`

      The master openMP thread of each MPI rank controls one GPU unit.
      Therefore, **you had better make the total number of MPI 
      processes be equal to that of the available GPU units** 
      (32 in the above case) in order to bring out the GPU's maximum 
      performance on your cluster machines.


6. PBI Job File

   Next, we demonstrate an example of
   the way OpenFMO is run on a PBS queuing system.
   When performing the same GPU-accelerated FMO-RHF 
   calculations described above,
   you need to write a PBS job file.
   The minimum example "job.sh" is as follows:

     code-block:: sh

        #!/bin/sh 
        #PBS -j oe
        #PBS -N JobName
        OFMOPATH="/OpenFMO/executables/install/directory"
        SCRDIR="/Pass/To/OpenFMO/Scratch/Files/Directory"
        LIBRARY_PATH=${LD_LIBRARY_PATH}
        cd ${PBS_O_WORKDIR}
        OMP_NUM_THREADS=4
        opt=""
        opt+=" -np 32"
        opt+=" -bindir ${OFMOPATH}"
        opt+=" -B 0"
        opt+=" -ng 15"
        opt+=" -scrdir ${SCRDIR}"
        date
        set -x
        mpiexec.hydra -np 1 -print-rank-map ${OFMOPATH}/ofmo-master $opt Input_Fine_Name
        set +x
        date

   Then, submit the PBS job:

   `$ qsub job.sh`


Input File Format
-----------------

1. "Skeleton-RHF"

   Here, we show an example of 
   a "skeleton-RHF" input file used for
   the RHF energy calculation
   of one water molecule with the 6-31G(d) basis sets:

     code-block:: input

       6-31G*
       0 3
       O  10.438  -22.339  1.220
       H  11.395  -22.339  1.220
       H  10.198  -21.412  1.220

   The input file format is schematically written by

     code-block:: input

       Name_of_Basis_Sets
       Molecular_Charge Number_of_Atoms
       Atomic_Name X Y Z
       Atomic_Name X Y Z
       Atomic_Name X Y Z
       ...

2. OpenFMO

   OpenFMO adopts the same input-file format
   as the FMO calculations implemented in
   GAMESS ( http://www.msg.ameslab.gov/gamess/ )
   *ab initio* quantum chemistry package.
   Since some of the input groups used in GAMESS are
   directory used in the OpenFMO,
   the GAMESS documentations,
   Input Description ( http://www.msg.ameslab.gov/gamess/GAMESS_Manual/input.pdf ) and Further Information ( http://www.msg.ameslab.gov/gamess/GAMESS_Manual/refs.pdf ) , 
   are also useful for the users of OpenFMO.


   See the details of input format, input groups, and their options 
   used in OpenFMO at OpenFMO official site ( http://www.openfmo.org/ ).

     file-block:: sample input

        $gddi niogroup=1 nioprocs=1 $end
        $fmo nfrag=3 ICHARG(1)= 0,0,0
        INDAT(1)=0, 1, -6, 0, 7, -13, 0, 14, -17,0  $end
        $basis gbasis=sto ngauss=3 $end
        $fmoxyz
        N  7.0   3.5584    0.0170    0.1638
        H  1.0   3.6446   -0.8687   -0.3332
        H  1.0   3.4912   -0.2124    1.1546
        C  6.0   2.3540    0.7121   -0.2674
        H  1.0   2.2350    1.6486    0.2858
        H  1.0   2.4304    0.9444   -1.3339
        C  6.0   1.1558   -0.1725    0.0097
        O  8.0   1.1192   -0.9807    0.9350
        N  7.0   0.1194    0.0665   -0.8809
        H  1.0   0.2322    0.7805   -1.5946
        C  6.0  -1.1505   -0.6217   -0.8231
        H  1.0  -0.9953   -1.6290   -0.4254
        H  1.0  -1.5383   -0.6729   -1.8442
        C  6.0  -2.1620    0.0850    0.0422
        O  8.0  -3.2962   -0.3376    0.2304
        O  8.0  -1.6980    1.2320    0.5903
        H  1.0  -2.3743    1.6726    1.1478
        $end
        $FMOLMO
        STO-3G 5 5
        1 0  -0.117784    0.542251    0.000000    0.000000    0.850774
        0 1  -0.117787    0.542269    0.802107    0.000000   -0.283586
        0 1  -0.117787    0.542269   -0.401054   -0.694646   -0.283586
        0 1  -0.117787    0.542269   -0.401054    0.694646   -0.283586
        0 1   1.003621   -0.015003    0.000000    0.000000    0.000000
        $end
        $FMOBND
        -4    7 STO-3G
        -11   14 STO-3G
        $end


Samples
-------

You can download some input/output files for "skeleton-RHF" and OpenFMO
at OpenFMO official site ( http://www.openfmo.org/ ).




