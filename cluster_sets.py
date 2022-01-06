def get_template(cluster):
    if cluster == 'cesup_fermi' or cluster == 'cesup_fermi_gpu' or 'cesup_fermi2' :
        TEMPLATE= """#!/bin/sh
#PBS -S /bin/sh
#PBS -l select={nodes}:ncpus={ncpus}:mpiprocs={ncpus}:ompthreads=1:ngpus=2:mem=92GB
#PBS -N {job_name}
#PBS -V
#PBS -j oe
cd $PBS_O_WORKDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/acml/ifort64/lib
# export MPI_BUFS_PER_HOST=1024
# export MPI_BUFS_PER_PROC=128
# export MPI_GROUP_MAX=1024
export OMP_NUM_THREADS=1
"""
    if cluster == 'cesup_gauss':
        TEMPLATE="""
#!/bin/sh
#PBS -S /bin/sh
#PBS -l nodes={nodes}:ppn={ncpus}
#PBS -N {job_name}
#PBS -V
#PBS -j oe
INPUT=CsSbBrUdj
cd $PBS_O_WORKDIR
# export ESPRESSO_TMPDIR=$SCRATCH/$USER
# if [ ! -d $ESPRESSO_TMPDIR ]; then mkdir -p $ESPRESSO_TMPDIR; fi
# export EXEC=$EXEC
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/acml/ifort64/lib
# export MPI_BUFS_PER_HOST=1024
# export MPI_BUFS_PER_PROC=128
# export MPI_GROUP_MAX=1024
export OMP_NUM_THREADS=1
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------" 
echo "____ {job_name}______"
        """
    if cluster == 'sdumont':
        TEMPLATE="""
#!/bin/bash
#SBATCH --nodes={nodes}                      # here the number of nodes
#SBATCH --ntasks={ntasks}                    # here total number of mpi tasks
## #SBATCH --ntasks-per-node=1            # here ppn = number of process per nodes
#SBATCH -p {queue}                 # target partition
#SBATCH --exclusive                    # to have exclusvie use of your nodes
#SBATCH --job-name=$INPUT
#SBATCH --output=log.txt
## #SBATCH --time={settime}
echo $SLURM_JOB_NODELIST
cd $SLURM_SUBMIT_DIR
"""
    if cluster == 'furg' or cluster == 'aws':
        TEMPLATE="#!/bin/bash \n "
    if cluster == 'local':
        TEMPLATE="#!/bin/bash \n"
    return TEMPLATE

def getQEexec(calc,cluster,user):
    """ returns full path of QEexecutable for given inputfile"""
    QEexec_dict={'scf':'pw.x','nscf':'pw.x','relax':'pw.x','vc-relax':'pw.x','bands':'pw.x','dos':'dos.x',
            'chargexsf':'pp.x','chargecube':'pp.x','spin-polarization':'pp.x','V11':'pp.x','Vtotal':'pp.x',
            'spol':'pp.x',
            'homo0':'pp.x','homo1':'pp.x','homo2':'pp.x','homo3':'pp.x',
            'lumo0':'pp.x','lumo1':'pp.x','lumo2':'pp.x','lumo3':'pp.x',
            'ppbands':'bands.x','pdos':'projwfc.x','pdos_bands':'projwfc.x','avg':'average.x',
            'nscf_seq':'pw.x',  # exclusive of AFLOWpi function
            'simple':'simple.x',  'simple_ip':'simple_ip.x',
            }
    QEdir=""
    if cluster =='cesup_fermi':
        QEdir="/home/u/{user}/qe-6.4.1/bin/".format(user=user)
    elif cluster =='cesup_fermi2' or cluster =='cesup_fermiNC' or cluster=='cesup_fermiDJ':
        QEdir="/home/u/{user}/qe-6.4.1_2/bin/".format(user=user)
    elif cluster =='cesup_gaussNC':
        QEdir="/home/u/{user}/q-e-qe-6.6/bin/".format(user=user)
    elif cluster =='cesup_fermi_gpu':
        QEdir="/opt/qe/bin/"
    elif cluster =='cesup_gauss':
        QEdir="/home/u/{user}/qe-6.4.1/bin/".format(user=user)
    elif cluster =='sdumont':
        QEdir="/scratch/lamai/rogerio.gouvea/qe-6.4.1/bin/"
    elif cluster =='furg':
        QEdir=""
    elif cluster == 'local':
        QEdir=""
    else:
        raise("cluster specified not implemented, edit cluster_sets.py on MinFlow directory.")
    #### special cases
    if cluster == 'cesup_fermi_gpu' and (calc == 'simple' or calc == 'simple_ip'):
        QEdir="/home/u/{user}/qe-6.4.1/bin/".format(user=user) ## simple not compiled in gpu version
    return QEdir+QEexec_dict[calc]

def getPseudoDir(cluster="local",user="user"):
    if (cluster =='cesup_fermi' or cluster =='cesup_fermi2' or 
        cluster =='cesup_fermi_gpu' or cluster =='cesup_gauss' ) :
        pseudodir="/home/u/{user}/pseudo".format(user=user)
    elif cluster=='cesup_fermiNC' or cluster=='cesup_gaussNC':
#        pseudodir="/home/u/raglamai/pseudo/SG15/"
        pseudodir="/home/u/{user}/pseudo/AFLOWPI_PSEUDOS".format(user=user)
    elif cluster=='cesup_fermiDJ':
        pseudodir="/home/u/raglamai/pseudo/djsr"

    elif cluster =='furg' :
        pseudodir="/home/mestrado/rogeriog/QE-calculos/pseudos"
    elif cluster == 'aws':
        pseudodir="/home/ubuntu/pseudo"
    elif cluster =='sdumont':
        pseudodir="/scratch/lamai/{user}/pseudo".format(user=user)
    elif cluster == 'local':
        pseudodir="./"
    else:
        raise("cluster specified not implemented, edit cluster_sets.py on MinFlow directory.")
    return pseudodir
