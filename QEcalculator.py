# coding: utf-8
from ase.io import read
from ase.build import make_supercell
import numpy as np
from ase.io.espresso import write_espresso_in
from ase.io.espresso import read_espresso_out
### there was an issue with version of ase espresso.py, I had to 
### modify according to followin issue: https://gitlab.com/ase/ase/-/issues/984
from ase.dft.kpoints import bandpath
from copy import deepcopy
from io import StringIO
import re
import subprocess
import os, shlex, sys,shutil, stat, signal
import MinFlow.cluster_sets as cluster_sets
from MinFlow.special_functions import set_element_magnetization
from ase.dft.kpoints import BandPath

class QEcalc:
    """
    prefix: prefix of the calculation
    calculation: 'relax','vc-relax','bands','scf','nscf'
    ecutwfc: cutoff in Ry
    scale_cutrho: to produce ecutrho as scale_cutrho*ecutwfc
    pseudodict: contains dir entry with pseudo_dir and potentials entry with each pseudopotential
    kspacing: finds kspacing according to https://github.com/rosswhitfield/ase/blob/a652830b9821bfde233463616a4249758b7ffb00/ase/io/espresso.py#L1599
    ktuple: follows quantum espresso automatic like (6,6,6,0,0,0) for 6x6x6 mokhorstpack centered in gamma.
    """
    def __init__(self,prefix=None,calculation=None,ecutwfc=60,pseudodict=None,structure=None,scale_cutrho=4, kpath=None, kspacing=None,ktuple=(1,1,1,0,0,0),nbnd="",**kwargs):
        if ( prefix is None ) or (calculation is None) or ( pseudodict is None) or (structure is None):
            raise Exception("You must provide prefix, calculation, pseudodict and structure to declare the calculator!!")
        self.prefix=prefix
        if "." in prefix:
            raise Exception("prefix cannot contain \".\" !")
        self.calculation=calculation
        self.ecutwfc=ecutwfc
        self.scale_cutrho=scale_cutrho
        self.pseudodict=pseudodict
        self.structure=structure
        self.kpath=kpath
        self.kspacing=kspacing
        self.ktuple=ktuple
        self.input_data = {     ### this is the default parameters for input
        "CONTROL" : {
	  'etot_conv_thr' :   1.0E-04,
	  'forc_conv_thr' :   1.0E-03,
	  'outdir' : './OUTPUT',
	  'verbosity ': 'high'
	}, \
	'SYSTEM': {
	  'occupations' : 'smearing',
	  'smearing' : 'gauss',
	  'degauss' :   0.0002,
	}, \
	'ELECTRONS': {
	  'conv_thr' :   1.0E-06,
	  'electron_maxstep' : 300,
	  'mixing_beta' :   0.4
	}
	}
        ## if a input_data is passed to the function its taken as reference
        ## in substitution of the default
        if kwargs.get('CUSTOM',None) is not None:
            CUSTOMdict=kwargs.get('CUSTOM')
            for subdicts in CUSTOMdict.keys():
                for flag in CUSTOMdict[subdicts].keys():
                    self.input_data[subdicts][flag]=CUSTOMdict[subdicts][flag]
        if kwargs.get('input_data',None) is not None:
            self.input_data=kwargs.get('input_data')
        ## write basic input_data, overwrite input_data
        self.input_data['calculation']=calculation
        self.input_data['prefix']=prefix
        self.input_data['pseudo_dir']=pseudodict['dir']
        self.input_data['ecutwfc']=ecutwfc
        self.input_data['ecutrho']=ecutwfc*scale_cutrho

        if nbnd is None:
            self.nbnd=nbnd
        else:
            self.input_data['SYSTEM']['nbnd']=int(nbnd)
        
        if self.input_data['SYSTEM']['nspin'] == 2:
            element_1=self.structure.get_chemical_symbols()[0]
            self.magmoms=set_element_magnetization(self.structure,element_1)
            self.structure.set_initial_magnetic_moments(self.magmoms)
        else:
            self.magmoms=""


        ## these are filled when methods are applied
        self.input_file=""
        self.output_file=""
        self.relaxed_structure=""
        self.nelectron=""
        self.band_from_str=False
        self.symbols=""
        


    def write_input(self):
        """ Writes the input file and returns its name"""
        self.input_file='{}.{}.in'.format(self.prefix,self.calculation)
        self.output_file="out."+self.input_file.rsplit(".",1)[0]
        ## for the case you read last structure with get_relaxed
        ## updates structure to print input w last coordinates
        if self.relaxed_structure != "":
            self.structure=self.relaxed_structure
        with open(self.input_file,'w') as f:
            if (self.kspacing is not None) : ## case kspacing is specified
                write_espresso_in(f,self.structure,
                        input_data=self.input_data,
                        pseudopotentials=self.pseudodict['potentials'],
                        kspacing=self.kspacing,
                        crystal_coordinates=True)
            elif (self.kpath is not None):  ## case for band structure calculations
                if isinstance(self.kpath,BandPath):
                    write_espresso_in(f,self.structure,
                            input_data=self.input_data,
                            pseudopotentials=self.pseudodict['potentials'],
                            kpts=self.kpath,
                            crystal_coordinates=True)
                if isinstance(self.kpath,str):
                    write_espresso_in(f,self.structure,
                        input_data=self.input_data,
                        pseudopotentials=self.pseudodict['potentials'],
                        kpts=None,  ## gamma point algorithm
                        crystal_coordinates=True)
                    self.band_from_str=True
                    ## gamma KPOINTS will be replaced in next code

            elif (self.ktuple[:3] == (1,1,1)) and (self.ktuple[3:] == (0,0,0)):  ## case for gamma only
                write_espresso_in(f,self.structure,
                        input_data=self.input_data,
                        pseudopotentials=self.pseudodict['potentials'],
                        kpts=None,  ## gamma point algorithm
                        crystal_coordinates=True)
            else: ## general case of MP grid
                write_espresso_in(f,self.structure,
                        input_data=self.input_data,
                        pseudopotentials=self.pseudodict['potentials'],
                        kpts=self.ktuple[:3],
                        koffset=self.ktuple[3:],
                        crystal_coordinates=True)
        if self.band_from_str:   ## for the special case of banduppy Kpoints that are passed as string
            new_data=""
            with open(self.input_file,'r') as inputf:
                data=inputf.read()
                kpath_str="K_POINTS (crystal) \n"+self.kpath
                new_data=re.sub("K_POINTS.*\n",kpath_str,data)
            with open(self.input_file,'w') as inputf:
                inputf.write(new_data)
            ## modify previous input to include band path passed as a string.
    def updateCalc(self,outfile):
        with open(outfile) as f:
            pwo_lines=f.readlines()
            initialize_magn_valence=True
            readvalence=False
            readmagn=False
            got_ntyp=False
            self.symbols=[]
            magn=[]
            i_magn=0
            valence=[]
            i_val=0
            for line in pwo_lines:
                if "number of Kohn-Sham states" in line:
                    self.nbnd=int(line.split()[4])
                if "number of electrons" in line:
                    self.nelectron=int(float(line.split()[4]))
                if "number of atomic types" in line and initialize_magn_valence:
                    self.ntyp=int(float(line.split()[5]))
                    got_ntyp=True
                    magn=np.zeros(self.ntyp)
                    valence=np.zeros(self.ntyp)
                    initialize_magn_valence=False ## otherwise may initialize again
                if "atomic species   valence" in line:
                    readvalence=True
                    continue ## go to next line
                if got_ntyp and i_val == self.ntyp:
                    readvalence=False ## stop reading new types
                if readvalence:
                    valence[i_val]=float(line.split()[1])
                    self.symbols.append(line.split()[0])
                    i_val+=1
                if "atomic species   magnetization" in line:
                    readmagn=True
                    continue ## go to next line
                if got_ntyp and i_magn == self.ntyp:
                    readmagn=False ## stop reading
                if readmagn:
                    magn[i_magn]=float(line.split()[1])
                    i_magn+=1
            if not np.all(magn) : ## if all magn is 0 no magmom is updated. 
                self.magmoms=np.round(valence*magn,2)


    ##############################################
    def get_relaxed_structure(self,outfile,index=-1,vcrelax=True):  ## by default it takes the last optimized
        """ get ase.Atoms corresponding to relaxed structure from QE output """
    ### modify ase espresso.py to include nbnd later ###
        self.updateCalc(outfile)
#        with open(outfile) as f:
#            pwo_lines=f.readlines()
#            for line in pwo_lines:
#                if "number of Kohn-Sham states" in line:
#                    self.nbnd=int(line.split()[4])
#                if "number of electrons" in line:
#                    self.nelectron=int(float(line.split()[4]))
    ##############################################
        if vcrelax: ## otherwise doesnt take last coordinate in vcrelax
            try:
                with open(outfile) as f:
                    r=read_espresso_out(f,index=index,results_required=False)
                    for relaxed in r:
                        ### codigo para ler magmons e colocar na structure
                        if isinstance(self.magmoms,np.ndarray):
                            for idx, symbol in enumerate(self.symbols):
                                magmoms=set_element_magnetization(relaxed, symbol,
                                                          starting_mag=self.magmoms[idx])
                                relaxed.set_initial_magnetic_moments(magmoms)
                        self.relaxed_structure=relaxed

            except FileNotFoundError as error:
                print("Could not load outfile.")
        else:
            try:
                with open(outfile) as f:
                    r=read_espresso_out(f,index=index,
                                        results_required=True)
                    for relaxed in r:
                        ### codigo para ler magmons e colocar na structure
                        if isinstance(self.magmoms,np.ndarray):
                            from MinFlow.special_functions import set_element_magnetization
                            for idx, symbol in enumerate(self.symbols):
                                magmoms=set_element_magnetization(relaxed, symbol,
                                                          starting_mag=self.magmoms[idx])
                                relaxed.set_initial_magnetic_moments(magmoms)
                        self.relaxed_structure=relaxed
            except FileNotFoundError as error:
                print("Could not load outfile.")

    def seekpath(self):
        """function takes a calc and gives seekpath path, it has to read from espresso input, ase atoms doesnt
        produce correct data for seekpath."""
        inputfile=self.input_file
        ## reads cell parameters
        ## read atomic positions, converts atomic symbol to number, and produces
        ## atomic positions in crystal
        ## atomic numbers in same order of atomic positions

def write_sh(calcs=[],user="raglamai",job_name="job",cluster="local", nodes="", ncpus="", 
         queue="",settime="",QE_parflags=""):
            inputfiles=[]
            if len(calcs) == 0:
                raise Exception("You must specify some calcs to submit.")
            template=cluster_sets.get_template(cluster)
            ## writes the executables to be run 
            with open(user+"."+job_name+".sh", 'w') as f:
                ## writes initial part of sh file with cluster settings
                if cluster == "sdumont": 
                    ntasks=ncpus*nodes
                    f.write(template.format(job_name=job_name, nodes=nodes, ntasks=ntasks, queue=queue,settime=settime))
                else:
                    f.write(template.format(job_name=job_name, nodes=nodes, ncpus=ncpus))
                for i in range(len(calcs)):
                    QE_exec=cluster_sets.getQEexec(calcs[i].calculation,cluster,user)
                    output_file=calcs[i].output_file  ## output file determined on write_input
                    input_file=calcs[i].input_file   ## input file determined on write_input
                    inputfiles.append(input_file)
                    if QE_exec=="average.x":  ## average.x demands run on single processsor
                        f.write("{QE_exec} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file))
                    elif cluster=="local" or cluster=="aws" or cluster=="furg": ## those which execute locally
                        f.write("mpirun --oversubscribe -np {ncpus} {QE_exec} {QE_parflags} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file,QE_parflags=QE_parflags,ncpus=ncpus))
                    elif cluster=="sdumont": ## clusters that use srun
                        f.write("srun {QE_exec} {QE_parflags} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file,QE_parflags=QE_parflags))
                    else: ## clusters that use mpiexec command
                        f.write("mpiexec {QE_exec} {QE_parflags} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file,QE_parflags=QE_parflags))
            print("Input files written:")
            [print(x) for x in inputfiles]
            ## for cesup due to sequential submission this modification is necessary 
            if 'cesup' in cluster:
                edited_sh=""
                with open(user+"."+job_name+".sh", 'r') as f:
                    data=f.read()
                    data=data.split('<',1)
                    datainp=re.sub('<','-inp',data[1])
                    data=data[0]+'<'+datainp
                    edited_sh=data
                with open(user+"."+job_name+".sh", 'w') as f:
                    f.write(edited_sh)
            return f.name


def wait_job(job, waittime=5): ## 5s for every check
    while len(job) != 0:
        for jobid in job:
            x = subprocess.Popen('qstat -f '+jobid, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std_out, std_err = x.communicate()
            if std_err :
                job.remove(jobid)
                break
        os.system("sleep " + str(waittime))

def execute(cmd):   ## execute 
    popen = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, universal_newlines=True)
    for jobid in iter(popen.stdout.readline, ""):
        return [jobid] 

def relax_with_gpu(sh_file, **kwargs):
    import os, shlex, sys,shutil
    prefix=sh_file.split(".")[1]
    OUTPUT="out."+prefix+".vc-relax"
    INPUT=prefix+".vc-relax.in"
    ibrav0=True ### looks for CELL_PARAMETERS
    try:
        with open(OUTPUT) as f:
            pass
    except:
        OUTPUT="out."+prefix+".relax"
        INPUT=prefix+".relax.in"
        ibrav0=False
    end_relaxation=False
    iteration=0
    while not end_relaxation and iteration < kwargs.get('iterations',50):
        iteration += 1
        nat=0
        data="" ##holds structure information
        read=False
        with open(OUTPUT) as f:
           lines=f.readlines()
           for line in lines:
               m=re.search("number of atoms/cell",line)
               if m != None:
                   nat=int(line.split()[4])
       #            print(nat)
               if ibrav0 :
                   m=re.search("CELL_PARAMETERS",line)
               else :
                   m=re.search("ATOMIC_POSITIONS",line)
               if m != None:
                   data=""  ## only last occurance will be taken
                   read=True
                   n=1
               if read : 
                   data=data+line
                   n=n+1
                   if ibrav0:
                       if n == nat+7:
                           read=False
                   else :
                       if n == nat+2:
                           read=False
        if data == "":
            end_relaxation=True

        with open(INPUT) as f:
          text=f.read()
          if ibrav0:
              pattern = re.compile('CELL_PARAMETERS.*', re.DOTALL|re.I)
          else:
              pattern = re.compile('ATOMIC_POSITIONS.*', re.DOTALL|re.I)
          newtext=re.sub(pattern, data, text)

        with open(INPUT,"w") as f:
           f.write(newtext)

       # copy last out
        shutil.copyfile(OUTPUT, OUTPUT+str(iteration))
        # Execute and print out, get jobid but wont use it
        jobid=execute("qsub {sh_file}".format(sh_file=sh_file))
        wait_job(jobid)

       # find jobdone
        with open(OUTPUT) as f:
           lines=f.readlines()
           for line in lines:
               done=re.search("JOB DONE.",line)
               if done != None :
                   end_relaxation=True





def run(sh_file,local=False,initialize_only=False,fermi_gpu=False,**kwargs):
        ## make sh file executable
        os.chmod(sh_file, os.stat(sh_file).st_mode | stat.S_IEXEC)
        ## qsub or execute file
        jobid=0
        print("running... {sh_file}".format(sh_file=sh_file))
        if not local:
            if initialize_only:
                jobid=execute("qsub {sh_file}".format(sh_file=sh_file))
            elif fermi_gpu:
                jobid=execute("qsub {sh_file}".format(sh_file=sh_file))
                wait_job(jobid)
                relax_with_gpu(sh_file)
            else:
                jobid=execute("qsub {sh_file}".format(sh_file=sh_file))
                wait_job(jobid)
        else:
            if initialize_only:
                P = subprocess.Popen("./{}".format(sh_file), 
                stdout=subprocess.PIPE, universal_newlines=True)
                os.system("sleep 5")
                P.kill()
                os.system("sleep 5")  
                ## have to be insistent otherwise file remains
                ## open when try to read it
                os.kill(int(P.pid), signal.SIGKILL)
                os.system("sleep 5")
                if kwargs["prefix"]:
                    execute("touch {}.EXIT".format(kwargs['prefix']))
            else:
                jobid=execute("./{sh_file}".format(sh_file=sh_file))
        print("done.")


class QEcalc2():
    def __init__(self,prefix, calculation):
        self.prefix=prefix
        self.calculation=calculation
        self.input_file="{}.{}.in".format(self.prefix,self.calculation)
        self.output_file="out.{}.{}".format(self.prefix,self.calculation)

def gencalc_absorption(calcdict,kgrid_for_absorption,simple_dict=None,
            simpleip_dict=None, from_previous_scf=False, from_previous_nscf=False,
            cluster="local"):
    ##calcdict may be the calcdict currently used to setup a calculation,
    ## kgrid_for_absorption, tuple like (k1,k2,k3) with kgrid dimensions for
    ## absorption calculations.
    ## simpledict and simpleip_dict will have non-default values 
    ## to substitute
    ##to generate required simple inputs.
    ## if from previous scf skip scf calculation, be sure they use same prefix
    ## given in calcdict otherwise nscf wont run.
    ## cluster, have to give cluster parameters to do tmp calculation that gets
    ## complete nscf grid.
    ##--set simple scf calculation to be run later--
    calcdict['calculation']='scf' ### have to be scf to do the test.
    calc=QEcalc(**calcdict)
    if not from_previous_scf or not from_previous_nscf:
        calc.write_input()
    ##read scf dict write scf input for given ktuple and execute a local pw.x
    ##just to get kpoints.--
    calcdict['prefix']=calcdict['prefix']+"_tmp"
    calcdict['CUSTOM']['SYSTEM']['noinv']=True
    calcdict['CUSTOM']['SYSTEM']['nosym']=True
    calcdict['ktuple']=kgrid_for_absorption+(0,0,0)
    calctmp=QEcalc(**calcdict)
    calctmp.write_input()
    shdict={'user':"raglamai",'job_name':calc.prefix,
    'calcs':[calctmp], 'cluster':"local", 'nodes':1,'ncpus':1}
    sh_file=write_sh(**shdict)
    run(sh_file, local=True,initialize_only=True,prefix=calcdict['prefix'])
    ##-- then, read out.scf to get coordinates of kpoint--
    outfile_tmp="out."+calctmp.prefix+".scf"
    kpts=""
    nks=0
    ## we will need nbnd and nelectron in simple
    nbnd=0
    nelectron=0
    with open(outfile_tmp,"r") as outtmp:
        lines=outtmp.readlines()
        read_kcryst_coords=False
        for line in lines:
            if "number of Kohn-Sham states" in line:
                nbnd=int(line.split()[4])
            if "number of electrons" in line:
                nelectron=int(float(line.split()[4]))
            if re.search("\s*\scryst. coord.\s\s*",line):
                read_kcryst_coords=True
                continue
            if re.search("^\s*\n",line):
                read_kcryst_coords=False
            if read_kcryst_coords:
                kpt=re.findall(".*\((.*)\).*",line)
                kpt=kpt[0].lstrip()
                kpts=kpts+kpt+' 1.0 \n'
                nks+=1 
    ##--write nscf input with noinv nosym and full kpoint grid.
    calcdict['prefix']=calc.prefix
    calcdict['calculation']='nscf'
    calcdict['kpath']=str(nks)+'\n'+kpts
    ### bands have to be substantially increased to get higher energies
    calcdict['CUSTOM']['SYSTEM']['nbnd']=nbnd+100
    calc_nscf=QEcalc(**calcdict)

    if not from_previous_nscf:
        calc_nscf.write_input()

    def_simple_dict={ 'prefix':calc.prefix,
    'calc_mode':1, 'num_nbndv': nbnd,
    'num_val':nbnd, 'num_cond':100,
    's_bands':0.1, 's_product':1.0,
    'nkpoints(1)':calctmp.ktuple[0],
    'nkpoints(2)':calctmp.ktuple[1],
    'nkpoints(3)':calctmp.ktuple[2],
    'nonlocal_commutator':'.true.',
    }
    ##--overwrite def simple_dict with the one passed
    if simple_dict is not None:
        for key,value in simple_dict.items():
            def_simple_dict[key]=value
    ## more keys may be necessary in the template
    inputsimple_str="""&inputsimple
prefix='{prefix}'
outdir='./OUTPUT'
calc_mode={calc_mode}
num_nbndv={num_nbndv}
num_val={num_val}
num_cond={num_cond}
s_bands={s_bands}
s_product={s_product}
nkpoints(1)={nkpoints(1)}
nkpoints(2)={nkpoints(2)}
nkpoints(3)={nkpoints(3)}
nonlocal_commutator={nonlocal_commutator}
/
""".format(**def_simple_dict)
    with open("{prefix}.simple.in".format(prefix=calc.prefix),'w') as f:
        f.write(inputsimple_str)
    calc_simple=QEcalc2(calc.prefix,'simple')
    
    ### now for simpleip input
    #write simpleip input from dict
    def_simpleip_dict={ 'prefix':calc.prefix,
    'interp_grid(1)':calctmp.ktuple[0],
    'interp_grid(2)':calctmp.ktuple[1],
    'interp_grid(3)':calctmp.ktuple[2],
    'fermi_degauss' :  0.002,
    'fermi_ngauss' :  -99,
    'drude_degauss':  0.002,
    'wmin': 0.0, 'wmax': 0.735, 'nw': 1000,
    'inter_broadening': 0.002,
    'intra_broadening': 0.002,
    'nonlocal_commutator':'.true.',
    'nonlocal_interpolation':'.false.',
    }
    ##--overwrite def simple_dict with the one passed
    if simpleip_dict is not None:
        for key,value in simpleip_dict.items():
            def_simpleip_dict[key]=value
    ## more keys may be necessary in the template
    inputsimpleip_str="""&inputsimpleip
simpleip_in%prefix = '{prefix}'
simpleip_in%outdir = './OUTPUT'
simpleip_in%interp_grid(1) = {interp_grid(1)}
simpleip_in%interp_grid(2) = {interp_grid(2)}
simpleip_in%interp_grid(3) = {interp_grid(3)}
simpleip_in%fermi_degauss = {fermi_degauss}
simpleip_in%fermi_ngauss = {fermi_ngauss}
simpleip_in%drude_degauss = {drude_degauss}
simpleip_in%wmin = {wmin}
simpleip_in%wmax = {wmax}
simpleip_in%nw = {nw}
simpleip_in%inter_broadening = {inter_broadening}
simpleip_in%intra_broadening = {intra_broadening}
simpleip_in%nonlocal_commutator= {nonlocal_commutator}
simpleip_in%nonlocal_interpolation= {nonlocal_interpolation}
/
""".format(**def_simpleip_dict)
    
    with open("{prefix}.simple_ip.in".format(prefix=calc.prefix),'w') as f:
        f.write(inputsimpleip_str)
    calc_simpleip=QEcalc2(calc.prefix,'simple_ip')
    
    if from_previous_scf:
        return [calc_nscf,calc_simple,calc_simpleip]
    elif from_previous_nscf:
        return [calc_simple,calc_simpleip]
    else:
        return [calc,calc_nscf,calc_simple,calc_simpleip]


