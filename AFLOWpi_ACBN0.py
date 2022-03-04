import AFLOWpi
import subprocess
import MinFlow.cluster_sets as cluster_sets
from MinFlow.QEcalculator import run
import re
import glob
import fileinput
def modify_AFLOWpi_pyfile(pyfile,cluster,local=False):
    # opening the file in read mode
    replacement = ""
    replacement += "import MinFlow \n"
    replacement += "from MinFlow.QEcalculator import run \n"
    replacement += "from MinFlow.AFLOWpi_ACBN0 import _run \n"
    special_scf_pdos=False
    with open(pyfile, "r") as pyf:
        for line in pyf:
            if re.search("\s*oneCalc\['_AFLOWPI_WORKFLOW_'\]=\['scf', 'dos'\]\s*",line):
                special_scf_pdos=True
        # using the for loop
    if not special_scf_pdos:
        with open(pyfile, "r") as pyf:
            for line in pyf:
                if re.search("oneCalc['_AFLOWPI_WORKFLOW_']=['scf', 'dos']",line):
                    special_scf_pdos=True
                line = line.rstrip()
                if re.search("AFLOWpi.run._oneRun",line):
                    if local:
                        line="        run('_'+ID+'.sh',local=True)"
                    else:
                        line="        run('_'+ID+'.sh')"
                if re.search("AFLOWpi.run._submitJob",line):
                        line=""
                if re.search("AFLOWpi.scfuj._run",line):
                            if local:
                                line="oneCalc, newUvals = _run(__submitNodeName__,oneCalc,ID,cluster='{cluster}',local=True,config=configFile,mixing=0.0,kp_mult=2.0,U_eff=True)".format(cluster=cluster)
                            else:
                                line="oneCalc, newUvals = _run(__submitNodeName__,oneCalc,ID,cluster='{cluster}',local=False,config=configFile,mixing=0.0,kp_mult=2.0,U_eff=True)".format(cluster=cluster)
                replacement = replacement + line + "\n"
    else:  ## special case
        with open(pyfile, "r") as pyf:
            data=pyf.read()
            data=re.sub("AFLOWpi.run._oneRun.*\n", "\n", data, 1)
            data=re.sub("AFLOWpi.run._oneRun.*calcType='dos'.*\n","AFLOWpi.run._qe__pre_run(oneCalc,ID,'dos',__submitNodeName__,'espresso')\n",data,1)
            data=re.sub("AFLOWpi.run._oneRun.*calcType='pdos'.*\n","AFLOWpi.run._qe__pre_run(oneCalc,ID,'pdos',__submitNodeName__,'espresso')\n",data,1)
            data=re.sub("#END_RUN_BLOCK.*\n","run('_'+ID+'.sh',local=True) \n#END_RUN_BLOCK \n",data,1)
            data=re.sub("AFLOWpi.run._submitJob.*\n","\n",data,1)
            replacement+=data
    with open(pyfile,"w") as pyf:
        pyf.write(replacement)


def writeAFLOWpi_sh(calcs=[],user="raglamai",job_name="job",cluster="local", nodes="", ncpus="", 
         queue="",settime="",QE_parflags="",AFLOWpiID="",AFLOWpiPATH=""):
            inputfiles=[]
            if len(calcs) == 0:
                raise Exception("You must specify some calcs to submit.")
            template=cluster_sets.get_template(cluster)
            ## writes the executables to be run 
            filename=AFLOWpiPATH+"_"+AFLOWpiID+".sh"
            with open(filename, 'w') as f:
                if cluster == "sdumont":
                    ntasks=ncpus*nodes
                    f.write(template.format(job_name=job_name, nodes=nodes, ntasks=ntasks, queue=queue,settime=settime))
                else:
                    f.write(template.format(job_name=job_name, nodes=nodes, ncpus=ncpus))
                for i in range(len(calcs)):
                    QE_exec=cluster_sets.getQEexec(calcs[i],cluster,user)
			## specific calculations have an additional sufix
                    sufix=""
                    if calcs[i]=="nscf_seq":
                        sufix="_nscf"
                    if calcs[i]=="dos":
                        sufix="_dos"
                    if calcs[i]=="pdos":
                        sufix="_pdos"
                ####
                    input_file=AFLOWpiID+sufix+".in"
                    output_file=AFLOWpiID+sufix+".out"  ## output file determined on write_input
                    if cluster=="local" or cluster=="aws" or cluster=="furg": ## those which execute locally
                            f.write("mpirun --oversubscribe -np {ncpus} {QE_exec} {QE_parflags} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file,QE_parflags=QE_parflags,ncpus=ncpus))
                    elif cluster=="sdumont": ## clusters that use srun
                            f.write("srun {QE_exec} {QE_parflags} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file,QE_parflags=QE_parflags))
                    else: ## clusters that use mpiexec command
                            f.write("mpiexec {QE_exec} {QE_parflags} < {input_file} > {output_file} \n".format(QE_exec=QE_exec,input_file=input_file,output_file=output_file,QE_parflags=QE_parflags))
            if 'cesup' in cluster:
                edited_sh=""
                with open(filename, 'r') as f:
                    data=f.read()
                    data=data.split('<',1)
                    datainp=re.sub('<','-inp',data[1])
                    data=data[0]+'<'+datainp
                    edited_sh=data
                with open(filename, 'w') as f:
                    f.write(edited_sh)


############# MODIFIED _run for ACBN0 ####################3

import os
import re
import logging
import copy 
import numpy as np
import atexit
import math 
import __main__
import itertools as it
import shutil 
import string
import itertools
import collections
import csv
import subprocess,stat
from .QEcalculator import run

def _run(__submitNodeName__,oneCalc,ID,cluster="cesup_fermiNC",local=False,config=None,mixing=0.10,kp_mult=1.6,U_eff=True):
        """there is additional variable cluster"""
        execPrefix = ''
        execPostfix = ''
        oneCalcID = ID


        def abortIFRuntimeError(subdir, ID):
            outfile = file(os.path.join(subdir, "%s.out"%ID)).read()
            errorList = re.findall(r'from (.*) : error #.*\n',outfile)
            if len(errorList) > 0:        
                logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
                print(("Error in %s.out -- ABORTING ACBN0 LOOP"%ID))                    
                raise SystemExit



        if '__runList__' not in list(oneCalc.keys()):
            oneCalc['__runList__']=[]

            
        if config is not None:
                AFLOWpi.prep._forceGlobalConfigFile(config)
                logging.debug('forced config %s' % config)
        else:
                try:
                        config = AFLOWpi.prep._getConfigFile()
                        AFLOWpi.prep._forceGlobalConfigFile(config)
                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)


        if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
            execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

        else:
            execPrefix=''


        if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
                execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
        else:
                execPostfix=''


        if AFLOWpi.prep._ConfigSectionMap('run','engine') == '':
                engine = AFLOWpi.prep._ConfigSectionMap('run','engine')
        else:
                engine = 'espresso'


        subdir = oneCalc['_AFLOWPI_FOLDER_']
        oneCalc['_AFLOWPI_CONFIG_']=config

        nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_mult,band_factor=1.0,wsyminv=False)    
##################################################################################################################
        
##################################################################################################################
        pdos_calc,pdos_ID = AFLOWpi.scfuj.projwfc(oneCalc,ID,ovp=True)

        if not re.match('northo',execPostfix) or not re.match('no',execPostfix):
            execPostfix+=' -northo 1'
        if local:
            run("_"+ID+".sh",local=True)
        else:
            run("_"+ID+".sh")
        splitInput = AFLOWpi.retr._splitInput(nscf_calc['_AFLOWPI_INPUT_'])
        modified_run_tb_ham_prep(__submitNodeName__,oneCalc,ID,kp_factor=kp_mult,cond=0,ovp=True,band_factor=1.0,wsyminv=True)

        AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.save'])
        AFLOWpi.scfuj._add_paopy_header(oneCalc,ID,shift_type=1,shift=1.0,thresh=0.90,tb_kp_mult=1.0,acbn0=True,ovp=True,smearing='gauss')
        AFLOWpi.scfuj._run_paopy(oneCalc,ID,acbn0=True)

        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        '''
            Will need to be filled in for the executable name with whatever wanT executable is called 
            and that executable needs to be moved to the calculation directory tree before this is called
        '''
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham_up.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham_down.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp_up.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp_down.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp.txt")

############
##################################################################################################################

        #Get new U values from acbn0.py 
        try:
            AFLOWpi.scfuj.acbn0(oneCalc, pdos_ID,execPrefix)
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
            raise SystemExit



        newUvals,newJvals = AFLOWpi.scfuj.getU_frmACBN0out(oneCalc,ID,U_eff=U_eff)
        Uvals=newUvals
        Jvals=newJvals

        #slight to help with oscillation
        try:
            old_U = oneCalc['_AFLOWPI_UVALS_']
        except:
            old_U = newUvals


        for spec,val in list(old_U.items()):
            newUvals[spec]=mixing*old_U[spec]+(1.0-mixing)*newUvals[spec]


        oneCalc['_AFLOWPI_UVALS_']=newUvals
        oneCalc['_AFLOWPI_JVALS_']=newJvals

        #Update Uvals
        oneCalc = AFLOWpi.scfuj.updateUvals(oneCalc, newUvals,newJvals,ID=ID,U_eff=U_eff)
        #update Uvals in _<ID>.py
        AFLOWpi.prep._modifyVarVal(oneCalc,ID,varName='uValue',value=newUvals)

        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        a = ID + '.in'
        
        inputfile = oneCalc['_AFLOWPI_INPUT_']
        with  open(os.path.join(subdir,a),'w') as new_inputfile:
                new_inputfile.write(inputfile)

        '''make sure the input is set to 'from_scratch'''
        AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','restart_mode',"'from_scratch'")

        '''
        return the new uVals to the main script where the 
        conditional statement will decide to rerun with 
        the new U vals as input or keep going on with the 
        rest of the script.
        '''

        '''save one calc with the full list so next iteration it runs everything fresh'''
        oneCalc['__runList__']=[]


        try:
            oneCalc['_SCFUJ_LoopCount_'] += 1
        except:
            oneCalc['_SCFUJ_LoopCount_']  = 1

        oneCalc['__status__']['SCFUJ Iteration']=oneCalc['_SCFUJ_LoopCount_']

#        try:
#            '''check for uval oscillation and if it's happening end the loop'''
#            if checkOscillation(ID,new_oneCalc)==True:
#                new_oneCalc['__status__']['Error']='U Value Oscillation'
#                AFLOWpi.prep._saveOneCalc(new_OneCalc,ID)
#                raise SystemExit
#        except Exception as e:
#            AFLOWpi.run._fancy_error_log(e)
            


        '''remove all the extra .oneCalc and input files files generated to start fresh'''
        for stage in ['nscf','pdos','WanT_bands']:
            try:
                os.remove('./_%s_%s.oneCalc'%(ID,stage))
#                os.remove('./%s_%s.in'%(ID,stage))
            except:
                pass
            

        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        try:
            pass
#            oneCalc,ID = AFLOWpi.prep._oneUpdateStructs(oneCalc,ID,override_lock=True)

        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


        try:
#            AFLOWpi.prep._clean_want_bands(oneCalc,ID)
            AFLOWpi.plot._bands(oneCalc,ID,yLim=[-15,5],postfix='acbn0_TB',tight_banding=True)
        except:
            pass


        return oneCalc, Uvals


def modified_run_tb_ham_prep(__submitNodeName__,oneCalc,ID,config=None,kp_factor=2.0,cond=1,ovp=False,band_factor=1.25,tetra_nscf=False,wsyminv=False):
        execPrefix = ''
        execPostfix = ''
        oneCalcID = ID

        try:
            pw_path = os.path.join(AFLOWpi.prep._ConfigSectionMap('prep','engine_dir'),'pw.x')
            shutil.copy(pw_path,oneCalc['_AFLOWPI_FOLDER_'])
        except: pass

        def abortIFRuntimeError(subdir, ID):
            with open(os.path.join(subdir, "%s.out"%ID)) as ifo:
                outfile =  ifo.read()
 
            errorList = re.findall(r'from (.*) : error #.*\n',outfile)
            if len(errorList) > 0:        
                logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
                print(("Error in %s.out -- ABORTING ACBN0 LOOP"%ID))                    
                raise SystemExit



        if '__runList__' not in list(oneCalc.keys()):
            oneCalc['__runList__']=[]

            
        if config is not None:
                AFLOWpi.prep._forceGlobalConfigFile(config)
                logging.debug('forced config %s' % config)
        else:
                try:
                        config = AFLOWpi.prep._getConfigFile()
                        AFLOWpi.prep._forceGlobalConfigFile(config)
                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)


        if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
            execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

        else:
            execPrefix=''


        if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
                execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
        else:
                execPostfix=''


        if AFLOWpi.prep._ConfigSectionMap('run','engine') == '':
                engine = AFLOWpi.prep._ConfigSectionMap('run','engine')
        else:
                engine = 'espresso'


        subdir = oneCalc['_AFLOWPI_FOLDER_']
        oneCalc['_AFLOWPI_CONFIG_']=config

        if 'scf' not in oneCalc['__runList__']:

            try:
                npool=AFLOWpi.retr._get_pool_num(oneCalc,ID)        

                if npool!=1:
                    if len(re.findall(r'npool[s]*\s*(?:\d*)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool[s]*\s*(?:\d*)','npool %s'%npool,execPostfix)
                        logging.debug(execPostfixPrime)

            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)


##################################################################################################################
#######            AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)


            oneCalc['__runList__'].append('scf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            



            nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_factor,unoccupied_states=cond,band_factor=band_factor,tetra_nscf=tetra_nscf,wsyminv=wsyminv)  



        else:
            '''if we are restarting from a job killed from going walltime 
            try to load ID_nscf and if we can't then just make a new one'''
            try:
                nscf_ID='%s_nscf' % ID
                nscf_calc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],nscf_ID)                
                '''we have to make sure nscf step has the correct walltime and start time if it's a restart'''
                nscf_calc['__walltime_dict__']=oneCalc['__walltime_dict__']
            except Exception as e:
                try:
                    nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_factor,
                                                                      band_factor=band_factor)  

                except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)

        nscf_exec_postfix = execPostfix_LOCAL = AFLOWpi.prep._ConfigSectionMap('TB','exec_postfix_nscf')            
        if nscf_exec_postfix != "":
            execPostfix = nscf_exec_postfix 
        
##################################################################################################################
        if 'nscf' not in oneCalc['__runList__']:


            try:
                npool=AFLOWpi.retr._get_pool_num(nscf_calc,nscf_ID)        

                if npool!=1:
                    if len(re.findall(r'npool\s*(?:\d+)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool\s*(?:\d+)','npool %s'%npool,execPostfix)
                        logging.debug(execPostfixPrime)

            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

#####            AFLOWpi.run._oneRun(__submitNodeName__,nscf_calc,nscf_ID,execPrefix=execPrefix,
    ###                            execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)
            AFLOWpi.retr._writeEfermi(nscf_calc,nscf_ID)

            abortIFRuntimeError(subdir, nscf_ID)

            oneCalc['__runList__'].append('nscf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)       
##################################################################################################################
        pdos_calc,pdos_ID = AFLOWpi.scfuj.projwfc(oneCalc,ID,paw=False,ovp=ovp)


        if not re.match('northo',execPostfix) or not re.match('no',execPostfix):
            execPostfix+=' -northo 1'

        execPrefix_LOCAL = AFLOWpi.prep._ConfigSectionMap('run','exec_prefix')
        execPostfix_LOCAL = AFLOWpi.prep._ConfigSectionMap('run','exec_postfix')            
        
        if 'pdos' not in oneCalc['__runList__']:
            pdosPath = os.path.join(AFLOWpi.prep._ConfigSectionMap('prep','engine_dir'),'projwfc.x')

            shutil.copy(pdosPath,oneCalc['_AFLOWPI_FOLDER_'])

            pdosPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'projwfc.x')

###            AFLOWpi.run._oneRun(__submitNodeName__,pdos_calc,pdos_ID,execPrefix=execPrefix,
  ###                              execPostfix=execPostfix,engine='espresso',calcType='custom',
     ###                           executable='projwfc.x',execPath=pdosPath)
#############
            oneCalc['__runList__'].append('pdos')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            abortIFRuntimeError(subdir, pdos_ID)



        eFermi=0.0


        AFLOWpi.prep._form_TB_dir(oneCalc,ID)
        eFermi=10.0

        splitInput = AFLOWpi.retr._splitInput(nscf_calc['_AFLOWPI_INPUT_'])
        del oneCalc['__runList__']

        dos_fermi = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_PAOFLOW_dos.efermi'%ID)

        with open(dos_fermi,'w') as ifo:
                ifo.write(str(0.0))

        return oneCalc,ID


def initAFLOW_ACBN0(prefix="",user="",cluster="cesup_fermiNC",local=False,nodes=1,ncpus=1):
    """
this will initialize an aflowpi calculation to determine U values
from ACBN0 implementation. It is going to look for a prefix.scf.in
file to generate inputs. Additionally a pseudopotential dir containing
NC pseudopotentials is required. Through the cluster command it will
find the appropriate QE executables.
Produces a folder ACBN0/prefix/ACBN0_prefix_0001 containing .py scripts
    """
	### most things are commented, maybe include further vc relax and bands, this is minimal to get U vals.
    ACBN0py="""
import AFLOWpi
# start the AFLOWpirame session
session = AFLOWpi.prep.init('ACBN0', '{prefix}',
                            config='./ACBN0.config')


# form the calculation set from ref input and allvars dict
calcs = session.from_file('{prefix}.scf.in')
# relax the structure
calcs.scf()
# calculate the the DOS and PDOS for Si2
calcs.dos()
calcs.plot.opdos(en_range=[-10,10],postfix='without_acbn0')
# calculate the bands for Si2
####calcs.bands(nk=200)
# do the plot the Electronic Band Structure
# and atom projected DOS for Si2
####calcs.plot.bands(en_range=[-10,10],DOSPlot='APDOS',
####                 postfix='without_acbn0')
# run the ACBN0 pseudo-hybrid functional to
# self-consistently get Hubbard U
calcs.acbn0(thresh=0.1,relax='scf',kp_factor=2.0)            
# calculate the the DOS and PDOS for PBE+U Si2
##calcs.vcrelax()
##calcs.dos()
# do the plot the Oribital Proj. DOS for PBE+U Si2
##calcs.plot.opdos(en_range=[-10,10],postfix='with_acbn0')
# calculate the bands for PBE+U Si2
##calcs.bands(nk=200)
# do the plot the Electronic Band Structure
# and atom projected DOS for PBE+U Si2
##calcs.plot.bands(en_range=[-10,10],DOSPlot='APDOS',
##                postfix='with_acbn0')
# run the calculation workflow
#calcs.submit()
with open("ID","w") as f:
    f.write(calcs.items()[0][0])
"""
    ACBN0config="""
[provenance]
title=acbn
author={user}
affiliation=ufrgs	
[prep]
pseudo_dir = {pseudo_dir}
engine_dir = {exec_dir}
work_dir = ./

[run]
python_command=python3
engine = espresso
	"""
    with open("ACBN0.py", "w") as f:
        f.write(ACBN0py.format(prefix=prefix))
    with open("ACBN0.config", "w") as f:
        pseudo_dir=cluster_sets.getPseudoDir(cluster=cluster,user=user)
        exec_dir=cluster_sets.getQEexec('scf',cluster,user).rsplit("/",1)[0] ## only exec dir
        f.write(ACBN0config.format(user=user, pseudo_dir=pseudo_dir, exec_dir=exec_dir))
#    if cluster == "local":  ## delete string matching pattern in place
#        for line in fileinput.input("ACBN0.config", inplace = True):
#            if not re.search(r'\bengine_dir\b', line):
#                print(line, end="")

    with open("ID","w") as f:
        f.write("")
    subprocess.run(['python3', 'ACBN0.py'])
    prefix_aflowpi=""
    with open("ID","r") as f:
        prefix_aflowpi=f.readline().split("_")[0]
    print('ID',prefix_aflowpi)

    filedir="ACBN0/{prefix}/ACBN0_{prefix}_0001/".format(prefix=prefix)

    ACBN0files=["_"+prefix_aflowpi+"_01.sh","_"+prefix_aflowpi+"_02.sh","_"+prefix_aflowpi+"_03.sh"]
    ACBN0calcs=[["scf"],["nscf","dos","pdos"],["scf","nscf_seq","pdos"]]

    for i,AFLOWpiID in enumerate(ACBN0files):
        if local:
            writeAFLOWpi_sh(calcs=ACBN0calcs[i],user=user,cluster="local",job_name="{prefix}_ACBN_{i}".format(prefix=prefix,i=i), nodes=nodes, ncpus=ncpus, QE_parflags="",AFLOWpiID=AFLOWpiID.lstrip("_").rstrip(".sh"),AFLOWpiPATH=filedir)
        else:
            writeAFLOWpi_sh(calcs=ACBN0calcs[i],user=user,cluster=cluster,job_name="{prefix}_ACBN_{i}".format(prefix=prefix,i=i), nodes=nodes, ncpus=ncpus, QE_parflags="",AFLOWpiID=AFLOWpiID.lstrip("_").rstrip(".sh"),AFLOWpiPATH=filedir)
    for i,AFLOWpiID in enumerate(ACBN0files):
        AFLOWpython_file=filedir+AFLOWpiID.rstrip(".sh")+".py"
        if local:
            modify_AFLOWpi_pyfile(AFLOWpython_file,cluster,local=True)
        else:
            modify_AFLOWpi_pyfile(AFLOWpython_file,cluster)



def runAFLOW_ACBN0(prefix="",local=False,mode="default"):
    import ast
    def readUvals():
        lastU_data=""
        with open("Uvals.log","r") as Ulog:
            for line in Ulog:
                pass
            lastU_data = line.rstrip()
        Uvals = ast.literal_eval(lastU_data)
        return Uvals

    if mode == "read_Uvals":
        try:
            Uvals=readUvals()
        except:
            with open("ID","r") as f:
                prefix_aflowpi=f.readline().split("_")[0]
            AFLOWdir="ACBN0/{prefix}/ACBN0_{prefix}_0001/".format(prefix=prefix)
            Uvals_log=prefix_aflowpi+"_03_uValLog.log"
            shutil.copyfile(AFLOWdir+Uvals_log, "Uvals.log")
            Uvals=readUvals
        return Uvals

    with open("ID","r") as f:
        prefix_aflowpi=f.readline().split("_")[0]
    AFLOWdir="ACBN0/{prefix}/ACBN0_{prefix}_0001/".format(prefix=prefix)
    os.chdir(AFLOWdir)
    ACBN0pyfiles=["_"+prefix_aflowpi+"_01.py","_"+prefix_aflowpi+"_02.py","_"+prefix_aflowpi+"_03.py"]
    for py_file in ACBN0pyfiles:
        subprocess.run(['python3', py_file])
    finish_acbn0=False
    notconverged=False
    LOGfile="../AFLOWpi/LOG.log"
    while not finish_acbn0:
        subprocess.run(['python3', ACBN0pyfiles[-1]])  ## keep sending last scfuj script if not converged.
        with open(LOGfile, "r") as log:
            for line in log:
                if re.search("Completed scfuj convergence for final U values",line):
                    finish_acbn0=True
                if re.search("Maximum no. of iterations reached. scfuj did not converge",line):
                    finish_acbn0=True
                    notconverged=True
    if finish_acbn0 and not notconverged:
        Uvals_log=prefix_aflowpi+"_03_uValLog.log"
        shutil.copyfile(Uvals_log, "../../../Uvals.log") ## copyfile to main folder
        os.chdir("../../../")
        Uvals=readUvals()
        return Uvals
    else:
        raise("ACBN0 Uvals not converged or incomplete, check what happened.")
   


