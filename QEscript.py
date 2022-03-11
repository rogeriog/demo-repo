
#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
Give in structure_filepath the path of the structure you want to run calculation,
may be a cif file, a xsf, POSCAR (VASP), Quantum espresso input or Quantum espresso
output file. If you specify as a list, structures will run in sequence.
"""
#########################################################################################################

structures_filepath=["structures/structure1.scf.in","structures/structure2.cif"]

#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
process_structure function determines how we read strucures in the files of structures_filepath,
and in there we can also specify alterations on the structures we are reading, such as change one
type of atom.
"""
#########################################################################################################

def process_structure(struct_filepath):
    ## never change these two lines 
    from MinFlow.special_functions import read_structure
    structure=read_structure(struct_filepath)

    ## extra modifications depending on what you need
    from MinFlow.special_functions import substitute_atom
#    structure=substitute_atom(structure,'Ba','Sr')
    return structure


#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
If automatic_prefix is chosen, prefix for calculation will be obtained from the
structure composition itself. Beware of possible prefix repetition.
"""
#########################################################################################################

automatic_prefix=False

#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
Give the prefix for the structure to be calculated, you should give a list if 
structure_filepath is also a list, prefix will be used to produce individual
folder for the calculation. Use sameFolderRun=True to not create a new folder.
"""
#########################################################################################################

prefixes=["str1","str2"]
sameFolderRun=False


#########################################################################################################
#########################################################################################################
#########################################################################################################
"""

cluster_tag : Specify cluster settings through tag as specified previously on MinFlow/cluster_sets.py

user : specify username in the cluster, that is necessary sometimes to find PP folder for calculation

nodes : nodes you want to allocate for calculation

ncpus : number of cpus per node to be allocated

PPsufix : Every type of pseudopotential has an specific termination, use this to specify
what source of pseudopotential you want to use, they have to be present in the
pseudopotential folder of the cluster that is specified in MinFlow/cluster_sets.py

kgrid : if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 

kgrid_shift : Specify kgrid shift

ecutwfc : wavefunction cutoff for calculation, as determined from convergence calculation
in the structure of interest.

scale_cutrho : number to multiply ecutwfc to obtain ecutrho, the cutoff on charge density.

kpath_density : if bands calculation specify kpath density, default is 32.

run_specs : leave as None, other options are "local", to run everything locally, usual for very
fast runs, and also "fermi_gpu_relax", specific due to problems in relaxation in fermi_gpu.

calculation : all traditional calculations available in QE such as relax, vc-relax, scf, nscf and bands,
there are also two special cases implemented:
   - 'acbn0' to determine hubbard values using acbn0 method, specify cluster tag which gives the PPs used
     for AFLOWpi preferentially, a previous scf run must be present because it read prefix.scf.in for the 
     structure.
   - 'absorption' using simple package to obtain optical data, it will calculate scf, nscf with high density
   of kpoints and then calculate optical parameters with simple. It requires NC PP, therefore make sure
   to specify the right cluster tag and PPsufix. You can specify kgrid for absorption with the key "kgrid_for_absorption"
   and it must be a tuple of (k1,k2,k3), if you dont specify it will just duplicate the kgrid given in calc_set.
   You can also pass custom parameters like "simple_dict", "simpleip_dict" and "simple_mode", which will specify
   the parameters of simple and simpleip calculations and if you want to start from a previous scf or nscf run.

   - 'bands', in the case of band calculations, the kpath we want to integrate over sometimes is not exactly
   of the structure we are calculating, that is the case in doped structure in which deformation slightly breaks
   symmetry. In these cases its useful to specify 'structure_for_bands' and give a structure filename to read
   a pristine structure to take the integration path from.

read_Uvals : if acbn0 calculation was performed and you want to read Uvals from that calculation, set it to True.
"""
#########################################################################################################


## few steps of a relax calculation with coarse kgrid and lower cutoff
calc_relax1={
'cluster_tag':'cesup_fermi',
'user':'raglamai',
'nodes':1,
'ncpus':24,
'PPsufix':'gbrv_us.upf',
'calculation':'vc-relax',
'kgrid':5,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
'kgrid_shift':(0,0,0),  ## kgrid shift 
'ecutwfc':50,
'scale_cutrho':6,
#'run_specs':'local', ## if you want to run locally
##'run_specs':'fermi_gpu_relax', ## specific for gpu relaxation

'CUSTOM':{'CONTROL':{'nstep':10},
          'SYSTEM':{'nspin':2,'input_dft':'pbe', 'degauss':0.002},
          'ELECTRONS':{'mixing_beta':0.3},
          'IONS':{},
          'CELL':{}}
}

### a relax calculation with converged kgrid and cutoff
calc_relax2={
'cluster_tag':'cesup_fermi',
'user':'raglamai',
'nodes':1,
'ncpus':24,
'PPsufix':'gbrv_us.upf',
'calculation':'vc-relax',
'kgrid':10,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
'kgrid_shift':(0,0,0),  ## kgrid shift 
'ecutwfc':50,
'scale_cutrho':6,
#'run_specs':'local', ## if you want to run locally
###'run_specs':'fermi_gpu_relax', ## specific for gpu relaxation
#
'CUSTOM':{'CONTROL':{'nstep':100},
          'SYSTEM':{'nspin':2,'input_dft':'pbe', 'degauss':0.002},
          'ELECTRONS':{'mixing_beta':0.3},
          'IONS':{},
          'CELL':{}}
}
#

## scf example calculation, previous relaxed structure will be read automatically if they
## are both in the calc_sets variable at the end. If you want to restart from scf calculation
## though, you should specify relaxed structure file (out.prefix.relax, for example), in structure_filepath.
calc_scf={
'cluster_tag':'cesup_fermiNC',
'user':'raglamai',
'nodes':1,
'ncpus':24,
'PPsufix':'gbrv_us.upf',
'calculation':'scf',
'kgrid':15,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
'kgrid_shift':(0,0,0),  ## kgrid shift 
'ecutwfc':50,
'scale_cutrho':6,
#'run_specs':'local', ## if you want to run locally
#'run_specs':'dryrun', ## dryrun for acbn0

'CUSTOM':{'CONTROL':{'nstep':100},
          'SYSTEM':{'nspin':2,'input_dft':'pbe', 'degauss':0.002},
          'ELECTRONS':{'mixing_beta':0.3},
          'IONS':{},
          'CELL':{}}
}
## example for ACBN0 calculation 
#calc_set_acbn0={
#'cluster_tag':'cesup_fermiNC',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'calculation':'acbn0',
#
#'kgrid':10,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift 
#'ecutwfc':70,
#'scale_cutrho':4,
##'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':100},
#          'SYSTEM':{'input_dft':'pbesol', 'degauss':0.02},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}
#
#calc_relax_U={
##'prefix':'NBNOCo_cU',
#'read_Uvals':True,
#
#'cluster_tag':'cesup_fermi',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'PPsufix':'gbrv_us.upf',
#'calculation':'vc-relax',
#'kgrid':10,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift 
#'ecutwfc':50,
#'scale_cutrho':6,
##'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':50},
#          'SYSTEM':{'nspin':2,'input_dft':'pbesol', 'degauss':0.02},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}
#
#calc_scf_U={
##'prefix':'NBNOCo_cU',
#'read_Uvals':True,
#
#'cluster_tag':'cesup_fermi2',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'PPsufix':'gbrv_us.upf',
#'calculation':'scf',
#'kgrid':15,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift 
#'ecutwfc':50,
#'scale_cutrho':6,
##'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':50},          
#          'SYSTEM':{'nspin':2,'input_dft':'pbesol', 'degauss':0.02},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}
#
#
#
#calc_nscf_U={
##'prefix':'NBNOCo_cU',
#'read_Uvals':True,
#
#'cluster_tag':'cesup_fermi2',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'PPsufix':'gbrv_us.upf',
#'calculation':'nscf',
#'kgrid':(8,8,8),  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift 
#'ecutwfc':50,
#'scale_cutrho':6,
##'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':50},          
#          'SYSTEM':{'nspin':2,'input_dft':'pbesol', 'degauss':0.02},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}


## nscf example calculation, only kgrid duplicated, nbnd  will be taken from previous
## relax or scf calculation, if you want to start from nscf calculation, then give
## also nbnd.
calc_nscf={
'cluster_tag':'cesup_fermi',
'user':'raglamai',
'nodes':1,
'ncpus':24,
'PPsufix':'gbrv_us.upf',
'calculation':'nscf',
'kgrid':15,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
'kgrid_shift':(0,0,0),  ## kgrid shift 
#'nbnd':235, ##only necessary if not continuing from previous scf or relax calculation
'ecutwfc':50,
'scale_cutrho':6,
#'run_specs':'local', ## if you want to run locally

'CUSTOM':{'CONTROL':{'nstep':100},
          'SYSTEM':{'nspin':2,'input_dft':'pbe', 'degauss':0.002},
          'ELECTRONS':{'mixing_beta':0.3},
          'IONS':{},
          'CELL':{}}
}

## bands calculation also gets nbnd from previous calculation, if starting directly from
## this nbnd should be specified to include some conduction band states. kgrid variables
## are not relevant for this calculation
calc_bands={
#'read_Uvals':True,

'cluster_tag':'cesup_fermi',
'user':'raglamai',
'nodes':1,
'ncpus':24,
'PPsufix':'gbrv_us.upf',
'calculation':'bands',
#'structure_for_bands':'../KNbO3.cif',## for bands sometimes we want another structure to get integration path.
'kgrid':15,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
'kgrid_shift':(0,0,0),  ## kgrid shift 
#'nbnd':235, ##only necessary if not continuing from previous scf or relax calculation
'ecutwfc':50,
'scale_cutrho':6,
#'run_specs':'local', ## if you want to run locally

'CUSTOM':{'CONTROL':{'nstep':100},          
          'SYSTEM':{'nspin':2,'input_dft':'pbe', 'degauss':0.002},
          'ELECTRONS':{'mixing_beta':0.3},
          'IONS':{},
          'CELL':{}}
}



## example for absorption calculation 
#calc_set_absorption={
#'cluster_tag':'cesup_fermi',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'PPsufix':'ONCV_PBE_sr.upf',
#'calculation':'absorption',
#'kgrid':15,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift
#'kgrid_for_absorption':(10,10,10)
#'ecutwfc':70,
#'scale_cutrho':4,
##'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':100},
#          'SYSTEM':{'input_dft':'pbe'},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}

## example for ACBN0 calculation 
#calc_set_acbn0={
#'cluster_tag':'cesup_fermiNC',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'calculation':'acbn0',
#
#'kgrid':10,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift 
#'ecutwfc':70,
#'scale_cutrho':4,
##'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':100},
#          'SYSTEM':{'input_dft':'pbe'},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}

## example of set to read hubbard values from acbn0 to initialize scf
#calc_set_scf_hubbard={
#'prefix':'Si2U',
#'read_Uvals':True,
#
#'cluster_tag':'cesup_fermi',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'PPsufix':'gbrv_us.upf',
#'calculation':'scf',
#'kgrid':10,  ## if integer given it is considered k_cutoff_length, otherwise use (x,x,x) 
#'kgrid_shift':(0,0,0),  ## kgrid shift 
#'ecutwfc':40,
#'scale_cutrho':6,
###'run_specs':'local', ## if you want to run locally
#
#'CUSTOM':{'CONTROL':{'nstep':100},
#          'SYSTEM':{'input_dft':'pbe'},
#          'ELECTRONS':{'mixing_beta':0.3},
#          'IONS':{},
#          'CELL':{}}
#}


#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
PP_sets will required less information, basically cluster related information
and pp_runs in which you should specify "pp", "bands" or combinations as "pp+bands"
If you want to perform bands first run a bands calculation.

"""
#########################################################################################################

## specify pp_runs
PP_sets1={
'cluster_tag':'cesup_fermi',
'user':'raglamai',
'nodes':1,
'ncpus':24,
'pp_runs':"pp+bands",
#'run_specs':'local',
}

## if you have more than one type of calculation in same structure folder
## you should specify the prefix of the one you want to calculate.
#PP_sets2={
#'cluster_tag':'cesup_fermi',
#'user':'raglamai',
#'nodes':1,
#'ncpus':24,
#'prefix':'Si2U',
#'pp_runs':"pp",
#'run_specs':'local',
#}

#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
Once you declare the dictionaries of calc_set, you should organize them in the following variable, calc_sets;
the variable must have two levels, that will determine which calculation will run together in same submission.
For example:
calc_sets=[[calc_set1],[calc_set2,calc_set3,PP_sets1,PP_sets2]]
will run calc_set1 in first submission, then send another jobs for calc_set2, calc_set3 and PP_sets1, PP_sets2 in
the second submission.

If you wanna ignore all calc_sets and go directly to post-processing, declare
calc_sets=[]

"""
#########################################################################################################

# for this calculation, relax1 and relax2 will be submitted in separated jobs, scf onwards will
# all be executed in the same job.
#calc_sets=[[calc_relax1],[calc_scf_dry,calc_set_acbn0,calc_relax_U],[calc_scf_U,calc_nscf_U,calc_bands_U,PP_sets1]]
calc_sets=[[calc_relax1],[calc_relax2],[calc_scf,calc_nscf,calc_bands,PP_sets1]]
# if you don't want to run any calculation, for example, if only doing postprocessing, keep it as:
# calc_sets=[]


#########################################################################################################
#########################################################################################################
#########################################################################################################
"""
PP_processing dict contains posprocessing which is usually made locally to generate, geometry and charge
information, pdos and band plots.
What is performed is defined by the string in 'pp' key, each is separated by '+' sign. Only perform bands
processing if you have run bands calculation followed by bands pp. 
There are also some options to customize plots of bands or PDOS.
"""
#########################################################################################################

## in this example, informations like charge geometry, gap, pdos and band plots will be generated,
## if there is more than one prefix you should specify through 'prefix' key.
PP_processing={
'pp':"charge_geometry+gap+pdos+bands",
#'pp':"bands",
'structure_for_bands':'../KNbO3.cif',
#'prefix':['NBNOCo_cU'],
'CUSTOMIZE_PLOT': {'fermi':0.0,
               'label':"",
               'color':"black",
               'RANGE_VB':6,
               'RANGE_CB':4,
             'redstates':[],  ## states colored red, you can always give range(92,96)
             'greenstates':[],  ## states colored green
             'bluestates':[], #range(188,192),  ## states colored blue, no states with []
             'labelsprojs':["","",""],  ## first red, second green, third blue
## notice if there is overlapping colors will blend
             'maxcontrast':True,  ## notice if maxcontrast=0, color will be as strong as the state proportion, 
             #usually this contribution is very low compared to total, therefore maxcontrast is default.
             'contrast_ratios':[1,1,1],  ## only if maxcontrast true # values should be in 0-1 range ## can use to
             ## obtain custom colors
             'plot_emass_fitting':False,
             'emass_calculation':False,
             }
}


###################################### SCRIPT, NOT TO CHANGE ############################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
from MinFlow.AFLOWpi_ACBN0 import *
from MinFlow.QEcalculator import QEcalc
import MinFlow.QEcalculator as QEcalculator
import MinFlow.QEpp as QEpp
import MinFlow.cluster_sets as cluster_sets
from MinFlow.special_functions import read_structure,substitute_atom
from MinFlow.special_functions import getPseudoDict,get_name_from_structure,get_data_from_outQE
import MinFlow.kgrid as kgrid_package
import numpy as np
import glob,sys,os

if isinstance(structures_filepath,list) or isinstance(structures_filepath,tuple):
   structures=structures_filepath
   prefixes=prefixes
else:
   structures=[structures_filepath]
   prefixes=[prefixes]
for i, i_structure in enumerate(structures):
    try:
        structure=process_structure(i_structure)

    except:
        print(f"{i_structure} not readable.")
        continue ## structure is not readable case, go to next
    if automatic_prefix:
        prefix=get_name_from_structure(structure)
    else:
        prefix=prefixes[i]
    if not sameFolderRun:
        ### make a prefix folder and cd to it ###
        try:
            os.mkdir(prefix) ## make prefix folder
        except OSError as error: 
            print(error)
        try:
            os.chdir(prefix)  ## enter prefix folder
        except FileNotFoundError as error:
            print(error)
    
    run_first=False
    relaxed_structure=None
    nbnd=None
    for i in range(len(calc_sets)): ## sets of calculation that will run individually
        calcs=[]
        for j in range(len(calc_sets[i])): ## sets of calculation that will run in the same batch
            print(i,j,calc_sets[i][j].keys())
            if 'pp_runs' in calc_sets[i][j].keys(): ## case for PP
                if relaxed_structure is not None:
                    structure=read_structure(relaxed_structure)
                if "prefix" in calc_sets[i][j].keys():
                    prefix=calc_sets[i][j]['prefix']
                if 'pp' in calc_sets[i][j]['pp_runs'].split("+"):
                    try:
                        PPcalcs=QEpp.preparePP(prefix=prefix,
                        nelectron=get_data_from_outQE('out.'+prefix+'.vc-relax','nelectron'),
                        structure=structure)
                    except FileNotFoundError:
                        try:
                            PPcalcs=QEpp.preparePP(prefix=prefix,
                            nelectron=get_data_from_outQE('out.'+prefix+'.relax','nelectron'),
                            structure=structure)
                        except FileNotFoundError:
                            try:
                                PPcalcs=QEpp.preparePP(prefix=prefix,
                                nelectron=get_data_from_outQE('out.'+prefix+'.scf','nelectron'),
                                structure=structure)
                            except FileNotFoundError:
                                print("It's necessary a previous relax or scf calculation to obtain nelectrons\
                                for HOMO LUMO PP calculation")

                    calcs.extend(PPcalcs)
                if 'bands' in calc_sets[i][j]['pp_runs'].split("+"):
                    PPcalcs_bands=QEpp.preparePP_bands(prefix=prefix)
                    calcs.extend(PPcalcs_bands)



            else: ## case for all other calculations 
                cluster=calc_sets[i][j]['cluster_tag']
                user=calc_sets[i][j]['user']
                nodes=calc_sets[i][j]['nodes']
                ncpus=calc_sets[i][j]['ncpus']
                if "prefix" in calc_sets[i][j].keys():
                    prefix=calc_sets[i][j]['prefix']
                if "PPsufix" in calc_sets[i][j].keys():
                    PPsufix=calc_sets[i][j]['PPsufix']
                else:
                    PPsufix=None
                calculation=calc_sets[i][j]['calculation']
                if "kgrid" in calc_sets[i][j].keys():
                    kgrid=calc_sets[i][j]['kgrid']
                if 'run_specs' in calc_sets[i][j].keys():
                    run_specs=calc_sets[i][j]['run_specs']
                else:
                    run_specs=None
                if 'kpath_density' in calc_sets[i][j].keys():
                    kpath_density=calc_sets[i][j]['kpath_density']
                else:
                    kpath_density=32
                if "kgrid_shift" in calc_sets[i][j].keys():
                    kgrid_shift=calc_sets[i][j]['kgrid_shift']
                if "ecutwfc" in calc_sets[i][j].keys():
                    ecutwfc=calc_sets[i][j]['ecutwfc']
                if "scale_cutrho" in calc_sets[i][j].keys():
                    scale_cutrho=calc_sets[i][j]['scale_cutrho']
                if "nbnd" in calc_sets[i][j].keys():
                    nbnd=calc_sets[i][j]['nbnd']
                
                if isinstance(kgrid,list) or isinstance(kgrid,tuple):
                    ktuple=kgrid+kgrid_shift
                else:
                    ktuple=kgrid_package.calc_kpt_tuple(structure,cutoff_length=kgrid)+(0,0,0)
                #########################

                ### get nbnd for nscf or bands calculation
                if calculation in ["nscf","bands"]:
                    if nbnd is None:
                        for i_calc in ["relax","vc-relax","scf"]:
                            try:
                                print(f"Reading nbnd from {'out.'+prefix+'.'+i_calc}.")
                                nbnd=get_data_from_outQE('out.'+prefix+'.'+i_calc,'nbnd')
                                print(f"nbnd is {nbnd}")
                                break
                            except Exception as e:
                                print("couldnt read nbnd!",e)
                        nbnd=nbnd+0.2*nbnd
                    else:
                        print("nbnd not found and not given, specify nbnd beforehand, otherwise,\
                        no conduction band will be present in nscf or bands calculation.")
                    


                ## to read structure from previous relaxation if available
                if relaxed_structure is not None:
                    structure=read_structure(relaxed_structure)
                if calculation not in ["acbn0"]:
                    pseudodict=getPseudoDict(structure,cluster=cluster, user=user,
                    enforce_PPsufix=PPsufix)
                    calcdict={'prefix': prefix,
                    'calculation':calculation,
                    'ecutwfc':ecutwfc,
                    'pseudodict':pseudodict,
                    'structure':structure,
                    'scale_cutrho':scale_cutrho, 
                    'ktuple':ktuple,
                    'nbnd':nbnd,
                    'CUSTOM':{'CONTROL':calc_sets[i][j]['CUSTOM']['CONTROL'],
                              'SYSTEM':calc_sets[i][j]['CUSTOM']['SYSTEM'],
                              'ELECTRONS':calc_sets[i][j]['CUSTOM']['ELECTRONS'],
                              'IONS':calc_sets[i][j]['CUSTOM']['IONS'],
                              'CELL':calc_sets[i][j]['CUSTOM']['CELL'],}
                    }
                    if 'read_Uvals' in calc_sets[i][j].keys():
                        if calc_sets[i][j]['read_Uvals'] :
                            Uvals=runAFLOW_ACBN0(prefix=prefix,mode="read_Uvals")
                            calcdict['CUSTOM']['SYSTEM']['lda_plus_u']=True
                    ## implement hubbard values for each element, check if ase version allows.
                            for idx, element in enumerate(pseudodict['potentials'].keys()):
                                calcdict['CUSTOM']['SYSTEM']['Hubbard_U('+str(idx+1)+')']=Uvals.get(element, 0.0)
                #########################
                    if calcdict['calculation'] == 'bands':
                        if 'structure_for_bands' in calc_sets[i][j].keys():
                            structure_bands=read_structure(calc_sets[i][j]['structure_for_bands'])
                            print('diditrun',structure_bands)
                        else:
                            structure_bands=structure
                        path=structure_bands.cell.bandpath().interpolate(density=kpath_density) 
                        print(path)
                        calcdict['kpath']=path
                    if calculation == "absorption":
                        if "kgrid_for_absorption" in calc_sets[i][j].keys(): 
                            kgrid_for_absorption=calc_sets[i][j]["kgrid_for_absorption"]
                        else:
                            print("You should give kgrid for absorption for absorption calculation.")
                            print("kgrid_for_absorption will be 2*calcdict kgrid")
                            kgrid_for_absorption=tuple(np.array(ktuple[:3])*2)
                        if "simple_dict" in calc_sets[i][j].keys(): 
                            simple_dict=calc_sets[i][j]['simple_dict']
                        else:
                            simple_dict=None

                        if "simpleip_dict" in calc_sets[i][j].keys(): 
                            simpleip_dict=calc_sets[i][j]['simpleip_dict']
                        else:
                            simpleip_dict=None

                        if "simple_mode" in calc_sets[i][j].keys(): 
                            simple_mode=calc_sets[i][j]['simple_mode']
                        else:
                            simple_mode=None
                        
                        if simple_mode is None:
                            abscalcs=QEcalculator.gencalc_absorption(calcdict,kgrid_for_absorption,
                                                                    simple_dict=simple_dict,                          
                                                                   simpleip_dict=simpleip_dict, 
                                                                   from_previous_scf=False, from_previous_nscf=False,
                                                                   cluster=cluster)
                        elif simple_mode == 'from_previous_scf':
                            abscalcs=QEcalculator.gencalc_absorption(calcdict,kgrid_for_absorption,
                                                                simple_dict=simple_dict,
                                                                simpleip_dict=simpleip_dict, 
                                                                from_previous_scf=True, from_previous_nscf=False,
                                                                cluster=cluster)
                        elif simple_mode == 'from_previous_nscf':
                            abscalcs=QEcalculator.gencalc_absorption(calcdict,kgrid_for_absorption,
                                                                simple_dict=simple_dict,
                                                                simpleip_dict=simpleip_dict, 
                                                                from_previous_scf=False, from_previous_nscf=True,
                                                                cluster=cluster)
                    if calculation not in ["absorption"]:
                        calc=QEcalc(**calcdict) ## just to print scf.in
                        inp=calc.write_input()
                        calcs.append(calc)
                    else:
                        print('do we get here?',abscalcs,len(abscalcs))
                        calcs.extend(abscalcs)
                    if calculation == "relax" or calculation == 'vc-relax':
                        relaxed_structure=calc.output_file
                else:
                    if calculation == "acbn0":
                        if 'run_specs' in calc_sets[i][j].keys():
                            print('test',user,cluster)
                            if calc_sets[i][j]['run_specs']=='local':
                                initAFLOW_ACBN0(prefix=prefix,user="raglamai",
                                cluster=cluster, local=True,nodes=nodes,ncpus=ncpus) ### fermiNC takes the PPs from acbn0
                                Uvals=runAFLOW_ACBN0(prefix=prefix,local=True)
                                continue ## shdict stuff should not run in this
                        else:
                            print('test2',user,cluster)
                            initAFLOW_ACBN0(prefix=prefix,user="raglamai",
                            cluster=cluster, nodes=nodes,ncpus=ncpus) ### fermiNC takes the PPs from acbn0
                            Uvals=runAFLOW_ACBN0(prefix=prefix)
                            continue ## shdict stuff should not run in this 

        shdict={'user':user,'job_name':calcs[0].prefix,'calcs':calcs, 'cluster':cluster,
                'nodes':nodes,'ncpus':ncpus}
        sh_file=QEcalculator.write_sh(**shdict)
        #### run sh file ############
        if run_specs is None:
            QEcalculator.run(sh_file)
        elif run_specs == "local":
            shdict={'user':user,'job_name':calcs[0].prefix,'calcs':calcs, 'cluster':'local',
                'nodes':nodes,'ncpus':ncpus}
            sh_file=QEcalculator.write_sh(**shdict)
            QEcalculator.run(sh_file,local=True)
        elif run_specs == "dryrun":
            shdict={'user':user,'job_name':calcs[0].prefix,'calcs':calcs, 'cluster':'local',
                'nodes':nodes,'ncpus':ncpus}
            sh_file=QEcalculator.write_sh(**shdict)
#            QEcalculator.run(sh_file,local=True)
        elif run_specs == "fermi_gpu_relax":
            QEcalculator.run(sh_file,fermi_gpu=True)
        run_first=True
                

    #### ppnscf ###
    ###############
    if PP_processing:
        directory_path = os.getcwd()
        print("My current directory is : " + directory_path)
        mode=PP_processing['pp']
        if 'structure_for_bands' in PP_processing.keys():
            structure=read_structure(PP_processing['structure_for_bands'])
        band_indexes=[x for x in structure.cell.bandpath().path if x != ',']
        if 'prefix' in PP_processing.keys():
            if isinstance(PP_processing['prefix'],list):
                for item in PP_processing['prefix']:
                    QEpp.PP_processing(item, mode=mode, indexes=band_indexes,**PP_processing['CUSTOMIZE_PLOT'])
            else:
                QEpp.PP_processing(PP_processing['prefix'], mode=mode, indexes=band_indexes, **PP_processing['CUSTOMIZE_PLOT'])
        else:
            QEpp.PP_processing(prefix, mode=mode, indexes=band_indexes, **PP_processing['CUSTOMIZE_PLOT'])



"""

## run acbn0 in each structure
## relax structure with hubbard
## scf,nscf,bands, pp on relaxed structure
## band unfolding (is it necessary? Ba doped new unit cell)
## absorption
for structure_file in glob.glob("structures/*"):
    try:
        structure=read_structure(structure_file)
        ## structure=substitute_atom(structure,'K','Na')
    except:
        continue
    if automatic_prefix:
        prefix=get_name_from_structure(structure)
    else:
        prefix=prefixes[i]
    prefix=structure_file.split("/")[-1].split(".")[0].replace("K","N")[:-1]
    ### make a prefix folder and cd to it ###
    try:
        os.mkdir(prefix)
    except OSError as error: 
        print(error)
    try:
        os.chdir(prefix)
    except FileNotFoundError as error:
        print(error)
    ##########################################

#    if prefix == 'KBNOCo_c2':
#        print('later on')
#        sys.exit()

    pseudodict=getPseudoDict(structure,cluster="cesup_fermiNC",
                        user="raglamai") #,enforce_PPsufix="gbrv_us.upf")
    ktuple=kgrid.calc_kpt_tuple(structure,cutoff_length=5)+(0,0,0)
   #########################
    calcdict={'prefix': prefix,
'calculation':'vc-relax',
'ecutwfc':70,
'pseudodict':pseudodict,
'structure':structure,
'scale_cutrho':4, 
'ktuple':ktuple,
'CUSTOM':{'CONTROL':{'nstep':5,},
          'SYSTEM':{'input_dft':'pbesol'},
          'ELECTRONS':{},
          'IONS':{},
          'CELL':{}}
}
   #########################
    calc= QEcalc(**calcdict) ## just to print scf.in
    inp=calc.write_input()
    shdict={'user':"raglamai",'job_name':calc.prefix,
    'calcs':[calc], 'cluster':"cesup_fermiNC", 'nodes':1,'ncpus':24}
    sh_file2=QEcalculator.write_sh(**shdict)
    #### run sh file ############
    if prefix != 'NBNO_cU':
       QEcalculator.run(sh_file2)#, local=True)
    calc.get_relaxed_structure(calc.output_file) ## to get relax
    calcdict['calculation']='scf'
    calc.relaxed_structure.cell=np.round(calc.relaxed_structure.cell,4)
    calcdict['structure']=calc.relaxed_structure
    calc=QEcalc(**calcdict)
    calc.write_input()
    ############

    ###############################################################
    #########  GET Uvals from ACBN0 ###############################
    ###############################################################
    Uvals=None
    if prefix != 'specific':
        from MinFlow.AFLOWpi_ACBN0 import *
        ### using pp from ACBN0 seems to converge faster, but pseudodojo
        ## scalar relativistic is more trustable
        ## it will look for prefix.scf.in, it must be written
        ## PPs are determined from cluster option
        initAFLOW_ACBN0(prefix=calc.prefix,user="raglamai",
        cluster="cesup_fermiNC", nodes=1,ncpus=24) ### fermiNC takes the PPs from acbn0
        Uvals=runAFLOW_ACBN0(prefix=calc.prefix)
        print(Uvals) ## we should use these uvals with pseudodojo and gbrv to compare
#    if prefix == '':
#        Uvals={'Nb': 0.287554, 'K': 0.007203, 'O': 6.530089, 'N': 0.693253}
    ###############################################################
    if prefix == 'specific':
        from MinFlow.AFLOWpi_ACBN0 import *
        Uvals=runAFLOW_ACBN0(prefix=calc.prefix,mode="read_Uvals" )
    print(Uvals)
    #sys.exit()
    ########## calculate bands with obtained U vals #########
    ### vc-relax with Uvals ##
    pseudodict=getPseudoDict(structure,cluster="cesup_fermi",
                        user="raglamai" ,enforce_PPsufix="gbrv_us.upf")
    ktuple=kgrid.calc_kpt_tuple(structure,cutoff_length=20)+(0,0,0)
#    if prefix == '':
#        structure = read_structure('KBNOCo_cU.structure')
    calcdict={'prefix': prefix,
'calculation':'scf',
'ecutwfc':50,
'pseudodict':pseudodict,
'structure':structure,
'scale_cutrho':6, 
'ktuple':ktuple,
'CUSTOM':{'CONTROL':{'nstep':50},
          'SYSTEM':{'degauss':0.1},
          'ELECTRONS':{'mixing_beta':0.08},
          'IONS':{},
          'CELL':{}}
}
    ktuple=kgrid.calc_kpt_tuple(structure,cutoff_length=15)+(0,0,0)
    calcdict['ktuple']=ktuple ## return to old ktuple
    calcdict['prefix']=prefix+'U'
    calcdict['calculation']='vc-relax'
    calcdict['CUSTOM']['SYSTEM']['lda_plus_u']=True
    ## implement hubbard values for each element, check if ase version allows.
    for idx, element in enumerate(pseudodict['potentials'].keys()):
        calcdict['CUSTOM']['SYSTEM']['Hubbard_U('+str(idx+1)+')']=Uvals.get(element, 0.0)

    calc_relU=QEcalc(**calcdict)
    calc_relU.write_input()
    ############################
#    sys.exit()
    #### set sh file for relax ###
    shdict={'user':"raglamai",'job_name':calc_relU.prefix,
    'calcs':[calc_relU], 'cluster':"cesup_fermiNC", 'nodes':1,'ncpus':24}
    sh_file2=QEcalculator.write_sh(**shdict)
    #### run sh file ############
    if prefix != '':
       QEcalculator.run(sh_file2)#, local=True)
    ##############################
    
    calcdict['CUSTOM']['SYSTEM']['degauss']=0.002
    calcdict['CUSTOM']['ELECTRONS']['mixing_beta']=0.3
    
    ### scf ###
    calc_relU.get_relaxed_structure(calc_relU.output_file) ## to get relax
    calcdict['calculation']='scf'
    calc_relU.relaxed_structure.cell=np.round(calc_relU.relaxed_structure.cell,4)
    calcdict['structure']=calc_relU.relaxed_structure
    calc_scfU=QEcalc(**calcdict)
    calc_scfU.write_input()
    ############

    ### nscf ###
    #calc_relU.get_relaxed_structure(calc_scfU.output_file) ## to get nbnd 
    calcdict['calculation']='nscf'
    nscf_ktuple=tuple(2*x for x in ktuple[:3])+(0,0,0) ## duplicate kgrid
    calcdict['CUSTOM']['SYSTEM']['nbnd']=calc_relU.nbnd+int(0.2*calc_relU.nbnd)
    calcdict['ktuple']=nscf_ktuple
    calc_nscfU=QEcalc(**calcdict)
    calc_nscfU.write_input()
    #############

    #### ppnscf ###
    PPcalcs=QEpp.preparePP(prefix=calc_scfU.prefix,nelectron=calc_relU.nelectron,
    structure=calc_relU.structure)
    ###############

    #### set sh file for relax ###
#    shdict={'user':"raglamai",'job_name':calc_relU.prefix,
#    'calcs':[calc_relU], 'cluster':"cesup_fermi", 'nodes':1,'ncpus':24}
#    sh_file=QEcalculator.write_sh(**shdict)
    ##############################

    #### bands calculation ####
    calcdict['calculation']='bands'
    path=calc_scfU.structure.cell.bandpath().interpolate(density=32) 
    calcdict['kpath']=path
    calc_bandsU= QEcalc(**calcdict) ## new calculator for bandstructure with Uvals
    calc_bandsU.write_input()
    
    ### PP bands ###
    PPcalcs_bands=QEpp.preparePP_bands(calc_scfU.prefix)
    ################
    
    ### setup sh file and run ##
    shdict['job_name']=calc_scfU.prefix
    shdict['calcs']=[calc_scfU,calc_nscfU,calc_bandsU]+PPcalcs+PPcalcs_bands
    #shdict['calcs']=[calc_nscfU,calc_bandsU]+PPcalcs+PPcalcs_bands
    sh_file3=QEcalculator.write_sh(**shdict)
    if prefix != '':
        QEcalculator.run(sh_file3)
    
    # post processing everything
    band_indexes=[x for x in calc_scfU.structure.cell.bandpath().path if x != ',']
    QEpp.PP_processing(calc_scfU.prefix, mode="charge_geometry+gap+pdos+bands",
                       indexes=band_indexes)
    #QEpp.PP_processing(calc_scfU.prefix, mode="bands",
    #                   indexes=band_indexes)
    ########################################################

    os.chdir("../") ## have to go back to main folder to do next structure














### make supercell
from ase.build import make_supercell
P=[ [2,0,0],
    [0,2,0],
    [0,0,2] ]
supercell=make_supercell(structure,P)
# equiavlent to:
#supercell2=structure*(2,2,2)

#from MinFlow.special_functions import get_positions_symbol
#get_positions_symbol(supercell2,"Nb")
supercell.symbols[1]="Co"  ## index obtained by checking Nbs previously
pseudodict=getPseudoDict(supercell,cluster="cesup_fermi_gpu",user="raglamai"
                        ,enforce_PPsufix="gbrv_us.upf")
ktuple=kgrid.calc_kpt_tuple(supercell,cutoff_length=15)+(0,0,0)
prefix=str(structure.symbols)+"CoSP"
#prefix=str(structure.symbols)+"Co"

#### IF SPIN POLARIZED ###########
from MinFlow.special_functions import set_element_magnetization
## break symmetry by applying starting magnetization in some atom
magmoms=set_element_magnetization(supercell,"Nb")
supercell.set_initial_magnetic_moments(magmoms)
##################################

### relax supercell ###

###########################
calcSC= QEcalc(**calcdict)
## change nstep/convergence_thr or k-resolution
inputSC=calcSC.write_input()
## now setting to write sh file to use in qsub
shdict={'user':"raglamai",'job_name':symbol+"_relax",
'calcs':[calcSC], 'cluster':"local", 'nodes':1,'ncpus':12}
sh_file=QEcalculator.write_sh(**shdict)
QEcalculator.run(sh_file,local=True)
calcSC.get_relaxed_structure(calcSC.output_file)
#############################

### scf supercell ###
calcdict['calculation']='scf'
###########################
calcSC= QEcalc(**calcdict)
inputSC=calcSC.write_input()

from MinFlow.special_functions import get_pathPBZ
from MinFlow.special_functions import get_nk
import MinFlow.banduppy as banduppy

labels,pathPBZ=get_pathPBZ(structure)
##############################
unfold=banduppy.UnfoldingPath(
            supercell=P,   # How the SC latticevectors are expressed in the PC basis (should be a 3x3 array of integers)
            pathPBZ=pathPBZ,
            # Path nodes in reduced coordinates in the primitive BZ. if the segmant is skipped, put a None between nodes
            nk=get_nk(pathPBZ,64),
            # number of k-points in each non-skipped segment.  Or just give one number, if they are equal
             labels=labels )   # GMKGALHALMHK for Γ—M—K—Γ—A—L—H—A|L—M|H—K
###############################

kpointsPBZ=unfold.kpoints_SBZ_str()   # as tring  containing the k-points to beentered into the PWSCF input file  after 'K_POINTS crystal ' line. maybe transformed to formats of other codes  if needed


calcdict['kpath']=kpointsPBZ
## nspin, mag and others used in previous scf
#calcdict['input_data']=calcSC.input_data  
calcdict['calculation']='bands'

calcSC.get_relaxed_structure(calcSC.output_file)
##
print(calcSC.nbnd+int(0.2*calcSC.nbnd))
calcdict['CUSTOM']['SYSTEM']['nbnd']=calcSC.nbnd+int(0.2*calcSC.nbnd)
calcSC2= QEcalc(**calcdict)
inputSC2=calcSC2.write_input()

## now setting to write sh file to use in qsub
#shdict={'user':"raglamai",'job_name':prefix,
#'calcs':[calcSC,calcSC2], 'cluster':"cesup_fermi_gpu", 'nodes':1,'ncpus':24}
shdict={'user':"raglamai",'job_name':prefix,
'calcs':[calcSC2], 'cluster':"cesup_fermi_gpu", 'nodes':1,'ncpus':24}
sh_file=QEcalculator.write_sh(**shdict)
#QEcalculator.run(sh_file)#, local=True)
sys.exit() ####################
import pickle ##pickle necessary, probably have to run more than once
#fermi=QEpp.get_Efermi(prefix)
fermi=7.8407
spin_pol=True
#print("fermi",fermi)
#bands=""
#try:
#    bands=pickle.load(open("bandstructure.pickle","rb"))  
#except Exception as err:
bands=unfold.getBandstructure(prefix,"espresso",spin_polarized=spin_pol)
#    pickle.dump(bands,open("bandstructure.pickle","wb"))
#    print(err)
#sys.exit()
"""

