import math
import re, os
import glob
import numpy as np
import subprocess
import itertools
import matplotlib.pyplot as plt
from scipy import integrate
import MinFlow.banduppy as banduppy
from MinFlow.special_functions import generate_xsf_file, get_geometry_data 
## QEpp.py
class PPcalc():
    def __init__(self,prefix, calculation, inputfile, inputdata):
        self.prefix=prefix
        self.calculation=calculation
        self.inputdata=inputdata
        self.input_file=inputfile
        self.output_file=self.input_file.split("/")[0]+"/out."+self.input_file.split("/")[1].rsplit(".",1)[0]




def preparePP(prefix, nelectron=0,structure=0):
    PPcalcs=[]
    HOMO_band=int(math.ceil(nelectron/2))
    LUMO_band=HOMO_band+1
    try: 
        os.mkdir("PP")
    except OSError as error: 
        print(error)

    ## to obtain the xsf file for the relaxed structure and
    ## to calculate all relevant geometry information
#    try:
#        generate_xsf_file('out.'+prefix+'.scf')
#    except:
#        try:
#            generate_xsf_file('out.'+prefix+'.vc-relax')
#        except:
#            try:
#                generate_xsf_file('out.'+prefix+'.relax')
#            except Exception as e:
#                print(e)
#    try:
#        get_geometry_data(structure)
#    except Exception as e:
#        print(e)

## DOS
    inp="PP/{prefix}.dos.in".format(prefix=prefix)
    with open(inp,"w") as f:
        DOSpp="""&dos  
prefix = '{prefix}' 
outdir = './OUTPUT'  
ngauss = 0
Emin=-24,  Emax=24, DeltaE=0.001
fildos='./PP/{prefix}.dos'
/
"""
        f.write(DOSpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'dos',inp,DOSpp))
## CHARGE XSF
    inp="PP/{prefix}.charge3d.in".format(prefix=prefix)
    with open(inp,"w") as f:
        CHARGEpp="""
&inputpp
prefix = '{prefix}' 
outdir = './OUTPUT'  
filplot = './PP/{prefix}.charge'
plot_num = 0
/
&plot
nfile = 1
filepp(1) = './PP/{prefix}.charge'
weight(1) = 1.0
iflag = 3
output_format = 5
fileout = './PP/{prefix}.charge.xsf'
/
"""
        f.write(CHARGEpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'chargexsf',inp,CHARGEpp))
## CHARGE CUBE
    inp="PP/{prefix}.chargecube.in".format(prefix=prefix)
    with open(inp,"w") as f:
        CHARGECUBEpp="""&inputpp
/
&plot
nfile = 1
filepp(1) = './PP/{prefix}.charge'
weight(1) = 1.0
iflag = 3
output_format = 6
fileout = './PP/{prefix}.charge.cube'
/
"""
        f.write(CHARGECUBEpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'chargecube',inp,CHARGECUBEpp))
## SPIN
    inp="PP/{prefix}.spol.in".format(prefix=prefix)
    with open(inp,"w") as f:
        SPpp="""&inputpp
prefix = '{prefix}' 
outdir = './OUTPUT'  
filplot = './PP/{prefix}.spol'
plot_num = 6
/
&plot
nfile = 1
filepp(1) = './PP/{prefix}.spol'
weight(1) = 1.0
iflag = 3
output_format = 6
fileout = './PP/{prefix}.spol.cube'
/
"""
        f.write(SPpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'spin-polarization',inp,SPpp))
## VTOTAL num=11
    inp="PP/{prefix}.V11.in".format(prefix=prefix)
    with open(inp,"w") as f:
        V11pp="""&inputpp
prefix = '{prefix}' 
outdir = './OUTPUT'  
filplot = './PP/{prefix}.vtot11'
plot_num = 11
/
&plot
nfile = 1
filepp(1) = './PP/{prefix}.vtot11'
weight(1) = 1.0
iflag = 3
output_format = 6
fileout = './PP/{prefix}.vtot11.cube'
/
"""
        f.write(V11pp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'V11',inp,V11pp))
## VTOTAL num=1
    inp="PP/{prefix}.Vtotal.in".format(prefix=prefix)
    with open(inp,"w") as f:
        Vtotalpp="""&inputpp
prefix = '{prefix}' 
outdir = './OUTPUT'  
filplot = './PP/{prefix}.vtot'
plot_num = 1
/
&plot
nfile = 1
filepp(1) = './PP/{prefix}.vtot'
weight(1) = 1.0
iflag = 3
output_format = 6
fileout = './PP/{prefix}.vtot.cube'
/
"""
        f.write(Vtotalpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'Vtotal',inp,Vtotalpp))
    ## have to automatize avg.in
    if structure != 0 :  ## structure has to be given to get points to integrate
        inp="PP/avg.in"
        with open(inp,"w") as f:
            AVGpp="""1
PP/{prefix}.vtot11
1.D0
{points}
3
{windowsize}
"""
            f.write(AVGpp.format(prefix=prefix,points=int(structure.cell.cellpar()[2]*100),
                windowsize=structure.cell.cellpar()[2]))
            PPcalcs.append(PPcalc(prefix,'avg',inp,AVGpp))

##### HOMO LUMO ####
    for homoindex in range(4):  ## print 4 below
        inp="PP/{prefix}.HOMO{homoindex}.in".format(prefix=prefix,homoindex=homoindex)
        with open(inp,"w") as f:
            HOMOpp="""&inputpp
 prefix ='{prefix}',
 outdir = './OUTPUT'  
filplot = './PP/{prefix}.HOMO' ! this will be rewritten
plot_num = 7
 kpoint = 1,
 kband = {HOMO_band},   
/
&plot
 nfile= 1,
 filepp(1)='./PP/{prefix}.HOMO',
 weight(1)= 1.0,
 iflag= 3,
 output_format=6, 
 fileout='./PP/{prefix}.HOMO{homoindex}.cube'
/
"""
            f.write(HOMOpp.format(prefix=prefix,HOMO_band=HOMO_band-homoindex,homoindex=homoindex))
            PPcalcs.append(PPcalc(prefix,'homo{homoindex}'.format(homoindex=homoindex),inp,HOMOpp))
    for lumoindex in range(4):  ## print 4 above
        inp="PP/{prefix}.LUMO{lumoindex}.in".format(prefix=prefix,lumoindex=lumoindex)
        with open(inp,"w") as f:
            LUMOpp="""&inputpp
 prefix ='{prefix}',
 outdir = './OUTPUT'  
filplot = './PP/{prefix}.LUMO' ! this will be rewritten
plot_num = 7
 kpoint = 1,
 kband = {LUMO_band},   
/
&plot
 nfile= 1,
 filepp(1)='./PP/{prefix}.LUMO',
 weight(1)= 1.0,
 iflag= 3,
 output_format=6, 
 fileout='./PP/{prefix}.LUMO{lumoindex}.cube'
/
"""
            f.write(LUMOpp.format(prefix=prefix,LUMO_band=LUMO_band+lumoindex,lumoindex=lumoindex))
            PPcalcs.append(PPcalc(prefix,'lumo{lumoindex}'.format(lumoindex=lumoindex),inp,LUMOpp))
    try: 
        os.mkdir("PP/PDOS")
    except OSError as error: 
        print(error)

#### PDOS ############
    inp="PP/{prefix}.pdos.in".format(prefix=prefix)
    with open(inp,"w") as f:
        PROJpp="""&projwfc
prefix = '{prefix}' 
outdir = './OUTPUT'  
ngauss = 0
Emin=-24,  Emax=24, DeltaE=0.001
filpdos='./PP/PDOS/{prefix}'
/
"""
        f.write(PROJpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'pdos',inp,PROJpp))
    return PPcalcs


def preparePP_bands(prefix):
    PPcalcs=[]
#### BANDS  ############
    inp="PP/{prefix}.bands.in".format(prefix=prefix)
    with open(inp,"w") as f:
        BANDSpp="""&bands
prefix = '{prefix}' 
outdir = './OUTPUT'  
filband = './PP/band_{prefix}.dat'
/
"""
        f.write(BANDSpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'ppbands',inp,BANDSpp))
    try: 
        os.mkdir("PP/PDOS_BANDS")
    except OSError as error: 
        print(error)
##### PDOS ############
    inp="PP/{prefix}.pdos_bands.in".format(prefix=prefix)
    with open(inp,"w") as f:
        PROJpp="""&projwfc
prefix = '{prefix}' 
outdir = './OUTPUT'  
ngauss = 0
Emin=-24,  Emax=24, DeltaE=0.001
filpdos='./PP/PDOS_BANDS/{prefix}'
/
    """
        f.write(PROJpp.format(prefix=prefix))
        PPcalcs.append(PPcalc(prefix,'pdos_bands',inp,PROJpp))
    return PPcalcs


############################# PYTHON FUNCTIONS OVER PP DATA #####################################
def PP_processing(prefix,mode="charge_geometry+gap+pdos",**kwargs):
    mode=mode.split("+")
    print(kwargs,kwargs.get('structure',None))
    if "charge_geometry" in mode:
        try:
            generate_xsf_file('out.'+prefix+'.scf')
        except:
            try:
                generate_xsf_file('out.'+prefix+'.vc-relax')
            except:
                try:
                    generate_xsf_file('out.'+prefix+'.relax')
                except Exception as e:
                    print(e)
        if kwargs.get('structure',None) is not None:
            try:
                get_geometry_data(kwargs['structure'])
            except Exception as e:
                print(e)
        directory_path = os.getcwd()
        print("My current directory is : " + directory_path)
        try:
            os.chdir("PP")
        except FileNotFoundError as error:
            print(error)
        baderPP(prefix)
        bader_average_and_bonds(prefix)
        os.chdir("../")
    if "potential" in mode:
        try:
            os.chdir("PP")
        except FileNotFoundError as error:
            print(error)
        try:
            os.rename("../avg.dat", "avg."+prefix+".dat")  ## avg.dat is written on folder of execution, have to mv it
            potential_average(filename="avg."+prefix+".dat")
        except Exception as error:
            print(error)
        os.chdir("../")
    if "gap" in mode:
        try:
            os.chdir("PP")
        except FileNotFoundError as error:
            print(error)
        getGap()
        os.chdir("../")
    if "pdos" in mode:
        try:
            os.chdir("PP/PDOS")
        except FileNotFoundError as error:
            print(error)
        sumPDOS(prefix)
        plot_PDOS(prefix,fermi=kwargs.get('fermi',0.0),
        RANGE_VB=kwargs.get('RANGE_VB',5),RANGE_CB=kwargs.get('RANGE_CB',5))
        os.chdir("../../")
    if "bands" in mode:
        try:
            os.chdir("PP")
        except FileNotFoundError as error:
            print(error)
        plot_bands(prefix,**kwargs)
        
#        indexes=kwargs.get('indexes'),fermi=kwargs.get('fermi',0.0),label=kwargs.get('label',''),
#        color=kwargs.get('color','black'),RANGE_VB=kwargs.get('RANGE_VB',1.5),RANGE_CB=kwargs.get('RANGE_CB',5),
#        redstates=kwargs.get('redstates',[]),greenstates=kwargs.get('greenstates',[]),bluestates=kwargs.get('bluestates',[]), labelsprojs=kwargs.get('labelsprojs',["","",""]),maxcontrast=kwargs.get("maxcontrast",True), 
#        contrast_ratios=kwargs.get("contrast_ratios",[1,1,1]), )
        #except Exception as error:
        #    print(error)
        os.chdir("../")


def PP_unfold(prefix,unfold,fermi=0.0,spin_polarized=False):
    if spin_polarized:
        os.chdir("OUTPUT")
        bands_dw=banduppy.BandStructure(code="espresso", prefix=prefix,spin_channel="dw")
        bands_up=banduppy.BandStructure(code="espresso", prefix=prefix,spin_channel="up")
        #pickle.dump(bands,open("bandstructure.pickle","wb"))
        result_up,result_down=unfold.unfold((bands_up,bands_dw),break_thresh=0.1,spin_polarized=True)
        #now plot the result as fat band
        unfold.plot(save_file="unfold_fatbandsp.png",plotSC=True,Emin=-5,Emax=5,Ef=3.988,fatfactor=5,mode='fatband',spin_polarized=True)
        #or as a colormap
        unfold.plot(save_file="unfold_densitysp.png",plotSC=True,Emin=-5,Emax=5,Ef=3.988,mode='density',smear=0.2,nE=200,spin_polarized=True)
        os.chdir("../")
    else:
        os.chdir("OUTPUT")
        bands=banduppy.BandStructure(code="espresso", prefix=prefix)
        #pickle.dump(bands,open("bandstructure.pickle","wb"))
        result=unfold.unfold(bands,break_thresh=0.1,spin_polarized=True)
        #now plot the result as fat band
        unfold.plot(save_file="unfold_fatband.png",plotSC=False,Emin=-5,Emax=5,Ef=fermi,fatfactor=5,mode='fatband')
        #or as a colormap
        unfold.plot(save_file="unfold_density.png",plotSC=False,Emin=-5,Emax=5,Ef=fermi,mode='density',smear=0.02,nE=200)
        os.chdir("../")


def baderPP(prefix):
        if os.path.isfile('bader_lnx_64.tar.gz'):
            print ("bader already downloaded.")
        else:
            subprocess.run(['wget', 'http://theory.cm.utexas.edu/henkelman/code/bader/download/bader_lnx_64.tar.gz'])
            subprocess.call(['tar', '-xzvf', 'bader_lnx_64.tar.gz'])
        subprocess.call(['./bader', '{prefix}.charge.cube'.format(prefix=prefix)])
        PREFIX=prefix
        os.rename("ACF.dat",prefix+"-ACF.dat")
        os.rename("AVF.dat",prefix+"-AVF.dat")
        os.rename("BCF.dat",prefix+"-BCF.dat")

	### funcao para fazer um grep em arquivos
        def grep(filename,regex):
            for file in glob.iglob(filename):
                for line in open(file, 'r'):
                    if re.search(regex, line):
                        return line
	
	
	### pegar nat e ntyp do input do QE
        subject=grep('../'+PREFIX+'.scf.in','nat')
        match = re.search("nat\s{0,}.\s{0,}(\d{1,})", subject)  ## . para ignorar o = entre paranteses o que quer dar match, \d{1,} significa uma ou mais ocorrencias de digito, igualmente \s{0,} vai considerar nenhum ou algum espaco
        if match:
            nat = match.group(1)
        else:
            print('Error: could not find number of atoms in file')
	#print(nat)
	
        subject=grep('../'+PREFIX+'.scf.in','ntyp')
        match = re.search("ntyp\s{0,}.\s{0,}(\d{1,})", subject)  ## . para ignorar o = entre paranteses o que quer dar match, \d{1,} significa uma ou mais ocorrencias de digito, igualmente \s{0,} vai considerar nenhum ou algum espaco
        if match:
            ntyp = match.group(1)
        else:
            print('Error: could not find ntyp in file')
	#print(ntyp)	
	
	### ler arquivo de carga do posprocessamento
        filecharge=PREFIX+'.charge'
        with open(filecharge) as fhand:
            lines = [line.rstrip("\n") for line in fhand]
            i=-1   ## i sera a linha da primeira ocorrencia da identificacao de carga
            for line in lines:
#      print(line)
                i=i+1
## iniciando com espacos depois digitos espacos letras espacos digitos, para pegar linhas de carga '    1   XX   20'
                if re.match(r"^\s+\d+\s+[a-zA-Z]+\s+\d+", line): 
                    break
        chargeatoms=lines[i:i+int(ntyp)]  ## pseudopotential charge for every atom
        chargeatomslist=[ line.split() for line in chargeatoms] ## dados para elementos
        chargedata=lines[i+int(ntyp):i+int(ntyp)+int(nat)] ## get elements and their coordinates in material
   #     print(chargeatoms, chargeatomslist, chargedata)
	### ler arquivo da analise de bader
        filebader=PREFIX+'-ACF.dat'
        with open(filebader) as fhand:
            lines = [line.rstrip("\n") for line in fhand]
        baderdata=lines[2:2+int(nat)]
	#transforma arquivos em nparrays
        bader = np.genfromtxt(baderdata)
        charge = np.genfromtxt(chargedata)
        chargeatoms = np.genfromtxt(chargeatoms)
        if chargeatoms.ndim == 1:
            chargeatoms=chargeatoms[2]
        else:
            chargeatoms = chargeatoms[:,2] ## PP atoms charge
	#for linha in bader verifica linha in chargedata, pega  da charge data, encontra o valor correspondente da  in chargeatoms , subtrai esse valor da  do bader
        atomlist=[] ### lista dos elementos em ordem
        for i in range(int(nat)):  ##  varre de 0 a nat
            if nat == 1:
                idatom=int(charge[4])-1
                chgatom=chargeatoms ## usa esse id para pegar a carga nominal do atomo
            else:
                idatom=int(charge[i,4])-1 ## pega o id do atomo no arquivo de carga ##pq posicao 0 e 1 aqui
                if chargeatoms.ndim == 0: ## ntyp=1
                    chgatom=chargeatoms ## usa esse id para pegar a carga nominal do atomo
                else:
                    chgatom=chargeatoms[idatom] ## usa esse id para pegar a carga nominal do atomo
            bader[i,4]=bader[i,4]-chgatom  ## pega a carga do atomo por bader e subtrai a carga nominal
            bader[0,0]=np.array(bader[0,0].astype(str))
            atomslist=atomlist.append(chargeatomslist[idatom][1])
	
        dt=np.dtype([('col0','U4'),('col1',float), ('col2',float),('col3',float), ('col4',float),('col5',float), ('col6',float)])
        baderfinal = np.empty((int(nat),), dtype=dt)
        baderfinal['col0'] = atomlist
        baderfinal['col1'] = bader[:,1]
        baderfinal['col2'] = bader[:,2]
        baderfinal['col3'] = bader[:,3]
        baderfinal['col4'] = bader[:,4]
        baderfinal['col5'] = bader[:,5]
        baderfinal['col6'] = bader[:,6]
	#print(baderfinal)
	
        np.savetxt('BaderChgAnalysis_'+PREFIX+'.out', baderfinal, delimiter=' ', header="Element, x, y, z, ChargeTransfered, MinDist, AtomVol", fmt="%3s %9f %9f %9f %9f %9f %9f")

def getGap():
	
	nspin=False
	RANGE_VB=3
	RANGE_CB=3
	CUTOFF_DOS=1 ## value above is considered 1
	fermi_idx_tol=500  ## difference in index from fermi index to consider gap close to fermi level
	
	
	datatxt=""
	for file in glob.glob("**/*.dos",recursive=True):
	    if len(file.split('\\')[-1].split('.')) > 2:
	        continue  ## case of out.*.dos files skipped
	    name=str(file)
	#    print("DOS FILE GIVEN:", name)
	    datatxt=datatxt+"DOS FILE GIVEN: "+name+"\n"
	    datatxt=datatxt+"----------------------------------------------------------- \n"
	    with open(file) as f:
	        line0=f.readlines()[0].split()[-2]  ## to get efermi on top of .dos file
	        fermi=float(line0)
	    data = np.genfromtxt(file)
	    #print(data[data[:,0]==fermi][0][1])
	    print("fermi",fermi)
	
	
	    ## test spin polarized
	    if len(data[0,:]) > 3:
	        nspin=True
	 #       print("SPIN POLARIZED CONSIDERED!")
	
	    ## only data in specified range
	    data_in_range=data[data[:,0]>(fermi-RANGE_VB)]
	    data_in_range=data_in_range[data_in_range[:,0]<(fermi+RANGE_CB)]
	    energy=data_in_range[:,0]
	    fermiidx=np.argwhere(energy==fermi)[0][0]
	    if nspin:
	        DOSup=np.where(data_in_range[:,1] > CUTOFF_DOS, 1, 0)
	        DOSdown=np.where(data_in_range[:,2] > CUTOFF_DOS, 1, 0)
	        DOSint=np.where(np.diff(data_in_range[:,3]) > 0.1, 1, 0)
	    else:
	        DOS=np.where(data_in_range[:,1] > CUTOFF_DOS, 1, 0)
	        DOSint=np.where(np.diff(data_in_range[:,2]) > 0.1, 1, 0)        
	        
	        ### detect all gaps in DOS, those closer to efermi are indicated.
	    if nspin:
	        counting=False
	        idx0=0
	        idx1=0
	        for idx, entry in enumerate(DOSup):
	            if entry == 0 and counting==False:
	                counting=True
	                idx0=idx
	                difftofermi=idx0-fermiidx 
	            if entry == 1 and counting==True:
	                counting=False
	                idx1=idx
	                gap=energy[idx1]-energy[idx0]
	                if gap > 0.1:
	#                    print("gap-SpinUP:",energy[idx1]-energy[idx0])
	                    datatxt=datatxt+"gap-SpinUP: "+str(energy[idx1]-energy[idx0])
	                    if abs(difftofermi) < fermi_idx_tol:
	 #                       print("Probably real gap, close to Efermi")
	                        datatxt=datatxt+"  Probably real gap, close to Efermi \n"
	                    else:
	                        datatxt=datatxt+"\n"
	
	                        
	        counting=False
	        idx0=0
	        idx1=0
	        for idx, entry in enumerate(DOSdown):
	            if entry == 0 and counting==False:
	                counting=True
	                idx0=idx
	                difftofermi=idx0-fermiidx 
	            if entry == 1 and counting==True:
	                counting=False
	                idx1=idx
	                gap=energy[idx1]-energy[idx0]
	                if gap > 0.1:
	  #                  print("gap-SpinDOWN:",energy[idx1]-energy[idx0])
	                    datatxt=datatxt+"gap-SpinDOWN: "+str(energy[idx1]-energy[idx0])
	                    if abs(difftofermi) < fermi_idx_tol:
	   #                     print("Probably real gap, close to Efermi")
	                        datatxt=datatxt+"  Probably real gap, close to Efermi \n"
	                    else:
	                        datatxt=datatxt+"\n"                 
	    
	        counting=False
	        idx0=0
	        idx1=0
	        for idx, entry in enumerate(DOSint):
	            if entry == 0 and counting==False:
	                counting=True
	                idx0=idx
	                difftofermi=idx0-fermiidx 
	            if entry == 1 and counting==True:
	                counting=False
	                idx1=idx
	                gap=energy[idx1]-energy[idx0]
	                if gap > 0.1:
	    #                print("gap-anySPIN:",energy[idx1]-energy[idx0])
	                    datatxt=datatxt+"gap-anySPIN: "+str(energy[idx1]-energy[idx0])
	                    if abs(difftofermi) < fermi_idx_tol:
	     #                   print("Probably real gap, close to Efermi")
	                        datatxt=datatxt+"  Probably real gap, close to Efermi \n"
	                    else:
	                        datatxt=datatxt+"\n"
	
	    else: ## no nspin
	        counting=False
	        idx0=0
	        idx1=0
	        for idx, entry in enumerate(DOS):
	            if entry == 0 and counting==False:
	                counting=True
	                idx0=idx
	                difftofermi=idx0-fermiidx 
	            if entry == 1 and counting==True:
	                counting=False
	                idx1=idx
	                gap=energy[idx1]-energy[idx0]
	                if gap > 0.1:
	      #              print("gap:",energy[idx1]-energy[idx0])
	                    datatxt=datatxt+"gap: "+str(energy[idx1]-energy[idx0])
	                    if abs(difftofermi) < fermi_idx_tol:
	       #                 print("Probably real gap, close to Efermi")
	                        datatxt=datatxt+"  Probably real gap, close to Efermi \n"
	                    else:
	                        datatxt=datatxt+"\n"
	    # plt.plot(data_in_range[:,0], DOSdown)
	    datatxt=datatxt+"----------------------------------------------------------- \n"
	with open("gaps.txt", "w") as text_file:
	        text_file.write(datatxt)

def bader_average_and_bonds(prefix):

	filename="BaderChgAnalysis_{prefix}.out".format(prefix=prefix)
	fileout="{prefix}_baderbonds.out".format(prefix=prefix)
	baderdata = np.loadtxt(filename,dtype={'names': ['element', 'x', 'y','z','charge'],
	                     'formats': ['U10', 'f4', 'f4', 'f4','f8']},skiprows=1)
	written_data=""
	written_data=written_data+str(baderdata)+"\n"
	
	baderdata=[list(baderdata[i]) for i in range(len(baderdata))] ## transform list makes easier to edit
	
	## get list of elements
	elements=[]
	for row in baderdata:
	    if row[0] not in elements:
	        elements.append(row[0]) # if row[0] not in elements
	
	
	## calculate average charge of each element
	avgcharge=[]
	for element in elements:
	    count=0
	    charge=0
	    for row in baderdata:
	        if element == row[0]:
	            count+=1 ##increment number of the element
	            charge+=row[4]
	    avgcharge.append(charge/float(count))
	#print('------------------------------------------------------- \n')
	#print('------------------------------------------------------- \n')
	written_data=written_data+'------------------------------------------------------- \n'+'------------------------------------------------------- \n'
	
	#print('Elements: ',elements)
	#print('Average charge of elements: ',avgcharge)   
	written_data=written_data+'Elements: '+str(elements)+"\n"+'Average charge of elements: '+ \
	             str(avgcharge)+"\n"
	
	#print('------------------------------------------------------- \n')
	#print('------------------------------------------------------- \n')
	written_data=written_data+'------------------------------------------------------- \n'+'------------------------------------------------------- \n'
	
	
	## labelling by element
	labelelements=np.zeros(len(elements))
	for row in baderdata:
	    indexelement=elements.index(row[0])
	    labelelements[indexelement]+=1 ## increase label of the element
	    row.append(elements[indexelement]+str(int(labelelements[indexelement])))
	
	
	## gets distance between different atoms within given cutoff
	cutoff_bond=5  ## in angstrom
	#print('-------------CUTOFF FOR BOND LENGTH :  '+str(cutoff_bond)+'  angstrom ----------------------- \n')
	written_data=written_data+'-------------CUTOFF FOR BOND LENGTH :  '+str(cutoff_bond)+'  angstrom ----------------------- \n'
	
	distbond_table=[]
	numberofbonds=0
	for row in baderdata:
	    for row2 in baderdata:
	       if row2 != row:
	             distxyz=np.array([(row2[i]-row[i])/1.88973 for i in range(1,4)])
	             distbond = np.linalg.norm(distxyz)  ## get modulus, distance between atoms
	             # print(distxyz,row[0],row2[0],distbond )
	             
	             if distbond < cutoff_bond:
	                 distbond_table.append([row[0],row2[0],distbond,row[4],row2[4]])
	                 print(str(row[5])+"-"+str(row2[5])+": "+str(distbond)+" ang" +" Charges: "+
	                       str(row[4])+")-("+str(row2[4]))
	                 written_data=written_data+str(row[5])+"-"+str(row2[5])+": "+str(distbond)+" ang" +" Charges: "+ \
	                       str(row[4])+")-("+str(row2[4])+"\n"
	                 numberofbonds+=1
	    #         labelelements2[indexelement]+=1 ## increase label of the element
	# print("DISTBONDTABLE",distbond_table)
	possible_element_bonds=[list(comb) for comb in itertools.combinations(elements,2)]
	
	avg_bond_length=np.zeros(len(possible_element_bonds))
	bond_lengths = [[] for _ in range(len(possible_element_bonds))] ## list of lists with each bond
	counter_bond=np.zeros(len(possible_element_bonds))
	for bond in distbond_table:
	    for index,typebond in enumerate(possible_element_bonds):
	        if (typebond[0] == bond[0] and typebond[1] == bond[1]) or (typebond[1] == bond[0] and typebond[0] == bond[1]):
	            avg_bond_length[index]+=bond[2]
	           # if bond[2] not in bond_lengths[index]: ## no duplicates
	            bond_lengths[index].append(bond[2])
	            counter_bond[index]+=1
	std_devs=[ np.std(bond_lengths[i]) for i in range(len(bond_lengths)) ]
	variances=[ np.var(bond_lengths[i]) for i in range(len(bond_lengths)) ]
	avs=[ np.average(bond_lengths[i]) for i in range(len(bond_lengths)) ]
	
	# print(std_devs)
	# print(variances)
	# print(avs)
	## how to ignore zeros in division
	avg_bond_length = np.divide(avg_bond_length, counter_bond, out=np.zeros_like(avg_bond_length), where=counter_bond!=0)
	# print(possible_element_bonds)
	# print(avg_bond_length)
	#print('------------------------------------------------------- \n')
	#print('------------------------------------------------------- \n')
	written_data=written_data+'------------------------------------------------------- \n'+'------------------------------------------------------- \n'
	
	#print("Average bond lengths, if 0 the atomic distance is larger than the cutoff to be considered a bond:")
	#[ print(possible_element_bonds[i][0]+"-"+possible_element_bonds[i][1]+":\n avg:"+str(avg_bond_length[i])+", std:"+ \
	#    str(std_devs[i])+", variance: "+str(variances[i])+" \n")  for i in range(len(possible_element_bonds)) ]
	written_data=written_data+"Average bond lengths, if 0 the atomic distance is larger than the cutoff to be considered a bond: \n"
	for i in range(len(possible_element_bonds)):
	    written_data=written_data+possible_element_bonds[i][0]+"-"+possible_element_bonds[i][1]+": \n avg:"+str(avg_bond_length[i])+", std:"+ \
	    str(std_devs[i])+", variance: "+str(variances[i])+" \n"
	# print(written_data)
	written_data=written_data+'------------------------------------------------------- \n'+'------------------------------------------------------- \n'
	#print('------------------------------------------------------- \n')
	#print('------------------------------------------------------- \n')
	
	with open(fileout,"w") as f:
	    f.write(written_data) 


def potential_average(filename="avg.dat",calculate_macavg=False,plot_potentials=False):
    import numpy as np
    if len(filename.split("."))>2:
        prefix=filename.split(".")[1]
    else:
        prefix=""
    def get_period_for_mavg(filename):  ## file avg.dat from quantum espresso
        data=np.genfromtxt(filename)
        x=data[:,0]
        pot=data[:,1]

        # print(data)
        ft = np.fft.rfft(pot)
        freqs = np.fft.rfftfreq(len(pot), x[1]-x[0]) # Get frequency axis from the time axis
        mags=abs(ft)
        periods=[1/x for x in freqs if x != 0 ]
        idxmax=np.where(mags[1:]==max(mags[1:]))[0][0]
        print("Higher frequency component is found at period: ",periods[idxmax])
        # print(periods)
        # print(mags)
        if plot_potentials:
            fig, ax = plt.subplots(1,1)
            plt.plot(periods,mags[1:])
            plt.show()
        # a=np.stack((periods,mags[1:]),axis=-1)
        # print(a)
        return periods[idxmax]
    
    
    def getMacroAverage(distance, potential, value_of_period):  ## distance and potential mustve same dimensions
    
        R=len(distance)
        index_of_period = (np.abs(distance - value_of_period)).argmin() ## gets index of period
        index_of_period2 = (np.abs(distance - value_of_period/2)).argmin() ## gets index of half period
        
        d_extend=np.concatenate((distance,distance[1:index_of_period+1]+max(distance)))  ## extend distance to (0,len(distance)+oneperiod)
        V_extend=np.concatenate((potential[-index_of_period2:],potential,potential[:index_of_period2+1]))
        ## extend potential to ((half period in the end),normal distance,(halfperiod from beginning)), results in -T/2 -- T/2 integral
        # V_extend2=np.concatenate((potential,potential[1:index_of_period+1]))  ## extend potential ## equivalent to 0 -- T integral
        result=np.array([integrate.simps(V_extend[x:index_of_period+x], d_extend[x:index_of_period+x])/value_of_period for x in range(0,R) ])
        return result
        
        
    if calculate_macavg:
        period_mavg=get_period_for_mavg(filename)
        return period_mavg
    
    data=np.genfromtxt(filename)
    distance=data[:,0]
    potential=data[:,1]
    mac_avg_qe=data[:,2]
    if plot_potentials:
        plt.plot(data[:,0]*0.52918,data[:,1]*13.6056622,color="blue") ## potential Angstrom X eV
        plt.plot(distance*0.52918,mac_avg_qe*13.6056622,color="green")  ## macavg from QE ## should be identical to the calculated one
        plt.axhline(y=max(mac_avg_qe*13.6056622),linestyle=':',color='blue',alpha=0.5)
        plt.axhline(y=min(mac_avg_qe*13.6056622),linestyle=':',color='blue',alpha=0.5)
        if calculate_macavg:	
            mac_avg_calc=getMacroAverage(distance,potential,period_mavg)
            plt.plot(distance*0.52918,mac_avg_calc*13.6056622,color="purple") ## potential Angstrom X eV
            plt.axhline(y=max(mac_avg_calc*13.6056622),linestyle=':',color='blue',alpha=0.5)
            plt.axhline(y=min(mac_avg_calc*13.6056622),linestyle=':',color='blue',alpha=0.5)
        plt.show()
    
    writethis="Potential difference is: "+str(max(mac_avg_qe*13.6056622)-min(mac_avg_qe*13.6056622))+" eV \n"
    writethis+="Potential average is: "+str((max(mac_avg_qe*13.6056622)+min(mac_avg_qe*13.6056622))/2)+" eV"
    with open("averagepot."+prefix+".out","w") as f:
        f.write(writethis)
    print(writethis)
	


def sumPDOS(prefix):
	import sys
	import os
	import fnmatch, re
	import linecache
	import matplotlib
	import matplotlib.pyplot as plt
	import glob
	import numpy as np
	from cycler import cycler
	
	
	### funcao para fazer um grep em arquivos
	def grep(filename,regex):
	      for file in glob.iglob(filename):
	         for line in open(file, 'r'):
	            if re.search(regex, line):
	               return line
	
	PREFIX=prefix
	
	# Some default variables
	fermi=0
	graphtitle=""
	min_x,max_x=-2,7
	min_y,max_y="",""
	wfcnumber='X'
	wfc='X'
	showgraph=0
	nspin=1
	
	
	
	### pegar nspin do input do QE
	subject=grep('../../'+PREFIX+'.scf.in','nspin')
	#print('subject:', subject)
	#print(subject.strip().startswith('!'))
	if subject != None and not subject.strip().startswith('!') :
	    match = re.search("nspin\s{0,}.\s{0,}(\d{1,})", subject)  ## . para ignorar o = entre paranteses o que quer dar match, \d{1,} significa uma ou mais ocorrencias de digito, igualmente \s{0,} vai considerar nenhum ou algum espaco
	    if match:
	        nspin = int(match.group(1))
	    else:
	        print('Error: could not find nspin in file. Using nspin=1.')
	else:
	    print('Error: could not find nspin in file. Using nspin=1.')
	
	
	
	datapdos= [ ]
	dataelem= [ ]
	for file in glob.glob(PREFIX+'.pdos_atm*'):
	    print(file)
	    re_atmnumber = re.findall(PREFIX+".pdos_atm#([0-9]{1,}).{1,}", file)
	    print(re_atmnumber)
	    re_atm = re.findall(PREFIX+".pdos_atm#[0-9]{1,}\(([a-zA-Z]{1,})\).{1,}", file)
	    print(re_atm)
	    re_wfcnumber = re.findall(PREFIX+".pdos_atm#[0-9]{1,}\([a-zA-Z]{1,}\)_wfc#([0-9])\([a-z]\)", file)
	    print(re_wfcnumber)
	    re_wfc = re.findall(PREFIX+".pdos_atm#[0-9]{1,}\([a-zA-Z]{1,}\)_wfc#[0-9]\(([a-z])\)", file)
	    print(re_wfc)
	    datapdos.append([ re_atm[0], re_wfcnumber[0], re_wfc[0]])
	    dataelem.append([ re_atm[0] ])
	datapdos=np.array(datapdos)  ### conjunto de dados de projecao
	datapdos=np.unique(datapdos,axis=0)   ### dados de projecao sem repeticao
	print(datapdos)
	dataelem=np.array(dataelem)  ### conjunto de dados de projecao
	dataelem=np.unique(dataelem,axis=0)   ### dados de projecao sem repeticao
	print(dataelem)
	
	### usar as entradas do datapdos para somar os arquivos de forma apropriada
	
	
	########################## FUNCAO PARA SOMAR OS PDOS #############################
	def sum_PDOS(PREFIX,element, wfcnumber, wfc, fermi, graphtitle, min_x, max_x, min_y, max_y, showgraph, nspin ):
	  dosfiles=[]
	  if element == 'XX':
	    element='.{1,}'
	  if wfc == 'X':
	    wfc='.'
	  if wfcnumber == 'X':
	    wfcnumber='[0-9]'
	
	  selat='^'+PREFIX+'.pdos_atm#[0-9]{1,}\('+element+'\)_wfc#'+wfcnumber+'\('+wfc+'\)'
	  if wfc == '.':
	    wfc=''
	  if wfcnumber == '[0-9]':
	    wfcnumber=''
	
	  for dfile in glob.glob(PREFIX+'*'):
	    print(dfile)
	    if re.match(selat, dfile):
	      dosfiles.append(dfile) 
	  print('DOS files matching regexp:')
	  print(dosfiles)
	  if len(dosfiles)==0:
	    print("ERROR: Provide a (list of) valid DOS file(s)")
	    sys.exit()
	
	  ## se fermi nao especificado procurar no arquivo dos na pasta PP
	  if fermi == 0 :
	    with open('../'+PREFIX+'.dos') as f:
	      first_line = f.readline()
	      result = re.findall(r"Fermi *[^\w ] *(.*) eV", first_line, re.IGNORECASE | re.MULTILINE)
	      result = result[0].split()
	      if len(result) > 1 :
	        if result[0] == result[1]:
	          fermi=round(float(result[0]),3)
	        else:
	          print("ERROR: spin up and down fermi are different")
	      else:
	        fermi=round(float(result[0]),3)
	
	  #  with open('../../out.'+PREFIX+'.nscf') as f:
	  #    for line in f :
	  #        if re.search("highest", line):
	  #          re_fermi = line.split()
	  #          fermi=round(float(re_fermi[6]),3)
	  #        if re.search("Fermi energy", line):
	  #          re_fermi = line.split()
	  #          fermi=round(float(re_fermi[4]),3)
	  print('fermi: ',fermi)
	    
	  
	  mat=[]  # matrix with total sum of ldos
	
	  if nspin == 1:
	      for i in range(len(dosfiles)):
	       mati=[] # temporal matrix for each DOS file "i"
	       k=0
	       for line in open(dosfiles[i],'r'):
	        if len(line) > 10 and line.split()[0] != "#":
	            if wfc == 's' or wfc == '':
	              mati.append([float(line.split()[0]),float(line.split()[1])])
	            if wfc == 'p' :
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])
	            if wfc == 'd' :
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5]),float(line.split()[6])])
	            if wfc == 'f' :
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),float(line.split()[8])])
	
	       if mat == []: # if it is the first dos file, copy total matrix (mat) = the first dos files's data
	          mat=mati[:]
	       else:
	          for j in range(len(mati)): # if it is not the first file, sum values
	              if wfc == 's' or wfc == '':
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1]]
	              if wfc == 'p':  
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2],mat[j][3]+mati[j][3],mat[j][4]+mati[j][4] ]  
	              if wfc == 'd':  
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2],mat[j][3]+mati[j][3],mat[j][4]+mati[j][4],mat[j][5]+mati[j][5],mat[j][6]+mati[j][6] ]  
	              if wfc == 'f':  
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2],mat[j][3]+mati[j][3],mat[j][4]+mati[j][4],mat[j][5]+mati[j][5],mat[j][6]+mati[j][6],mat[j][7]+mati[j][7],mat[j][8]+mati[j][8] ]  
	
	
	      print("...ploting...")
	      if wfc == 's' or wfc == '':
	        x,y=[],[]    
	      if wfc == 'p' :
	        x,y,y1,y2,y3=[],[],[],[],[]    
	      if wfc == 'd':
	        x,y,y1,y2,y3,y4,y5=[],[],[],[],[],[],[]
	      if wfc == 'f':
	        x,y,y1,y2,y3,y4,y5,y6,y7=[],[],[],[],[],[],[],[],[]
	
	      for i in mat:
	        if wfc == 's' or wfc == '':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	        if wfc == 'p':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          y1.append(i[2])
	          y2.append(i[3])
	          y3.append(i[4])
	        if wfc == 'd':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          y1.append(i[2])
	          y2.append(i[3])
	          y3.append(i[4])
	          y4.append(i[5])
	          y5.append(i[6])
	        if wfc == 'f':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          y1.append(i[2])
	          y2.append(i[3])
	          y3.append(i[4])
	          y4.append(i[5])
	          y5.append(i[6])
	          y6.append(i[7])
	          y7.append(i[8])
	      if wfc == 's' or wfc == '':
	        x=np.array(x)
	        y=np.array(y)
	        total=np.stack((x, y), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f")
	        plt.plot(x,y,linewidth=1.0)
	        plt.fill(x,y,color='0.8')
	
	      if wfc == 'p' :
	        x=np.array(x)
	        y=np.array(y)
	        y1=np.array(y1)
	        y2=np.array(y2)
	        y3=np.array(y3)
	        total=np.stack((x, y, y1,y2,y3), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0)
	        plt.plot(x,y1,linewidth=1.0, label='$\mathregular{pz}$')
	        plt.plot(x,y2,linewidth=1.0, label='$\mathregular{px}$')
	        plt.plot(x,y3,linewidth=1.0, label='$\mathregular{py}$')
	        plt.fill(x,y,color='0.8')
	        plt.fill(x,y1,color='0.9')
	        plt.fill(x,y2,color='0.9')
	        plt.fill(x,y3,color='0.9')
	
	      if wfc == 'd' :
	        x=np.array(x)
	        y=np.array(y)
	        y1=np.array(y1)
	        y2=np.array(y2)
	        y3=np.array(y3)
	        y4=np.array(y4)
	        y5=np.array(y5)
	        total=np.stack((x, y, y1,y2,y3,y4,y5), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f %.6f %.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0)
	        plt.plot(x,y1,linewidth=1.0, label='$\mathregular{dz^2}$')
	        plt.plot(x,y2,linewidth=1.0, label='$\mathregular{dzx}$')
	        plt.plot(x,y3,linewidth=1.0, label='$\mathregular{dzy}$')
	        plt.plot(x,y4,linewidth=1.0, label='$\mathregular{dx^2-y^2}$')
	        plt.plot(x,y5,linewidth=1.0, label='$\mathregular{dxy}$')
	        plt.fill(x,y,color='0.8')
	        plt.fill(x,y1,color='0.9')
	        plt.fill(x,y2,color='0.9')
	        plt.fill(x,y3,color='0.9')
	        plt.fill(x,y4,color='0.9')
	        plt.fill(x,y5,color='0.9')
	
	      if wfc == 'f' :
	        x=np.array(x)
	        y=np.array(y)
	        y1=np.array(y1)
	        y2=np.array(y2)
	        y3=np.array(y3)
	        y4=np.array(y4)
	        y5=np.array(y5)
	        y6=np.array(y6)
	        y7=np.array(y7)
	        total=np.stack((x, y, y1,y2,y3,y4,y5,y6,y7), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0, label='f' )
	        plt.plot(x,y1,linewidth=1.0, label='f1')
	        plt.plot(x,y2,linewidth=1.0, label='f2')
	        plt.plot(x,y3,linewidth=1.0, label='f3')
	        plt.plot(x,y4,linewidth=1.0, label='f4')
	        plt.plot(x,y5,linewidth=1.0, label='f5')
	        plt.plot(x,y6,linewidth=1.0, label='f6')
	        plt.plot(x,y7,linewidth=1.0, label='f7')
	        plt.fill(x,y,color='0.8')
	        plt.fill(x,y1,color='0.9')
	        plt.fill(x,y2,color='0.9')
	        plt.fill(x,y3,color='0.9')
	        plt.fill(x,y4,color='0.9')
	        plt.fill(x,y5,color='0.9')
	        plt.fill(x,y6,color='0.9')
	        plt.fill(x,y7,color='0.9')
	      # if there is matplotlib, generate a plot with it
	      plt.title(graphtitle)
	      plt.xlabel('E (eV)')
	      plt.ylabel('States')
	      plt.legend(loc='upper right')
	    #  plt.grid(True)
	      # plt.rcParams.update({'font.size': 22})
	      if min_x and max_x:
	       fromx,tox=min_x,max_x 
	      margin=max(y)*0.3
	      plt.axis([fromx, tox, -margin, max(y)+margin])
	
	      fig = plt.gcf()
	      if showgraph == 1:
	        plt.show()   
	      fig.savefig(PREFIX+'_PDOS-'+element+wfcnumber+wfc+'.svg', format='svg', dpi=1000)
	      fig.clf()
	      plt.close()
	
	  if nspin == 2:
	      for i in range(len(dosfiles)):
	       mati=[] # temporal matrix for each DOS file "i"
	       k=0
	       for line in open(dosfiles[i],'r'):
	        if len(line) > 10 and line.split()[0] != "#":
	            if wfc == 's' or wfc == '':
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])
	            if wfc == 'p' :
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),float(line.split()[8])])
	            if wfc == 'd' :
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),float(line.split()[8]),float(line.split()[9]),float(line.split()[10]),float(line.split()[11]),float(line.split()[12])])
	            if wfc == 'f' :
	              mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4]),float(line.split()[5]),float(line.split()[6]),float(line.split()[7]),float(line.split()[8]),float(line.split()[9]),float(line.split()[10]),float(line.split()[11]),float(line.split()[12]),float(line.split()[13]),float(line.split()[14]),float(line.split()[15]),float(line.split()[16])])
	       if mat == []: # if it is the first dos file, copy total matrix (mat) = the first dos files's data
	          mat=mati[:]
	       else:
	          for j in range(len(mati)): # if it is not the first file, sum values
	              if wfc == 's' or wfc == '':
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2]]
	              if wfc == 'p':  
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2],mat[j][3]+mati[j][3],mat[j][4]+mati[j][4],mat[j][5]+mati[j][5],mat[j][6]+mati[j][6],mat[j][7]+mati[j][7],mat[j][8]+mati[j][8] ] 
	              if wfc == 'd':  
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2],mat[j][3]+mati[j][3],mat[j][4]+mati[j][4],mat[j][5]+mati[j][5],mat[j][6]+mati[j][6],mat[j][7]+mati[j][7],mat[j][8]+mati[j][8],mat[j][9]+mati[j][9],mat[j][10]+mati[j][10],mat[j][11]+mati[j][11],mat[j][12]+mati[j][12] ]  
	              if wfc == 'f':  
	                mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2],mat[j][3]+mati[j][3],mat[j][4]+mati[j][4],mat[j][5]+mati[j][5],mat[j][6]+mati[j][6],mat[j][7]+mati[j][7],mat[j][8]+mati[j][8],mat[j][9]+mati[j][9],mat[j][10]+mati[j][10],mat[j][11]+mati[j][11],mat[j][12]+mati[j][12],mat[j][13]+mati[j][13],mat[j][14]+mati[j][14],mat[j][15]+mati[j][15],mat[j][16]+mati[j][16] ]  
	
	
	      print("...ploting...")
	      if wfc == 's' or wfc == '':
	        x,y,yu=[],[],[]    
	      if wfc == 'p' :
	        x,y,yu,y1,y1u,y2,y2u,y3,y3u=[],[],[],[],[],[],[],[],[]    
	      if wfc == 'd':
	        x,y,yu,y1,y1u,y2,y2u,y3,y3u,y4,y4u,y5,y5u=[],[],[],[],[],[],[],[],[],[],[],[],[]
	      if wfc == 'f':
	        x,y,yu,y1,y1u,y2,y2u,y3,y3u,y4,y4u,y5,y5u,y6,y6u,y7,y7u=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	
	      for i in mat:
	        if wfc == 's' or wfc == '':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          yu.append(i[2])
	        if wfc == 'p':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          yu.append(i[2])
	          y1.append(i[3])
	          y1u.append(i[4])
	          y2.append(i[5])
	          y2u.append(i[6])
	          y3.append(i[7])
	          y3u.append(i[8])
	        if wfc == 'd':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          yu.append(i[2])
	          y1.append(i[3])
	          y1u.append(i[4])
	          y2.append(i[5])
	          y2u.append(i[6])
	          y3.append(i[7])
	          y3u.append(i[8])
	          y4.append(i[9])
	          y4u.append(i[10])
	          y5.append(i[11])
	          y5u.append(i[12])
	        if wfc == 'f':
	          x.append(i[0]-fermi)
	          y.append(i[1])
	          yu.append(i[2])
	          y1.append(i[3])
	          y1u.append(i[4])
	          y2.append(i[5])
	          y2u.append(i[6])
	          y3.append(i[7])
	          y3u.append(i[8])
	          y4.append(i[9])
	          y4u.append(i[10])
	          y5.append(i[11])
	          y5u.append(i[12])
	          y6.append(i[13])
	          y6u.append(i[14])
	          y7.append(i[15])
	          y7u.append(i[16])
	      if wfc == 's' or wfc == '':
	        x=np.array(x)
	        y=np.array(y)
	        yu=np.array(yu)
	        total=np.stack((x, y, yu), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0)
	        plt.fill(x,y,color='0.8')
	        plt.plot(x,-yu,linewidth=1.0)
	        plt.fill(x,-yu,color='0.7')
	
	      if wfc == 'p' :
	        x=np.array(x)
	        y=np.array(y)
	        yu=np.array(yu)
	        y1=np.array(y1)
	        y1u=np.array(y1u)
	        y2=np.array(y2)
	        y2u=np.array(y2u)
	        y3=np.array(y3)
	        y3u=np.array(y3u)
	        total=np.stack((x, y, yu, y1, y1u, y2, y2u, y3, y3u), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0)
	        plt.plot(x,y1,linewidth=1.0, label='$\mathregular{pz}$')
	        plt.plot(x,y2,linewidth=1.0, label='$\mathregular{px}$')
	        plt.plot(x,y3,linewidth=1.0, label='$\mathregular{py}$')
	        plt.fill(x,y,color='0.8')
	        plt.fill(x,y1,color='0.9')
	        plt.fill(x,y2,color='0.9')
	        plt.fill(x,y3,color='0.9')
	        plt.plot(x,-yu,linewidth=1.0)
	        plt.plot(x,-y1u,linewidth=1.0, label='$\mathregular{pz}$')
	        plt.plot(x,-y2u,linewidth=1.0, label='$\mathregular{px}$')
	        plt.plot(x,-y3u,linewidth=1.0, label='$\mathregular{py}$')
	        plt.fill(x,-yu,color='0.7')
	        plt.fill(x,-y1u,color='0.8')
	        plt.fill(x,-y2u,color='0.8')
	        plt.fill(x,-y3u,color='0.8')
	
	      if wfc == 'd' :
	        x=np.array(x)
	        y=np.array(y)
	        yu=np.array(yu)
	        y1=np.array(y1)
	        y1u=np.array(y1u)
	        y2=np.array(y2)
	        y2u=np.array(y2u)
	        y3=np.array(y3)
	        y3u=np.array(y3u)
	        y4=np.array(y4)
	        y4u=np.array(y4u)
	        y5=np.array(y5)
	        y5u=np.array(y5u)
	        total=np.stack((x, y, yu, y1, y1u, y2, y2u, y3, y3u, y4, y4u, y5, y5u), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0)
	        plt.plot(x,y1,linewidth=1.0, label='$\mathregular{dz^2}$')
	        plt.plot(x,y2,linewidth=1.0, label='$\mathregular{dzx}$')
	        plt.plot(x,y3,linewidth=1.0, label='$\mathregular{dzy}$')
	        plt.plot(x,y4,linewidth=1.0, label='$\mathregular{dx^2-y^2}$')
	        plt.plot(x,y5,linewidth=1.0, label='$\mathregular{dxy}$')
	        plt.fill(x,y,color='0.8')
	        plt.fill(x,y1,color='0.9')
	        plt.fill(x,y2,color='0.9')
	        plt.fill(x,y3,color='0.9')
	        plt.fill(x,y4,color='0.9')
	        plt.fill(x,y5,color='0.9')
	        plt.plot(x,-yu,linewidth=1.0)
	        plt.plot(x,-y1u,linewidth=1.0, label='$\mathregular{dz^2}$')
	        plt.plot(x,-y2u,linewidth=1.0, label='$\mathregular{dzx}$')
	        plt.plot(x,-y3u,linewidth=1.0, label='$\mathregular{dzy}$')
	        plt.plot(x,-y4u,linewidth=1.0, label='$\mathregular{dx^2-y^2}$')
	        plt.plot(x,-y5u,linewidth=1.0, label='$\mathregular{dxy}$')
	        plt.fill(x,-yu,color='0.7')
	        plt.fill(x,-y1u,color='0.8')
	        plt.fill(x,-y2u,color='0.8')
	        plt.fill(x,-y3u,color='0.8')
	        plt.fill(x,-y4u,color='0.8')
	        plt.fill(x,-y5u,color='0.8')
	
	      if wfc == 'f' :
	        x=np.array(x)
	        y=np.array(y)
	        yu=np.array(yu)
	        y1=np.array(y1)
	        y1u=np.array(y1u)
	        y2=np.array(y2)
	        y2u=np.array(y2u)
	        y3=np.array(y3)
	        y3u=np.array(y3u)
	        y4=np.array(y4)
	        y4u=np.array(y4u)
	        y5=np.array(y5)
	        y5u=np.array(y5u)
	        y6=np.array(y6)
	        y6u=np.array(y6u)
	        y7=np.array(y7)
	        y7u=np.array(y7u)
	        total=np.stack((x, y, yu, y1, y1u, y2, y2u, y3, y3u, y4, y4u, y5, y5u, y6, y6u, y7, y7u), axis=-1)
	        np.savetxt('sumpdos_'+PREFIX+'.'+element+wfcnumber+wfc, total, delimiter='', fmt="%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f")
	        plt.plot(x,y,linewidth=1.0, label='f' )
	        plt.plot(x,y1,linewidth=1.0, label='f1')
	        plt.plot(x,y2,linewidth=1.0, label='f2')
	        plt.plot(x,y3,linewidth=1.0, label='f3')
	        plt.plot(x,y4,linewidth=1.0, label='f4')
	        plt.plot(x,y5,linewidth=1.0, label='f5')
	        plt.plot(x,y6,linewidth=1.0, label='f6')
	        plt.plot(x,y7,linewidth=1.0, label='f7')
	        plt.fill(x,y,color='0.8')
	        plt.fill(x,y1,color='0.9')
	        plt.fill(x,y2,color='0.9')
	        plt.fill(x,y3,color='0.9')
	        plt.fill(x,y4,color='0.9')
	        plt.fill(x,y5,color='0.9')
	        plt.fill(x,y6,color='0.9')
	        plt.fill(x,y7,color='0.9')
	        plt.plot(x,-yu,linewidth=1.0, label='f' )
	        plt.plot(x,-y1u,linewidth=1.0, label='f1')
	        plt.plot(x,-y2u,linewidth=1.0, label='f2')
	        plt.plot(x,-y3u,linewidth=1.0, label='f3')
	        plt.plot(x,-y4u,linewidth=1.0, label='f4')
	        plt.plot(x,-y5u,linewidth=1.0, label='f5')
	        plt.plot(x,-y6u,linewidth=1.0, label='f6')
	        plt.plot(x,-y7u,linewidth=1.0, label='f7')
	        plt.fill(x,-yu,color='0.7')
	        plt.fill(x,-y1u,color='0.8')
	        plt.fill(x,-y2u,color='0.8')
	        plt.fill(x,-y3u,color='0.8')
	        plt.fill(x,-y4u,color='0.8')
	        plt.fill(x,-y5u,color='0.8')
	        plt.fill(x,-y6u,color='0.8')
	        plt.fill(x,-y7u,color='0.8')
	      # if there is matplotlib, generate a plot with it
	      plt.title(graphtitle)
	      plt.xlabel('E (eV)')
	      plt.ylabel('States')
	      plt.legend(loc='upper right')
	    #  plt.grid(True)
	      # plt.rcParams.update({'font.size': 22})
	      if min_x and max_x:
	       fromx,tox=min_x,max_x 
	      margin=max(y)*0.3
	      plt.axis([fromx, tox, min(-y)-margin, max(y)+margin])
	
	      fig = plt.gcf()
	      if showgraph == 1:
	        plt.show()   
	      fig.savefig(PREFIX+'_PDOS-'+element+wfcnumber+wfc+'.svg', format='svg', dpi=1000)
	      fig.clf()
	      plt.close()
	
	####################################################################################
	####################################################################################
	
	########### GERAR ARQUIVOS DE PROJECAO ORBITAIS E ATOMICOS #############
	
	for projection in datapdos:
	  sum_PDOS(PREFIX, projection[0], projection[1], projection[2], fermi, graphtitle, min_x, max_x, min_y, max_y, showgraph, nspin)
	print(dataelem)
	for projection in dataelem:
	  sum_PDOS(PREFIX, projection[0], 'X', 'X', fermi, graphtitle, min_x, max_x, min_y, max_y, showgraph, nspin)  
	####################################################################################

def get_Efermi(prefix):
    fermi=0.0
    try:
        with open('PP/'+prefix+'.dos') as f:
           first_line = f.readline()
           result = re.findall(r"Fermi *[^\w ] *(.*) eV", first_line, re.MULTILINE)
           result = result[0].split()
           if len(result) > 1 :
             if result[0] == result[1]:
               fermi=round(float(result[0]),3)
             else:
               print("ERROR: spin up and down fermi are different")
           else:
             fermi=round(float(result[0]),3)
        return fermi
    except FileNotFoundError:
        try:    
            with open('out.'+prefix+'.scf') as f:
                lines=f.readlines()
                for line in lines:
                    result=re.match(
                    ".*the Fermi energy is\s*(\d+\.\d+)\s*",line)
                    if result:
                        fermi=float(result.group(1))
                        return fermi
        except:
            print("Could not find fermi energy.")

def plot_PDOS(prefix,RANGE_VB=2,RANGE_CB=7,fermi=0):
	import sys
	import os
	import fnmatch, re
	import linecache
	import matplotlib
	import matplotlib.pyplot as plt
	import glob
	import numpy as np
	from cycler import cycler
	
	### se alinha do nspin comeca com ! deveria ser ignorado
	### implementar para qualquer variavel
	
	### funcao para fazer um grep em arquivos
	def grep(filename,regex):
	      for file in glob.iglob(filename):
	         for line in open(file, 'r'):
	            if re.search(regex, line):
	               return line
	
	PREFIX=prefix
	
	# Some default variables
	nspin=1
	graphtitle=""
	min_x,max_x=-RANGE_VB,RANGE_CB
	min_y,max_y="",""
	wfcnumber='X'
	wfc='X'
	showgraph=0
	
	######### PLOTTING ONLY  ##################
	### plotar DOStotal com proj de cada orbital
	  ## se fermi nao especificado procurar no arquivo dos na pasta PP
	if fermi == 0 :
	    with open('../'+PREFIX+'.dos') as f:
	      first_line = f.readline()
	      result = re.findall(r"Fermi *[^\w ] *(.*) eV", first_line, re.IGNORECASE | re.MULTILINE)
	      result = result[0].split()
	      if len(result) > 1 :
	        if result[0] == result[1]:
	          fermi=round(float(result[0]),3)
	        else:
	          print("ERROR: spin up and down fermi are different")
	      else:
	        fermi=round(float(result[0]),3)
	
	#    with open('../../out.'+PREFIX+'.nscf') as f:
	#      for line in f :
	#          if re.search("highest", line):
	#            re_fermi = line.split()
	#            fermi=round(float(re_fermi[6]),3)
	#          if re.search("Fermi energy", line):
	#            re_fermi = line.split()
	#            fermi=round(float(re_fermi[4]),3)
	print('fermi: ',fermi)
	
	## LER DADOS DO DOS
	dostotdata =  np.genfromtxt(PREFIX+'.pdos_tot')
	energy = dostotdata[:,0]-fermi
	dostot = dostotdata[:,1]
	
	### pegar nspin do input do QE
	subject=grep('../../'+PREFIX+'.scf.in','nspin')
	#print('subject:', subject)
	#print(subject.strip().startswith('!'))
	if subject != None and not subject.strip().startswith('!') :
	    match = re.search("nspin\s{0,}.\s{0,}(\d{1,})", subject)  ## . para ignorar o = entre paranteses o que quer dar match, \d{1,} significa uma ou mais ocorrencias de digito, igualmente \s{0,} vai considerar nenhum ou algum espaco
	    if match:
	        nspin = int(match.group(1))
	    else:
	        print('Error: could not find nspin in file. Using nspin=1.')
	else:
	    print('Error: could not find nspin in file. Using nspin=1.')
	
	
	######## CALCULO DO GAP ################################
	# indice do nivel de fermi no DOS
	## energy == -0.05 para conseguir encontrar o gap em sistemas que a energia
	## de fermi tem um dos muito pequeno.
	index_fermi = np.nonzero(energy==0)[0]
	# pega o DOS acima do nivel de Fermi
	DOS_overfermi = dostot[energy>-0.5] 
	# pega os indices do DOS com valor consideravel
	index_DOS=np.nonzero(DOS_overfermi>0.1)[0]
	for i in range(len(index_DOS)-1):
	#se os indices apresentam uma descontinuidade e a regiao de gap
	  if (index_DOS[i] != index_DOS[i+1]-1) : #and (index_DOS[i] > 500 ):   
	    idxVB=index_DOS[i]  # indice da VB
	    idxCB=index_DOS[i+1] # indice da CB
	    break  # so pega a primeira descontinuidade
	try:
	  idxVB
	except:
	  print("Nao encontrou GAP!")
	  GAP=0
	else:
	  print("Existe gap, calculando...")
	  index_topVB = index_fermi  # indice da VB no DOS inteiro
	  index_botCB = idxCB-idxVB+index_fermi  # indice da CB no DOS inteiro
	# index_topVB = idxVB+index_fermi  # indice da VB no DOS inteiro
	# index_botCB = idxCB+index_fermi  # indice da CB no DOS inteiro
	  #print(index_topVB,index_botCB)
	  GAP=energy[index_botCB[0]]-energy[index_topVB[0]]
	   #gap e a subtraçao das energias
	  print('Valor do GAP pelo DOS:',GAP,' eV.')
	# ax.axvline(x=0, color='black',linestyle='--')  ## linha do nivel de FERMI
	
	if not os.path.isfile('../../gap_nscf'):
	  gap_nscf = GAP
	  print('There is no NSCF gap, using gap calculated from DOS.')
	else :
	  gap_nscf = float(np.genfromtxt('../../gap_nscf'))
	  print('Gap given by NSCF calculation will be used.')
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	ax.plot(energy,dostot,linewidth=1.0, label='Total DOS')
	ax.fill(energy,dostot,color='0.8')
	##################################################################
	
	
	# funcao para pegar os sumpdos_ de interesse
	# agrupa os arquivos em pdosfiles
	pdosfiles= [ ]
	for file in glob.glob('sumpdos_'+PREFIX+'.*'):
	#    print(file)
	    re_file = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}[0-9][a-z]', file)
	 #   print(re_file)
	    if re_file != [] :
	      pdosfiles.append(re_file[0])
	#print(pdosfiles)
	
	for file in pdosfiles :
	  ### encontrar os dados do orbital do arquivo em quetao
	  element = re.findall('sumpdos_'+PREFIX+'.([a-zA-Z]{1,})', file)
	  wfcnumber = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}([0-9])', file)
	  wfc = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}[0-9]([a-z])', file)
	#  print(file)
	  ### importa os dados para plotar
	  pdosdata = np.genfromtxt(file)
	  pdos = pdosdata[:,1]
	  if nspin == 2:
	      pdosup = pdosdata[:,2]
	  ax.plot(energy,pdos,linewidth=1.0, label=element[0]+wfcnumber[0]+wfc[0])
	  if nspin == 2:
	      ax.plot(energy,-pdosup,linewidth=1.0, label=element[0]+wfcnumber[0]+wfc[0])
	plt.title(graphtitle)
	plt.xlabel('E (eV)')
	plt.ylabel('States')
	plt.legend(loc='upper right')
	# plt.grid(True)
	# plt.rcParams.update({'font.size': 22})
	if min_x and max_x:
	 fromx,tox=min_x,max_x 
	DOSMAX=max(dostot[(energy>min_x) & (energy<max_x)])
	margin=DOSMAX*0.4
	ax.axis([fromx, tox, -margin, DOSMAX + margin])  
	if nspin == 2:
	  ax.axis([fromx, tox, -DOSMAX-margin, DOSMAX + margin])  
	### FLECHA DUPLA INDICANDO BAND GAP
	if (GAP != 0) and (GAP > 1.3 ):
	  ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	  ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	
	  ax.text(energy[index_topVB[0]]+GAP/2, 0.14*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	if (GAP != 0) and (GAP < 1.3 ):
	  ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	  ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	
	  ax.text(3, 0.92*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	fig = plt.gcf()
	if showgraph == 1:
	  plt.show()   
	fig.savefig(PREFIX+'_PDOS'+'.svg', format='svg', dpi=1000)  ## vai plotar PDOS de cada orbital, sem discriminar o m, up e down em cima e embaixo caso spinpolarizado
	fig.clf()
	plt.close()
	
	
	
	
	#####################################################
	### plotar DOStotal proj cada orbital e cada m embaixo
	#####################################################
	if nspin == 1:
	    NUM_COLORS=30
	    cm = plt.get_cmap('gist_rainbow')
	    fig = plt.figure()
	    ax = fig.add_subplot(111)
	    ax.plot(energy,dostot,linewidth=1.0, label='Total DOS')
	    ax.fill(energy,dostot,color='0.8')
	    ax.set_prop_cycle(cycler(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]))
	    # funcao para pegar os sumpdos_ de interesse
	    # agrupa os arquivos em pdosfiles
	    pdosfiles= [ ]
	    for file in glob.glob('sumpdos_'+PREFIX+'.*'):
	    #  print(file)
	        re_file = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}[0-9][a-z]', file)
	    #   print(re_file)
	        if re_file != [] :
	          pdosfiles.append(re_file[0])
	    # print(pdosfiles)
	    # funcao para pegar os sumpdos_ de interesse
	    # agrupa os arquivos em pdosfiles
	
	    for file in pdosfiles :
	      ### encontrar os dados do orbital do arquivo em quetao
	      element = re.findall('sumpdos_'+PREFIX+'.([a-zA-Z]{1,})', file)[0]
	      wfcnumber = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}([0-9])', file)[0]
	      wfc = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}[0-9]([a-z])', file)[0]
	
	      pdosdata = np.genfromtxt(file)
	      pdos = pdosdata[:,1]
	      try:
	        pdos_p1 = -pdosdata[:,2]
	        pdos_p2 = -pdosdata[:,3]
	        pdos_p3 = -pdosdata[:,4]
	        try:
	          pdos_p4 = -pdosdata[:,5]
	          pdos_p5 = -pdosdata[:,6]
	          try:
	            pdos_p6 = -pdosdata[:,7]
	            pdos_p7 = -pdosdata[:,8]
	          except:
	            pass
	        except:
	          pass
	      except:
	        pass
	      ax.plot(energy,pdos,linewidth=1.0, label=element+wfcnumber+wfc)
	      if wfc == 'p' :
	        ax.plot(energy,pdos_p1,linewidth=1.0, linestyle=':',label=element+wfcnumber+'$\mathregular{pz}$')
	        ax.plot(energy,pdos_p2,linewidth=1.0, linestyle='--',label=element+wfcnumber+'$\mathregular{px}$')
	        ax.plot(energy,pdos_p3,linewidth=1.0, linestyle='-.',label=element+wfcnumber+'$\mathregular{py}$')
	      if wfc == 'd' :
	        ax.plot(energy,pdos_p1,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{dz^2}$')
	        ax.plot(energy,pdos_p2,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{dzx}$')
	        ax.plot(energy,pdos_p3,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{dzy}$')
	        ax.plot(energy,pdos_p4,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{dx^2-y^2}$')
	        ax.plot(energy,pdos_p5,linewidth=1.0,linestyle='-',label=element+wfcnumber+'$\mathregular{dxy}$')
	      if wfc == 'f' :
	        ax.plot(energy,pdos_p1,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{f1}$')
	        ax.plot(energy,pdos_p2,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{f2}$')
	        ax.plot(energy,pdos_p3,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{f3}$')
	        ax.plot(energy,pdos_p4,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{f4}$')
	        ax.plot(energy,pdos_p5,linewidth=1.0,linestyle='-',label=element+wfcnumber+'$\mathregular{f5}$')
	        ax.plot(energy,pdos_p6,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{f6}$')
	        ax.plot(energy,pdos_p7,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{f7}$')
	    plt.title(graphtitle)
	    plt.xlabel('E (eV)')
	    plt.ylabel('States')
	
	
	    # Shrink current axis by 20%
	    box = ax.get_position()
	    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	    # Put a legend to the right of the current axis
	    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	    #plt.legend(loc='upper right')
	    # plt.grid(True)
	    # plt.rcParams.update({'font.size': 22})
	    if min_x and max_x:
	     fromx,tox=min_x,max_x 
	    margin=max(dostot[(energy>min_x) & (energy<max_x)])*0.4
	    ax.axis([fromx, tox, min(dostot[(energy>min_x) & (energy<max_x)])-margin, max(dostot[(energy>min_x) & (energy<max_x)]) + margin])  
	
	
	    DOSMAX=max(dostot[(energy>min_x) & (energy<max_x)])
	    ### FLECHA DUPLA INDICANDO BAND GAP
	    if (GAP != 0) and (GAP > 1.3 ):
	      ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	      ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	
	      ax.text(energy[index_topVB[0]]+GAP/2, 0.14*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	    if (GAP != 0) and (GAP < 1.3 ):
	      ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	      ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	
	      ax.text(3, 0.92*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	    fig = plt.gcf()
	    if showgraph == 1:
	      ax.show()  
	    fig.savefig(PREFIX+'_PDOS_orbitais_m'+'.svg', format='svg', dpi=1000)
	    fig.clf()
	    plt.close()
	
	if nspin == 2:
	    NUM_COLORS=30
	    cm = plt.get_cmap('gist_rainbow')
	    fig = plt.figure()
	    ax = fig.add_subplot(111)
	    ax.plot(energy,dostot,linewidth=1.0, label='Total DOS')
	    ax.fill(energy,dostot,color='0.8')
	    ax.set_prop_cycle(cycler(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]))
	    # funcao para pegar os sumpdos_ de interesse
	    # agrupa os arquivos em pdosfiles
	    pdosfiles= [ ]
	    for file in glob.glob('sumpdos_'+PREFIX+'.*'):
	    #  print(file)
	        re_file = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}[0-9][a-z]', file)
	    #   print(re_file)
	        if re_file != [] :
	          pdosfiles.append(re_file[0])
	    # print(pdosfiles)
	    # funcao para pegar os sumpdos_ de interesse
	    # agrupa os arquivos em pdosfiles
	
	    for file in pdosfiles :
	      ### encontrar os dados do orbital do arquivo em quetao
	      element = re.findall('sumpdos_'+PREFIX+'.([a-zA-Z]{1,})', file)[0]
	      wfcnumber = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}([0-9])', file)[0]
	      wfc = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}[0-9]([a-z])', file)[0]
	
	      pdosdata = np.genfromtxt(file)
	      pdos = pdosdata[:,1]
	      pdosup = pdosdata[:,2]
	      try:
	        pdos_p1 = -pdosdata[:,3]
	        pdosup_p1 = -pdosdata[:,4]
	        pdos_p2 = -pdosdata[:,5]
	        pdosup_p2 = -pdosdata[:,6]
	        pdos_p3 = -pdosdata[:,7]
	        pdosup_p3 = -pdosdata[:,8]
	        try:
	          pdos_p4 = -pdosdata[:,9]
	          pdosup_p4 = -pdosdata[:,10]
	          pdos_p5 = -pdosdata[:,11]
	          pdosup_p5 = -pdosdata[:,12]
	          try:
	            pdos_p6 = -pdosdata[:,13]
	            pdosup_p6 = -pdosdata[:,14]
	            pdos_p7 = -pdosdata[:,15]
	            pdosup_p7 = -pdosdata[:,16]
	          except:
	            pass
	        except:
	          pass
	      except:
	        pass
	    ### plot pdosdown
	      ax.plot(energy,pdos,linewidth=1.0, label=element+wfcnumber+wfc)
	      if wfc == 'p' :
	        ax.plot(energy,pdos_p1,linewidth=1.0, linestyle=':',label=element+wfcnumber+'$\mathregular{pz}$')
	        ax.plot(energy,pdos_p2,linewidth=1.0, linestyle='--',label=element+wfcnumber+'$\mathregular{px}$')
	        ax.plot(energy,pdos_p3,linewidth=1.0, linestyle='-.',label=element+wfcnumber+'$\mathregular{py}$')
	      if wfc == 'd' :
	        ax.plot(energy,pdos_p1,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{dz^2}$')
	        ax.plot(energy,pdos_p2,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{dzx}$')
	        ax.plot(energy,pdos_p3,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{dzy}$')
	        ax.plot(energy,pdos_p4,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{dx^2-y^2}$')
	        ax.plot(energy,pdos_p5,linewidth=1.0,linestyle='-',label=element+wfcnumber+'$\mathregular{dxy}$')
	      if wfc == 'f' :
	        ax.plot(energy,pdos_p1,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{f1}$')
	        ax.plot(energy,pdos_p2,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{f2}$')
	        ax.plot(energy,pdos_p3,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{f3}$')
	        ax.plot(energy,pdos_p4,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{f4}$')
	        ax.plot(energy,pdos_p5,linewidth=1.0,linestyle='-',label=element+wfcnumber+'$\mathregular{f5}$')
	        ax.plot(energy,pdos_p6,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{f6}$')
	        ax.plot(energy,pdos_p7,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{f7}$')
	    plt.title(graphtitle)
	    plt.xlabel('E (eV)')
	    plt.ylabel('States')
	
	
	    # Shrink current axis by 20%
	    box = ax.get_position()
	    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	    # Put a legend to the right of the current axis
	    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	    #plt.legend(loc='upper right')
	    # plt.grid(True)
	    # plt.rcParams.update({'font.size': 22})
	    if min_x and max_x:
	     fromx,tox=min_x,max_x 
	    margin=max(dostot[(energy>min_x) & (energy<max_x)])*0.4
	    ax.axis([fromx, tox, min(dostot[(energy>min_x) & (energy<max_x)])-margin, max(dostot[(energy>min_x) & (energy<max_x)]) + margin])  
	
	
	    DOSMAX=max(dostot[(energy>min_x) & (energy<max_x)])
	    ### FLECHA DUPLA INDICANDO BAND GAP
	    if (GAP != 0) and (GAP > 1.3 ):
	      ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	      ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	
	      ax.text(energy[index_topVB[0]]+GAP/2, 0.14*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	    if (GAP != 0) and (GAP < 1.3 ):
	      ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	      ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	
	      ax.text(3, 0.92*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	    fig = plt.gcf()
	    if showgraph == 1:
	      ax.show()  
	    fig.savefig(PREFIX+'_PDOSdown_orbitais_m'+'.svg', format='svg', dpi=1000)
	    fig.clf()
	    plt.close()
	
	    ### plot pdosup
	    ax.plot(energy,pdosup,linewidth=1.0, label=element+wfcnumber+wfc)
	    if wfc == 'p' :
	        ax.plot(energy,pdosup_p1,linewidth=1.0, linestyle=':',label=element+wfcnumber+'$\mathregular{pz}$')
	        ax.plot(energy,pdosup_p2,linewidth=1.0, linestyle='--',label=element+wfcnumber+'$\mathregular{px}$')
	        ax.plot(energy,pdosup_p3,linewidth=1.0, linestyle='-.',label=element+wfcnumber+'$\mathregular{py}$')
	    if wfc == 'd' :
	        ax.plot(energy,pdosup_p1,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{dz^2}$')
	        ax.plot(energy,pdosup_p2,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{dzx}$')
	        ax.plot(energy,pdosup_p3,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{dzy}$')
	        ax.plot(energy,pdosup_p4,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{dx^2-y^2}$')
	        ax.plot(energy,pdosup_p5,linewidth=1.0,linestyle='-',label=element+wfcnumber+'$\mathregular{dxy}$')
	    if wfc == 'f' :
	        ax.plot(energy,pdosup_p1,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{f1}$')
	        ax.plot(energy,pdosup_p2,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{f2}$')
	        ax.plot(energy,pdosup_p3,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{f3}$')
	        ax.plot(energy,pdosup_p4,linewidth=1.0,linestyle=':',label=element+wfcnumber+'$\mathregular{f4}$')
	        ax.plot(energy,pdosup_p5,linewidth=1.0,linestyle='-',label=element+wfcnumber+'$\mathregular{f5}$')
	        ax.plot(energy,pdosup_p6,linewidth=1.0,linestyle='-.',label=element+wfcnumber+'$\mathregular{f6}$')
	        ax.plot(energy,pdosup_p7,linewidth=1.0,linestyle='--',label=element+wfcnumber+'$\mathregular{f7}$')
	    plt.title(graphtitle)
	    plt.xlabel('E (eV)')
	    plt.ylabel('States')
	
	
	    # Shrink current axis by 20%
	    box = ax.get_position()
	    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	    # Put a legend to the right of the current axis
	    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	    #plt.legend(loc='upper right')
	    # plt.grid(True)
	    # plt.rcParams.update({'font.size': 22})
	    if min_x and max_x:
	     fromx,tox=min_x,max_x 
	    margin=max(dostot[(energy>min_x) & (energy<max_x)])*0.4
	    ax.axis([fromx, tox, min(dostot[(energy>min_x) & (energy<max_x)])-margin, max(dostot[(energy>min_x) & (energy<max_x)]) + margin])  
	
	
	    DOSMAX=max(dostot[(energy>min_x) & (energy<max_x)])
	    ### FLECHA DUPLA INDICANDO BAND GAP
	    if (GAP != 0) and (GAP > 1.3 ):
	      ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	      ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	
	      ax.text(energy[index_topVB[0]]+GAP/2, 0.14*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	    if (GAP != 0) and (GAP < 1.3 ):
	      ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	      ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	
	      ax.text(3, 0.92*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	    fig = plt.gcf()
	    if showgraph == 1:
	      ax.show()  
	    fig.savefig(PREFIX+'_PDOSup_orbitais_m'+'.svg', format='svg', dpi=1000)
	    fig.clf()
	    plt.close()
	
	
	### plotar DOStotal cada atomo
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(energy,dostot,linewidth=1.0, label='Total DOS')
	ax.fill(energy,dostot,color='0.8')
	
	# funcao para pegar os sumpdos_ de interesse, nesse caso somente 
	# os do elementos individuais
	# agrupa os arquivos em pdosfiles
	pdosfiles= [ ]
	
	for file in glob.glob('sumpdos_'+PREFIX+'.*'):
	    print(file)
	    re_file = re.findall('sumpdos_'+PREFIX+'.[a-zA-Z]{1,}$', file)
	#    print(re_file)
	    if re_file != [] :
	      pdosfiles.append(re_file[0])
	# print(pdosfiles)
	
	
	for file in pdosfiles :
	  element = re.findall('sumpdos_'+PREFIX+'.([a-zA-Z]{1,})', file)[0]
	  pdosdata = np.genfromtxt(file)
	  pdos = pdosdata[:,1]
	  if nspin == 2:
	      pdosup = pdosdata[:,2]
	  ax.plot(energy,pdos,linewidth=1.0, label=element[0]+wfcnumber[0]+wfc[0])
	  if nspin == 2:
	      ax.plot(energy,-pdosup,linewidth=1.0, label=element[0]+wfcnumber[0]+wfc[0])
	plt.title(graphtitle)
	plt.xlabel('E (eV)')
	plt.ylabel('States')
	ax.legend(loc='upper right')
	# plt.grid(True)
	# plt.rcParams.update({'font.size': 22})
	
	if min_x and max_x:
	 fromx,tox=min_x,max_x 
	DOSMAX=max(dostot[(energy>min_x) & (energy<max_x)])
	margin=DOSMAX*0.4
	ax.axis([fromx, tox, -margin, DOSMAX + margin])  
	if nspin == 2:
	  ax.axis([fromx, tox, -DOSMAX-margin, DOSMAX + margin])  
	### FLECHA DUPLA INDICANDO BAND GAP
	if (GAP != 0) and (GAP > 1.3 ):
	  ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	  ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.6, head_length=0.5, color = 'black',  length_includes_head = True)
	
	  ax.text(energy[index_topVB[0]]+GAP/2, 0.14*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	if (GAP != 0) and (GAP < 1.3 ):
	  ax.arrow(energy[index_topVB[0]], 0.04*DOSMAX, GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	  ax.arrow(energy[index_botCB[0]], 0.04*DOSMAX, -GAP, 0, head_width=0.3, head_length=0.2, color = 'black',  length_includes_head = True)
	
	  ax.text(3, 0.92*DOSMAX, r'$E_g= $'+str(round(gap_nscf,2))+' eV', horizontalalignment='center')
	
	fig = plt.gcf()
	if showgraph == 1:
	  plt.show()
	fig.savefig(PREFIX+'_PDOS_elements'+'.svg', format='svg', dpi=1000)
	fig.clf()
	plt.close()

################################### PLOT BANDS ######################################################
################################### PLOT BANDS ######################################################
################################### PLOT BANDS ######################################################
################################### PLOT BANDS ######################################################
################################### PLOT BANDS ######################################################


def plot_bands(prefix,indexes,fermi=0.0,label="",color="black",RANGE_VB=1.5,RANGE_CB=5,
             redstates=[],  ## states colored red, you can always give range(92,96)
             greenstates=[],  ## states colored green
             bluestates=[],  ## states colored blue, no states with []
             labelsprojs=["","",""],  ## first red, second green, third blue
## notice if there is overlapping colors will blend
             maxcontrast=True,  ## notice if maxcontrast=0, color will be as strong as the state proportion, 
             #usually this contribution is very low compared to total, therefore maxcontrast is default.
             contrast_ratios=[1,1,1],  ## only if maxcontrast true # values should be in 0-1 range ## can use to
             ## obtain custom colors
             plot_emass_fitting=False,
             emass_calculation=False,
             print_svg=False,
             kpoint_section_thresold=0.2,  ### this MAY REQUIRE FINE TUNING, check the number of sections on output to see it it is correct.
             **kwargs
      ):
###################### SETTING VARIABLES ######################################
## files needed  ## give in list format
    PREFIX=[prefix]
    filesdos=[prefix+'.dos']  ## dos file from QE ## will be used only to find fermi if not given
    filespdos=['out.'+prefix+".pdos_bands"]  ## pdos output of a bandstructure calculation  ## contains most of data we need
    symfiles=['out.'+prefix+".bands"]  ## bands output of bandstructure calculation ## contains high symmetry k points
    fermis = [fermi]   ## keeping 0 the script will look for fermi value in filedos file.
    labels=[label] ## if label "" no label will be plotted
    colors=[color]
    
    #indexes=['$\Gamma$','M','K','$\Gamma$','A','L','H','A','L','M',
    #         'H','K']      ## these are the indexes of the bandstructure
    # Γ—M—K—Γ—A—L—H—A|L—M|H—K ## band path P-3m1 trigonal
    #   Γ—Y—F0|Δ0—Γ—Z—B0|G0—T—Y|Γ—S—R—Z—T  ## band path hexagonal CsSbI 63/mmc
    #   Γ—X—S—Y—Γ—Z—U—R—T—Z|X—U|Y—T|S—R   ## band path ortorhombic CsSbCl 
    
    #kpoint_section_thresold=0.1  ### this MAY REQUIRE FINE TUNING, check the number of sections on output to see it it is correct.
    ### less section than expected, reduce this value, more sections than expected, increase this value. 
    legendposition_tuple_lines=(1,1)
    legendposition_tuple_patches=(1,1)
    
    #RANGE_VB=1.5 ## ev to plot below efermi
    #RANGE_CB=5 ## ev to plot above efermi
    ## legend positioning in the end of script as well as output options.
    
    ## sometimes we want to shift conduction band applying a scissor operation, 
    # define its value in here in eV ## default is 0
    scissor_shift=0.0
    
    ## COLOR ASSIGNMENT
    projected_bands=True ## bands will be projected with different colors in different states.
    ## to fill appropriatedly run beforehand this script with "checkstates" argument.
    ## once you write down the states youre interested to color in the bandstructure proceed.
    ## the numbers given are the number after # in the states list given by checkstates run.
    #redstates=[]  ## states colored red, you can always give range(92,96)
    #greenstates=[]  ## states colored green
    #bluestates=[]  ## states colored blue, no states with []
    #labelsprojs=["","","Mo"]  ## first red, second green, third blue
    
    ## notice if there is overlapping colors will blend
    
    #maxcontrast=True  ## notice if maxcontrast=0, color will be as strong as the state proportion, 
    #usually this contribution is very low compared to total, therefore maxcontrast is default.
    #contrast_ratios=[1,1,1]  ## only if maxcontrast true # values should be in 0-1 range
    #if [1,0.5,0.5], red will be 2 times more intense than highest green and blue.
    
    ########### MANUAL INDEXING ##################################
    # maybe the script does not position the labels of high symmetry appropriatedly, this section is to correct
    # this problem manually, as well as few changes on plot settings
    manual_indexing=False
    indexes_pos=[]  ## ought to have same size of indexes list
    
    ########### EFFECTIVE MASS SETUP ###f############################
    #emass_calculation=False ## will print effective mass on high symmetry points, use manual indexing for arbitrary points
    #plot_emass_fitting=False ## will show the fitted curve to obtain effective masses
    
    if plot_emass_fitting :
        print(f"K point coordinate values must be scalated according to QE alat. Beware when calculating effective mass!!")

    #alat=16.365  ## alat used by quantum espresso have to know to change coords of k points to calculate eff mass
    #k_width=0.15 ## range above and below hspoint in default k-coordinates to consider in band 0.15 usually enough, 
    # can always check turning plot_emass_fitting to True.
    ylevel_refpoint=0.5  ## ref point to calculate effective mass remember fermi level y=0. Important this level between gap!
    ## shift ylevel_refpoint to a value between valence and conduction band to calculate bands.
    ##########################################################################################
    ########### LEAVE EVERYTHING ELSE AS IT IS, EXCEPT IF YOU KNOW WHAT YOURE DOING ##########
    ##########################################################################################
    import re
    import numpy as np
    from math import sqrt
    import matplotlib.pyplot as plt
    import os,sys
    import matplotlib as mpl
    from matplotlib.collections import LineCollection
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from numpy.lib.function_base import diff
    import colorsys
    
    
    def modulo(x,y):
        """Module of the difference between two vectors"""
        return np.linalg.norm(np.subtract(x,y))
    
    
    
    ### this will allow to plot specific states as different colors
    def rgbline(x, y, red, green, blue, alpha=0, linestyles="solid",
                linewidth=2.5):
        """Get a RGB coloured line for plotting.
        Args:
            x (list): x-axis data.
            y (list): y-axis data (can be multidimensional array).
            red (list): Red data (must have same shape as ``y``).
            green (list): Green data (must have same shape as ``y``).
            blue (list): blue data (must have same shape as ``y``).
            alpha (:obj:`list` or :obj:`int`, optional): Alpha (transparency)
                data (must have same shape as ``y`` or be an :obj:`int`).
            linestyles (:obj:`str`, optional): Linestyle for plot. Options are
                ``"solid"`` or ``"dotted"``.
        """
        y = np.array(y)
        if len(y.shape) == 1:
            y = np.array([y])
            red = np.array([red])
            green = np.array([green])
            blue = np.array([blue])
            alpha = np.array([alpha])
        elif isinstance(alpha, int):
            alpha = [alpha] * len(y)
    
        seg = []
        colours = []
        for yy, rr, gg, bb, aa in zip(y, red, green, blue, alpha):
            pts = np.array([x, yy]).T.reshape(-1, 1, 2)
            seg.extend(np.concatenate([pts[:-1], pts[1:]], axis=1))
    
            nseg = len(x) - 1
            r = [0.5 * (rr[i] + rr[i + 1]) for i in range(nseg)]
            g = [0.5 * (gg[i] + gg[i + 1]) for i in range(nseg)]
            b = [0.5 * (bb[i] + bb[i + 1]) for i in range(nseg)]
            a = np.ones(nseg, np.float) * aa
            colours.extend(list(zip(r, g, b, a)))
    
        lc = LineCollection(seg, colors=colours, rasterized=True,
                            linewidth=linewidth, linestyles=linestyles)
        return lc
    
    ## cubic interpolation is necessary for our bands
    def cubic_interp1d(x0, x, y):
        """
        Interpolate a 1-D function using cubic splines.
        x0 : a float or an 1d-array
        x : (N,) array_like
            A 1-D array of real/complex values.
        y : (N,) array_like
            A 1-D array of real values. The length of y along the
            interpolation axis must be equal to the length of x.
    
        Implement a trick to generate at first step the cholesky matrice L of
        the tridiagonal matrice A (thus L is a bidiagonal matrice that
        can be solved in two distinct loops).
    
        additional ref: www.math.uh.edu/~jingqiu/math4364/spline.pdf 
        """
        x = np.asfarray(x)
        y = np.asfarray(y)
    
        # remove non finite values
        # indexes = np.isfinite(x)
        # x = x[indexes]
        # y = y[indexes]
    
        # check if sorted
        if np.any(np.diff(x) < 0):
            indexes = np.argsort(x)
            x = x[indexes]
            y = y[indexes]
    
        size = len(x)
    
        xdiff = np.diff(x)
        ydiff = np.diff(y)
    
        
        # allocate buffer matrices
        Li = np.empty(size)
        Li_1 = np.empty(size-1)
        z = np.empty(size)
    
        # fill diagonals Li and Li-1 and solve [L][y] = [B]
        Li[0] = sqrt(2*xdiff[0])
        Li_1[0] = 0.0
        B0 = 0.0 # natural boundary
        z[0] = B0 / Li[0]
    
        for i in range(1, size-1, 1):
            Li_1[i] = xdiff[i-1] / Li[i-1]
        # print(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
        # print(i,ydiff)
            Li[i] = sqrt(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
    
            Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
            z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]
    
        i = size - 1
        Li_1[i-1] = xdiff[-1] / Li[i-1]
        Li[i] = sqrt(2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])
        Bi = 0.0 # natural boundary
        z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]
    
        # solve [L.T][x] = [y]
        i = size-1
        z[i] = z[i] / Li[i]
        for i in range(size-2, -1, -1):
            z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]
    
        # find index
        index = x.searchsorted(x0)
        np.clip(index, 1, size-1, index)
    
        xi1, xi0 = x[index], x[index-1]
        yi1, yi0 = y[index], y[index-1]
        zi1, zi0 = z[index], z[index-1]
        hi1 = xi1 - xi0
    
        # calculate cubic
        f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
            zi1/(6*hi1)*(x0-xi0)**3 + \
            (yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
            (yi0/hi1 - zi0*hi1/6)*(xi1-x0)
        return f0
    
    
    def Symmetries(fstring):
        f = open(fstring,'r')
        x = np.zeros(0)
        for i in f:
            if "high-symmetry" in i:
                x = np.append(x,float(i.split()[-1]))
        f.close()
        return x
    ###################### FUNCTIONS AND PROCEDURES ##########################################
    ##########################################################################################
    
    n_bandsections=-1
    for id_structure in range(len(PREFIX)):
        filedos=filesdos[id_structure] 
        filepdos=filespdos[id_structure] 
        symfile=symfiles[id_structure]
        fermi=fermis[id_structure]  
        label_structure=labels[id_structure]
        color=colors[id_structure]
    
        ## following snips determine fermi energy, if spin polarized and pdos file type.
        if fermi == 0 :
            try:
                with open(filedos) as f:
                    first_line = f.readline()
                    result = re.findall(r"Fermi *[^\w ] *(.*) eV", first_line, re.IGNORECASE | re.MULTILINE)
                    fermi=round(float(result[0]),3)
            except IOError:
                print("Dos file not accessible, cannot proceed.")
        print("Fermi ",fermi)
    
        nspin=-1 ## have to declare first otherwise is local
        try:
            with open(filedos) as f:
                lines = f.readlines()
                if len(lines[2].split()) > 3:
                    nspin=True
                    print("Spin polarized detected. Up and down spins will be plotted in the same sections. ",fermi)
                else:
                    nspin=False
                    print("This calculation is not spin polarized",fermi)
        except IOError:
            print("Dos file not accessible, cannot proceed.")
    
    
        bandtype=-1 ## energies dont have ==== in pdos file  ## two types of bandtype
        try:
            with open(filepdos) as f:
                pdosdata = f.read()  ## gets all data in single string
                if pdosdata.find("==== e") == -1:
                    bandtype=1
                    print("Bandtype 1. ")
                else:
                    bandtype=0
                    print("Bandtype 0. ")
        except IOError:
            print("PDOS file not accessible, cannot proceed.")
    
    
        ##########################################################################################
        ##########################################################################################
    
    
    
        ## the following loop will look on filepdos file and find in the text, first the different states contained in the structure,
        ## each orbital of each atom will produce a different state. Its important to know these states so as to find the
        ## contributions we want colored in the final plot.
        ## the code will look then for each kpoint and add to the kpoints array and look for each energy of this kpoint
        ## this will be the raw bandplot, without the compositions
        kpoints=[]
        fullenergykpts=[]
        energykpt=[]
        states=[]
        with open(filepdos, "r") as f:
            lines=f.readlines()
            for line in lines:
                state=re.findall(r"^\s*state(.*)",line)
                if state != []: 
                    states.append([re.findall(r".*(#\s*[0-9]*):.*",state[0]), state]) 
        #            states.append(state)
                if line.startswith(" k"):
                    kpt = re.findall(r"([0-9]\.[0-9]+)",line)
                    kpoints.append(list(map(float,kpt)))
                    if energykpt != [] :  ## to close energykpt if new kpoint is found
                        fullenergykpts.append(list(map(float,energykpt)))
                        energykpt=[]
        #        energyk=re.findall(r"^====.*=.*(-?[0-9]\.[0-9]+).*",line)
                if bandtype == 0 : ## there are two band types in .pdos in QE, have to check beforehand
                    if line.startswith("===="):
                    #energyk=re.findall(r"^====*.*(.[0-9]+\.[0-9]+).*",line)[0]
                        energyk=re.findall(r"^====*.* ([-+]?\d+\.\d+).*",line)[0]
                        energykpt.append(energyk)
                if bandtype == 1 :
                    if re.match(r'\s+e =', line):
                        energyk=re.findall(r"^\s*e =*.* ([-+]?\d+\.\d+).*",line)[0]
                        energykpt.append(energyk)
                if line.startswith("   JOB DONE."):  ## check job done and appends last energykpt
                    fullenergykpts.append(list(map(float,energykpt)))
                    energykpt=[]
        fullenergykpts=np.array(fullenergykpts)
    
        ################ CHECKSTATES RUN #############################
        if (len(sys.argv)>1) and (sys.argv[1] == "checkstates"):
            for i in range(len(states)):
                print(states[i][1][0])
            sys.exit()  ### check states will only print the states to check which ones to plot colored
        #######################################################################################
    
        ## states executed before gave state found in the format '# N' and the line found, we just want the '# N' match.
        statesidx=[ states[i][0][0] for i in range(len(states))]
    
        ## once raw band with number of kpoints and energies for each was found, we declare an zeroed
        ## numpy array that will contain the value of contribution of each state for each energy of each kpoint
        projections=np.zeros((len(kpoints),len(fullenergykpts[0]),len(states)))
    
        ## before filling projections array, we need to estabilish the kpointspath coordinate,
        ##  this will be the x of our bandstructure. There are complications when the band is divided in sections.
        mod=0
        mod_section_list=[0]
        fullkpointspath=[0]
        ## creates kpoint list direct from file
        for i in range(len(kpoints)-1): 
            mod=modulo(kpoints[i+1],kpoints[i])+mod
            fullkpointspath.append(mod) 
        ## the kpoints in reciprocal space have their distance calculated, there will be huge differences when
        ## the band structure is discontinuous.
        ##find sections of kpointspath where there is band discontinuity, a 0.2 value of change is taken as enough
        ## for a discontinuity to be defined, usually the difference in kpoints should be less than 0.1 in a good
        ## band structure.
        print('kst',kpoint_section_thresold)
        bandsections=np.where(np.diff(fullkpointspath)>kpoint_section_thresold)[0]+1 ## add 1 because ndiff takes the previous value
        print(" A total of %2d sections were found. Is that what you expected? " % (len(bandsections)+1) )
        ### bandsections holds indexes of the discontinuity the fullkpointspath 
        bandsections=np.insert(bandsections, 0, 0, axis=0)  ## add 0 in beginning
        bandsections=np.insert(bandsections, len(bandsections), -1, axis=0) 
    
    
        ## add last term as -1, that is, end of kpointspath
    
        ## calculate the path differences to correct the discontinuities and have continuous x along the kpath
        endpath=[fullkpointspath[i-1] for i in bandsections[1:-1]] 
        beginnextpath=[fullkpointspath[i] for i in bandsections][1:-1]
        pathdifferences=np.array(beginnextpath)-np.array(endpath)
        pathdifferences=np.insert(pathdifferences, 0, 0, axis=0)  ## add 0 in beginning
        pathdifferences=np.cumsum(pathdifferences)  ## pathdifferences is cummulative 
        ## pathdifference will be subtracted on each beggining of subsection.
    
        ### kpoints path subdivided
        ### there is a calculation of deltas which is the normalized size of each subsection
        ### this comes handy for the plotting of the final bandstructure to be proportional
        ### to the size of each section.
        bandsections=list(bandsections[:-1])+[None]  
        ## had to change bandsections again otherwise kpointspath is cut short in last section
        deltas=[]
        kpointspaths=[]
        lastkpointcoord=0
        fullkpointspath=np.array(fullkpointspath)
        for m in range(len(bandsections)-1):
            kpointspath=fullkpointspath[bandsections[m]:bandsections[m+1]]-pathdifferences[m]
            kpointspaths.append(kpointspath)
            deltas.append(max(kpointspath)-min(kpointspath))
        deltas=np.array(deltas)  
        sumdelta=np.sum(np.array(deltas))
        deltas=deltas/sumdelta
    
        ## if spin polarized there will be double of sections detected, have to convolute these in the same graph
        
        if id_structure == 0 :  ## only does it for the first plot.
            if nspin:
                n_bandsections=int((len(bandsections)-1)/2)  ## bandsections will be half
                print(bandsections,n_bandsections)
                fig, ax = plt.subplots(1,n_bandsections,sharey=True,
                                squeeze=True, gridspec_kw={'wspace': 0.3,'width_ratios': deltas[:-n_bandsections]})
    
            else:
                n_bandsections=(len(bandsections)-1)
                print(n_bandsections)
                fig, ax = plt.subplots(1,n_bandsections,sharey=True,
                                squeeze=True, gridspec_kw={'wspace': 0.3,'width_ratios': deltas})    
    
#        if nspin and n_bandsections == 2: ## no subsections on band plot
#            ax=[ax]  ## no subscript problem
    
        if n_bandsections == 1: ## no subsections on band plot
            ax=[ax]  ## no subscript problem
        print('axes',ax,ax[0],'nspin',nspin) 
        ax[0].set_ylabel("Energy (eV)")
    
        hspoints = Symmetries(symfile)
        print("A total of %2d highsymmetry points were found by QE. Is that expected?" % (len(hspoints)))
        ## hspoints contains the points found in symfile the out.prefix.bands file in which 
        # QE finds symmetry lines of band structure
    
        ### finally we proceed to fill the projections array. 
        ind=0 ## index of letters of high symmetry 
        ## has to be outside sectioning loop to proceed lettering otherwise restarts every section
        for m in range(n_bandsections):  ## this loop is repeated for each band section
            energykpts=fullenergykpts[bandsections[m]:bandsections[m+1]]
            i=-1
            j=-1
            with open(filepdos, "r") as f:
                lines=f.readlines()
                for line in lines:
                    if line.startswith(" k"):
                        i+=1
                        j=-1 ##resets others
                    if bandtype == 0 : ## there are two band types in .pdos in QE, have to check beforehand
                        if line.startswith("===="):
                            j+=1
                            energyk=re.findall(r"^====*.* ([-+]?\d+\.\d+).*",line)[0]
                            energykpt.append(energyk)
                    if bandtype == 1 :
                        if re.match(r'\s+e =', line):
                            j+=1
                    ## this will look for the projection and for the identification of the state
                    if re.match(r"^\s* psi =", line): 
                        result=re.findall(r"(\d+\.\d+)",line)
                        label=re.findall(r"\[([^]]*)\]", line)
                        newlabel = [ statesidx.index(label[i]) for i in range(len(label)) ]
                        results= [[x,y] for x,y in zip(result,newlabel)]
                        for jj in range(len(results)):
                            projections[i][j][results[jj][1]]=results[jj][0]
    
                    if re.match(r"^\s* \+",line):
                        result=re.findall(r"(\d+\.\d+)",line)
                        label=re.findall(r"\[([^]]*)\]", line)
                        newlabel = [ statesidx.index(label[i]) for i in range(len(label)) ]
                        results= [[x,y] for x,y in zip(result,newlabel)]
                        for jj in range(len(results)):
                            projections[i][j][results[jj][1]]=results[jj][0]
            ## remember the band plot has to be interpolated, so after the states are filled we proceed to this
            kinterp = np.linspace(min(kpointspaths[m]),  max(kpointspaths[m]), num=10*len(kpointspaths[m]), endpoint=True)
            #interpolation generates 10 times the previous kpoints in each section
    
            if not plot_emass_fitting: ## useful in plot_emass
                ax[m].set_xticklabels([]) ## no xticklabels we want the special characters
            axis = [min(kinterp),max(kinterp),-RANGE_VB, RANGE_CB]  ## valence band and conduction band ranges are set in here
            ax[m].set_ylim([axis[2],axis[3]])
            ax[m].set_xlim([axis[0],axis[1]])
            ### SETTING TICKS IN Y AXIS
            majoryLocator   = MultipleLocator(1)
            majoryFormatter = FormatStrFormatter('%d')
            minoryLocator   = MultipleLocator(0.5)
            ax[m].yaxis.set_major_locator(majoryLocator)
            ax[m].yaxis.set_major_formatter(majoryFormatter)
            ax[m].yaxis.set_minor_locator(minoryLocator)
            
            for i in range(len(energykpts[0])):  ## for every band 
                ## here interpolation function is called
                if min(energykpts[:,i]) > fermi+0.1:   ##### creating an artificial shift on the graph to match NC PP
                    bandinter= cubic_interp1d(kinterp, kpointspaths[m], energykpts[:,i]+scissor_shift-fermi)
                else:
                    bandinter= cubic_interp1d(kinterp, kpointspaths[m], energykpts[:,i]-fermi)
                
                ## now the color vectors same size of kinterp are created
                red=np.absolute(cubic_interp1d(kinterp, kpointspaths[m],np.zeros(len(kpointspaths[m]))))
                for r in np.array(redstates)-1: ## each color states are user defined after checkstates run
                    red+=np.absolute(cubic_interp1d(kinterp, kpointspaths[m],projections[:,i,r]))
                if maxcontrast: ## only if maxcontrast option is true
                    red=(contrast_ratios[0]/max(red))*red if max(red) != 0 else red ## get maxs contrast
                    ## contrast ratios values are inserted in this expression
    
                green=np.absolute(cubic_interp1d(kinterp, kpointspaths[m],np.zeros(len(kpointspaths[m]))))
                for g in np.array(greenstates)-1:
                    green+=np.absolute(cubic_interp1d(kinterp, kpointspaths[m],projections[:,i,g]))
                if maxcontrast:
                    green=(contrast_ratios[1]/max(green))*green if max(green) != 0 else green ## get maxs contrast
    
                blue=np.absolute(cubic_interp1d(kinterp, kpointspaths[m],np.zeros(len(kpointspaths[m]))))               
                for b in np.array(bluestates)-1:
                    blue+=np.absolute(cubic_interp1d(kinterp, kpointspaths[m],projections[:,i,b]))
                if maxcontrast:  
                    blue=(contrast_ratios[2]/max(blue))*blue   if max(blue) != 0 else blue ## get maxs contrast
    
                ## need to find hue of rgb 
                hue=np.zeros(len(red)) 
                for i in range(len(hue)):
                    hue[i]=colorsys.rgb_to_hsv(red[i],green[i],blue[i])[0] ## get hue
                ## then saturation based on weights of each of the considered states
                satur=np.zeros(len(red)) 
                for i in range(len(satur)):
                    satur[i]=red[i]+green[i]+blue[i]
                if maxcontrast:
                    satur=(1/max(satur))*satur   if max(satur) != 0 else satur 
                
                ## now satur hue and value=1 will lead to the new red green blue values
                for i in range(len(hue)):
                    red[i],green[i],blue[i]=colorsys.hsv_to_rgb(hue[i], satur[i], 1)
                ## now we should colors that vary from white to the mixed color depending on weights on band states
            
    
                ## here colored lines are finally ploted
                if projected_bands:
                    ax[m].plot(kinterp, bandinter, color=color,label=label_structure,lw=1)
                    lc = rgbline(kinterp, bandinter, red, green, blue, alpha=1, linestyles='solid',
                                         linewidth=(mpl.rcParams['lines.linewidth'] * 2.75))
                    ax[m].add_collection(lc)
                else:
                    ax[m].plot(kinterp, bandinter, color=color,label=label_structure,lw=1.5)
    
                
                    
    
            ### now last touch is to have the labels of high symmetry points placed appropriatedly on the band structure
            ## the symbols are given by indexes function, the number should naturally match number of highsym points in 
            ## the structure
    
            if id_structure == 0 :
                if not manual_indexing :   ## only first plot needs this setting.
                    lastp=-1 ## last point coordinate to not put two labels at same place. Changes after first run.
                    #print(hspoints)
                    xticks=[] ## xticklabels to be set later
                    for p in hspoints: #This is the high symmetry lines
                        ## checks if p is within range of the max and min kinterp, therefore placeble and if his distance from last
                        ## placed point is larger than 0.1
                        if p <= max(kinterp)+0.1 and p >= min(kinterp)-0.1 and abs(lastp-p)>0.1:
                            lastp=p
                            x1 = [p,p]
                            x2 = [axis[2],axis[3]]
                            if not emass_calculation: ## these lines affect emass calculation
                                ax[m].plot(x1,x2,'--',lw=0.55,color='black',alpha=0.75)
                            xticks.append(p)
                            ax[m].text(p-0.07,axis[2]-(axis[3]-axis[2])/10, indexes[ind],size=14)  ## to adjust labels position
                            ind+=1
                    ax[m].set_xticks(xticks)
                else: ## manual_indexing set to True
                    for i in range(len(indexes)):
                        x1 = [indexes_pos[i],indexes_pos[i]]
                        # print([j,j])
                        x2 = [axis[2],axis[3]]
                        ax[m].plot(x1,x2,'--',lw=0.55,color='black',alpha=0.75)
                        ax[m].text(indexes_pos[i],axis[2]-(axis[3]-axis[2])/10, indexes[i],size=14)
                        #ax[m].set_xticklabels(indexes)  ## works but labels are not uniform
    
        if nspin:  ## this will plot other spin on top of first spin sections
            
            ## the .pdos file is organized such that kpoints of other spin just repeat after first spin is over
            
            for m in range(n_bandsections):  ## this loop is repeated for each band section
                # print('section',m)
                energykpts=fullenergykpts[bandsections[m+n_bandsections]:bandsections[m+n_bandsections+1]]
                
                ## now bandsections are further to get other spin
                i=-1
                j=-1
                with open(filepdos, "r") as f:
                    lines=f.readlines()
                    for line in lines:
                        if line.startswith(" k"):
                            i+=1
                            j=-1 ##resets others
                        if bandtype == 0 : ## there are two band types in .pdos in QE, have to check beforehand
                            if line.startswith("===="):
                                j+=1
                                energyk=re.findall(r"^====*.* ([-+]?\d+\.\d+).*",line)[0]
                                energykpt.append(energyk)
                        if bandtype == 1 :
                            if re.match(r'\s+e =', line):
                                j+=1
                        ## this will look for the projection and for the identification of the state
                        if re.match(r"^\s* psi =", line): 
                            result=re.findall(r"(\d+\.\d+)",line)
                            label=re.findall(r"\[([^]]*)\]", line)
                            newlabel = [ statesidx.index(label[i]) for i in range(len(label)) ]
                            results= [[x,y] for x,y in zip(result,newlabel)]
                            for jj in range(len(results)):
                                projections[i][j][results[jj][1]]=results[jj][0]
    
                        if re.match(r"^\s* \+",line):
                            result=re.findall(r"(\d+\.\d+)",line)
                            label=re.findall(r"\[([^]]*)\]", line)
                            newlabel = [ statesidx.index(label[i]) for i in range(len(label)) ]
                            results= [[x,y] for x,y in zip(result,newlabel)]
                            for jj in range(len(results)):
                                projections[i][j][results[jj][1]]=results[jj][0]
                ## remember the band plot has to be interpolated, so after the states are filled we proceed to this
                kinterp = np.linspace(min(kpointspaths[m+n_bandsections]),  max(kpointspaths[m+n_bandsections]), 
                                    num=10*len(kpointspaths[m+n_bandsections]), endpoint=True)
                #interpolation generates 10 times the previous kpoints in each section
                ### have to shift back the second spin k to overlap first spin k point region
                previous_kinterp = np.linspace(min(kpointspaths[m]),  max(kpointspaths[m]), 
                                    num=10*len(kpointspaths[m]), endpoint=True)
                diffk_2spin=min(kinterp)-min(previous_kinterp)  ## get difference of nspin 1 and nspin 2 x of bands
                kinterp=kinterp-diffk_2spin  ## make 2spin overlap previous 1spin.
                for i in range(len(energykpts[0])):  ## for every band 
                    ## here interpolation function is called
                    if min(energykpts[:,i]) > fermi+0.1:   ##### creating an artificial shift on the graph to match NC PP
                        bandinter= cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin, energykpts[:,i]+scissor_shift-fermi)
                    else:
                        bandinter= cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin, energykpts[:,i]-fermi)
                        
                    
                    ## now the color vectors same size of kinterp are created
                    red=np.absolute(cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin,np.zeros(len(kpointspaths[m+n_bandsections]))))
                    for r in np.array(redstates)-1: ## each color states are user defined after checkstates run
                        red+=np.absolute(cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin,projections[:,i,r]))
                    if maxcontrast: ## only if maxcontrast option is true
                        red=(contrast_ratios[0]/max(red))*red if max(red) != 0 else red ## get maxs contrast
                        ## contrast ratios values are inserted in this expression
    
                    green=np.absolute(cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin,np.zeros(len(kpointspaths[m+n_bandsections]))))
                    for g in np.array(greenstates)-1:
                        green+=np.absolute(cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin,projections[:,i,g]))
                    if maxcontrast:
                        green=(contrast_ratios[1]/max(green))*green if max(green) != 0 else green ## get maxs contrast
    
                    blue=np.absolute(cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin,np.zeros(len(kpointspaths[m+n_bandsections]))))          
                    for b in np.array(bluestates)-1:
                        blue+=np.absolute(cubic_interp1d(kinterp, kpointspaths[m+n_bandsections]-diffk_2spin,projections[:,i,b]))
                    if maxcontrast:  
                        blue=(contrast_ratios[2]/max(blue))*blue   if max(blue) != 0 else blue ## get maxs contrast
    
                    ## need to find hue of rgb 
                    hue=np.zeros(len(red)) 
                    for i in range(len(hue)):
                        hue[i]=colorsys.rgb_to_hsv(red[i],green[i],blue[i])[0]
                    ## then saturation based on weights of each of the considered states
                    satur=np.zeros(len(red)) 
                    for i in range(len(satur)):
                        satur[i]=red[i]+green[i]+blue[i]
                    if maxcontrast:
                        satur=(1/max(satur))*satur   if max(satur) != 0 else satur 
                    
                    ## now satur hue and value=1 will lead to the new red green blue values
                    for i in range(len(hue)):
                        red[i],green[i],blue[i]=colorsys.hsv_to_rgb(hue[i], satur[i], 1)
                    ## now we should colors that vary from white to the mixed color depending on weights on band states
    
                    ## here colored lines are finally ploted
                    if projected_bands:
                        ax[m].plot(kinterp, bandinter, linestyle='--',color=color,label=label_structure,lw=1)
                        lc = rgbline(kinterp, bandinter, red, green, blue, alpha=1, linestyles='solid',
                                            linewidth=(mpl.rcParams['lines.linewidth'] * 2.75))
                        ax[m].add_collection(lc)
                    else:
                        ax[m].plot(kinterp, bandinter,linestyle='--', color=color,label=label_structure,lw=1.5)
    
    
        ######################## START CALCULATION OF EFFECTIVE MASS ##########################
        #######################################################################################
    
        if emass_calculation and id_structure == 0:  ## only runs on first one, needs implementation for more at same time.
            pairs=list(zip(indexes,hspoints))  ## tuples of hspoints and their k value
            ## First we need all points contained in the plot in a huge array
            for m in range(n_bandsections):
                all_data=np.array([[0,0]])
                for i in range(len(ax[m].lines)):
                    all_data=np.concatenate((all_data,ax[m].lines[i].get_xydata()),axis=0)
                all_data=np.delete(all_data, 0, axis=0)  # delete first initialized element
                all_data=np.round(all_data,decimals=3)
                mink=min(all_data[:,0])  ## important to define edges
                maxk=max(all_data[:,0])
                
                ## function to calculate effective mass
                def get_emass(band_section,alat,k_hspoint):
                    #### PREPARE DATA WE NEED 0,0 on MAX OR MIN OF SEGMENT ####
                    ind=np.argsort(band_section[:,0])  ## get indices sorted by lowest x
                    band_section=band_section[ind] ## reorganize 
    
                    ### check type of segment 
                    hspoint_type=None  ##  will be "max" or "min"
                    signcheck=np.sign(np.sum(np.diff(band_section[:,1])))  ## does it increase or decrease as k increase
                    diffmax=k_hspoint-np.max(band_section[:,0])
                    diffmin=k_hspoint-np.max(band_section[:,0])
                    diffcheck=diffmax<diffmin  ## true if closer to max
                    if signcheck==-1 and diffcheck:
                        hspoint_type == "min"
                    elif signcheck==-1 and not diffcheck:
                        hspoint_type == "max"
                    elif signcheck==1 and diffcheck:
                        hspoint_type == "min"
                    elif signcheck==1 and not diffcheck:
                        hspoint_type == "min"
    
                    ## max or min has to become 0,0 point
                    band_section[:,0]=band_section[:,0]-k_hspoint
                    if hspoint_type == "min":
                        band_section[:,1]=band_section[:,1]-min(band_section[:,1])
                    if hspoint_type == "max":
                        band_section[:,1]=band_section[:,1]-max(band_section[:,1])
                    ###### DATA READY ##########
    
                    ### k coordinates are 2pi/alat, have to convert to 2pi/bohr
                    ### multiply 1/alat
                    c=np.polynomial.polynomial.polyfit(band_section[:,0]*(2*np.pi)/alat,band_section[:,1]/27.2114,[0,2])
                    #print("coeficients from lower to higher",c)
                    ### if units of k are bohr-1 and energy are hartree, effmass will be in units of electron mass.
                    effmass=1/(2*c[2])  ## have to multiply by 2 because is the second derivative 
                    
                    # ### if you want to check data in a plot
                    # x=np.linspace(min(band_section[:,0]*(2*np.pi)/alat),max(band_section[:,0]*(1/alat)),num=10)
                    # p=np.polynomial.polynomial.polyval(x,c)
                    # fig, ax = plt.subplots(1,1)
                    # ax.plot(x,p)
                    # ax.plot(band_section[:,0]*(2*np.pi)/alat,band_section[:,1]/27.2114)
                    return effmass
    
                def plotfit(band_section,ax,**kwargs):
                    ## plot points
                    ax.plot(band_section[:,0],band_section[:,1],'o-',color=kwargs['color'])
                    ## parabolas are wrong commented
                    # c_plot=np.polynomial.polynomial.polyfit(band_section[:,0],band_section[:,1],[0,2])
                    # x=np.linspace(0,k_width,num=100)
                    # p=np.polynomial.polynomial.polyval(x,c_plot)
                    # x_hspoint=kwargs['x_hs']
                    # ax.plot(x+x_hspoint,p,color="black")
    
    
                ref_point=[0.000, ylevel_refpoint]
                ## now we run every hspoint, check if in edges or middle and within data
                for idx_pair, hspoint in enumerate(pairs):
                    left_edge,right_edge=(False,False)
                    if hspoint[1] > maxk+0.001 or hspoint[1] < mink-0.001: ## is the point valid for the axis
                        ## print("Invalid hspoint, maybe it is in another axis. HSP: {0} MINK: {1} MAXK: {2}".format(hspoint[1],mink,maxk))
                        continue ## if it isnt go to next
                    if hspoint[1]-0.1 < mink:  ## left edge point
                        print("left edge")
                        left_edge=True
                    if hspoint[1]+0.1 > maxk:  ## right edge point
                        print("right edge")
                        right_edge=True
    
                    ref_point[0]=hspoint[1]  ## ref point k equals the high symmetry point k
                    print("hspoint:",hspoint[0])
                    ### takes a subset containing all points in the k_width range to the hspoint
                    subset=all_data[np.where(np.isclose(all_data[:,0],ref_point[0],atol=k_width))]
                    
                    
                    ### conduction band ####
                    cbsubset=subset[subset[:,1]>ref_point[1]]
                    ind=np.argsort(cbsubset[:,1])  ## get indices sorted by lowest y
                    cbsubset=cbsubset[ind] ## reorganize subset
                    lower_cbsubset=np.array([[-1,-1]])
                    for point in cbsubset:
                        ## check if point x is already in the list
                        # if any point already added has x close to the new point ## these are points in superior bands
                        if np.any(np.isclose(lower_cbsubset[:,0],point[0],atol=1e-5)): 
                            continue
                        else:
                            lower_cbsubset=np.append(lower_cbsubset,[point],axis=0)
                    lower_cbsubset=np.delete(lower_cbsubset,0,axis=0)
    
                    if not left_edge and not right_edge:  ## if its not an edge find and add median.
                        median_cbsubset=np.median(lower_cbsubset,axis=0)
                        idx = (np.abs(lower_cbsubset[:,0] - median_cbsubset[0])).argmin()
                        median_cbsubset[1]=lower_cbsubset[idx,1]  ## median point have same y as its neighbour, not y average
                        if np.any(np.isin(median_cbsubset,lower_cbsubset)): ## if median not in set add
                            lower_cbsubset=np.append(lower_cbsubset,[median_cbsubset],axis=0)
                            ind=np.argsort(lower_cbsubset[:,0])  ## get indices sorted by lowest x
                            lower_cbsubset=lower_cbsubset[ind] ## reorganize subset
                        median_idx=np.where(lower_cbsubset==median_cbsubset)[0][0]
                        lefthand_cband=lower_cbsubset[:median_idx+1]
                        righthand_cband=lower_cbsubset[median_idx:]
                    if left_edge:
                        righthand_cband=lower_cbsubset
                    if right_edge: 
                        lefthand_cband=lower_cbsubset
    
    
                ### valence band ####
                    vbsubset=subset[subset[:,1]<ref_point[1]]
                    ind=np.argsort(vbsubset[:,1])  ## get indices sorted by lowest y
                    vbsubset=vbsubset[ind] ## reorganize subset
                    vbsubset=vbsubset[::-1] ## reverse
                    higher_vbsubset=np.array([[-1,-1]])
                    for point in vbsubset:
                        ## check if point x is already in the list
                        # if any point already added has x close to the new point ## these are points in superior bands
                        if np.any(np.isclose(higher_vbsubset[:,0],point[0],atol=1e-5)): 
                            continue
                        else:
                            higher_vbsubset=np.append(higher_vbsubset,[point],axis=0)
                    higher_vbsubset=np.delete(higher_vbsubset,0,axis=0)
    
                    if not left_edge and not right_edge:  ## if its not an edge find and add median.
                        median_vbsubset=np.median(higher_vbsubset,axis=0)
                        idx = (np.abs(higher_vbsubset[:,0] - median_cbsubset[0])).argmin()
                        median_vbsubset[1]=higher_vbsubset[idx,1]  ## median point have same y as its neighbour, not y average
                        if np.any(np.isin(median_vbsubset,higher_vbsubset)): ## if median not in set add
                            higher_vbsubset=np.append(higher_vbsubset,[median_vbsubset],axis=0)
                            ind=np.argsort(higher_vbsubset[:,0])  ## get indices sorted by lowest x
                            higher_vbsubset=higher_vbsubset[ind] ## reorganize subset
                        median_idx=np.where(higher_vbsubset==median_vbsubset)[0][0]
                        lefthand_vband=higher_vbsubset[:median_idx+1]
                        righthand_vband=higher_vbsubset[median_idx:]
                    if left_edge:
                        righthand_vband=higher_vbsubset
                    if right_edge: 
                        lefthand_vband=higher_vbsubset
    
                    if left_edge:
                        if plot_emass_fitting:
                            plotfit(righthand_cband,ax[m],color="purple",k_hspoint=hspoint[1])
                            plotfit(righthand_vband,ax[m],color="purple",k_hspoint=hspoint[1])
                        print('cbandmass '+pairs[idx_pair][0]+'-'+pairs[idx_pair+1][0],get_emass(righthand_cband,alat,hspoint[1]))
                        print('vbandmass '+pairs[idx_pair][0]+'-'+pairs[idx_pair+1][0],get_emass(righthand_vband,alat,hspoint[1]))
    
    
                
                    if right_edge:
                        if plot_emass_fitting:
                            plotfit(lefthand_cband,ax[m],color="orange",x_hspoint=hspoint[1]) 
                            plotfit(lefthand_vband,ax[m],color="orange",x_hspoint=hspoint[1]) 
                        print('cbandmass '+pairs[idx_pair-1][0]+'-'+pairs[idx_pair][0],get_emass(lefthand_cband,alat,hspoint[1]))
                        print('vbandmass '+pairs[idx_pair-1][0]+'-'+pairs[idx_pair][0],get_emass(lefthand_vband,alat,hspoint[1]))
    
    
                    if not left_edge and not right_edge:  ## if its not an edge find and add median.
                        if plot_emass_fitting:
                            plotfit(lefthand_cband,ax[m],color="orange",x_hspoint=hspoint[1]) 
                            plotfit(lefthand_vband,ax[m],color="orange",x_hspoint=hspoint[1])
                            plotfit(righthand_cband,ax[m],color="purple",x_hspoint=hspoint[1])
                            plotfit(righthand_vband,ax[m],color="purple",x_hspoint=hspoint[1])
                        print('cbandmass '+pairs[idx_pair-1][0]+'-'+pairs[idx_pair][0],get_emass(lefthand_cband,alat,hspoint[1]))
                        print('vbandmass '+pairs[idx_pair-1][0]+'-'+pairs[idx_pair][0],get_emass(lefthand_vband,alat,hspoint[1]))
                        print('cbandmass '+pairs[idx_pair][0]+'-'+pairs[idx_pair+1][0],get_emass(righthand_cband,alat,hspoint[1]))
                        print('vbandmass '+pairs[idx_pair][0]+'-'+pairs[idx_pair+1][0],get_emass(righthand_vband,alat,hspoint[1]))
    
#    if len(bluestates) != 0 and len(redstates) != 0 :
#        def print_colormap():
#            import matplotlib as mpl
#            import numpy as np
#            import matplotlib.pyplot as plt
#            from matplotlib.colors import LinearSegmentedColormap
#
#            #thresh = 0.2
#            nodes = [0,1.0]
#            colors = ["red", "blue"]
#            cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
#            cmap.set_under("white")
#            fig_cb, ax_cb = plt.subplots()
#
#            data = np.zeros((10,10))-1
#            print(data)
#
#            im = ax_cb.imshow(data, cmap=cmap, vmin=0, vmax=1)
#            cbar = fig_cb.colorbar(im, ticks=[0, 1],shrink=0.25,aspect=5)
#            cbar.ax_cb.set_yticklabels([labelsprojs[0],labelsprojs[2]])
#            cbar.set_label('Element contribution', rotation=270,  labelpad = 10)
#
#            fig_cb.savefig('colormap.svg', format='svg', dpi=1000,bbox_inches='tight')
#        print_colormap()
#
#        import matplotlib
#        from matplotlib.colors import LinearSegmentedColormap
#        #thresh = 0.2
#        nodes = [0,1.0]
#        colors = ["red", "blue"]
#        cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
#        cmap.set_under("white")
#        #figc, axc = plt.subplots()
#
#        data = np.arange(1, 0, -0.01).reshape(10, 10)
#
#        fig, ax = plt.subplots()
#        cax = fig.add_axes([0.27, 0.8, 0.5, 0.05])
#
#        im = ax.imshow(data, cmap=cmap)
#
#        cbar=fig.colorbar(im, cax=cax, ticks=[0, 1],shrink=0.25,aspect=5)
#        #matplotlib.colorbar.ColorbarBase(fig, cmap=cmap)
#
#        cbar.ax.set_yticklabels([labelsprojs[0],labelsprojs[2]])
#        cbar.set_label('Element contribution', rotation=270,  labelpad = 10)

    ## this will coalesce legends into the different color ones with their corresponding color.
    final_handles = [] 
    handlescolors=[]
    final_labels = []
    handles,labels = ax[0].get_legend_handles_labels()
    for idx, handle in enumerate(handles):
        color=handle.get_color()
        if color not in handlescolors:
            final_handles.append(handle)
            handlescolors.append(color)
            final_labels.append(labels[idx])  ## add corresponding label 
    
    linelegend=plt.legend(final_handles,final_labels,bbox_to_anchor=legendposition_tuple_lines,labelspacing=0.3,
                      loc='upper left', frameon=False, prop={"size":12},
                      handletextpad=0.1)  ## place upper left corner of legend on 0.61,0.5 axes coordinates.
    # Add the legend manually to the current Axes.
    ax = plt.gca().add_artist(linelegend) ## manual insertion to make two separate legends
    
    import matplotlib.patches as mpatches
    
    patches=[]
    red_patch = mpatches.Patch(color='red', label=labelsprojs[0])
    green_patch = mpatches.Patch(color='green', label=labelsprojs[1])
    blue_patch = mpatches.Patch(color='blue', label=labelsprojs[2])
    if len(redstates) != 0 :
        patches.append(red_patch)
        plt.legend(handles=patches,bbox_to_anchor=legendposition_tuple_patches, frameon=False, prop={"size":12},)
    elif len(greenstates) != 0 :
        patches.append(green_patch)
        plt.legend(handles=patches,bbox_to_anchor=legendposition_tuple_patches, frameon=False, prop={"size":12},)
    elif len(bluestates) != 0 :
        patches.append(blue_patch)
        plt.legend(handles=patches,bbox_to_anchor=legendposition_tuple_patches, frameon=False, prop={"size":12},)


#plt.savefig("cm.png")
   
    
    
    #plt.show()
    fig.savefig(PREFIX[0]+'_bandproj'+'.png', format='png', dpi=1000,bbox_inches='tight')
    if print_svg:
        fig.savefig(PREFIX[0]+'_bandproj'+'.svg', format='svg', dpi=1000,bbox_inches='tight')


