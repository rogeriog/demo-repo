import numpy as np
import itertools,math
import sys,glob,os,re

def get_nk(pathPBZ, kdensity):
    nk=[]
    def pairwise(iterable):
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b) 
    for k1, k2 in pairwise(pathPBZ):
#        print(k1,k2)
        if k1 is None or k2 is None:
            continue
        else:
            dist=np.linalg.norm(np.array(k2)-np.array(k1))
        nk.append(math.ceil(dist*kdensity))
    return tuple(nk)

from ase.io import read
from ase.io.vasp import read_vasp
from ase.io.espresso import read_espresso_in
from ase.io.espresso import read_espresso_out

def read_structure(filename):
    try:
        structure = read(filename)  ## cif
    except:
        try:
            structure = read_espresso_in(filename) ##espresso input
        except:
            try:
                with open(filename) as f:
                    r=read_espresso_out(f,index=-1,results_required=False) ##espresso output
                    for relaxed in r:
                        structure = relaxed
                ### results_required False so really takes last coordinate
            except:
                try:
                    with open(filename) as f:
                        r=read_espresso_out(f,index=-2,results_required=False)
                        print("Warning: last coordinate probably didnt converge!")
                        ##espresso output gets in case of non convergence
                        for relaxed in r:
                            structure = relaxed
                except:
                    try:
                        structure = read_vasp(filename) ## vasp input
                    except:
                        print("The structure is not readable")
                    raise AssertionError
                else:
                    return structure
            else:
                return structure
        else:
            return structure
    else:
        return structure


def substitute_atom(structure, atom,new_atom):
    """ Given an ASE structure, substitute atom by new_atom, atom and new_atom are 
      both atomic symbols, strings. """
    for i in range(len(structure)):
        if structure[i].symbol == atom:
            structure[i].symbol=new_atom
    return structure

def get_data_from_outQE(outfile,data):
    with open(outfile) as f:
        pwo_lines=f.readlines()
        initialize_magn_valence=True
        readvalence=False
        readmagn=False
        got_ntyp=False
        symbols=[]
        magn=[]
        i_magn=0
        valence=[]
        i_val=0
        for line in pwo_lines:
            if "number of Kohn-Sham states" in line:
                nbnd=int(line.split()[4])
            if "number of electrons" in line:
                nelectron=int(float(line.split()[4]))
            if "number of atomic types" in line and initialize_magn_valence:
                ntyp=int(float(line.split()[5]))
                got_ntyp=True
                magn=np.zeros(ntyp)
                valence=np.zeros(ntyp)
                initialize_magn_valence=False ## otherwise may initialize again
            if "atomic species   valence" in line:
                readvalence=True
                continue ## go to next line
            if got_ntyp and i_val == ntyp:
                readvalence=False ## stop reading new types
            if readvalence:
                valence[i_val]=float(line.split()[1])
                symbols.append(line.split()[0])
                i_val+=1
            if "atomic species   magnetization" in line:
                readmagn=True
                continue ## go to next line
            if got_ntyp and i_magn == ntyp:
                readmagn=False ## stop reading
            if readmagn:
                magn[i_magn]=float(line.split()[1])
                i_magn+=1
        if not np.all(magn) : ## if all magn is 0 no magmom is updated. 
            magmoms=np.round(valence*magn,2)
        if data == "nelectron":
            return nelectron
        if data == "ntyp":
            return ntyp
        if data == "nbnd":
            return nbnd
        if data == "magn":
            return magmoms

def get_total_energy(filename):
    """Search for the energy in file, if there are multiple gets last one"""
    string_to_search="!" ## identifier of energy
    total_energy = 0
    # Open the file in read only mode
    with open(filename, 'r') as read_obj:
        # Read all lines in the file one by one
        lines=read_obj.readlines()
        for line in lines:
            if string_to_search in line:
                #list_of_results.append(line.rstrip())  ## if youd like all results
                total_energy=line.rstrip().split()[4]
    return float(total_energy)

import MinFlow.cluster_sets as cluster_sets
def getPseudoDict(structure,cluster="local",user="user",enforce_PPsufix=None):
    """
    Pseudopotentials as previously determined by the pseudodir of the cluster/user, 
    must be in format XX.some_identifier.upf where XX is element symbol
    for this function to work, wont work well on XX_any_thing_here.upf
     enforce_PPsufix is useful if more than one type of PP in same folder, use like this:
     enforce_PPsufix="some_identifier.upf"
    """
    lowElements=list(map(str.lower,structure.get_chemical_symbols()))
    Elements=structure.get_chemical_symbols()
    pseudodir=cluster_sets.getPseudoDir(cluster=cluster,user=user)
    ppmatches=[]
    ppfile_list=[ppfile for ppfile in glob.glob(pseudodir+"/*") if not os.path.isdir(ppfile) and ppfile.split(".")[-1].lower() == "upf" ]
    if enforce_PPsufix is not None and isinstance(enforce_PPsufix,str):
        ppfile_list=[ppfile for ppfile in ppfile_list if ppfile.split("/")[-1].split(".",1)[1]==enforce_PPsufix]
        #print("PP list reduced to "+str(ppfile_list))
            ## only if matches the sufix specified
    for target_element in lowElements:
        for ppfile in ppfile_list:
            with open(ppfile) as ppf:
                for x in range(100):
                    line = ppf.readline() ## only first 100 lines
                    element_match=re.search('\s*element\s*=\s*"\s*(\w+)\s*"\s*',line)
                    if element_match is None:
                        element_match=re.search('\s*(\w+)\s*Element\s*',line)  ## format in gbrv_us
                    if element_match is not None:
                        element=element_match.groups()[0]
                        if target_element==element.lower():
                            ppmatches.append(ppfile.split("/")[-1])
                            break ## suspends search element found
            if element_match is not None:
                break  ## stops looking for match in other pp files
        else:
            print("PP not found for this element: "+target_element) 
            sys.exit(0)
    zipobj=zip(Elements, ppmatches )
    pseudodict={'potentials': dict(zipobj), 'dir':pseudodir }
    print("Check pseudo dict: \n "+str(pseudodict)+"\n")
    return pseudodict



def get_positions_symbol(structure,symbol):
    from ase.symbols import symbols2numbers
    z=symbols2numbers(symbol)
    array=[]
    for i,t in enumerate(structure.numbers==z):
        if t:
            array.append(i)
    positions=list(zip(structure.positions[array],array))
    #print(positions)
    return positions

def set_element_magnetization(structure,symbol,starting_mag=3.0):
    from ase.symbols import symbols2numbers
    z=symbols2numbers(symbol)
    mask=(np.array(structure.numbers)==z)
    magmom=structure.get_initial_magnetic_moments()
    magmom[mask]=starting_mag
    return magmom
def get_pathPBZ(structure):
    specialpoints_dict=structure.cell.bandpath().special_points
    ## we need lists not np.arrays for the coordinates
    specialpoints_dict = {k: list(v) for k, v in specialpoints_dict.items()}
    specialpoints_dict[","]=None  ## None necessary for breaks in banduppy
    path=[specialpoints_dict[x] for x in structure.cell.bandpath().path]
    labels=structure.cell.bandpath().path.replace(",","")
    return (labels,path)

def generate_xsf_file(input_file,**kwargs):
    import ase.io.xsf 
    structure=read_structure(input_file)
    xsffile=kwargs.get("filename",None)
    if xsffile is None:
        if input_file.split(".")[0] == "out":
            xsffile=input_file.split(".")[1]+".out.xsf"
        else:
            xsffile=input_file.split(".")[0]+".xsf"
    with open(xsffile, "w") as f:
        ase.io.xsf.write_xsf(f,[structure])
    print("XSF file generated.")

def get_name_from_structure(structure):
    symbols=structure.get_chemical_symbols()
    from collections import OrderedDict  ## OrderedDict is predictable set is not.
    different_symbols=list(OrderedDict.fromkeys(symbols))
    name=""
    for ele in different_symbols:
        count=symbols.count(ele)
        name+=ele+str(count)
    return name

def get_atom_color(atom):
    """
    Returns RGB color tuple corresponding to atom

    :param str atom: symbol for atom
    """
    import pandas as pd
    import MinFlow
    folder=os.path.dirname(MinFlow.__file__)
    df_colors = pd.read_csv(folder+"/jmolcolors.csv")
    r = df_colors[df_colors['atom'] == atom]['R'].values[0] / 255.0
    g = df_colors[df_colors['atom'] == atom]['G'].values[0] / 255.0
    b = df_colors[df_colors['atom'] == atom]['B'].values[0] / 255.0
    return (r, g, b)


def plot_cluster_geometry(structure,center,measurement='rxy_z',identifier="",showplot=False,
                          label="",**kwargs):
    if isinstance(structure,list):
        mode="comparative"
    else:
        mode="default"

    if mode == "default":
        prefix=get_name_from_structure(structure)
        positions=structure.get_positions()
        symbols=structure.get_chemical_symbols()

        from collections import OrderedDict  ## OrderedDict is predictable set is not.
        different_symbols=list(OrderedDict.fromkeys(symbols))
        colors=[get_atom_color(element) for element in different_symbols]
         
        
        if measurement=='rxy_z': 

            positions_rxy=positions[:,:2] ## get x,y coords
            positions_rxy=positions_rxy-np.array(center)[:2]
            positions_rxy=np.array([np.linalg.norm(r) for r in positions_rxy])
            symbol_positions_rxy=np.column_stack((symbols,positions_rxy))
            
            positions_rxy_by_symbol=[]
            for element in different_symbols:
                matches=np.where(symbol_positions_rxy[:,0]==element)
                positions_rxy_by_symbol.append(symbol_positions_rxy[matches[0]])
            positions_rxy=np.array([ pos[:,1].astype(np.float) for pos in positions_rxy_by_symbol],dtype=object)
            range_rxy=max([max(a) for a in positions_rxy])-min([min(a) for a in positions_rxy])

            positions_z=positions[:,2] ## get z coords
            positions_z=positions_z-np.array(center)[2]
            symbol_positions_z=np.column_stack((symbols,positions_z))
            
            positions_z_by_symbol=[]
            for element in different_symbols:
                matches=np.where(symbol_positions_z[:,0]==element)
                positions_z_by_symbol.append(symbol_positions_z[matches[0]])
            positions_z=np.array([ pos[:,1].astype(np.float) for pos in positions_z_by_symbol],dtype=object)
            range_z=max([max(a) for a in positions_z])-min([min(a) for a in positions_z])
#            print(get_data_range(positions_z))
#            print(get_data_range(positions_rxy))
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1,2)
            
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
            plt.rcParams['mathtext.default'] = 'regular'
            # https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
    
            ## plot RXY
            y,x,_=ax[0].hist(positions_rxy,bins=int(range_rxy*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
            ymax=np.amax(y)
            from matplotlib.ticker import MultipleLocator, FormatStrFormatter
            ax[0].set_xlabel(r'Distance from center ($\AA$)',fontsize=15)
            ax[0].set_ylabel(r'Element counts',fontsize=15)
    #        ax.set_xlim([0.0,angstrom_cutoff+0.1])
            ax[0].set_ylim([0,ymax+1])
            ### SETTING TICKS IN Y AXIS
            majoryLocator   = MultipleLocator(2)
            majoryFormatter = FormatStrFormatter('%d')
            minoryLocator   = MultipleLocator(1)
            ax[0].yaxis.set_major_locator(majoryLocator)
            ax[0].yaxis.set_major_formatter(majoryFormatter)
            ax[0].yaxis.set_minor_locator(minoryLocator)
            majorxLocator   = MultipleLocator(0.5)
            minorxLocator   = MultipleLocator(0.25)
            majorxFormatter = FormatStrFormatter('%.1f')
            ax[0].xaxis.set_major_locator(majorxLocator)
            ax[0].xaxis.set_major_formatter(majorxFormatter)
            ax[0].xaxis.set_minor_locator(minorxLocator)
            ax[0].text(0.2,0.85,"XY-plane\n"+label,transform=ax[0].transAxes, horizontalalignment='center', 
            verticalalignment='center', size=15)
    
            ## plot Z
            y,x,_=ax[1].hist(positions_z,bins=int(range_z*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
            ymax=np.amax(y)
            from matplotlib.ticker import MultipleLocator, FormatStrFormatter
            ax[1].set_xlabel(r'Distance from center ($\AA$)',fontsize=15)
            ax[1].set_ylabel(r'Element counts',fontsize=15)
    #        ax.set_xlim([0.0,angstrom_cutoff+0.1])
            ax[1].set_ylim([0,ymax+1])
            ### SETTING TICKS IN Y AXIS
            majoryLocator   = MultipleLocator(2)
            majoryFormatter = FormatStrFormatter('%d')
            minoryLocator   = MultipleLocator(1)
            ax[1].yaxis.set_major_locator(majoryLocator)
            ax[1].yaxis.set_major_formatter(majoryFormatter)
            ax[1].yaxis.set_minor_locator(minoryLocator)
            majorxLocator   = MultipleLocator(0.5)
            minorxLocator   = MultipleLocator(0.25)
            majorxFormatter = FormatStrFormatter('%.1f')
            ax[1].xaxis.set_major_locator(majorxLocator)
            ax[1].xaxis.set_major_formatter(majorxFormatter)
            ax[1].xaxis.set_minor_locator(minorxLocator)
    
    
            ax[1].text(0.2,0.85,"Z-direction\n"+label,transform=ax[1].transAxes, horizontalalignment='center', 
            verticalalignment='center', size=15)
    
            plt.legend(fontsize=15)
            import matplotlib
            matplotlib.use("Qt5Agg") ## otherwise doesnt plot after importing minflow.
            figure_height_inches=kwargs.get('height_inches',5)
            figure_width_inches=kwargs.get('width_inches',15)
            fig.set_figheight(figure_height_inches)
            fig.set_figwidth(figure_width_inches)
    
            fig.savefig(prefix+'-cluster_geometry'+identifier+'.png')
            fig.savefig(prefix+'-cluster_geometry'+identifier+'.svg')
    
            if showplot:
                plt.show() #.savefig('anything.png')
        
    if mode == "comparative":
        prefix=get_name_from_structure(structure[0])
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2,2)
        max_rxy=0
        min_rxy=0
        max_z=0
        min_z=0
        ymax=np.ones((2,2))
        for idx in range(len(structure)):
            positions=structure[idx].get_positions()
            symbols=structure[idx].get_chemical_symbols()
         
            from collections import OrderedDict  ## OrderedDict is predictable set is not.
            different_symbols=list(OrderedDict.fromkeys(symbols))
            colors=[get_atom_color(element) for element in different_symbols]
            if measurement=='rxy_z':
                positions_rxy=positions[:,:2] ## get x,y coords
                positions_rxy=positions_rxy-np.array(center[idx])[:2]
                positions_rxy=np.array([np.linalg.norm(r) for r in positions_rxy])
                symbol_positions_rxy=np.column_stack((symbols,positions_rxy))
                
                positions_rxy_by_symbol=[]
                for element in different_symbols:
                    matches=np.where(symbol_positions_rxy[:,0]==element)
                    positions_rxy_by_symbol.append(symbol_positions_rxy[matches[0]])
                positions_rxy=np.array([ pos[:,1].astype(np.float) for pos in positions_rxy_by_symbol],dtype=object)
                range_rxy=max([max(a) for a in positions_rxy])-min([min(a) for a in positions_rxy])
                max_rxy=max([max(a) for a in positions_rxy]+[max_rxy])
                min_rxy=min([min(a) for a in positions_rxy]+[min_rxy])
                
                positions_z=positions[:,2] ## get z coords
                positions_z=positions_z-np.array(center[idx])[2]
                symbol_positions_z=np.column_stack((symbols,positions_z))
                
                positions_z_by_symbol=[]
                for element in different_symbols:
                    matches=np.where(symbol_positions_z[:,0]==element)
                    positions_z_by_symbol.append(symbol_positions_z[matches[0]])
                positions_z=np.array([ pos[:,1].astype(np.float) for pos in positions_z_by_symbol],dtype=object)
                range_z=max([max(a) for a in positions_z])-min([min(a) for a in positions_z])
                max_z=max([max(a) for a in positions_z]+[max_z])
                min_z=min([min(a) for a in positions_z]+[min_z])
    #            print(get_data_range(positions_z))
    #            print(get_data_range(positions_rxy))
                import matplotlib.pyplot as plt
                
                plt.rcParams['font.family'] = 'serif'
                plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
                plt.rcParams['mathtext.default'] = 'regular'
                # https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
        
                ## plot RXY
                if idx == 0:
                    y,x,_=ax[idx][0].hist(positions_rxy,bins=int(range_rxy*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
                    ymax[idx][0]=np.amax(y)
                elif idx == 1:
                    y,x,_=ax[idx][0].hist(positions_rxy,bins=int(range_rxy*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols ,alpha=0.5,linestyle="--",edgecolor="black") 
                    ymax[idx][0]=np.amax(y)
                ## plot Z
                if idx == 0:
                    y,x,_=ax[idx][1].hist(positions_z,bins=int(range_z*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
                    ymax[idx][1]=np.amax(y)
                elif idx == 1:
                    y,x,_=ax[idx][1].hist(positions_z,bins=int(range_z*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols ,alpha=0.5,linestyle="--",edgecolor="black")
                    ymax[idx][1]=np.amax(y)


        for idx in range(len(structure)):
            if measurement=='rxy_z':
                from matplotlib.ticker import MultipleLocator, FormatStrFormatter
                if idx == 1:
                    ax[idx][0].set_xlabel(r'Distance from center ($\AA$)',fontsize=15)
                ax[idx][0].set_ylabel(r'Element counts',fontsize=15)
                ax[idx][0].set_xlim([min_rxy,max_rxy+0.1])
                ax[idx][0].set_ylim([0,max(ymax[0][0],ymax[1][0])+1])
                ### SETTING TICKS IN Y AXIS
                majoryLocator   = MultipleLocator(2)
                majoryFormatter = FormatStrFormatter('%d')
                minoryLocator   = MultipleLocator(1)
                ax[idx][0].yaxis.set_major_locator(majoryLocator)
                ax[idx][0].yaxis.set_major_formatter(majoryFormatter)
                ax[idx][0].yaxis.set_minor_locator(minoryLocator)
                majorxLocator   = MultipleLocator(0.5)
                minorxLocator   = MultipleLocator(0.25)
                majorxFormatter = FormatStrFormatter('%.1f')
                ax[idx][0].xaxis.set_major_locator(majorxLocator)
                ax[idx][0].xaxis.set_major_formatter(majorxFormatter)
                ax[idx][0].xaxis.set_minor_locator(minorxLocator)
                ax[idx][0].text(0.2,0.85,"XY-plane\n"+label[idx],transform=ax[idx][0].transAxes, horizontalalignment='center', 
                verticalalignment='center', size=15)
        
                from matplotlib.ticker import MultipleLocator, FormatStrFormatter
                if idx == 1:
                    ax[idx][1].set_xlabel(r'Distance from center ($\AA$)',fontsize=15)
                ax[idx][1].set_ylabel(r'Element counts',fontsize=15)
                ax[idx][1].set_xlim([min_z-0.1,max_z+0.1])
                ax[idx][1].set_ylim([0,max(ymax[0][1],ymax[1][1])+1])
                ### SETTING TICKS IN Y AXIS
                majoryLocator   = MultipleLocator(2)
                majoryFormatter = FormatStrFormatter('%d')
                minoryLocator   = MultipleLocator(1)
                ax[idx][1].yaxis.set_major_locator(majoryLocator)
                ax[idx][1].yaxis.set_major_formatter(majoryFormatter)
                ax[idx][1].yaxis.set_minor_locator(minoryLocator)
                majorxLocator   = MultipleLocator(0.5)
                minorxLocator   = MultipleLocator(0.25)
                majorxFormatter = FormatStrFormatter('%.1f')
                ax[idx][1].xaxis.set_major_locator(majorxLocator)
                ax[idx][1].xaxis.set_major_formatter(majorxFormatter)
                ax[idx][1].xaxis.set_minor_locator(minorxLocator)
        
                ax[idx][1].legend(fontsize=15,loc="upper left", bbox_to_anchor= (1.01, 1), 
                                   borderaxespad=0, frameon=False)
        
                ax[idx][1].text(0.2,0.85,"Z-direction\n"+label[idx],transform=ax[idx][1].transAxes, horizontalalignment='center', 
                verticalalignment='center', size=15)
        
        import matplotlib
        matplotlib.use("Qt5Agg") ## otherwise doesnt plot after importing minflow.
        figure_height_inches=kwargs.get('height_inches',5)
        figure_width_inches=kwargs.get('width_inches',15)
        fig.set_figheight(figure_height_inches)
        fig.set_figwidth(figure_width_inches)
        
        fig.savefig(prefix+'-cluster_geometry'+identifier+'_comp.png')
        fig.savefig(prefix+'-cluster_geometry'+identifier+'_comp.svg')
   
        if showplot:
            plt.show()




def plot_surface_geometry(structure,surface_direction='z',identifier="",showplot=False,
                          angstrom_cutoff=7,label=""):
    if isinstance(structure,list):
        mode="comparative"
    else:
        mode="default"

    if surface_direction == 'x':
        direction_idx=0
    elif surface_direction == 'y':
        direction_idx=1
    elif surface_direction == 'z':
        direction_idx=2

    if mode == "default":
        prefix=get_name_from_structure(structure)
        positions=structure.get_positions()
        symbols=structure.get_chemical_symbols()
        
        if surface_direction == 'x':
            direction_idx=0
        elif surface_direction == 'y':
            direction_idx=1
        elif surface_direction == 'z':
            direction_idx=2
        
        surfacemax=max(positions[:,direction_idx])
        symbol_positions=np.column_stack((symbols,positions))
            
        from collections import OrderedDict  ## OrderedDict is predictable set is not.
        different_symbols=list(OrderedDict.fromkeys(symbols))
        colors=[get_atom_color(element) for element in different_symbols]
        positions_by_symbol=[]
        for element in different_symbols:
            matches=np.where(symbol_positions[:,0]==element)
            positions_by_symbol.append(symbol_positions[matches[0]])
            
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        surface_zcoords_all=[]
        for i,positions in enumerate(positions_by_symbol):
        ## trim coordinates to last 5 angstrom
            surface_zcoords=positions[:,direction_idx+1].astype(np.float)
            surface_zcoords=surface_zcoords-surfacemax
            surface_zcoords=surface_zcoords[surface_zcoords>-angstrom_cutoff]
            surface_zcoords_all.append(-1*surface_zcoords)
            # https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
        y,x,_=ax.hist(surface_zcoords_all,bins=int(angstrom_cutoff*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
        ymax=np.amax(y)
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter
        ax.set_xlabel(r'Distance from surface ($\AA$)')
        ax.set_ylabel(r'Element counts')
        ax.set_xlim([0.0,angstrom_cutoff+0.1])
        ax.set_ylim([0,ymax+1])
        ### SETTING TICKS IN Y AXIS
        majoryLocator   = MultipleLocator(1)
        majoryFormatter = FormatStrFormatter('%d')
        #minoryLocator   = MultipleLocator(0.5)
        ax.yaxis.set_major_locator(majoryLocator)
        ax.yaxis.set_major_formatter(majoryFormatter)
        #ax[m].yaxis.set_minor_locator(minoryLocator)
        majorxLocator   = MultipleLocator(0.5)
        minorxLocator   = MultipleLocator(0.25)
        majorxFormatter = FormatStrFormatter('%.1f')
        ax.xaxis.set_major_locator(majorxLocator)
        ax.xaxis.set_major_formatter(majorxFormatter)
        ax.xaxis.set_minor_locator(minorxLocator)

        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        plt.rcParams['mathtext.default'] = 'regular'
        ax.text(0.1,0.9,label,transform=ax.transAxes, horizontalalignment='center', 
        verticalalignment='center', size=12)

        plt.legend()
        import matplotlib
        matplotlib.use("Qt5Agg") ## otherwise doesnt plot after importing minflow.
        fig.savefig(prefix+'-surf_geometry'+identifier+'.png')
        fig.savefig(prefix+'-surf_geometry'+identifier+'.svg')
        if showplot:
            plt.show() #.savefig('anything.png')
        
    if mode == "comparative":
        print("comparative")
        ymax=0
        prefix=get_name_from_structure(structure[0])
        import matplotlib.pyplot as plt
        """
        fig, ax = plt.subplots(1,1)
        for idx in reversed(range(len(structure))):
            positions=structure[idx].get_positions()
            symbols=structure[idx].get_chemical_symbols()
         
            surfacemax=max(positions[:,direction_idx])
            symbol_positions=np.column_stack((symbols,positions))
        	
            from collections import OrderedDict  ## OrderedDict is predictable set is not.
            different_symbols=list(OrderedDict.fromkeys(symbols))
            colors=[get_atom_color(element) for element in different_symbols]
        
            positions_by_symbol=[]
            for element in different_symbols:
                matches=np.where(symbol_positions[:,0]==element)
                positions_by_symbol.append(symbol_positions[matches[0]])
        
            surface_zcoords_all=[]
            for i,positions in enumerate(positions_by_symbol):
            ## trim coordinates to last 5 angstrom
                surface_zcoords=positions[:,direction_idx+1].astype(np.float)
                surface_zcoords=surface_zcoords-surfacemax
                surface_zcoords=surface_zcoords[surface_zcoords>-angstrom_cutoff]
            #    print(surface_zcoords)
                surface_zcoords_all.append(-1*surface_zcoords)
            #    print(surface_zcoords_all)
        # https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
            if idx == 0:
                ax.hist(surface_zcoords_all,bins=int(angstrom_cutoff*10),color=colors,histtype='bar', density=False,
                    stacked=True,label=different_symbols)
            elif idx==1: ## ref colored more lightly
#                colors=[np.array(color)+(np.array([256,256,256])-np.array(color)*255)/(4*255) for color in colors]
                ax.hist(surface_zcoords_all,bins=50,color=colors,histtype='bar', density=False,
                    stacked=True, label=different_symbols,alpha=0.5,linestyle="--",edgecolor="black")
"""
        fig, (ax1,ax2) = plt.subplots(2,1)
        for idx in reversed(range(len(structure))):
            positions=structure[idx].get_positions()
            symbols=structure[idx].get_chemical_symbols()
         
            surfacemax=max(positions[:,direction_idx])
            symbol_positions=np.column_stack((symbols,positions))
        	
            from collections import OrderedDict  ## OrderedDict is predictable set is not.
            different_symbols=list(OrderedDict.fromkeys(symbols))
            colors=[get_atom_color(element) for element in different_symbols]
        
            positions_by_symbol=[]
            for element in different_symbols:
                matches=np.where(symbol_positions[:,0]==element)
                positions_by_symbol.append(symbol_positions[matches[0]])
        
            surface_zcoords_all=[]
            for i,positions in enumerate(positions_by_symbol):
            ## trim coordinates to last 5 angstrom
                surface_zcoords=positions[:,direction_idx+1].astype(np.float)
                surface_zcoords=surface_zcoords-surfacemax
                surface_zcoords=surface_zcoords[surface_zcoords>-angstrom_cutoff]
            #    print(surface_zcoords)
                surface_zcoords_all.append(-1*surface_zcoords)
            #    print(surface_zcoords_all)
        # https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
            if idx==0:
                y,x,_ = ax1.hist(surface_zcoords_all,bins=int(angstrom_cutoff*10),color=colors,histtype='bar', density=False, stacked=True,label=different_symbols)
                ymax=max(ymax,np.amax(y))
            elif idx==1: ## ref colored more lightly
#                colors=[np.array(color)+(np.array([256,256,256])-np.array(color)*255)/(4*255) for color in colors]
                y,x,_ = ax2.hist(surface_zcoords_all,bins=int(angstrom_cutoff*10),color=colors,histtype='bar', density=False, stacked=True, label=different_symbols,alpha=0.5,linestyle="--",edgecolor="black")
                ymax=max(ymax,np.amax(y))
#        surface_zcoords_all=[]
#        for i,positions in enumerate(positions_by_symbol):
        ## trim coordinates to last 5 angstrom
#            surface_zcoords=positions[:,direction_idx+1].astype(np.float)
#            surface_zcoords=surface_zcoords-surfacemax
#            surface_zcoords=surface_zcoords[surface_zcoords>-angstrom_cutoff]
#            surface_zcoords_all.append(-1*surface_zcoords)
            # https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
#            ax.hist(surface_zcoords_all,bins=int(angstrom_cutoff*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter
        plt.xlabel(r'Distance from surface ($\AA$)')
        ax1.set_ylabel(r'Element counts')
        ax2.set_ylabel(r'Element counts')
        ax1.set_xlim([0.0,angstrom_cutoff+0.1])
        ax2.set_xlim([0.0,angstrom_cutoff+0.1])
        ax1.set_ylim([0,ymax+1])
        ax2.set_ylim([0,ymax+1])
        ### SETTING TICKS IN Y AXIS
        majoryLocator   = MultipleLocator(1)
        majoryFormatter = FormatStrFormatter('%d')
        #minoryLocator   = MultipleLocator(0.5)
        ax1.yaxis.set_major_locator(majoryLocator)
        ax1.yaxis.set_major_formatter(majoryFormatter)
        ax2.yaxis.set_major_locator(majoryLocator)
        ax2.yaxis.set_major_formatter(majoryFormatter)
        #ax[m].yaxis.set_minor_locator(minoryLocator)
        majorxLocator   = MultipleLocator(0.5)
        minorxLocator   = MultipleLocator(0.25)
        majorxFormatter = FormatStrFormatter('%.1f')
        ax1.xaxis.set_major_locator(majorxLocator)
        ax1.xaxis.set_major_formatter(majorxFormatter)
        ax1.xaxis.set_minor_locator(minorxLocator)
        ax2.xaxis.set_major_locator(majorxLocator)
        ax2.xaxis.set_major_formatter(majorxFormatter)
        ax2.xaxis.set_minor_locator(minorxLocator)
        ax1.legend()
        ax2.legend()

        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        plt.rcParams['mathtext.default'] = 'regular'
        ax1.text(0.1,0.9,label[0],transform=ax1.transAxes, horizontalalignment='center', 
        verticalalignment='center', size=12)
        ax2.text(0.1,0.9,label[1],transform=ax2.transAxes, horizontalalignment='center', 
        verticalalignment='center', size=12)

        import matplotlib
        matplotlib.use("Qt5Agg") ## otherwise doesnt plot after importing minflow.
        fig.savefig(prefix+'-surf_geometry'+identifier+'.png')
        fig.savefig(prefix+'-surf_geometry'+identifier+'.svg')
        if showplot:
            plt.show() #.savefig('anything.png')

def plot_badercharges(structure,prefix_baderanalysis="",
                      mode="default",direction='z',identifier="",showplot=False,
                          label="",**kwargs):

    import pandas as pd
    df=pd.read_csv("BaderChgAnalysis_"+prefix_baderanalysis+".out",delimiter=r"\s+")
    badercharges=df.iloc[:,4].to_numpy()

    if direction == 'x':
        direction_idx=0
    elif direction == 'y':
        direction_idx=1
    elif direction == 'z':
        direction_idx=2

    prefix=get_name_from_structure(structure)
    if kwargs.get('prefix',False):
        prefix=kwargs['prefix']
    positions=structure.get_positions()
    symbols=structure.get_chemical_symbols()


    surfacemax=max(positions[:,direction_idx])
    symb_pos_bader=np.column_stack((symbols,positions,badercharges))
    print(symb_pos_bader)

    from collections import OrderedDict  ## OrderedDict is predictable set is not.
    different_symbols=list(OrderedDict.fromkeys(symbols))
    colors=[get_atom_color(element) for element in different_symbols]
    data_by_symbol=[]
    for element in different_symbols:
        matches=np.where(symb_pos_bader[:,0]==element)
        data_by_symbol.append(symb_pos_bader[matches[0]])
#    print(data_by_symbol)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,1)
    for i in range(len(data_by_symbol)):
        pos=data_by_symbol[i][:,direction_idx+1].astype(float)
        bader=data_by_symbol[i][:,4].astype(float)
        p = pos.argsort()
        pos=pos[p]
        bader=bader[p]
        ax.plot(pos,bader,color=colors[i],linestyle='-',label=different_symbols[i])
        ax.plot(pos,bader,'o',color=colors[i])

    if mode=="potential_interface":
        print("A avg.dat file is necessary for potential.")
        filename="avg.dat"
        data=np.genfromtxt(filename)
        distance=data[:,0]
        potential=data[:,1]
        mac_avg_qe=data[:,2]
        ax.plot(data[:,0]*0.52918,data[:,1]*13.6056622,color="dimgray",alpha=0.5,label=r"$\overline{V}$") ## potential Angstrom X eV
        ax.plot(distance*0.52918,mac_avg_qe*13.6056622,color="black",label=r"$\overline{\overline{V}}$")  ## macavg from QE ## should be identical to the calculated one

        ax.axhline(y=max(mac_avg_qe*13.6056622),linestyle='--',color='dimgray',alpha=0.7)
        ax.axhline(y=min(mac_avg_qe*13.6056622),linestyle='--',color='dimgray',alpha=0.7)


    leg=ax.legend(bbox_to_anchor=(1,0.9),labelspacing=0.3,
                loc='upper left', frameon=False, prop={"size":13})
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)

    if kwargs.get('reverse',False):
        print('reversed x axis')
        ax.invert_xaxis()
    if kwargs.get('xlim',False):
        ax.set_xlim(kwargs['xlim'])
    if kwargs.get('ylim',False):
        ax.set_ylim(kwargs['ylim'])
#    ax.set_ylim([0,ymax+1])
    #(surface_zcoords_all,bins=int(angstrom_cutoff*10),color=colors,histtype='bar', density=False,stacked=True,label=different_symbols)
#    ymax=np.amax(y)
#    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    ax.set_xlabel(direction+r' ($\AA$)')
    ax.set_ylabel(r'Atomic charge ($e$)')

    ### SETTING TICKS IN Y AXIS
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    majoryLocator   = MultipleLocator(0.5)
    majoryFormatter = FormatStrFormatter('%.1f')
    minoryLocator   = MultipleLocator(0.1)
    ax.yaxis.set_major_locator(majoryLocator)
    ax.yaxis.set_major_formatter(majoryFormatter)
    ax.yaxis.set_minor_locator(minoryLocator)

    majorxLocator   = MultipleLocator(5)
    minorxLocator   = MultipleLocator(1)
    majorxFormatter = FormatStrFormatter('%.1f')
    ax.xaxis.set_major_locator(majorxLocator)
    ax.xaxis.set_major_formatter(majorxFormatter)
    ax.xaxis.set_minor_locator(minorxLocator)

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['mathtext.default'] = 'regular'
#    ax.text(0.1,0.9,label,transform=ax.transAxes, horizontalalignment='center',
#    verticalalignment='center', size=12)

    import matplotlib
    matplotlib.use("Qt5Agg") ## otherwise doesnt plot after importing minflow.
    fig.savefig(prefix+'-bader_plot'+identifier+'.png',bbox_inches='tight', dpi=1000)
    fig.savefig(prefix+'-bader_plot'+identifier+'.svg',bbox_inches='tight', dpi=1000)
    if showplot:
        plt.show() #.savefig('anything.png')





def get_geometry_data(structure):
    import itertools
    from ase.geometry.analysis import Analysis
    symbols_list=np.unique(np.array([str(structure[i].symbol) for i in range(len(structure))]))
    elements_combinations = [p for p in itertools.combinations(symbols_list,2)]
    elements_combinations_angles = [p for p in itertools.product(symbols_list,symbols_list,symbols_list)]
    ana = Analysis(structure)
    text="-------------------------------------------------------------------------------------\n"
    text+="-------------------------------- GEOMETRY DATA --------------------------------------\n"
    text+="-------------------------------------------------------------------------------------\n"
    text+=f"\tCELL PARAMETERS:\n"
    text+=f" v1={structure.cell[0]},\n v2={structure.cell[1]},\n v3={structure.cell[2]}\n\n"
    text+=f" a={np.round(structure.cell.lengths()[0],4)},\n b={np.round(structure.cell.lengths()[1],4)},\n \
    c={np.round(structure.cell.lengths()[2],4)}\n\n"
    text+=f" alpha={np.round(structure.cell.angles()[0],4)}, \n beta={np.round(structure.cell.angles()[1],4)},\n\
    gamma={np.round(structure.cell.angles()[2],4)}\n"
    text+=f"""ab_Area={np.linalg.norm(np.cross(structure.cell[0],structure.cell[1]))},
    bc_Area={np.linalg.norm(np.cross(structure.cell[1],structure.cell[2]))},
    ac_Area={np.linalg.norm(np.cross(structure.cell[0],structure.cell[2]))},
    cell_volume={np.dot(structure.cell[0],np.cross(structure.cell[1],structure.cell[2]))}"""

    text+=f"\n\tNUMBER AND AVERAGE LENGTH OF BONDS:\n"
    text+=f"--------------------------------------------------------\n"
    for combination in elements_combinations:
        symbol1,symbol2=combination
        Bonds = ana.get_bonds(symbol1, symbol2, unique=True)
        text+="There are {} {}-{} bonds in the structure.\n".format(len(Bonds[0]),symbol1,symbol2) 
        if np.shape(Bonds)[1] == 0:
            continue
        BondValues = ana.get_values(Bonds)
        text+="The average {}-{} bond length is {}.\n".format(symbol1,symbol2,np.average(BondValues))
        text+="------------------------------------\n"
    
    text+="\n\tLENGTH OF EVERY BOND: \n"
    text+=f"--------------------------------------------------------\n"
    text+=" !!! Element indexes as shown in VESTA. Open structure in VESTA to check. !!! \n".upper()
    for i in range(len(structure)):
        text+=f"Atom: {structure[i].symbol}, Index:{structure[i].index},\n"
        bonded_atoms=ana.all_bonds[0][i]
        bond_lists=list(zip(list(map(lambda x: x*i,[1]*len(bonded_atoms))),bonded_atoms))
        bond_lengths=ana.get_values([bond_lists])
        for idx,bond in enumerate(bond_lists):
            idx1_bond,idx2_bond=bond
            text+=f"{structure[idx1_bond].symbol}{structure[idx1_bond].index+1}-{structure[idx2_bond].symbol}{structure[idx2_bond].index+1}: {np.round(bond_lengths[0][idx],3)} angstrom\n"
    
    text+=f"\n\tANGLES BETWEEN TRIPLE OF ELEMENTS:\n"
    text+=f"--------------------------------------------------------\n"
    for combination in elements_combinations_angles:
        symbol1,symbol2,symbol3=combination
        Angles = ana.get_angles(*combination, unique=True)
        if np.shape(Angles)[1] == 0:
            continue
        text+="There are {} {}-{}-{} angles in the structure.\n".format(len(Angles[0]),symbol1,symbol2,symbol3)
        AngleValues = ana.get_values(Angles)
        text+="The average {}-{}-{} angle is {}.\n".format(symbol1,symbol2,symbol3,np.average(AngleValues))
        for idx,angle in enumerate(Angles[0]):
            ele1,ele2,ele3=angle
            text+=f"{structure[ele1].symbol}{structure[ele1].index+1}-{structure[ele2].symbol}{structure[ele2].index+1}-{structure[ele3].symbol}{structure[ele3].index+1}: {np.round(AngleValues[0][idx],3)}\n"
        text+=f"--------------------------------------------------------\n" 
    
    print(text)
    with open("Geometry_DATA.txt","w") as f:
        f.write(text)
