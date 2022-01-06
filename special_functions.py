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

def generate_xsf_file(input_file):
    import ase.io.xsf 
    structure=read_structure(input_file)
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
    df_colors = pd.read_csv("jmolcolors.csv")
    r = df_colors[df_colors['atom'] == atom]['R'].values[0] / 255.0
    g = df_colors[df_colors['atom'] == atom]['G'].values[0] / 255.0
    b = df_colors[df_colors['atom'] == atom]['B'].values[0] / 255.0
    return (r, g, b)

def plot_surface_geometry(structure,surface_direction='z'):
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
        surface_zcoords=surface_zcoords[surface_zcoords>-5]
    #    print(surface_zcoords)
        surface_zcoords_all.append(-1*surface_zcoords)
    #    print(surface_zcoords_all)
# https://www.weirdgeek.com/2018/11/plotting-stacked-histogram
    ax.hist(surface_zcoords_all,bins=50,color=colors,histtype='bar', density=False,
            stacked=True,label=different_symbols)

    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    ax.set_xlabel(r'Distance from surface ($\AA$)')
    ax.set_ylabel(r'Element counts')
    ax.set_xlim([0.0,5.1])
    ### SETTING TICKS IN Y AXIS
    majoryLocator   = MultipleLocator(1)
    majoryFormatter = FormatStrFormatter('%d')
    #minoryLocator   = MultipleLocator(0.5)
    ax.yaxis.set_major_locator(majoryLocator)
    ax.yaxis.set_major_formatter(majoryFormatter)
    #ax[m].yaxis.set_minor_locator(minoryLocator)
    majorxLocator   = MultipleLocator(0.5)
    minorxLocator   = MultipleLocator(0.25)
    majorxFormatter = FormatStrFormatter('%.2f')
    ax.xaxis.set_major_locator(majorxLocator)
    ax.xaxis.set_major_formatter(majorxFormatter)
    ax.xaxis.set_minor_locator(minorxLocator)

    plt.legend()
    import matplotlib
    matplotlib.use("Qt5Agg") ## otherwise doesnt plot after importing minflow.
    fig.savefig(prefix+'-surf_geometry.png')
    fig.savefig(prefix+'-surf_geometry.svg')
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
    text+=f" ab_Area={np.linalg.norm(np.cross(structure.cell[0],structure.cell[1]))}, \n 
    bc_Area={np.linalg.norm(np.cross(structure.cell[1],structure.cell[2]))}, \n
    ac_Area={np.linalg.norm(np.cross(structure.cell[0],structure.cell[2]))}, \n
    cell_volume={np.dot(structure.cell[0],np.cross(structure.cell[1],structure.cell[2]))} \n" 

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
