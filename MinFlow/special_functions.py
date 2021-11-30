import MinFlow.cluster_sets as cluster_sets
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
                for relaxed in read_espresso_out(filename,results_required=False): ## espresso output
                    structure = relaxed
                ### results_required False so really takes last coordinate
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


