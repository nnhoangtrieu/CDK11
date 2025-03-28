import copy 
import rdkit 
from tqdm import tqdm 
from rdkit import Chem
from rdkit.Chem import Draw, SanitizeMol, EditableMol
from rdkit.Chem.rdmolops import CombineMols


def idx_annotate(x):
    """
    Taken from: https://iwatobipen.wordpress.com/2017/02/25/draw-molecule-with-atom-index-in-rdkit/
    
    Add label to each molecule's atom: "atom_name:atom_index"
    """
    mol = copy.deepcopy(x)
    for idx in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(idx).SetProp(
            'molAtomMapNumber', 
            str(mol.GetAtomWithIdx(idx).GetIdx())
        )
    return mol

def mol(x) : 
    return Chem.MolFromSmiles(x)

def smiles(x) : 
    return Chem.MolToSmiles(x)

def draw(x):
    if isinstance(x, str):
        x = Chem.MolFromSmiles(x)
    img = Draw.MolToImage(x, size=(800, 800))
    img.show() 


def select_atom_to_add(x) : 
    if type(x) == str : mol_x = mol(x)
    draw(idx_annotate(mol_x))
    add_atom = int(input(f"Select atom in {x} to add to: "))
    return add_atom

def auto_add(x, y, selected_atom=None) : 
    if type(x) == str : x = mol(x)
    if type(y) == str : y = mol(y)
    
    combo = CombineMols(x, y) 
    output = []

    for i in range(x.GetNumAtoms()) :
        if selected_atom is not None: i = selected_atom
        for j in range(x.GetNumAtoms(), combo.GetNumAtoms()) :
            for b in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE] :
                combo_editable = Chem.EditableMol(combo)
                combo_editable.AddBond(i, j, order=b)

                try : 
                    Chem.SanitizeMol(combo_editable.GetMol())
                    output.append(Chem.MolToSmiles(combo_editable.GetMol()))
                except : pass
        if selected_atom is not None: break
    return output



def read_smi(path, delimiter='\t', titleLine=False) : 
    result = [] 
    if path.endswith('.txt') : 
        with open(path, 'r') as f : 
            for smi in tqdm(f.readlines(), desc='Reading SMILES') : 
                if Chem.MolFromSmiles(smi) is not None : 
                    result.append(smi.strip())
    elif path.endswith('.sdf') : 
        supplier = Chem.SDMolSupplier(path)
        for mol in tqdm(supplier, desc='Reading SMILES') : 
            if mol is None : 
                continue 
            result.append(Chem.MolToSmiles(mol))
    elif path.endswith('.smi') : 
        supplier = Chem.SmilesMolSupplier(path, delimiter=delimiter, titleLine=titleLine)
        for mol in tqdm(supplier, desc='Reading SMILES') : 
            if mol is None : 
                continue 
            result.append(Chem.MolToSmiles(mol))
    return result

def extract_smiles(x) :
    smi_list = read_smi(x)
    if not smi_list : 
        smi = x 
        if not mol(smi) : print(f'{smi} is not a valid SMILES string'); exit() 
        else : return [smi] 
    return smi_list

def save(data, path, mode='w') : 
    if path.endswith('.txt') : 
        with open(path, mode) as f : 
            for line in data : 
                f.write(line+'\n')


