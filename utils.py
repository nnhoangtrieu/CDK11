import copy 
import rdkit 
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


def draw(x) : 
    if type(x) == str : x = mol(x)
    Draw.MolToImage(x).show()



def auto_add(x, y) : 
    if type(x) == str : x = mol(x)
    if type(y) == str : y = mol(y)
    
    combo = CombineMols(x, y) 
    output = []

    for i in range(x.GetNumAtoms()) :
        for j in range(x.GetNumAtoms(), combo.GetNumAtoms()) :
            for b in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE] :
                combo_editable = EditableMol(combo)
                combo_editable.AddBond(i, j, order=b)

                try : 
                    SanitizeMol(combo_editable.GetMol())
                    output.append(smiles(combo_editable.GetMol()))
                except :
                    print('Failed')
    return output


