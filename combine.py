import os 
import argparse
import datetime
from tqdm import tqdm 
from rdkit import RDLogger
from utils import * 
RDLogger.DisableLog('rdApp.*')


script_dir = os.path.dirname(os.path.realpath(__file__))
time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

parser = argparse.ArgumentParser(description='Combine molecules')
parser.add_argument('-b', '--base', type=str, help='Can be either path to a list of SMILES or a SMILES string', default='Cn1cnc2c(N[C@H]3CCCNC3)ncnc12')
parser.add_argument('-a', '--add', type=str, help='Can be either path to a list of SMILES or a SMILES string', default='CNC[C@@H](C)O')
parser.add_argument('-ms', '--manual_select', type=bool, help='Manually select which atom to combine to', default=False)
parser.add_argument('-o', '--output', type=str, help='Output file path (Either point to a folder or .txt file)', default=None)
args = parser.parse_args()

base = extract_smiles(args.base)
add = extract_smiles(args.add)
selected_atom = [] 

if len(base) > 1 and args.manual_select :
    print(f'You have {len(base)} base molecules in your list')
    print(f'As you have selected to mannually select the atom to combine, you will be prompted to select the atom for every base molecules in your list')
    agree = input('Do you want to continue? (y/n) ')
    if agree.lower() != 'y' : exit()   
    else : 
        for x in base : 
            selected_atom.append(select_atom_to_add(x))

if len(base) == 1 and args.manual_select :
    selected_atom.append(select_atom_to_add(base[0]))


for i, b in enumerate(tqdm(base, desc='Combining...')) : 
    for a in tqdm(add, desc=f'Combine base molecule {i+1}') : 
        combinable_mol = auto_add(b, a, selected_atom[i] if args.manual_select else None)

    if args.output : 
        if (args.output).endswith('.txt') : save(combinable_mol, args.output, mode='a') 
        else : 
            if not os.path.exists(args.output) : os.makedirs(args.output, exist_ok=True)
            save(combinable_mol, os.path.join(args.output, f'{time}.txt'))
    else : 
        save(combinable_mol, os.path.join(script_dir, 'output', f'{time}.txt'), mode='a')