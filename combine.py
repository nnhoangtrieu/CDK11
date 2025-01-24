import os 
import argparse
import datetime
from tqdm import tqdm 
from rdkit import RDLogger
from utils import * 

# purine : C1=C2C(=NC=N1)N=CN2
# test block 1 : CNC[C@@H](C)O

script_dir = os.path.dirname(os.path.realpath(__file__))
time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
parser = argparse.ArgumentParser(description='Combine molecules')
parser.add_argument('-b', '--base', type=str, help='Can be either path to a list of SMILES or a SMILES string')
parser.add_argument('-a', '--add', type=str, help='Can be either path to a list of SMILES or a SMILES string')
parser.add_argument('-o', '--output', type=str, help='Output file path (Either point to a folder or .txt file)', default=None)
args = parser.parse_args()



base = extract_smiles(args.base)
add = extract_smiles(args.add)



RDLogger.DisableLog('rdApp.*')
for i, b in enumerate(tqdm(base, desc='Combining...')) : 
    for a in tqdm(add, desc=f'Combine base molecule {i+1}') : 
        combinable_mol = auto_add(b, a)

    if args.output : 
        if (args.output).endswith('.txt') : save(combinable_mol, args.output, mode='a') 
        else : 
            if not os.path.exists(args.output) : os.makedirs(args.output, exist_ok=True)
            save(combinable_mol, os.path.join(args.output, f'{time}.txt'))
    else : 
        save(combinable_mol, os.path.join(script_dir, 'output', f'{time}.txt'), mode='a')


