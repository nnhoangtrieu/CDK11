import argparse
from rdkit import RDLogger
from utils import * 


parser = argparse.ArgumentParser(description='Combine molecules')

parser.add_argument('-b', '--base', type=str, help='Can be either path to a list of SMILES or a SMILES string')
parser.add_argument('-a', '--add', type=str)

args = parser.parse_args()


base = validate_arg(args.base)
add = validate_arg(args.add)


output = []


RDLogger.DisableLog('rdApp.*')
for b in base : 
    for a in add : 
        output += auto_add(b, a)


for o in output : 
    print(o)

