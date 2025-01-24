# CDK 11

A Python script to generate combinatorial molecular structures by connecting two molecules via single, double, or triple bonds between all possible atom pairs.

---

## Features
- Combine molecules from **SMILES strings** or files (`.txt`, `.smi`, `.sdf`).
- Generates all valid combinations with **single/double/triple bonds** between atom pairs.
- Outputs results to a timestamped file or custom directory.
- Automatically filters invalid structures (sanitization errors).

---

## Installation

### Prerequisites
- Python 3.6+
- Required packages: `rdkit` and `tqdm`.

```bash
pip install rdkit tqdm
```
--- 

## Usage 

### Basic Command 

```bash
python combine.py -b <base_molecule> -a <add_molecule> [-o OUTPUT_PATH]
```

#### Arguments: 

* `-b/--base`: Path to a SMILES file `(.txt, .smi, .sdf)` or a single SMILES string.
* `-a/--add`: Path to a SMILES file or a single SMILES string.
* `-o/--output`: (optional): Output file `(.txt)` or folder. Default: `./output/<timestamp>.txt`

#### Examples: 

1. Combine two SMILES strings: 

```bash
python combine.py -b "C1=C2C(=NC=N1)N=CN2" -a "CNC[C@@H](C)O" 
```

2. Use files as input: 

```bash
python combine.py -b base_molecules.sdf -a add_molecules.txt 
```


#### Output Behavior: 
* If `-o` points to a `.txt` file: Appends results to the file.
* If `-o` points to a folder: Creates a timestamped `.txt` file in the folder.
* If `-o` is omitted: Saves to `./output/<timestamp>.txt`.

---