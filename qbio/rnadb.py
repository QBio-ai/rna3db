import json
import os
from Bio.PDB import MMCIFParser, PDBIO
import glob

trainset = glob.glob("/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/train_set/*/*/*.cif")
print("Total CIF files in train_set:", len(trainset))
testset = glob.glob("/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/test_set/*/*/*.cif")
print("Total CIF files in test_set:", len(testset))

traindir = "/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/train_set/SSfiles"
os.makedirs(traindir, exist_ok=True)
testdir = "/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/test_set/SSfiles"
os.makedirs(testdir, exist_ok=True)

with open(f"{traindir}/processed_train_files.txt", "w") as f:
    f.write("# Non processed files\n")

for cif_file in trainset:
    filename = os.path.basename(cif_file).replace(".cif", "")
    pdb_file = f"{traindir}/output.pdb"
    try:
        parser = MMCIFParser()
        structure = parser.get_structure("structure", cif_file)
        parser = MMCIFParser()
        structure = parser.get_structure("structure", cif_file)
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname in ["U", "URA"]:
                        residue.resname = "  T" if resname == "U" else "THY"
        io = PDBIO()
        for model in structure:
            for chain in model:
                if len(chain.id) > 1:
                    chain.id = chain.id[0]
        io.set_structure(structure)
        io.save(pdb_file)
        cmd = f"/root/qbio/RNA-BRiQ/build/bin/BRiQ_AssignSS {pdb_file} {traindir}/ss_{filename}.txt"
        os.system(cmd)
        print(f"Secondary structure assignment completed for {filename}")
        with open(f"{traindir}/ss_{filename}.txt", "r") as ss_file:
            lines = ss_file.readlines()
        lines = [line for line in lines if not line.lstrip().startswith("break")]
        with open(f"{traindir}/ss_{filename}.txt", "w") as ss_file:
            ss_file.writelines(lines)
    except Exception as e:
        with open(f"{traindir}/processed_train_files.txt", "a") as f:
            f.write(f"{filename}\n")
            f.write(f"# Error during processing: {e}\n")
        
with open(f"{testdir}/processed_test_files.txt", "w") as f:
    f.write("# Non processed files\n")

for cif_file in testset:
    filename = os.path.basename(cif_file).replace(".cif", "")
    pdb_file = f"{testdir}/output.pdb"
    try:
        parser = MMCIFParser()
        structure = parser.get_structure("structure", cif_file)
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname in ["U", "URA"]:
                        residue.resname = "  T" if resname == "U" else "THY"
        io = PDBIO()
        for model in structure:
            for chain in model:
                if len(chain.id) > 1:
                    chain.id = chain.id[0]
        io.set_structure(structure)
        io.save(pdb_file)
        cmd = f"/root/qbio/RNA-BRiQ/build/bin/BRiQ_AssignSS {pdb_file} {testdir}/ss_{filename}.txt"
        os.system(cmd)
        print(f"Secondary structure assignment completed for {filename}")
        with open(f"{testdir}/ss_{filename}.txt", "r") as ss_file:
            lines = ss_file.readlines()
        lines = [line for line in lines if not line.lstrip().startswith("break")]
        with open(f"{testdir}/ss_{filename}.txt", "w") as ss_file:
            ss_file.writelines(lines)
    except Exception as e:
        with open(f"{testdir}/processed_test_files.txt", "a") as f:
            f.write(f"{filename}\n")
            f.write(f"# Error during processing: {e}\n")

output_test = f"{testdir}/output.pdb"
output_train = f"{traindir}/output.pdb"
os.remove(output_test)
os.remove(output_train)

import hashlib

def file_hash(filepath):
    hasher = hashlib.md5()
    with open(filepath, 'rb') as f:
        buf = f.read()
        hasher.update(buf)
    return hasher.hexdigest()

temp_dup_dir = os.path.join(traindir, "duplicates")
os.makedirs(temp_dup_dir, exist_ok=True)

seen_hashes = {}
for fname in os.listdir(traindir):
    fpath = os.path.join(traindir, fname)
    if os.path.isfile(fpath):
        h = file_hash(fpath)
        if h in seen_hashes:
            print(f"Duplicate found: {fname} (same as {seen_hashes[h]}), moving to duplicates folder.")
            os.rename(fpath, os.path.join(temp_dup_dir, fname))
        else:
            seen_hashes[h] = fname

temp_dup_dir = os.path.join(testdir, "duplicates")
os.makedirs(temp_dup_dir, exist_ok=True)

seen_hashes = {}
for fname in os.listdir(testdir):
    fpath = os.path.join(testdir, fname)
    if os.path.isfile(fpath):
        h = file_hash(fpath)
        if h in seen_hashes:
            print(f"Duplicate found: {fname} (same as {seen_hashes[h]}), moving to duplicates folder.")
            os.rename(fpath, os.path.join(temp_dup_dir, fname))
        else:
            seen_hashes[h] = fname
