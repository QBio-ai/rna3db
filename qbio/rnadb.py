import json
import os
from Bio.PDB import MMCIFParser, PDBIO
import glob

'''
tempdir = "/root/qbio/RNA-BRiQ/temp"
if not os.path.exists(tempdir):
    os.makedirs(tempdir)

jsonfile = "/root/qbio/RNA-BRiQ/RNAdataset/split.json"
with open(jsonfile, "r") as f:
    split_data = json.load(f)
    print("keys in split_data:", split_data.keys())
train_ids = split_data["train_set"]
valid_ids = split_data["valid_set"]
test_ids = split_data["test_set"]

print(f"Train IDs:", len(train_ids))
print(f"Valid IDs:", len(valid_ids))
print(f"Test IDs:", len(test_ids))

total_items = 0
for component in train_ids.values():
    total_items += len(component)
print("Train Total items inside all components:", total_items)

total_items2 = 0
for component in test_ids.values():
    total_items2 += len(component)
print("Test Total items inside all components:", total_items2)



def count_sequences(dataset):
    pdbchain = 0
    seq = 0
    for component in dataset.values():
        for group in component.values():
            pdbchain += 1
            for item in group.values():
                if isinstance(item, dict) and "sequence" in item:
                    seq += 1
    return seq, pdbchain

train_sequences, train_pdb_chain = count_sequences(train_ids)
test_sequences, test_pdb_chain = count_sequences(test_ids)
print("Number of PDB chains in train_set:", train_pdb_chain)
print("Number of PDB chains in test_set:", test_pdb_chain)
print("Number of sequences in train_set:", train_sequences)
print("Number of sequences in test_set:", test_sequences)
'''
'''
for component_key, component in train_ids.items():
    for group_key, group in component.items():
        for item_key, item in group.items():
            if isinstance(item, dict) and "sequence" in item:
                pathtofile = f"/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/train_set/{component_key}/{group_key}/{item_key}.cif"
                print("Sequence from train_set:", item["sequence"])
                if pathtofile:
                    pdb_file = f"{tempdir}/output.pdb"
                    parser = MMCIFParser()
                    structure = parser.get_structure("structure", pathtofile)
                    print(f"Processing structure for {item_key}")
                    # Convert U/URA residues to T/THY
                    for model in structure:
                        for chain in model:
                            for residue in chain:
                                resname = residue.get_resname().strip()
                                if resname in ["U", "URA"]:
                                    residue.resname = "  T" if resname == "U" else "THY"
                    io = PDBIO()
                    print("Saving PDB file to:", pdb_file)
                    # Ensure chain IDs are a single character
                    for model in structure:
                        for chain in model:
                            if len(chain.id) > 1:
                                print(f"Chain ID '{chain.id}' is too long, truncating to '{chain.id[0]}'")
                                chain.id = chain.id[0]
                    io.set_structure(structure)
                    io.save(pdb_file)
                    print(f"PDB file saved for {item_key} at {pdb_file}")
                    # Call external tool for secondary structure assignment
                    cmd = f"/root/qbio/RNA-BRiQ/build/bin/BRiQ_AssignSS {pdb_file} {tempdir}/ss_{item_key}.txt"
                    os.system(cmd)
                    print(f"Secondary structure assignment completed for {item_key}")
                    with open(f"{tempdir}/ss_{item_key}.txt", "r") as ss_file:
                        ss_data = ss_file.read()
                        print("Secondary structure assignment:\n", ss_data)
                else:
                    raise FileNotFoundError(f"CIF file not found for {item_key}")
                
'''

trainset = glob.glob("/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/train_set/*/*/*.cif")
print("Total CIF files in train_set:", len(trainset))
testset = glob.glob("/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/test_set/*/*/*.cif")
print("Total CIF files in test_set:", len(testset))

traindir = "/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/train_set/SSfiles"
os.makedirs(traindir, exist_ok=True)
testdir = "/root/qbio/RNA-BRiQ/RNAdataset/rna3db-mmcifs/test_set/SSfiles"
os.makedirs(testdir, exist_ok=True)

trainprocessed = glob.glob(f"{traindir}/ss_*.txt")
trainfiles = set(os.path.basename(f).replace("ss_", "").replace(".txt", "") for f in trainprocessed)
with open(f"{traindir}/processed_train_files.txt", "w") as f:
    f.write("# Non processed files\n")

for cif_file in trainset:
    filename = os.path.basename(cif_file).replace(".cif", "")
    if filename in trainfiles:
        print(f"Skipping already processed file: {filename}")
        continue
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