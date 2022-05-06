import os

with open("list_pdb.txt", "r") as filin:
    list_names = []
    for ligne in filin:
        list_names.append(ligne.strip())

#print(list_names)

os.chdir("/Users/MAEL/Documents/M1_BI/Stage/research/new_models/PDB")

dico_AA = {
    "ALA" : "A",
    "ARG" : "R",
    "ASN" : "N",
    "ASP" : "D",
    "CYS" : "C",
    "GLU" : "E",
    "GLN" : "Q",
    "GLY" : "G",
    "HIS" : "H",
    "ILE" : "I",
    "LEU" : "L",
    "LYS" : "K",
    "MET" : "M",
    "PHE" : "F",
    "PRO" : "P",
    "SER" : "S",
    "THR" : "T",
    "TRP" : "W",
    "TYR" : "Y",
    "VAL" : "V"
}

fileTest = "../../recherche_perso_cys/comparaison_sequences/synthese.txt"

try:
    os.remove(fileTest)
except OSError as e:
    print(e)
else:
    print("File is deleted successfully")

for name in list_names:
    
    with open(name, "r") as pdb_file, open("../../recherche_perso_cys/comparaison_sequences/dossier_seq/"+name+".txt", "w") as file_out:
        seq_AA = []
        ligne1 = pdb_file.readline()
        num_ligne1 = ligne1[22:26].strip()
        aa_eq = dico_AA[ligne1[17:20].strip()]
        seq_AA.append(aa_eq)
        for line in pdb_file:
            if (line[22:26].strip() != num_ligne1):
                num_ligne1 = line[22:26].strip()
                aa_eq = dico_AA[line[17:20].strip()]
                seq_AA.append(aa_eq)
                #print(line[17:20].strip())
        seq_AA = "".join(seq_AA)
        file_out.write(">"+name+"\n")
        file_out.write(seq_AA)
        with open("../../recherche_perso_cys/comparaison_sequences/synthese.txt", "a") as file_comp:
            file_comp.write(">"+name+"\n"+seq_AA+"\n")
            #file_comp.write(seq_AA)
            #file_comp.write("\n")
        #print(seq_AA)

