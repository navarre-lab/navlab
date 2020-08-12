#!/Users/sornnujahkathirgamanathan/miniconda/bin/python
import os
import nav
import parser

project_name = input("What is the name of the dict and blast db for this project?")
pct_id = input("What pct ID cutoff for the tab files?")

genbank_file_list = []
directory_list = []
for entry in os.scandir('.'):
    # check if entry is a ".gbff" file
    if entry.is_file() and (entry.name[0] != ".") and (entry.name.split(".")[-1] == "gbff"):
        with open(entry) as fp:
            line = str(fp.readline())
            if "LOCUS" in line:
                # print(str(i)+". "+entry.name + "\t{}".format(line.strip()))
                genbank_file_list.append(entry.name)
    elif entry.is_dir():
        directory_list.append(entry.name)

# if not already in the current directory make the directories gbk and fna
if "gbk" not in directory_list:
    os.mkdir("gbk")
if "fna" not in directory_list:
    os.mkdir("fna")
if "faa" not in directory_list:
    os.mkdir("faa")
if "xml" not in directory_list:
    os.mkdir("xml")
if "tab" not in directory_list:
    os.mkdir("tab")

# Now that we have a list of gbk like files we will make a dictionary out of
# them and rename them to something better than what NCBI does

project_dictionary = []

for gbff_file in genbank_file_list:
    dict_entry = nav.dictmaker(gbff_file)
    new_name = dict_entry.split(" | ")[0]
    os.rename(gbff_file, str("gbk/" + new_name + ".gbk"))
    project_dictionary.append(dict_entry)

dictionary_file = open(project_name + ".dict", "w")
dictionary_file.write("CLUSTER | | | \n")
dictionary_file.write("\n".join(project_dictionary))
dictionary_file.close()

# part 2 - now create fna files and put in the fna directory

for entry in os.scandir('gbk'):
    if entry.is_file():
        if str(entry.name)[0] != "." and str(entry.name.split(".")[-1]) == "gbk":
            pass_file = "gbk/" + str(entry.name)
            nav.convert_gbk_to_fna(pass_file)
            new_name = str(entry.name.split(".")[0]) + ".fna"
            os.rename("gbk/" + new_name, "fna/" + new_name)

# part 3 - now create faa files and put in the faa directory_list
for entry in os.scandir('fna'):
    if entry.is_file():
        if str(entry.name)[0] != "." and str(entry.name.split(".")[-1]) == "fna":
            pass_file = "fna/" + str(entry.name)
            nav.convert_fna_to_faa(pass_file)
            new_name = str(entry.name.split(".")[0]) + ".faa"
            os.rename("fna/" + new_name, "faa/" + new_name)

# make a blast db
cat_faa_cmd = "cat faa/*.faa > faa/temp.faa"  # now you have all proteins together
os.system(cat_faa_cmd)

parser.pseudoremover("faa/temp.faa")

# perform mmseq2 clustering
cluster_cmd1 = "mmseqs createdb faa/pseudofree.faa DB"
cluster_cmd2 = "mmseqs cluster DB clusteredDB tmp --min-seq-id " + str(float(pct_id) / 100) + " -c 0.9 --cov-mode 0"
cluster_cmd3 = "mmseqs createseqfiledb DB clusteredDB clu_seq"
cluster_cmd4 = "mmseqs result2flat DB DB clu_seq clu_seq.fasta"

os.system(cluster_cmd1)
os.system(cluster_cmd2)
os.system(cluster_cmd3)
os.system(cluster_cmd4)

# remove the concatenated fna files
remove_cmd = "rm faa/temp.faa && rm faa/pseudofree.faa"
os.system(remove_cmd)

parser.parse_cluster("clu_seq.fasta")

# concatenate faa file contents to temp.faa
cat_faa_cmd = "cat faa/*.faa > faa/temp.faa"
os.system(cat_faa_cmd)

# make a blast db
blast_cmd1 = "makeblastdb -parse_seqids -in faa/temp.faa -dbtype prot -out " + project_name
os.system(blast_cmd1)
blast_cmd2 = "diamond makedb -p N --in faa/temp.faa -d " + project_name
os.system(blast_cmd2)
# remove the concatenated fna files
remove_cmd = "rm faa/temp.faa && rm faa/temp1.faa"
os.system(remove_cmd)

faa_list = []
for entry in os.scandir('faa'):
    if entry.is_file() and (entry.name[0] != ".") and (entry.name.split(".")[-1]) == "faa":
        faa_list.append(entry.name)

for faa_file in faa_list:
    xml_name = "xml/" + faa_file.split(".")[0] + ".xml"
    tab_name = "tab/" + faa_file.split(".")[0] + ".tab"
    dict_name = project_name + ".dict"
    blast_command = "diamond blastp --max-target-seqs 0 --outfmt 5 --query-cover 90 --subject-cover 90 --query faa/" +\
                    faa_file + " --db " + project_name + " --out " + xml_name
    os.system(blast_command)
    parser.pxml(xml_name, tab_name, dict_name, pct_id)
    os.system("rm " + xml_name)
