import os
import pandas as pd
import glob
import copy
import yaml
import numpy as np

per_sample_config = {}
if os.path.exists(str(config["samples"])):
    print("Reading sample specific data from file '"+str(config["samples"]))
    sample_data = pd.read_table(config["samples"], sep=",", na_values=None).set_index("ID", drop=False)
    IDS = sample_data.index.values.tolist()
    ind_options = {'ID': "ID",
        'forward': "forward",
        'reverse': "reverse",
        'Assembler': 'Assembler',
        'seed': "seed",
        'SRA': "SRA", 
        'Circules': "Assembler_options:mitobim:circules",
        'Clade': "Assembler_options:mitoflex:clade",
        'Code': "Assembler_options:mitoflex:code",
        'novoplasty_kmer': "Assembler_options:novoplasty:kmer",
        'Read_length': "Assembler_options:novoplasty:readlength",
        'Adapter': "trimming:trimmomatic_options:adapter",
        'Type': "Assembler_options:getorganelle:type",
        'GO_Rounds': "Assembler_options:getorganelle:GO_rounds"}
    
    #remove irrelevant columns from list
    present = list(ind_options.keys())
    for c in present:
        if not c in sample_data.columns:
            print(c+" is not in there - remove from list")
            del(ind_options[c])
    present = ind_options.keys()
    
    for ID in IDS:
        print(ID)
        if "sample_config" in sample_data.columns and os.path.exists(str(sample_data.loc[ID]["sample_config"])):
            yamlfile = sample_data.loc[ID]["sample_config"]
            print("attempting to read in sample specific config file - expecting yaml: "+str(yamlfile))
            with open(yamlfile) as f:
                per_sample_config[ID] = yaml.safe_load(f)
            continue
        per_sample_config[ID] = copy.deepcopy(config)
        print(per_sample_config[ID])
        print("checking for individual options")
        for c in present:
            print(c)
            li = ind_options[c].split(":")
            print(li,len(li))
            if len(li) == 1:
                print("level 1:")
                if isinstance(per_sample_config[ID][li[0]], str):
                    in_config = str(per_sample_config[ID][li[0]])
                    in_data = str(sample_data.loc[ID][c])
                    print("in config: "+in_config)
                    print("in data: "+in_data)
                elif isinstance(per_sample_config[ID][li[0]], list):
                    in_config = sorted(per_sample_config[ID][li[0]])
                    print(sample_data.loc[ID][c])
                    print(type(sample_data.loc[ID][c]))
                    if str(sample_data.loc[ID][c]) == "nan":
                        print("Nonetype")
                        in_data = str(sample_data.loc[ID][c])
                    else:
                        print("splitting")
                        in_data = sorted(str(sample_data.loc[ID][c]).split("|"))
                    print("in config: "+str(in_config))
                    print("in data: "+str(in_data))
                elif isinstance(per_sample_config[ID][li[0]], type(None)):
                    print("Found NoneType: "+str(per_sample_config[ID][li[0]]))
                    in_config = str(per_sample_config[ID][li[0]])
                    in_data = str(sample_data.loc[ID][c])
                    print("in config: "+in_config)
                    print("in data: "+in_data)
                if in_config != in_data and in_data != "nan":
                    print("need to update with: "+str(in_data))
                    per_sample_config[ID][li[0]] = in_data
            elif len(li) == 2:
                print("level 2:")
                in_config = str(per_sample_config[ID][li[0]][li[1]])
                in_data = str(sample_data.loc[ID][c])
                print("in config: "+in_config)
                print("in data: "+in_data)
                if in_config != in_data:
                    print("need to update: "+str(in_data))
                    per_sample_config[ID][li[0]][li[1]] = in_data
            elif len(li) == 3:
                print("level 3:")
                in_config = str(per_sample_config[ID][li[0]][li[1]][li[2]])
                in_data = str(sample_data.loc[ID][c])
                print("in config: "+in_config)
                print("in data: "+in_data)
                if in_config != in_data and in_data != "nan":
                    print("need to update: "+str(in_data))
                    per_sample_config[ID][li[0]][li[1]][li[2]] = in_data
else:
    print("Reading from config file only")
    IDS = [str(config["ID"])]
    per_sample_config[str(config["ID"])] = config.copy()
    
           
print("\n## See what's what")
for ID in per_sample_config.keys():
        print(ID+":")
        print(per_sample_config[ID])

#Assembler = config["Assembler"] 
#sub = config["sub"]

def get_accession(wildcards):
        return per_sample_config[wildcards.id]["SRA"]

def trimin_forward(wildcards):
        if per_sample_config[wildcards.id]["SRA"]:
            return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_1.fastq.gz"
        else:
            return "output/{id}/reads/local_reads/"+wildcards.id+"_1.fastq.gz"
def trimin_reverse(wildcards):
        if per_sample_config[wildcards.id]["SRA"]:
            return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_2.fastq.gz"
        else:
            return "output/{id}/reads/local_reads/"+wildcards.id+"_2.fastq.gz"

def get_seed(wildcards):
        seedfile = per_sample_config[wildcards.id]["seed"]
        #now we make sure that seedfile string to avoid problems with nan
        seedfile = str(seedfile)
        if os.path.exists(seedfile):
            return seedfile
        else:
            print("WARNING:\nThe assemblers GetOrganelle, Novoplasty and Mitobim require a valid seed file. You don't seem to be providing one (currently: '"+seedfile+"') for sample: "+wildcards.id+". Please doublecheck.\n")
            os._exit(1)


def get_reads_for_assembly_fw(wildcards):
	if (wildcards.sub == "all"):
		if (per_sample_config[wildcards.id]["trimming"]["skip"] == "yes"):
			if per_sample_config[wildcards.id]["SRA"]:
				return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_1.fastq.gz"
			else:
				return "output/{id}/reads/local_reads/"+wildcards.id+"_1.fastq.gz"
		else:
			return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_1P_trim.fastq.gz"	
	else:
		return "output/{id}/reads/sub/{sub}/{id}_1.fastq.gz"

def get_reads_for_assembly_rv(wildcards):
	if (wildcards.sub == "all"):
		if (per_sample_config[wildcards.id]["trimming"]["skip"] == "yes"):
			if per_sample_config[wildcards.id]["SRA"]:
				return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_2.fastq.gz"
			else:
				return "output/{id}/reads/local_reads/"+wildcards.id+"_2.fastq.gz"
		else:
			return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_2P_trim.fastq.gz"	
	else:
		return "output/{id}/reads/sub/{sub}/{id}_2.fastq.gz"

def get_trimmed_reads_fw(wildcards):
	if (per_sample_config[wildcards.id]["trimming"]["skip"] == "yes"):
		if per_sample_config[wildcards.id]["SRA"]:
			return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_1.fastq.gz"
		else:
			return "output/{id}/reads/local_reads/"+wildcards.id+"_1.fastq.gz"
	else:
		return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_1P_trim.fastq.gz"	

def get_trimmed_reads_rv(wildcards):
	if (per_sample_config[wildcards.id]["trimming"]["skip"] == "yes"):
		if per_sample_config[wildcards.id]["SRA"]:
			return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_2.fastq.gz"
		else:
			return "output/{id}/reads/local_reads/"+wildcards.id+"_2.fastq.gz"
	else:
		return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_2P_trim.fastq.gz"	


def get_clade(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["mitoflex"]["clade"]

def get_code(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["mitoflex"]["code"]

def get_circules_mode(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["mitobim"]["circules"]

def get_forward(wildcards):
        if os.path.exists(str(per_sample_config[wildcards.id]["forward"])):
                return per_sample_config[wildcards.id]["forward"]
        else:
                return

def get_reverse(wildcards):
        if os.path.exists(str(per_sample_config[wildcards.id]["reverse"])):
                return per_sample_config[wildcards.id]["reverse"]
        else:
                return

def get_kmer(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["novoplasty"]["kmer"]

def get_readlength(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["novoplasty"]["readlength"]

def get_adapter(wildcards):
        return per_sample_config[wildcards.id]["trimming"]["trimmomatic_options"]["adapter"]

def get_type(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["getorganelle"]["type"]

def get_rounds(wildcards):
        return per_sample_config[wildcards.id]["Assembler_options"]["getorganelle"]["GO_rounds"]

def gather_assemblies(wildcards):
        return  glob.glob("output/gathered_assemblies/*.fasta")


def trigger_gather(wildcards):
    pull_list = []
    for i in [ str(i) for i in per_sample_config.keys()]:
        for s in [ str(s) for s in per_sample_config[i]["sub"]]:
            for a in [ str(a) for a in per_sample_config[i]["Assembler"]]:
                print("Looking for: output/"+i+"/assemblies/"+s+"/"+a+"/"+i+"."+s+"."+a+".fasta .. ", end="")
                if os.path.exists(str("output/"+i+"/assemblies/"+s+"/"+a+"/"+i+"."+s+"."+a+".fasta")):
                    print("found")
                    pull_list.append("output/"+i+"/assemblies/"+s+"/"+a+"/"+i+"."+s+"."+a+".fasta")
                else:
                    print("nothing")
    print("PULL_LIST: "+str(pull_list))
    return pull_list
        

# determine which conbinations to process
# if assembly or all - do all combinations that are possilbe
# if annotate - read in from gathered assembly and only do those
to_process = {"id": [], "sub": [], "assembler": []}
if os.environ["RUNMODE"] == "annotate":
    for f in glob.glob("output/gathered_assemblies/*.fasta"):
        (i,s,a) = os.path.basename(f).split(".")[:-1]
        if (i in IDS) and (s in per_sample_config[i]["sub"]) and (a in per_sample_config[i]["Assembler"]):
            to_process["id"].append(i)
            to_process["sub"].append(s)
            to_process["assembler"].append(a)
else:
    for i in IDS:
        for s in per_sample_config[i]["sub"]:
            for a in per_sample_config[i]["Assembler"]:    
                to_process["id"].append(i)
                to_process["sub"].append(s)
                to_process["assembler"].append(a)

print("to_process")
print(to_process)
def pick_assembly(wildcards):
    # this functions controls which assemblers are used for which sample - called in the Snakefile by the assembly_only rule
    pull_list = []
    for j in range(len(to_process["id"])):
        i = to_process["id"][j]
        s = to_process["sub"][j]
        a = to_process["assembler"][j]
        pull_list.append("output/"+i+"/assemblies/"+s+"/"+a+"/"+a+".ok")
    return pull_list

##use this instead of expand to input all mitos.done files for annotation_stats rule

def pick_stats(wildcards):
    pull_list = []
    for j in range(len(to_process["id"])):
        i = to_process["id"][j]
        s = to_process["sub"][j]
        a = to_process["assembler"][j]
        pull_list.append("output/"+i+"/annotation/mitos/"+i+"."+s+"."+a+".mitos.done")
    return pull_list

##use this instead of expand to input all second_mitos.done files for gene_positions rule

def pick_mitos2(wildcards):
    pull_list = []
    for j in range(len(to_process["id"])):
        i = to_process["id"][j]
        s = to_process["sub"][j]
        a = to_process["assembler"][j]
        pull_list.append("output/"+i+"/annotation/second_mitos/"+i+"."+s+"."+a+".second_mitos.done")
    return pull_list

##use this as the driver to pick runmode from report rule

def pick_mode(wildcards):
    pull_list = []
    for j in range(len(to_process["id"])):
        i = to_process["id"][j]
        s = to_process["sub"][j]
        a = to_process["assembler"][j]
        pull_list.append("output/"+i+"/annotation/compare/CCT/"+i+"."+s+"."+a+".CCT.done")

    if os.environ["RUNMODE"] == "annotate":
        print("Mode is 'annotate': ", len(pull_list), "target files.")
    else:
        print("Mode is 'all': ", len(pull_list), "target files.")
    
    for f in pull_list:
        print(f)
    return pull_list

def trigger_all(wildcards):
    if os.environ["RUNMODE"] == "annotate":
        if len(glob.glob("output/gathered_assemblies/*.fasta")) > 0:
            return["output/report/report.html"]
    else:
        return["output/report/report.html"]

