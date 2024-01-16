import os
import pandas as pd
import glob

IDS = sample_data.index.values.tolist()
Assembler = config["Assembler"] 
sub = config["sub"]

def get_accession(wildcards):
        return sample_data.loc[(wildcards.id), ["SRA"]].values[0]

def get_seed(wildcards):
        seedfile = sample_data.loc[(wildcards.id), ["seed"]].values[0]
        #now we make sure that seedfile string to avoid problems with nan
        seedfile = str(seedfile)
        if os.path.exists(seedfile):
            return seedfile
        else:
            print("The assemblers GetOrganelle, Novoplasty and Mitobim require a valid seed file. You don't seem to be providing one ('"+seedfile+"') for sample: "+wildcards.id+". Please doublecheck the column 'seed' in your "+config["samples"]+".")
            os._exit(1)


def get_reads_for_assembly_fw(wildcards):
	if (wildcards.sub == "all"):
		if (config["skip_trimming"] == "yes"):
			if sample_data.loc[(wildcards.id), ["SRA"]].any():
				return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_1.fastq.gz"
			else:
				return "output/{id}/reads/local_reads/"+wildcards.id+"_1.fastq.gz"
		else:
			return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_1P_trim.fastq.gz"	
	else:
		return "output/{id}/reads/sub/{sub}/{id}_1.fastq.gz"

def get_reads_for_assembly_rv(wildcards):
	if (wildcards.sub == "all"):
		if (config["skip_trimming"] == "yes"):
			if sample_data.loc[(wildcards.id), ["SRA"]].any():
				return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_2.fastq.gz"
			else:
				return "output/{id}/reads/local_reads/"+wildcards.id+"_2.fastq.gz"
		else:
			return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_2P_trim.fastq.gz"	
	else:
		return "output/{id}/reads/sub/{sub}/{id}_2.fastq.gz"

def get_trimmed_reads_fw(wildcards):
	if (config["skip_trimming"] == "yes"):
		if sample_data.loc[(wildcards.id), ["SRA"]].any():
			return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_1.fastq.gz"
		else:
			return "output/{id}/reads/local_reads/"+wildcards.id+"_1.fastq.gz"
	else:
		return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_1P_trim.fastq.gz"	

def get_trimmed_reads_rv(wildcards):
	if (config["skip_trimming"] == "yes"):
		if sample_data.loc[(wildcards.id), ["SRA"]].any():
			return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_2.fastq.gz"
		else:
			return "output/{id}/reads/local_reads/"+wildcards.id+"_2.fastq.gz"
	else:
		return "output/{id}/reads/trimmed/"+config["trimming"]["software"]+"/{id}_2P_trim.fastq.gz"	


def get_clade(wildcards):
        return sample_data.loc[(wildcards.id), ["Clade"]].dropna().values[0]

def get_code(wildcards):
        return sample_data.loc[(wildcards.id), ["Code"]].dropna().values[0]

def get_forward(wildcards):
        if len(sample_data.loc[(wildcards.id), ["forward"]].dropna()) == 0:
                return
        else:
                return sample_data.loc[(wildcards.id), ["forward"]].dropna().values[0]

def get_reverse(wildcards):
        if len(sample_data.loc[(wildcards.id), ["reverse"]].dropna()) == 0:
                return
        else:
                return sample_data.loc[(wildcards.id), ["reverse"]].dropna().values[0]

def get_kmer(wildcards):
        return sample_data.loc[(wildcards.id), ["novoplasty_kmer"]].dropna().values[0]

def get_readlength(wildcards):
        return sample_data.loc[(wildcards.id), ["Read_length"]].dropna().values[0]

def get_adapter(wildcards):
        return sample_data.loc[(wildcards.id), ["Adapter"]].dropna().values[0]

def get_type(wildcards):
        return sample_data.loc[(wildcards.id), ["Type"]].dropna().values[0]

def get_rounds(wildcards):
        return sample_data.loc[(wildcards.id), ["GO_Rounds"]].dropna().values[0]

def gather_assemblies(wildcards):
        return  glob.glob("output/gathered_assemblies/*.fasta")

##use this instead of expand to input all mitos.done files for annotation_stats rule

def pick_stats(wildcards):
    pull_list = []
    if os.environ["RUNMODE"] == "annotate":
        for f in glob.glob("output/gathered_assemblies/*.fasta"):
            (i,s,a) = os.path.basename(f).split(".")[:-1]
            pull_list.append("output/"+i+"/annotation/mitos/"+i+"."+s+"."+a+".mitos.done")
        return pull_list
    else:
        pull_list = []
        for i in IDS:
            for s in sub:
                for a in Assembler:
                    pull_list.append("output/"+i+"/annotation/mitos/"+i+"."+s+"."+a+".mitos.done")
        PULL_LIST = pull_list
        return PULL_LIST 

##use this instead of expand to input all second_mitos.done files for gene_positions rule

def pick_mitos2(wildcards):
    pull_list = []
    if os.environ["RUNMODE"] == "annotate":
        for f in glob.glob("output/gathered_assemblies/*.fasta"):
            (i,s,a) = os.path.basename(f).split(".")[:-1]
            pull_list.append("output/"+i+"/annotation/second_mitos/"+i+"."+s+"."+a+".second_mitos.done")
        return pull_list
    else:
        pull_list = []
        for i in IDS:
            for s in sub:
                for a in Assembler:
                    pull_list.append("output/"+i+"/annotation/second_mitos/"+i+"."+s+"."+a+".second_mitos.done")
        PULL_LIST = pull_list
        return PULL_LIST 

##use this as the driver to pick runmode from report rule

def pick_mode(wildcards):
    pull_list = []
    if os.environ["RUNMODE"] == "annotate":
        for f in glob.glob("output/gathered_assemblies/*.fasta"):
            (i,s,a) = os.path.basename(f).split(".")[:-1]
            pull_list.append("output/"+i+"/annotation/compare/CCT/"+i+"."+s+"."+a+".CCT.done")
        print("Mode is 'annotate': ", len(pull_list), "target files.")
        for f in pull_list:
            print(f)
        return pull_list
    elif os.environ["RUNMODE"] == "assembly":
        for i in IDS:
            for s in sub:
                for a in Assembler:    
                    pull_list.append("output/"+i+"/assemblies/"+s+"/"+a+"/"+a+".ok")
        print("Mode is 'assembly': ", len(pull_list), "target files.")
        for f in pull_list:
            print(f)
        return pull_list
    else:
        pull_list = []
        for i in IDS:
            for s in sub:
                for a in Assembler:
                    pull_list.append("output/"+i+"/annotation/compare/CCT/"+i+"."+s+"."+a+".CCT.done")
        PULL_LIST = pull_list
        print("Mode is 'all': ", len(PULL_LIST), "target files.")
        return PULL_LIST 

