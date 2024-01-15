def trimin_forward(wildcards):
        if sample_data.loc[(wildcards.id), ["SRA"]].any():
            return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_1.fastq.gz"
        else:
            return "output/{id}/reads/local_reads/"+wildcards.id+"_1.fastq.gz"
def trimin_reverse(wildcards):
        if sample_data.loc[(wildcards.id), ["SRA"]].any():
            return "output/{id}/reads/downloaded_reads/"+wildcards.id+"_2.fastq.gz"
        else:
            return "output/{id}/reads/local_reads/"+wildcards.id+"_2.fastq.gz"

rule trimgalore:
        input:
                f = trimin_forward,
                r = trimin_reverse,
        params:
                wd = os.getcwd(),
                sample = "{id}",
                trimgalore_params = config["trimming"]["trimgalore_options"]["parameters"]
        singularity:
                "docker://chrishah/trim_galore:0.6.0"
        log:
                stdout = "results/{id}/logs/trimgalore.stdout.txt",
                stderr = "results/{id}/logs/trimgalore.stderr.txt"
        output:
                f_trimmed = "output/{id}/reads/trimmed/trimgalore/{id}_1P_trim.fastq.gz",
                f_orphans = "output/{id}/reads/trimmed/trimgalore/{id}_1P_unpaired.fastq.gz",
                r_trimmed = "output/{id}/reads/trimmed/trimgalore/{id}_2P_trim.fastq.gz",
                r_orphans = "output/{id}/reads/trimmed/trimgalore/{id}_2P_unpaired.fastq.gz",
                ok = "output/{id}/reads/trimmed/trimgalore/trim_{id}.ok"
        shadow: "minimal"
        threads: config["threads"]["trimmomatic"] 
        shell:
                """
                trim_galore \
                {params.trimgalore_params} \
                {input.f} {input.r} 1> {log.stdout} 2> {log.stderr}

                mv $(find ./ -name "*_val_1.fq.gz") {output.f_trimmed}
                mv $(find ./ -name "*_val_2.fq.gz") {output.r_trimmed}
	
                if [[ -f $(find ./ -name "*_unpaired_1.fq.gz") ]]; then mv $(find ./ -name "*_unpaired_1.fq.gz") {output.f_orphans}; else touch {output.f_orphans}; fi
                if [[ -f $(find ./ -name "*_unpaired_2.fq.gz") ]]; then mv $(find ./ -name "*_unpaired_2.fq.gz") {output.r_orphans}; else touch {output.r_orphans}; fi

                touch {output.ok}

                """
        
rule trimmomatic:
    input:
        f = trimin_forward,
        r = trimin_reverse
    output:
        fout = "output/{id}/reads/trimmed/trimmomatic/{id}_1P_trim.fastq.gz",
        funp = "output/{id}/reads/trimmed/trimmomatic/{id}_1P_unpaired.fastq.gz",
        rout = "output/{id}/reads/trimmed/trimmomatic/{id}_2P_trim.fastq.gz",
        runp = "output/{id}/reads/trimmed/trimmomatic/{id}_2P_unpaired.fastq.gz",
        ok = "output/{id}/reads/trimmed/trimmomatic/trim_{id}.ok"
    params:
        adapter = get_adapter,
        minlength = config["trimming"]["trimmomatic_options"]["minlength"],
        windowsize = config["trimming"]["trimmomatic_options"]["windowsize"],
        stepsize = config["trimming"]["trimmomatic_options"]["stepsize"],
        quality = config["trimming"]["trimmomatic_options"]["quality"],
        required_quality = config["trimming"]["trimmomatic_options"]["required_quality"],
        seed_mismatches = config["trimming"]["trimmomatic_options"]["seed_mismatches"],
        palindrome_clip = config["trimming"]["trimmomatic_options"]["palindrome_clip"],
        simple_clip = config["trimming"]["trimmomatic_options"]["simple_clip"]
    threads: config["threads"]["trimmomatic"] 
    singularity: 
        "docker://reslp/trimmomatic:0.38"
    shell:
        """
        if [[ ! -f {params.adapter} ]]; then
            echo "Adpater file not found. Please check your config files." >&2
            exit 1
        fi
        trimmomatic PE -threads {threads} {input.f} {input.r} {output.fout} {output.funp} {output.rout} {output.runp} ILLUMINACLIP:{params.adapter}:{params.seed_mismatches}:{params.palindrome_clip}:{params.simple_clip} LEADING:{params.quality} TRAILING:{params.quality} SLIDINGWINDOW:{params.windowsize}:{params.required_quality} MINLEN:{params.minlength}
        touch {output.ok}
        """
