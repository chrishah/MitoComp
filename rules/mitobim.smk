rule interleave:
    input:
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv
    threads: config["threads"]["interleave"]
    output:
        "output/{id}/reads/interleave/{sub}/{id}_{sub}_interleaved.fastq.gz"
    singularity:
        "docker://reslp/bbmap:38.90"
    shell:
        """
        reformat.sh in1={input.f} in2={input.r} out={output}
        """

rule MITObim:
    input:
        rules.interleave.output
    output:
        ok = "output/{id}/assemblies/{sub}/mitobim/mitobim_preclip.ok"
    params:
        id = "{id}",
        seed = get_seed,
        wd = os.getcwd(),
	end = config["Assembler_options"]["mitobim"]["iterations"],
	non_default = config["Assembler_options"]["mitobim"]["non_default"],
        outdir = "output/{id}/assemblies/{sub}/mitobim"
    log: 
        stdout = "output/{id}/assemblies/{sub}/mitobim/stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/mitobim/stderr.txt"
    benchmark: "output/{id}/assemblies/{sub}/mitobim/{id}.{sub}.mitobim.benchmark.txt"
    singularity:
        "docker://chrishah/mitobim:v.1.9.1"
    threads: config["threads"]["mitobim"] 
    shell:
        """
        if [ -d {params.outdir}/run ]; then rm -rf {params.outdir}/run; fi
        mkdir -p {params.outdir}/run
        cd {params.outdir}/run
        
        # run mitobim - capture returncode, so if it fails, the pipeline won't stop
        MITObim.pl -sample {params.id} -ref seed -readpool {params.wd}/{input} --quick {params.wd}/{params.seed} -end {params.end} --NFS_warn_only {params.non_default} 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr} && returncode=$? || returncode=$?
        if [ $returncode -gt 0 ]
        then
            echo -e "\\n#### [$(date)]\\tmitobim exited with an error - moving on - for details see: {params.wd}/{log.stderr}" 1>> {params.wd}/{log.stdout}
            touch ../{wildcards.id}.{wildcards.sub}.mitobim.fasta.error
        fi

        #if the expected final assembly exists, get a copy
        final_fasta=$(find ./ -name "*noIUPAC.fasta")
        # check if the search returned only one file and copy if yes
        if [[ -z $final_fasta ]]
        then
            echo -e "\\n#### [$(date)]\\tthere was no error, but mitobim has not produced a final assembly - for details see: {params.wd}/{log.stdout} - moving on" 1>> {params.wd}/{log.stdout}
            touch ../{wildcards.id}.{wildcards.sub}.mitobim.fasta.missing
        elif [ "$(echo $final_fasta | tr ' ' '\\n' | wc -l)" -eq 1 ] && [ $(grep "^>" $final_fasta | wc -l) -eq 1 ]
        then
            echo -e "\\n#### [$(date)]\\tmitobim seems to have produced a final assembly - Kudos!" 1>> {params.wd}/{log.stdout}
            ln -sf {params.wd}/{params.outdir}/run/$final_fasta $(pwd)/../{wildcards.id}.{wildcards.sub}.mitobim-preclip.fasta
            echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitobim-preclip.fasta" 1>> {params.wd}/{log.stdout}
        else
            echo -e "\\n#### [$(date)]\\tmitobim seems to have produced multiple assemblies or assemblies containing multiple sequences - don't know which to pick - moving on" 1>> {params.wd}/{log.stdout}
            touch ../{wildcards.id}.{wildcards.sub}.mitobim.fasta.needs_attention
        fi
        touch {params.wd}/{output.ok}       
        """


rule MITObim_circules:
    input:
        rules.MITObim.output
    output:
        ok = "output/{id}/assemblies/{sub}/mitobim/mitobim.ok",
    params:
        id = "{id}",
        wd = os.getcwd(),
        outdir = "output/{id}/assemblies/{sub}/mitobim",
        readlength = get_readlength
    log: 
        stdout = "output/{id}/assemblies/{sub}/mitobim/mitobim_circules.stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/mitobim/mitobim_circules.stderr.txt"
    benchmark: "output/{id}/assemblies/{sub}/mitobim/{id}.{sub}.mitobim_circules.benchmark.txt"
    singularity:
        "docker://chrishah/mitobim:v.1.9.1"
    threads: 1
    shell:
        """
        cd {params.outdir}
        #make sure all remnants from previous attempts are gone
        rm {wildcards.id}.{wildcards.sub}.mitobim.*

        echo -e "\\n#### [$(date)]\\tChecking for mitobim assembly" 1> {params.wd}/{log.stdout}
        if [[ -f {wildcards.id}.{wildcards.sub}.mitobim-preclip.fasta ]]
        then
            echo -e "\\n#### [$(date)]\\tfound {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitobim-preclip.fasta" 1>> {params.wd}/{log.stdout}
            echo -e "#### [$(date)]\\tchecking for signs of circularity with circules.py\n" 1>> {params.wd}/{log.stdout}
            circules.py -f {wildcards.id}.{wildcards.sub}.mitobim-preclip.fasta 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

            echo -e "\\n#### [$(date)]\\tparsing output and clip assembly if criteria for circularity are fullfilled" 1>> {params.wd}/{log.stdout}
            {params.wd}/bin/parse_and_clip.sh {wildcards.id}.{wildcards.sub}.mitobim-preclip.fasta {params.wd}/{log.stdout} {params.readlength} {wildcards.id}.{wildcards.sub}.mitobim 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

            if [[ -f {wildcards.id}.{wildcards.sub}.mitobim.fasta ]]
            then            
                echo -e "\\n#### [$(date)]\\tmitobim seems to have produced a circular assembly - Kudos!" 1>> {params.wd}/{log.stdout}
                echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitobim.fasta" 1>> {params.wd}/{log.stdout}
            fi
        else
             echo -e "\\n#### [$(date)]\\tIt seems there is nothing worth checking - moving on" 1>> {params.wd}/{log.stdout}
        fi
        
        touch {params.wd}/{output.ok}
        """
