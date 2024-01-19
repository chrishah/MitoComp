rule get_organelle:
    input:
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv
    output:
        ok = "output/{id}/assemblies/{sub}/getorganelle/getorganelle.ok"
    params:
        wd = os.getcwd(),
        outdir = "output/{id}/assemblies/{sub}/getorganelle",
        seed = get_seed,
        type = get_type,
        rounds = get_rounds
    singularity:"docker://reslp/getorganelle:1.7.1"
    log:
        stdout = "output/{id}/assemblies/{sub}/getorganelle/stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/getorganelle/stderr.txt" 
    benchmark: "output/{id}/assemblies/{sub}/getorganelle/{id}.{sub}.getorganelle.benchmark.txt"
    threads: config["threads"]["getorganelle"] 
    shell:
        """
        cd {params.outdir}

        # run getorganelle - capture returncode, so if it fails, the pipeline won't stop
        get_organelle_from_reads.py -1 {params.wd}/{input.f} -2 {params.wd}/{input.r} -o {params.wd}/{params.outdir}/run -F {params.type} -t {threads} -R {params.rounds} -s {params.wd}/{params.seed} 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr} && returncode=$? || returncode=$?
        if [ $returncode -gt 0 ]
        then
            echo -e "\\n#### [$(date)]\\tgetorganelle exited with an error - moving on - for details see: {params.wd}/{log.stderr}" 1>> {params.wd}/{log.stdout}
            touch {wildcards.id}.{wildcards.sub}.getorganelle.fasta.error
        fi

        #if the expected final assembly exists, get a copy
        if [[ -d {params.wd}/{params.outdir} ]]
        then #check first of folder exists in case get_organelle exits with the wrong exit code.
            final_fasta=$(find {params.wd}/{params.outdir}/ -maxdepth 1 -name "*.fasta")
        fi
        # check if the search returned only one file and copy if yes -- also check only 1 sequence in final fasta
        if [[ -z $final_fasta ]] 
        then
            echo -e "\\n#### [$(date)]\\tgetorganelle has not produced the final assembly - moving on" 1>> {params.wd}/{log.stdout}
            touch {wildcards.id}.{wildcards.sub}.getorganelle.fasta.missing
        elif [ "$(echo $final_fasta | tr ' ' '\\n' | wc -l)" -eq 1 ] && [ $(grep "^>" $final_fasta | wc -l) -eq 1 ]
        then
            echo -e "\\n#### [$(date)]\\tgetorganelle seems to have produced a final assembly - Kudos!" 1>> {params.wd}/{log.stdout}
            cp $final_fasta {wildcards.id}.{wildcards.sub}.getorganelle.fasta
            echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.getorganelle.fasta" 1>> {params.wd}/{log.stdout}
            cp $final_fasta {params.wd}/output/gathered_assemblies/{wildcards.id}.{wildcards.sub}.getorganelle.fasta
        else
            echo -e "\\n#### [$(date)]\\tgetorganelle seems to have produced multiple assemblies or assemblies containing multiple sequences - don't know which to pick - moving on" 1>> {params.wd}/{log.stdout}
            touch {wildcards.id}.{wildcards.sub}.getorganelle.fasta.needs_attention
        fi
        touch {params.wd}/{output.ok}
        """
