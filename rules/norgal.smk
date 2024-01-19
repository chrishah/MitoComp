rule norgal:
    input:
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv
    output:
        ok = "output/{id}/assemblies/{sub}/norgal/norgal.ok"
    params:
        wd = os.getcwd(),
        outdir = "output/{id}/assemblies/{sub}/norgal"
    log:
        stdout = "output/{id}/assemblies/{sub}/norgal/stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/norgal/stderr.txt"
    benchmark: "output/{id}/assemblies/{sub}/norgal/{id}.{sub}.norgal.benchmark.txt"
    threads: config["threads"]["norgal"] 
    singularity: "docker://reslp/norgal:1.0"
    shell:
        """
        # add a directory from the container to the PATH, so norgal finds all necessary executables
        export PATH="/software/norgal/binaries/linux:$PATH"
        # if norgal was run before, remove the previous run
        if [ -d {params.outdir}/run ]; then rm -rf {params.outdir}/run; fi
	cd {params.outdir}

        # run norgal - capture returncode, so if it fails, the pipeline won't stop
        norgal.py -i {params.wd}/{input.f} {params.wd}/{input.r} -o {params.wd}/{params.outdir}/run --blast -t {threads} 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr} && returncode=$? || returncode=$? 
        if [ $returncode -gt 0 ]
        then
            echo -e "\\n#### [$(date)]\\tnorgal exited with an error - moving on - for details see: {params.wd}/{log.stderr}" 1>> {params.wd}/{log.stdout}
            touch {wildcards.id}.{wildcards.sub}.norgal.fasta.error
        fi

        #if the expected final assembly exists, get a copy
        final_fasta=$(find {params.wd}/{params.outdir}/run -name "circular.candidate.fa")
        # check if the search returned only one file and copy if yes
        if [[ -z $final_fasta ]]
        then
            echo -e "\\n#### [$(date)]\\tnorgal has not produced a circularized assembly - moving on" 1>> {params.wd}/{log.stdout}
            touch {wildcards.id}.{wildcards.sub}.norgal.fasta.missing
        elif [ "$(echo $final_fasta | tr ' ' '\\n' | wc -l)" -eq 1 ] && [ $(grep "^>" $final_fasta | wc -l) -eq 1 ]
        then
            echo -e "\\n#### [$(date)]\\tnorgal seems to have produced a final circularized assembly - Kudos!" 1>> {params.wd}/{log.stdout}
            cp $final_fasta {wildcards.id}.{wildcards.sub}.norgal.fasta
            echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.norgals.fasta" 1>> {params.wd}/{log.stdout}
            cp $final_fasta {params.wd}/output/gathered_assemblies/{wildcards.id}.{wildcards.sub}.norgal.fasta
        else
            echo -e "\\n#### [$(date)]\\tnorgal seems to have produced multiple circularized assemblies or assemblies containing multiple sequences - don't know which to pick - moving on" 1>> {params.wd}/{log.stdout}
            touch {wildcards.id}.{wildcards.sub}.norgal.fasta.needs_attention
        fi
        touch {params.wd}/{output.ok}
        """

