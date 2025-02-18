rule mitoz:
    input:
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv,
    output:
        ok = "output/{id}/assemblies/{sub}/mitoz/mitoz.ok",
    params:
        wd = os.getcwd(),
        outdir = "output/{id}/assemblies/{sub}/mitoz",
        id = "{id}",
	non_default = config["Assembler_options"]["mitoz"]["non_default"],
        clade = get_clade,
        genetic_code = get_code,
    log:
        stdout = "output/{id}/assemblies/{sub}/mitoz/stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/mitoz/stderr.txt"
    benchmark: "output/{id}/assemblies/{sub}/mitoz/{id}.{sub}.mitoz.benchmark.txt"
    threads: config["threads"]["mitoz"] 
    singularity: "docker://guanliangmeng/mitoz:3.4"
    shell:
        """
        cd {params.outdir}

        # run mitoz - capture returncode, so if it fails, the pipeline won't stop 
        mitoz assemble --outprefix {wildcards.id} --thread_number {threads} --fq1 {params.wd}/{input.f} --fq2 {params.wd}/{input.r} --genetic_code {params.genetic_code} --clade {params.clade} --requiring_taxa {params.clade} {params.non_default} 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr} && returncode=$? || returncode=$?
        if [ $returncode -gt 0 ]
        then
            echo -e "\\n#### [$(date)]\\tmitoz exited with an error - moving on - for details see: {params.wd}/{log.stderr}" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoz.fasta.error
        fi
 
        #if the expected final assembly exists, get a copy
        final_fasta=$(find -name "MitoFlex.picked.fa")
        # check if the search returned only one file and copy if yes
        if [[ -z $final_fasta ]]
        then
            echo -e "\\n#### [$(date)]\\tmitoz has not produced a final assembly - moving on" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoz.fasta.missing
        elif [ "$(echo $final_fasta | tr ' ' '\\n' | wc -l)" -eq 1 ] && [ $(grep "^>" $final_fasta | wc -l) -eq 1 ]
        then
            echo -e "\\n#### [$(date)]\\tmitoz seems to have produced a final assembly - Kudos!" 1>> {params.wd}/{log.stdout}
            ln -sf $final_fasta {wildcards.id}.{wildcards.sub}.mitoz.fasta
            echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoz.fasta" 1>> {params.wd}/{log.stdout}
        else
            echo -e "\\n#### [$(date)]\\tmitoz seems to have produced multiple assemblies or assemblies containing multiple sequences - don't know which to pick - moving on" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoz.fasta.needs_attention
        fi
        touch {params.wd}/{output.ok}
        """
