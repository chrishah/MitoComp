rule setup_mitoflex_db:
    output: 
        ok = "bin/MitoFlex/mitoflex.db.status.ok",
    params:
        wd = os.getcwd(),
    singularity: "docker://samlmg/mitoflex:v0.2.9"
    threads: 1
    shell:
        """
        cp -pfr /MitoFlex/* bin/MitoFlex/
        cp bin/ncbi_custom.py bin/MitoFlex/ncbi.py
        cd bin/MitoFlex/
        export HOME=$(pwd)
        echo $HOME
        #execute modified ncbi.py script with 'y' or 'n' as additional options (our modification allows for non-interactive use of the script) - 'y' -> download taxdump; 'n' -> use existing taxdump
        ./ncbi.py y 
        touch {params.wd}/{output.ok}
        """

rule mitoflex:
    input:
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv,
        db = rules.setup_mitoflex_db.output
    output:
        ok = "output/{id}/assemblies/{sub}/mitoflex/mitoflex.ok",
    params:
        wd = os.getcwd(),
        outdir = "output/{id}/assemblies/{sub}/mitoflex",
        id = "{id}",
	non_default = config["Assembler_options"]["mitoflex"]["non_default"],
        clade = get_clade,
        genetic_code = get_code,
    log:
        stdout = "output/{id}/assemblies/{sub}/mitoflex/stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/mitoflex/stderr.txt"
    benchmark: "output/{id}/assemblies/{sub}/mitoflex/{id}.{sub}.mitoflex.benchmark.txt"
    threads: config["threads"]["mitoflex"] 
    singularity: "docker://samlmg/mitoflex:v0.2.9"
    shell:
        """
        cd {params.outdir}
        export HOME="{params.wd}/bin/MitoFlex"

            echo -e "\\n#### [$(date)]\\tFind full mitoflex log here: {params.outdir}/MitoFlex/MitoFlex.log" 1> {params.wd}/{log.stdout}

        # run mitoflex - capture returncode, so if it fails, the pipeline won't stop 
        {params.wd}/bin/MitoFlex/MitoFlex.py assemble --workname MitoFlex --threads {threads} --fastq1 {params.wd}/{input.f} --fastq2 {params.wd}/{input.r} {params.non_default} 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr} && returncode=$? || returncode=$?

	#findmitoscaf --workname MitoFlex --threads {threads} --fastq1 {params.wd}/{input.f} --fastq2 {params.wd}/{input.r}--genetic-code {params.genetic_code} --clade {params.clade} {params.non_default}
        if [ $returncode -gt 0 ]
        then
            echo -e "\\n#### [$(date)]\\tmitoflex exited with an error - moving on - for details see: {params.wd}/{log.stderr}" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoflex.fasta.error
        fi
 
        #if the expected final assembly exists, get a copy
        final_fasta=$(find -name "MitoFlex.picked.fa")
        # check if the search returned only one file and copy if yes
        if [[ -z $final_fasta ]]
        then
            echo -e "\\n#### [$(date)]\\tmitoflex has not produced a final assembly - moving on" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoflex.fasta.missing
        elif [ "$(echo $final_fasta | tr ' ' '\\n' | wc -l)" -eq 1 ] && [ $(grep "^>" $final_fasta | wc -l) -eq 1 ]
        then
            echo -e "\\n#### [$(date)]\\tmitoflex seems to have produced a final assembly - Kudos!" 1>> {params.wd}/{log.stdout}
            ln -sf $final_fasta {wildcards.id}.{wildcards.sub}.mitoflex.fasta
            echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoflex.fasta" 1>> {params.wd}/{log.stdout}
        else
            echo -e "\\n#### [$(date)]\\tmitoflex seems to have produced multiple assemblies or assemblies containing multiple sequences - don't know which to pick - moving on" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.mitoflex.fasta.needs_attention
        fi
        touch {params.wd}/{output.ok}
        """
