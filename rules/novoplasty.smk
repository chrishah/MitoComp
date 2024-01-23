rule NOVOconfig:
    input:
        "bin/NOVOconfig.txt",
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv,
    output:
        "output/{id}/assemblies/{sub}/novoplasty/NOVOconfig_{id}_{sub}.txt"
    params:
        project_name = "{id}_{sub}",
        wd = os.getcwd(),
        seed = get_seed,
        log = "output/{id}/assemblies/{sub}/novoplasty/NOVOconfig_{id}_{sub}_log.txt",
        kmer = get_kmer,
        Read_length = get_readlength
    shell:
        """
        forward=$(realpath {params.wd}/{input.f})
        reverse=$(realpath {params.wd}/{input.r})

        cp {input[0]} {output}
        sed -i 's?^Project name.*?Project name = {params.project_name}?g' {output}
        sed -i 's?^Seed Input.*?Seed Input = {params.wd}/{params.seed}?g' {output}
        sed -i 's?^Extended log.*?Extended log = {params.wd}/{params.log}?g' {output}
        sed -i "s?^Forward reads.*?Forward reads = $forward?g" {output}
        sed -i "s?^Reverse reads.*?Reverse reads = $reverse?g" {output}
        sed -i 's?^K-mer.*?K-mer = {params.kmer}?g' {output}
        sed -i 's?^Read Length.*?Read Length = {params.Read_length}?g' {output}
        """

rule NOVOplasty:
    input:
        config = rules.NOVOconfig.output,
        f = get_reads_for_assembly_fw,
        r = get_reads_for_assembly_rv
    output: 
        ok = "output/{id}/assemblies/{sub}/novoplasty/novoplasty.ok"
    params:
        wd = os.getcwd(),
        outdir = "output/{id}/assemblies/{sub}/novoplasty/"
    log:
        stdout = "output/{id}/assemblies/{sub}/novoplasty/stdout.txt",
        stderr = "output/{id}/assemblies/{sub}/novoplasty/stderr.txt"
    benchmark: "output/{id}/assemblies/{sub}/novoplasty/{id}.{sub}.novoplasty.benchmark.txt"
    threads: config["threads"]["novoplasty"] 
#    threads: per_sample_config["Sy04"]["threads"]["novoplasty"] 
    singularity: "docker://reslp/novoplasty:4.2"
    shell:
        """
        # if novoplasty was run before, remove the previous run
        if [ -d {params.outdir}/run ]; then rm -rf {params.outdir}/run; fi
        mkdir -p {params.outdir}/run
        cd {params.outdir}/run

        NOVOPlasty.pl -c {params.wd}/{input.config} 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr} && returncode=$? || returncode=$?
        if [ $returncode -gt 0 ]
        then
            echo -e "\\n#### [$(date)]\\tnovoplasty exited with an error - moving on - for details see: {params.wd}/{log.stderr}" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.novoplasty.fasta.error
        fi

        cd ..
        # find the expected final assembly file
        final_fasta=$(find ./ -name "Circularized_assembly*")
        # check if the variable is empty
        if [[ -z $final_fasta ]]
        then
            echo -e "\\n#### [$(date)]\\tnovoplasty has not produced a circularized assembly - moving on" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.novoplasty.fasta.missing
        elif [ "$(echo $final_fasta | tr ' ' '\\n' | wc -l)" -eq 1 ] && [ $(grep "^>" $final_fasta | wc -l) -eq 1 ]
        then
            echo -e "\\n#### [$(date)]\\tnovoplasty seems to have produced a circularized assembly - Kudos!" 1>> {params.wd}/{log.stdout}
            ln -sf $final_fasta {wildcards.id}.{wildcards.sub}.novoplasty.fasta
            echo -e "\\n#### [$(date)]\\tfind a copy of the assembly at: {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.novoplasty.fasta" 1>> {params.wd}/{log.stdout}
        else
            echo -e "\\n#### [$(date)]\\tnovoplasty seems to have produced multiple circularized assemblies or assemblies containing multiple sequences - don't know which to pick - moving on" 1>> {params.wd}/{log.stdout}
            touch {params.wd}/{params.outdir}/{wildcards.id}.{wildcards.sub}.novoplasty.fasta.needs_attention
        fi
        touch {params.wd}/{output.ok}
        """
