rule fastqdump:
        params:
                accession = get_accession,
                wd = os.getcwd()
        output:
                f = "output/{id}/reads/downloaded_reads/{id}_1.fastq.gz",
                r = "output/{id}/reads/downloaded_reads/{id}_2.fastq.gz"
        log: 
                stdout = "output/{id}/reads/downloaded_reads/stdout.txt",
                stderr = "output/{id}/reads/downloaded_reads/stderr.txt"
        singularity:
                "docker://reslp/sra-tools:2.10.9"
        threads: config["threads"]["download"] 
        shadow: "minimal"
        shell:
                """
                # configuration of sra-tools is messed up in singularity. This is connected with these issues:
                # https://github.com/ncbi/sra-tools/issues/291
                # https://standage.github.io/that-darn-cache-configuring-the-sra-toolkit.html
                echo "no local files, downloading from SRA: {params.accession}" 1> {log.stdout} 2> {log.stderr}
                export HOME=$(pwd)
                mkdir -p $HOME/.ncbi
                printf '/LIBS/GUID = "%s"\\n' `uuidgen` > $HOME/.ncbi/user-settings.mkfg
                mkdir -p $HOME/tmp
                echo "/repository/user/main/public/root = \'$HOME/tmp\'" >> $HOME/.ncbi/user-settings.mkfg

                # download
                #prefetch --max-size 1024000000 {params.accession} 1>> {log.stdout} 2>> {log.stderr}
                fastq-dump --split-files --gzip --defline-seq '@$ac-$sn/$ri' --defline-qual '+' {params.accession} 1>> {log.stdout} 2>> {log.stderr}

                #rename to expected output files
                mv {params.accession}_1.fastq.gz {output.f} 1>> {log.stdout} 2>> {log.stderr}
                mv {params.accession}_2.fastq.gz {output.r} 1>> {log.stdout} 2>> {log.stderr}

                """

rule prep_local_reads:
        input:
                f = get_forward,
                r = get_reverse
        params:
                wd = os.getcwd()
        output:
                f = "output/{id}/reads/local_reads/{id}_1.fastq.gz",
                r = "output/{id}/reads/local_reads/{id}_2.fastq.gz"
        threads: 1
        shell:
                """
                if [[ -f "{params.wd}/{input.f}" ]] && [[ -f "{params.wd}/{input.r}" ]]; then 
                    echo "using local fastq.gz files"
                    ln -s {params.wd}/{input.f} {output.f}
                    ln -s {params.wd}/{input.r} {output.r}
                fi
                """
