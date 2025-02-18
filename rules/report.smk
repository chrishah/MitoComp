rule annotation_statsII:
    input:
        pick_mode
    params:
        ids = get_ids()
    output:
        done = "output/stats/annotation_statsII.done",
    shell:
        """
	if [[ -f output/stats/assembly_pathsII.txt ]]; then rm output/stats/assembly_pathsII.txt; fi
	if [[ -f output/stats/bed_pathsII.txt ]]; then rm output/stats/bed_pathsII.txt; fi
	if [[ -f output/stats/map_pathsII.txt ]]; then rm output/stats/map_pathsII.txt; fi
        for id in $(echo "{params.ids}")
        do
            find ./output/$id/annotation/alignment/ -name "$id*.final.fasta" >> output/stats/assembly_pathsII.txt
            find ./output/$id/annotation/second_mitos/ -name "result.bed" >> output/stats/bed_pathsII.txt
            find ./output/$id/annotation/compare/CCT/*.png >> output/stats/map_pathsII.txt
        done
        scripts/annotate.py output/stats/bed_pathsII.txt output/stats/assembly_pathsII.txt output/stats/GenesII.txt
        touch {output.done}
        """

rule report:
    input:
        rules.annotation_statsII.output
    output:
        "output/reports/"+config["report_prefix"]+"_report/report.html"
    params:
        prefix = config["report_prefix"],
        wd = os.getcwd()
    singularity:
        "docker://reslp/rmarkdown:4.0.3"
    shell:
        """     
        # gather bedfiles
        # bedfiles all have the same name, therfore this hack to find, rename and copy the files to a location in the report directory.
        mkdir -p output/reports/{params.prefix}_report/bedfiles
        for f in $(cat output/stats/bed_pathsII.txt)
        do
            name=$(echo $f | awk -F/ '{{print $6}}')
            cp $f output/reports/{params.prefix}_report/bedfiles/$name.bed
        done
        
        # gather assemblies assemblies and rename sequences for report
        mkdir -p output/reports/{params.prefix}_report/assemblies
        cp $(cat output/stats/assembly_pathsII.txt) output/reports/{params.prefix}_report/assemblies/
        for f in $(ls -1 output/reports/{params.prefix}_report/assemblies/*)
	do
            header=$(basename $f | cut -d "." -f 1-3)
            echo $header
            if [ ! -z "$header" ]
            then
                sed -i "s#>.*#>$header#g" output/reports/{params.prefix}_report/assemblies/$header.final.fasta
            else
                echo $header
                break 
            fi
        done
        
        # gather maps   
        mkdir -p output/reports/{params.prefix}_report/maps
        cp $(cat output/stats/map_pathsII.txt) output/reports/{params.prefix}_report/maps
                
        # copy genes file
        cp output/stats/GenesII.txt output/reports/{params.prefix}_report/GenesII.txt
                        
        # create report
        sed 's/PREFIX/{params.prefix}/' ./scripts/report.Rmd > output/reports/{params.prefix}_report/report.Rmd 
	cd output/reports/{params.prefix}_report/
        Rscript -e 'rmarkdown::render("./report.Rmd")'
        cd -

        # clean up
        tar -pcf {params.wd}/output/reports/{params.prefix}_report.tar -C {params.wd}/output/reports {params.prefix}_report
        """     
