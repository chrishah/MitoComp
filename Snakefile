import os
import pandas as pd

include: "rules/functions.smk"
include: "rules/download.smk"
include: "rules/trimming.smk"
include: "rules/subsample.smk"
include: "rules/norgal.smk"
include: "rules/getorganelle.smk"
include: "rules/mitoflex.smk"
include: "rules/mitoz.smk"
include: "rules/novoplasty.smk"
include: "rules/mitobim.smk"
include: "rules/annotation.smk"
include: "rules/alignment.smk"
include: "rules/annotationII.smk"
include: "rules/CCT.smk"
include: "rules/report.smk"


rule all:
    input:
        "output/reports/"+config["report_prefix"]+"_report/report.html"

rule assembly_only:
    input:
        pick_assembly

rule gather:
    input:
        trigger_gather
    output:
        "output/gathered_assemblies/gathered.done"
    log:
        "output/gathered_assemblies/gathered.log"
    singularity: "docker://chrishah/mitobim:v.1.9.1" #this is only to make sure that rsync is stable
    shell:
        """
        rsync -avpuzP -L {input} output/gathered_assemblies/ &> {log}
        touch {output}
        """
