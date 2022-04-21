rule all:
    input:
        "NA_db/uniprot-viruses_Euk_20211025.dmnd"

rule makedb:
    input:
        "uniprot-viruses_Euk_20211025.fas"
    output:
        "NA_db/uniprot-viruses_Euk_20211025.dmnd"
    log:
        "AA_db/db.log"
    shell:
        "diamond makedb --in {input} -d {output}"

