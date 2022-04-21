rule all:
    input:
        "NT_db/gbvrl_genome_only_euk_host.200.index"

rule makeblastdb:
    input:
        "gbvrl_genome_only_euk_host.200.fasta"
    output:
        "NT_db/gbvrl_genome_only_euk_host.200.index"
    params:
        tp="nucl"
    log:
        "NT_db/db.log"
    shell:
        "makeblastdb -in {input} -dbtype {params.tp} -out {output}"
