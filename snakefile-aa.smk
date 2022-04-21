#Program:Iterative scrutiny pipeline for EVRD-aa


########################################################################
##Note:This script consists of 4 parts in total                       ##
##Part I.Host genomes iterative scrutiny                              ##
##Part II.Vector sequence iterative scrutiny                          ##
##Part III.Annotation cross scrutiny                                  ##
##Part IV.Viral metagenome cross scrutiny                             ##
##Finally,you need to review the final documents produced by each part##
########################################################################

sra_id_lst = ["Avians","bat","bovine","human","pig","Rodents","Tick"]

rule all:
    input:
        "Host_scrutiny_aa/Uniprot-2-host-90-100.blastpout",
        "Host_scrutiny_aa/subject-seq-2-nr.blastpout",
        "Vector_scrutiny_aa/Uniprot-2-NVPC-90-100.blastpout",
        "Cross_validation_aa/Uniprot-2-uniprot-e50-100.blastpout",
        "Cross_validation_aa/error_seq_comment",
        expand("Viral_metagenome_aa/{SRA}.megahit/{SRA}-1000.fas",SRA = sra_id_lst),
        expand("Viral_metagenome_aa/{SRA}.blastxout",SRA = sra_id_lst),
        expand("Viral_metagenome_aa/{SRA}-uniq.txt",SRA = sra_id_lst),
        expand("Viral_metagenome_aa/{SRA}-all.txt",SRA = sra_id_lst),
        "Viral_metagenome_aa/Metagenome-subject-seq.fasta",
        "Viral_metagenome_aa/subject-seq-2-nr.blastpout"





diamond_db = 'NA_db/uniprot-viruses_Euk_20211025.dmnd'
diamond_nr = '/biostack/database/nr/nr_20191213'
diamond_NVPC = '/project/cjj/test/NVPC.dmnd'

##################################################################################################################
###########Part I.Host genomes iterative scrutiny                                                   ##############
###########Host genome were used to blastp search against UniProt database                          ##############
###########maximum of 1000 subjects to show alignments(length ≥ 100 and identity ≥ 90%)             ##############
##################################################################################################################

rule host_diamond:
    input:
        "host.faa"
    output:
        temp("Host_scrutiny_aa/Uniprot-2-host.blastpout")
    shell:
        "diamond blastp -q {input} -k 1000 -p 36 -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -d {diamond_db} -o {output}"


rule host_awk:
    input:
        "Host_scrutiny_aa/Uniprot-2-host.blastpout"
    output:
        "Host_scrutiny_aa/Uniprot-2-host-90-100.blastpout"
    shell:
        "awk '$3>=90' {input} | awk '$7>=100' > {output}"


rule host_subject:
    input:
        "Host_scrutiny_aa/Uniprot-2-host-90-100.blastpout"
    output:
        "Host_scrutiny_aa/subject-seq.fasta"
    shell:
        "awk '{{print $2 \"|\" $10 \"_\" $11 \"\\t\" $13}}' {input} | \
         seqkit tab2fx  > {output} " 

rule host_sed:
    input:
        "Host_scrutiny_aa/subject-seq.fasta"
    output:
        temp("Host_scrutiny_aa/subject-seq-sed.fasta")
    shell:
        "sed 's/-//g' {input} > {output}"


########################################################################################################################
######The aligned sequences of subject were extracted and subjected to blastp search against nr database ###############
########################################################################################################################

rule host_diamond2nr:
    input:
        "Host_scrutiny_aa/subject-seq-sed.fasta"
    output:
        "Host_scrutiny_aa/subject-seq-2-nr.blastpout"
    shell:
        "diamond blastp -q {input}  -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -k 100 -p 36 -d {diamond_nr} -o {output}"



##################################################################################################################
###########Part II.Vector sequence iterative scrutiny                                               ##############
###########non-viral protein core were used to blastp search against UniProt database               ##############
###########maximum of 1000 subjects to show alignments(length ≥ 100 and identity ≥ 90%)             ##############
##################################################################################################################

rule vector_diamond:
    input:
        "uniprot-viruses_Euk_20211025.fas"
    output:
        "Vector_scrutiny_aa/Uniprot-2-NVPC.blastpout"
    shell:
        "diamond blastp -q {input}  -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -k 1000 -p 36 -d {diamond_NVPC} -o {output}"

rule vector_awk:
    input:
        "Vector_scrutiny_aa/Uniprot-2-NVPC.blastpout"
    output:
        "Vector_scrutiny_aa/Uniprot-2-NVPC-90-100.blastpout"
    shell:
        "awk '$3>=90' {input} | awk '$7>=100' > {output}"


##################################################################################################################
###########Part III.Annotation cross scrutiny                                                       ##############
###########UniProt database were used to blastp search against UniProt database                     ##############
###########maximum of 1000 subjects to show alignments(evalue ≤ 1e-50 and length ≥ 100)             ##############
##################################################################################################################

rule cross_diamond:
    input:
        "uniprot-viruses_Euk_20211025.fas"
    output:
        temp("Cross_validation_aa/Uniprot-2-uniprot.blastpout")
    shell:
        "diamond blastp -q {input}  -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -k 1000 -p 36 -d {diamond_db} -o {output}"


rule cross_awk:
    input:
        "Cross_validation_aa/Uniprot-2-uniprot.blastpout"
    output:
        "Cross_validation_aa/Uniprot-2-uniprot-e50-100.blastpout"
    shell:
        "awk '$4<=1e-50' {input} | awk '$7>=100' > {output}"


rule cross_query:
    input:
        "Cross_validation_aa/Uniprot-2-uniprot-e50-100.blastpout"
    output:
        temp("Cross_validation_aa/query")
    shell:
        "awk '{{print $1}}' {input} | awk -F '|' '{{print $1\",\"$2}}' | \
        awk -F ',' '{{print $1\"\t\"$6}}' > {output}"



rule cross_subject:
    input:
        "Cross_validation_aa/Uniprot-2-uniprot-e50-100.blastpout"
    output:
        temp("Cross_validation_aa/subject")
    shell:
        "awk '{{print $2}}' {input} | awk -F '|' '{{print $1\",\"$2}}' |  \
        awk -F ',' '{{print $1\"\t\"$6}}' > {output}"


rule cross_paste:
    input:
        "Cross_validation_aa/query",
        "Cross_validation_aa/subject"
    output:
        "Cross_validation_aa/query-subject"
    shell:
        "paste {input[0]} {input[1]} > {output}"

##################################################################################################################
###########Find the different family annotated                                                      ##############
##################################################################################################################

rule error_seq:
    input:
        "Cross_validation_aa/query-subject"
    output:
        "Cross_validation_aa/error_seqid"
    shell:
        "awk '$2!=$4' {input} | awk '$4!=\"\"' | \
        awk '{{print $1}}' | sort -u > {output}"


rule corss_py:
    input:
        "Cross_validation_aa/query-subject",
        "Cross_validation_aa/error_seqid"
    output:
        "Cross_validation_aa/error_seq_comment"
    shell:
        "python while.py {input[0]} {input[1]} {output}"


##################################################################################################################
###########Part IV.Viral metagenome cross scrutiny                                                  ##############
###########SRA were used to blastx search against UniProt database                                  ##############
###########maximum of 1000 subjects to show alignments(length ≥ 100 and identity ≥ 80%)             ##############
##################################################################################################################

rule sra_seqkit:
    input:
        "Viral_metagenome_aa/{SRA}.megahit/final.contigs.fa"
    output:
        "Viral_metagenome_aa/{SRA}.megahit/{SRA}-1000.fas"
    shell:
        "seqkit seq -m 1000 {input} > {output}"

rule sra_diamond:
    input:
        "Viral_metagenome_aa/{SRA}.megahit/{SRA}-1000.fas"
    output:
        temp("Viral_metagenome_aa/{SRA}.blastxout")
    shell:
        "diamond blastx -q {input}  -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -k 1000 -p 36 -d {diamond_db} -o {output}"

rule sra_awk:
    input:
        "Viral_metagenome_aa/{SRA}.blastxout"
    output:
        temp("Viral_metagenome_aa/{SRA}-uniq.txt")
    shell:
        "awk '$3>=80' {input} | awk '$7 >= 100' \
        | awk '{{print $2}}' |  sort -u > {output}"


rule cat:
    input:
        expand("Viral_metagenome_aa/{SRA}-uniq.txt", SRA = sra_id_lst)
    output:
        "Viral_metagenome_aa/alluniq.txt"
    shell:
        "cat {input} | sort | uniq -c | awk '$1 >=2' | \
        awk '{{print $2}}' > {output}"


rule sra_sort:
    input:
        "Viral_metagenome_aa/{SRA}-uniq.txt",
        "Viral_metagenome_aa/alluniq.txt"
    output:
        "Viral_metagenome_aa/{SRA}-all.txt"
    shell:
        "sort {input[0]} {input[1]} | uniq -d > {output}"


rule sra_NF:
    input:
        "Viral_metagenome_aa/{SRA}.blastxout"
    output:
        "Viral_metagenome_aa/{SRA}-80-100.blastxout"
    shell:
        "awk '$3>=80' {input} | awk '$7 >= 100' | \
        awk -F ' ' '{{for(i=2;i<=NF;++i)printf $i \"\\t\";printf\"\\n\"}}' > {output}"


rule sra_py:
    input:
        "Viral_metagenome_aa/{SRA}-80-100.blastxout",
        "Viral_metagenome_aa/{SRA}-all.txt"
    output:
        "Viral_metagenome_aa/{SRA}-all.blastxout"
    shell:
        "python while.py {input[0]} {input[1]} {output}" 

rule sra_all:
    input:
        expand("Viral_metagenome_aa/{SRA}-all.blastxout", SRA = sra_id_lst)
    output:
        "Viral_metagenome_aa/all.blastxout"
    shell:
        "cat {input} > {output}"

rule sra_sed:
    input:
        "Viral_metagenome_aa/all.blastxout"
    output:
        "Viral_metagenome_aa/Metagenome-subject-seq.fasta"
    shell:
        "awk '{{print $1 \"|\" $9 \"_\" $10 \"\t\" $12}}' {input} | \
        seqkit tab2fx | sed 's/-//g' > {output}"

########################################################################################################################
######The aligned sequences of subject were extracted and subjected to blastp search against nr database ###############
########################################################################################################################

rule sra_diamond2nr:
    input:
        "Viral_metagenome_aa/Metagenome-subject-seq.fasta"
    output:
        "Viral_metagenome_aa/subject-seq-2-nr.blastpout"
    shell:
        "diamond blastp -q {input} -d {diamond_nr} -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -k 100 -p 36 -o {output}"

