#Program:Iterative scrutiny pipeline for EVRD-nt

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
        "Host_scrutiny_nt/gbvrl-2-host-85-150.blastnout",
        "Host_scrutiny_nt/subject-seq-2-nt.blastnout",
        "Vector_scrutiny_nt/gbvrl-2-NVPC-99-60.blastxout",
        "Vector_scrutiny_nt/gbvrl-2-UniVec-85-150.blastnout",
        "Vector_scrutiny_nt/subject-seq-2-nt.blastnout",
        "Cross_validation_nt/gbvrl-2-gbvrl-e50-500.blastnout",
        "Cross_validation_nt/error_seq_comment",
        expand("Viral_metagenome_nt/{SRA}.megahit/{SRA}-1000.fas",SRA = sra_id_lst),
        expand("Viral_metagenome_nt/{SRA}.blastnout",SRA = sra_id_lst),
        expand("Viral_metagenome_nt/{SRA}-uniq.txt",SRA = sra_id_lst),
        expand("Viral_metagenome_nt/{SRA}-all.txt",SRA = sra_id_lst),
        "Viral_metagenome_nt/Metagenome-subject-seq.fasta",
        "Viral_metagenome_nt/subject-seq-2-nt.blastnout"






Blast_db = 'NT_db/gbvrl_genome_only_euk_host.200.index'
Blast_nt = '/biostack/database/nt/nt'
diamond_NVPC = '/project/cjj/test/NVPC.dmnd'


##################################################################################################################
###########Part I.Host genomes iterative scrutiny                                                   ##############
###########Host genome were used to blastn search against Genbank database                          ##############
###########maximum of 1000 subjects to show alignments(length ≥ 150 and identity ≥ 85%)             ##############
##################################################################################################################

rule host_blastn:
    input:
        "host1.fasta"
    output:
        temp("Host_scrutiny_nt/gbvrl-2-host.blastnout")
    shell:
        "blastn -query {input} -db {Blast_db} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"


rule host_awk:
    input:
        "Host_scrutiny_nt/gbvrl-2-host.blastnout"
    output:
        "Host_scrutiny_nt/gbvrl-2-host-85-150.blastnout"
    shell:
        "awk '$3>=85' {input} | awk '$7>=150' | grep -v \'Retroviridae\' > {output}"


rule host_subject:
    input:
        "Host_scrutiny_nt/gbvrl-2-host-85-150.blastnout"
    output:
        "Host_scrutiny_nt/subject-seq.fasta"
    shell:
        "awk '{{print $2 \"|\" $10 \"_\" $11 \"\\t\" $13}}' {input} | \
         seqkit tab2fx  > {output} " 

rule host_sed:
    input:
        "Host_scrutiny_nt/subject-seq.fasta"
    output:
        temp("Host_scrutiny_nt/subject-seq-sed.fasta")
    shell:
        "sed 's/-//g' {input} > {output}"


########################################################################################################################
######The aligned sequences of subject were extracted and subjected to blastn search against nt database ###############
########################################################################################################################

rule host_blast2nt:
    input:
        "Host_scrutiny_nt/subject-seq-sed.fasta"
    output:
        "Host_scrutiny_nt/subject-seq-2-nt.blastnout"
    shell:
        "blastn -query {input} -db {Blast_nt} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"


##################################################################################################################
###########Part II.Vector sequence iterative scrutiny                                               ##############
###########non-viral protein core were used to blastx search against Genbank database               ##############
###########maximum of 1000 subjects to show alignments(length ≥ 60 and identity ≥ 99%)              ##############
##################################################################################################################

rule vector_diamond:
    input:
        "gbvrl_genome_only_euk_host.200.fasta"
    output:
        "Vector_scrutiny_nt/gbvrl-2-NVPC.blastxout"
    shell:
        "diamond blastx -q {input}  -f 6 \
        qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq \
        -k 1000 -p 36 -d {diamond_NVPC} -o {output}"

rule vector_awk:
    input:
        "Vector_scrutiny_nt/gbvrl-2-NVPC.blastxout"
    output:
        "Vector_scrutiny_nt/gbvrl-2-NVPC-99-60.blastxout"
    shell:
        "awk '$3>=99' {input} | awk '$7>=60' > {output}"

rule vector_blastn:
    input:
        "UniVec.fasta"
    output:
        temp("Vector_scrutiny_nt/gbvrl-2-UniVec.blastnout")
    shell:
        "blastn -query {input} -db {Blast_db} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"

rule UniVec_awk:
    input:
        "Vector_scrutiny_nt/gbvrl-2-UniVec.blastnout"
    output:
        "Vector_scrutiny_nt/gbvrl-2-UniVec-85-150.blastnout"
    shell:
        "awk '$3>=85' {input} | awk '$7>=150' > {output}"

rule vector_subject:
    input:
        "Vector_scrutiny_nt/gbvrl-2-UniVec-85-150.blastnout"
    output:
        "Vector_scrutiny_nt/subject-seq.fasta"
    shell:
        "awk '{{print $2 \"|\" $10 \"_\" $11 \"\\t\" $13}}' {input} | \
         seqkit tab2fx  > {output}"

rule vector_sed:
    input:
        "Vector_scrutiny_nt/subject-seq.fasta"
    output:
        temp("Vector_scrutiny_nt/subject-seq-sed.fasta")
    shell:
        "sed 's/-//g' {input} > {output}"


########################################################################################################################
######The aligned sequences of subject were extracted and subjected to blastn search against nt database ###############
########################################################################################################################


rule blast2nt:
    input:
        "Vector_scrutiny_nt/subject-seq-sed.fasta"
    output:
        "Vector_scrutiny_nt/subject-seq-2-nt.blastnout"
    shell:
        "blastn -query {input} -db {Blast_nt} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"


##################################################################################################################
###########Part III.Annotation cross scrutiny                                                       ##############
###########Genbank database were used to blastn search against Genbank database                     ##############
###########maximum of 1000 subjects to show alignments(evalue ≤ 1e-50 and length ≥ 500)             ##############
##################################################################################################################

rule cross_blastn:
    input:
        "gbvrl_genome_only_euk_host.200.fasta"
    output:
        temp("Cross_validation_nt/gbvrl-2-gbvrl.blastnout")
    shell:
        "blastn -query {input} -db {Blast_db} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"


rule cross_awk:
    input:
        "Cross_validation_nt/gbvrl-2-gbvrl.blastnout"
    output:
        "Cross_validation_nt/gbvrl-2-gbvrl-e50-500.blastnout"
    shell:
        "awk '$4<=1e-50' {input} | awk '$7>=500' > {output}"


rule cross_query:
    input:
        "Cross_validation_nt/gbvrl-2-gbvrl-e50-500.blastnout"
    output:
        temp("Cross_validation_nt/query")
    shell:
        "awk '{{print $1}}' {input} | awk -F '|' '{{print $1\",\"$2}}' | \
        awk -F ',' '{{print $1\"\t\"$6}}' > {output}"


rule cross_subject:
    input:
        "Cross_validation_nt/gbvrl-2-gbvrl-e50-500.blastnout"
    output:
        temp("Cross_validation_nt/subject")
    shell:
        "awk '{{print $2}}' {input} | awk -F '|' '{{print $1\",\"$2}}' |  \
        awk -F ',' '{{print $1\"\t\"$6}}' > {output}"

rule cross_paste:
    input:
        "Cross_validation_nt/query",
        "Cross_validation_nt/subject"
    output:
        "Cross_validation_nt/query-subject"
    shell:
        "paste {input[0]} {input[1]} > {output}"

rule error_seq:
    input:
        "Cross_validation_nt/query-subject"
    output:
        "Cross_validation_nt/error_seqid"
    shell:
        "awk '$2!=$4' {input} | awk '$4!=\"\"' | \
        awk '{{print $1}}' | sort -u > {output}"


rule corss_py:
    input:
        "Cross_validation_nt/query-subject",
        "Cross_validation_nt/error_seqid"
    output:
        "Cross_validation_nt/error_seq_comment"
    shell:
        "python while.py {input[0]} {input[1]} {output}"


##################################################################################################################
###########Part IV.Viral metagenome cross scrutiny                                                  ##############
###########SRA were used to blastn search against Genbank database                                  ##############
###########maximum of 1000 subjects to show alignments(length ≥ 150 and identity ≥ 80%)             ##############
##################################################################################################################

rule SRA_seqkit:
    input:
        "Viral_metagenome_nt/{SRA}.megahit/final.contigs.fa"
    output:
        "Viral_metagenome_nt/{SRA}.megahit/{SRA}-1000.fas"
    shell:
        "seqkit seq -m 1000 {input} > {output}"

rule SRA_blast:
    input:
        "Viral_metagenome_nt/{SRA}.megahit/{SRA}-1000.fas"
    output:
        temp("Viral_metagenome_nt/{SRA}.blastnout")
    shell:
        "blastn -query {input} -db {Blast_db} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"

rule sra_awk:
    input:
        "Viral_metagenome_nt/{SRA}.blastnout"
    output:
        temp("Viral_metagenome_nt/{SRA}-uniq.txt")
    shell:
        "awk '$3>=80' {input} | awk '$7 >= 150' \
        | awk '{{print $2}}' |  sort -u > {output}"


rule sra_cat:
    input:
        expand("Viral_metagenome_nt/{SRA}-uniq.txt", SRA = sra_id_lst)
    output:
        "Viral_metagenome_nt/alluniq.txt"
    shell:
        "cat {input} | sort | uniq -c | awk '$1 >=2' | \
        awk '{{print $2}}' > {output}"


rule sra_sort:
    input:
        "Viral_metagenome_nt/{SRA}-uniq.txt",
        "Viral_metagenome_nt/alluniq.txt"
    output:
        "Viral_metagenome_nt/{SRA}-all.txt"
    shell:
        "sort {input[0]} {input[1]} | uniq -d > {output}"


rule sra_NF:
    input:
        "Viral_metagenome_nt/{SRA}.blastnout"
    output:
        "Viral_metagenome_nt/{SRA}-80-150.blastnout"
    shell:
        "awk '$3>=80' {input} | awk '$7 >= 150' | \
        awk -F ' ' '{{for(i=2;i<=NF;++i)printf $i \"\\t\";printf\"\\n\"}}' > {output}"


rule sra_py:
    input:
        "Viral_metagenome_nt/{SRA}-80-150.blastnout",
        "Viral_metagenome_nt/{SRA}-all.txt"
    output:
        "Viral_metagenome_nt/{SRA}-all.blastnout"
    shell:
        "python while.py {input[0]} {input[1]} {output}" 

rule sra_all:
    input:
        expand("Viral_metagenome_nt/{SRA}-all.blastnout", SRA = sra_id_lst)
    output:
        "Viral_metagenome_nt/all.blastnout"
    shell:
        "cat {input} > {output}"

rule sra_seqkit:
    input:
        "Viral_metagenome_nt/all.blastnout"
    output:
        "Viral_metagenome_nt/Metagenome-subject-seq.fasta"
    shell:
        "awk '{{print $1 \"|\" $9 \"_\" $10 \"\t\" $12}}' {input} | \
        seqkit tab2fx | sed 's/-//g' > {output}"

########################################################################################################################
######The aligned sequences of subject were extracted and subjected to blastn search against nt database ###############
########################################################################################################################

rule sra_blast_nt:
    input:
        "Viral_metagenome_nt/Metagenome-subject-seq.fasta"
    output:
        "Viral_metagenome_nt/subject-seq-2-nt.blastnout"
    shell:
        "blastn -query {input} -db {Blast_nt} -outfmt \
        \"6 qseqid stitle pident evalue qlen slen length qstart qend sstart send qseq sseq\" \
         -num_alignments 1000 -num_threads 36 -out {output}"

