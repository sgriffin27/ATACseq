# ATACseq
modified from AD_MapATACseq
#Loading files and making constructs
(1) follow instructions in lewis lab repository to load files (you will need you ATAC sequences and their reference genomes (located in config_bwt))
(2) transfer files onto your computer following the instructions in the lewis lab repository 
(3) configure bowtie2 files of the genomes using:
    module load Bowtie2/2.5.2-GCC-11.3.0

    bowtie2-build [options]* <reference_in> <bt2_base>

    #where reference_in is your reference genome FASTA file and bt2_base is the basename for the genome.
(4) transfer onto computer from the cluster

#Analyzing with IGV
(1) get ready to analyze on igv
(2) load fai files
(3) load gff files
(4) load big wig files
(5) when done with igv, you can save the session 

#Mapping ATAC reads
(1)
