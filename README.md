# ATACseq
modified from AD_MapATACseq

(1) follow instructions in lewis lab repository to load files (you will need you ATAC sequences and their reference genomes (located in config_bwt))
(2) transfer files onto your computer following the instructions in the lewis lab repository 
(3) configure bowtie2 files of the genomes using:
    module load Bowtie2/2.5.2-GCC-11.3.0

    bowtie2-build [options]* <reference_in> <bt2_base>

    #where reference_in is your reference genome FASTA file and bt2_base is the basename for the genome.
(4) transfer onto computer from the cluster
(5) get ready to analyze on igv
(6) load fai files
(7) load gff files
(8) load big wig files
(9) when done with igv, you can save the session 
