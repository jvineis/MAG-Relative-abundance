# MAG-Relative-abundance
Steps to calculate relative abundance of MAGs derived from separate assemblies across samples
This is just one of the ways to do this... There are many others.. some are easier, some are harder, but you will get an idea of how it all works.

The process begins with a collection of MAGs that you have dereplicaed (or not) depending on what you see fit.
Create a list of the fastq files that you would like to map back to your collection of MAGs.  Maybe something like this, which you would execute in the directory containing the fasta files for each MAG.
    
       ls *.fa | sed 's/\.fa//g' > genomes.txt

You will also need a list of the fastq files that you would like to map back to the MAGs.  Should look something like this

    sample1.fastq
    sample2.fastq
    sample3.fastq

Then we need to fix the names of the contigs to ensure that they are unique among all contigs in your collection of genomes.  like thus

    or i in `cat genomes.txt`; do sed -i "s/c_/${i}_c/g" $i.fa; done
    
Now we concatenate all fasta files and build a bowtie2 database

    cat *.fa > ALL-NON-REDUNDANT-MAGS.fa
    bowtie2-build ALL-NON-REDUNDANT-MAGS.fa ALL-NON-REDUNDANT-MAGS
    
If you have a cluster handy... Great!  Lets do this.  If not... I think you might need one.  This is the command to map each separate fastq onto your concatenated fasta of all the MAGs

    for i in `cat samples.txt`; do clusterize bowtie2 --very-sensitive -x ALL-NON-REDUNDANT-MAGS -U ../../../BOWEN-METAGENOMICS/"$i" -S "$i".sam; done
    
In the world of Slurm management, the batch script might look like something like this if your fasta files are contained in a directory called QUALITY_READS.  The name of this script is "00_mapping_master.shx" which you can see called in the bash script used to execute this command.  Don't forget to specify the "-f" flag in your bowtie2 command if you are working with merged read fasta instead of paired fastq.

    #!/bin/bash
    #
    #SBATCH --nodes=2
    #SBATCH --tasks-per-node=8
    #SBATCH --mem=80Gb
    #SBATCH --partition=general

    bowtie2 --very-sensitive -x ALL-NON-REDUNDANT-MAGS -U QUALITY_READS/${READS} -S MAPPING/${READS}.sam
    
Which you would execute using a bash script that looks something like this.

    #!bin/bash

    for READS in `cat samples.txt`; do echo "${READS}"; export READS; sbatch 00_mapping_master.shx; sleep 1; done

    
Then you will want to filter the sam files and convert them to bam files, sort, index, and generate a report of the contig coverages.  Here is how to do this.

    for i in `cat samples.txt`; do samtools view -b -q 10 "$i".sam > "$i".bam; done
    for i in `cat samples.txt`; do samtools sort "$i".bam > "$i"-sorted.bam; done
    for i in `cat samples.txt`; do samtools index "$i"-sorted.bam;done
    for i in `cat samples.txt`; do samtools idxstats "$i"-sorted.bam > "$i"-contig-coverages.txt; done
    

