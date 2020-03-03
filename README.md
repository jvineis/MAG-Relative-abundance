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

Then we need to fix the names of the contigs to ensure that they are unique among all contigs in your collection of genomes.  like thus. but there are many other ways to do this.  

    for i in `cat genomes.txt`; do sed -i "s/c_/${i}_c/g" $i.fa; done
    
Now we concatenate all fasta files and build a bowtie2 database

    cat *.fa > ALL-NON-REDUNDANT-MAGS.fa
    bowtie2-build ALL-NON-REDUNDANT-MAGS.fa ALL-NON-REDUNDANT-MAGS
    
If you have a cluster handy... Great!  Lets do this.  If not... I think you might need one.  This is the command to map each separate fastq onto your concatenated fasta of all the MAGs

    for i in `cat samples.txt`; do clusterize bowtie2 --very-sensitive -f -x ALL-NON-REDUNDANT-MAGS -U ../../../BOWEN-METAGENOMICS/"$i" -S "$i".sam; done
    
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
    
#### This is another iteration of the same analysis.. with a few exceptions.  I did not use the --very-sensitive flag in bowtie2 or the -q flag. These flags were a little too stringent. I belive this because when mapping to my collection of non-redundant MAGs, some of them did not rec coverage or there was very low coverage due to the inability to map properly to more than one location and the likely presence of naturally occurring SNPs in the population. 

##### 1. Calculate ANI - assumes that you have already constructed anvio contig dbs for each of your MAGs and you are using an external_genomes.txt file. 
    clusterize anvi-compute-genome-similarity -e external_genomes.txt -o x_ANI -T 40

##### 2. Dereplicate the genomes based on the ANI calculation
    clusterize anvi-dereplicate-genomes --ani-dir x_ANI/ -o x_ANI_dereplication --program pyANI --method ANIb --use-full-percent-identity --min-full-percent-identity 0.90 --similarity-threshold 0.95

##### 3. Collect the names of the unique mags, copy the fasta files to a new location and then fix the deflines - this step will be very personalized based on the names of you fasta files, but you can probably make the append-name-to-fasta-deflines.py work for you 
    cut -f 3 CLUSTER_REPORT.txt > ../x_unique_mags.tx
    for i in `cat x_unique_mags.txt`; do cp $i'.fa' x_ANI_dereplication; done
    for i in `cat ../x_unique_mags.txt`; do python ~/scripts/append-name-to-fasta-deflines.py $i'.fa' $i'-modified.fa'; done
    
##### 4. Concatenate the fresh fasta files and build a bowtie2 database
    cat *-modified.fa > x_all-uniqe-concatenated.fa
    bowtie2-build x_all-uniqe-concatenated.fa x_all-uniqe-concatenated
    
##### 5. Map each of your datasets to the concatenated fasta
    for i in `cat /workspace/jvineis/BOWEN_METAGENOMICS/samples.txt`; do clusterize bowtie2 -f -x x_all-uniqe-concatenated -U '/workspace/jvineis/BOWEN_METAGENOMICS/'$i'.gz' -S $i'.sam'; done

##### 6. Use samtools to convert, filter, sort and index
    for i in `cat /workspace/jvineis/BOWEN_METAGENOMICS/samples.txt`; do clusterize samtools view $i'.sam' -b -o $i'.bam' -F 4 ;done
    for i in `cat /workspace/jvineis/BOWEN_METAGENOMICS/samples.txt`; do clusterize samtools sort $i'.bam' -o $i'-sorted.bam'; done
    for i in `cat /workspace/jvineis/BOWEN_METAGENOMICS/samples.txt`; do clusterize samtools index $i'-sorted.bam'; done
    
##### 7.  Generate the number of reads that were mapped to each of the sequences in the concatenated fasta
    for i in `cat /workspace/jvineis/BOWEN_METAGENOMICS/samples.txt`; do samtools idxstats $i'-sorted.bam' > $i'-contig-coverages.txt';done
    
##### 8.  Report the number of reads that mapped to each of your bins using the script below.  It might require some tweaking for it to work properly, but its a really simple script

    for i in `cat /workspace/jvineis/BOWEN_METAGENOMICS/samples.txt`; do python ~/scripts/calculate-mag-coverage.py -cov $i'-contig-coverages.txt' -mags x_MAGS-to-splits.txt -out $i'-mag-coverages.txt';done
    
##### 9.  At this point you are pretty close to wrapping this up.  I paste together all of the *mag-coverages.txt* files and then fix things up in excel (I know thats lame).  One handy way to correct your relative abundance for the number of reads per sample and the size of the MAG

    RelAbund = Reads in sample/Minimum reads in collection * (reads mapped*average read length)/genome-size
