# Xome-Blender
Xome-Blender is a collection of bash and R scripts based on SAMtools and GATK functions, useful to generate synthetic cancer genomes with different contamination level and intra-tumor heterogeneity and devoid of any synthetic element.
It is composed of two modules: InXalizer and Xome-Blender. The first prepare the data for the second module by calculating
the mean coverage and by detecting the germline variants of each input file. The second is main module that takes as input the sequencing data of two or more samples in BAM format and generates mixed samples with user-defined proportions and coverages.
###### Supported on Linux and Mac OS X.

## Requirements 
* R (https://www.r-project.org)
* SAMtools 0.1.16 or above (http://www.htslib.org/download/)
* GATK (https://www.broadinstitute.org/gatk/)

## Installation
#### Clone the Xome-Blender repository
    git clone git://github.com/rsemeraro/Xome-Blender.git
You can execute it as a simple bash script by typing ./Xome-Blender, or you can create a symbolic link to the files.

## Usage
* First run InXalizer. It requires the BAM files, their ID (to avoid the use of long file names), a reference genome file and a label for the germline variant file. Moreover, it is possible to load a target file, in BED format, to work on exomic data.

        InXalizer -f file1,file2,file3 -i Id1,Id2,Id3 -l My-favurite-label -r MyRef.fasta -b MyTarget.bed

  *Multiple elements must be comma separated.
* The data generated by InXalizer are then used to configure the Xome-Blender analysis.
The main module takes as input the BAM files, their ID, the coverages previously calculated, the percentage of each BAM in order to create the mixed sample, the desired coverage of the output file and the variant file generated in the previous step. <br /> It can be run in **single mode** by typing all the options in a shell, as in the example below.

        Xome-Blender -f file1,file2,file3 -i Id1,Id2,Id3 -c 127,138,90 -p 30,40,30 -tc 130 -v My-favurite-label.vcf
  \*Multiple elements must be comma separated. <br />
  \**The order of samples must be respected in each option!
  
  Alternatively, it can be run in **automated mode** by using the --list \(-l) option. The activation of this parameter requires a list file, by assigning this it's possible to run multiple consecutive analyses.

      Xome-Blender -l list_file.txt
   The list_file is a tab separated file containing different anlayses (one per row). Each row must contain all the options above.
  ##### List_file example
      NA18501.bam,NA12889.bam	NA18501,NA12889	145,193	20,80	200	My_Favourite_label_1.vcf
      NA18501.bam,NA12889.bam	NA18501,NA12889	145,193	30,70	150	My_Favourite_label_2.vcf
      NA18501.bam,NA12889.bam	NA18501,NA12889	145,193	40,60	100	My_Favourite_label_3.vcf
      NA18501.bam,NA12889.bam	NA18501,NA12889	145,193	50,50	180	My_Favourite_label_4.vcf 
