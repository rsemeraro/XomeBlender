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
First run InXalizer. It requires the BAM files, their ID (to avoid the use of long file names), a reference genome file and a label for the germline variant file. Moreover, it is possible to load a target file, in BED format, to work on exomic data.
##### Example
    InXalizer -f file1,file2,file3 -i Id1,Id2,Id3 -l My-favurite-label -r MyRef.fasta -b MyTarget.bed
*The elements must be comma separated    
