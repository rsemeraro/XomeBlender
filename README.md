# Xome-Blender
Xome-Blender is a collection of bash, R and C++ scripts based on SAMtools, GATK, Picard and VcfTools functions that allows to generate synthetic cancer genomes with user defined features such as the number of subclones, the number of somatic variants and the presence of CNV, without the addition of any synthetic element. It is composed of two modules: InXalizer and Xome-Blender. The first module is devoted to the blending process initialization. It takes as input a single BAM file, a set of user-defined parameters and returns the coverage of the sample and the input-files for the second module (Xome-Blender). Optionally, it creates a file containing the coordinates to insert CNV in the final product.
The second module generates the synthetic heterogeneous sample.
###### Supported on Linux and Mac OS X.

## Requirements 
* R (https://www.r-project.org)
* SAMtools 1.1 or above (http://www.htslib.org/download/)
* GATK 3.3 or above (https://www.broadinstitute.org/gatk/)
* VcfTools (https://sourceforge.net/projects/vcftools/)
* Picard (https://broadinstitute.github.io/picard/)

## Installation

* ### Clone the Xome-Blender repository
    
        git clone git://github.com/rsemeraro/Xome-Blender.git
* ### Compile
    The project can be compiled by calling make in the top-level directory:    

        make

    ###### Check the deps.txt file after make. If any link to a dependency is wrong, substitute it with the right one.
## Usage
First run InXalizer. It requires a BAM file, a label for it and a reference genome. By means of four parameters it is possible to tune the initialization process: 
 1. Subclone number = the number of subclones that will compose the final product (```-scn```).
 2. Variants number = the number of somatic variants that will appear in the final product (```-vn```).
 3. Subclonal architecture = the evolution model for the sample synthesis, it can be Linear or Branched (```-sa```).
 4. CNV = it's an option that allows for the generation of a CNV file, defining their number and length (```-c```).

Optionally, it is possible to use a target file, in bed format, to edit only defined portions of the BAM file (whole-exome or target sequencing experiments).

Examples:
* Running InXalizer with minimum requirements:

        InXalizer -f file.bam -l my_label -r MyRef.fasta -scn 2 -vn 50 -sa Branched
        
* Running InXalizer with target file:

        InXalizer -f file.bam -l my_label -r MyRef.fasta -scn 2 -vn 50 -sa Branched -b MyTargetfile.bed

* Running InXalizer with CNVs:

        InXalizer -f file.bam -l my_label -r MyRef.fasta -scn 2 -vn 50 -sa Branched -c 3,1000000

InXalizer can be run to generate CNV files only, by using the *"CNV list"*  function, invoked as ```-cl```, and a reference sequence.

        InXalizer -r ref.fa -cl cnv_list.txt
        
Each line of the ```cnv_list``` file must contain a label for the cnv file to be generated, the number and the size of CNV events and, optionally, the path to a target file.
  ##### List_file example:
      My_label_1	23,10000000
      My_label_2    3,100000000 path/to/target.bed

The data generated by InXalizer are then used to proceed with the blending process, performed by Xome-Blender.
The main module takes as input the BAM files produced by InXalizer, the sample's label (the same used in InXalizer), the starting coverage (```-sc```) (stored into the *.cov*  file), the percentage (```-p```) of each subclone in order to create the mixed sample, the desired coverage (```-fc```) of the output file and the subclone variant file (```-v```) generated in the previous step. If provided, it can use a cnv file generated with InXalizer to add CNV events to the final product. <br /> Xome-Blender can be run in *_"single mode"_* by typing all the options in a shell, as in the examples below.

Examples:
* Running Xome-Blender with minimum requirements:

        Xome-Blender -f Control.bam,Subclone1.bam,Subclone2.bam -la my_label -sc 127 -p 30,40,30 -fc 110 -v Subclone1.vcf,Subclone2.vcf
* Running Xome-Blender with CNVs:

        Xome-Blender -f Control.bam,Subclone1.bam,Subclone2.bam -la my_label -sc 127 -p 30,40,30 -fc 110 -v Subclone1.vcf,Subclone2.vcf -cnv my_label_CNV.txt
    ######  \*Multiple elements must be comma separated. <br />
  
Alternatively, it can be run in *_"automated mode"_* by using the ```--list (-l)``` option. The activation of this parameter requires a list file, containing the info to run multiple consecutive analyses.

      Xome-Blender -l list_file.txt
   The list_file is a tab separated file containing different anlayses (one per row). Each row must contain all the parameters above.
  ##### List_file example
      NA18501_Control.bam,NA18501_Subclone1.bam	NA18501	145	20,80   120 -v Subclone1.vcf
      NA18501_Control.bam,NA18501_Subclone1.bam	NA18501	145	30,70	90  -v Subclone1.vcf
      NA18501_Control.bam,NA18501_Subclone1.bam	NA18501	145	40,60	140 -v Subclone1.vcf
      NA18501_Control.bam,NA18501_Subclone1.bam	NA18501	145	50,50	50  -v Subclone1.vcf

## Contacts

This program has been developed by Roberto Semeraro, Department of Experimental and Clinical Medicine, University of Florence
