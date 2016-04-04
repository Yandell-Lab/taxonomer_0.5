#taxonomer
=========

##metagenomics toolkit 

###Building the code:

This repository is written in c and python, with cython bringing the two together.  

python dependencies -- we HIGHLY recommend using the anaconda distribution of python:

1. Python 2.7+
2. cython
3. biopython

To build, cd into taxonomer/ -- directory where all the code lives, type the following commands:

1. `rm scripts/cython/*.so`
2. `rm scripts/cython/*.c`
3. `python setup.py build_ext --inplace`

####Building a nucleotide database:

To build a database, you must use kanalyze to create the kmer file. You do this with kanalyze, which
is available here: http://sourceforge.net/projects/kanalyze/files/

The following is a command line example of kanalyze with the necessary paramaters (-k is the kmer size):      
  
`java -jar kanalyze.jar count -k 31 -d 30 -l 30 -o \<output file\> -m hex -rcanonical -f fasta \<input fasta\>`

Then use build_db.py to construct the database for classification

To see help, use python build_db.py -h

####Building a protein database:

To build a protein database, you must follow these steps:

1.  Use convert_protein_db.py to do necessary conversion before using kanalyze
2.  Run kanalyze to create kmer file, but with the following arguments (-k is the kmer size, no -rcanonical):  
  `java -jar kanalyze.jar count -k 30 -d 30 -l 30 -o \<output file\> -m hex -f fasta \<input fasta\>`
3.  Use python build_db.py, but specify --protein 1 to build protein database

####Classifying reads

Use classify_reads.py to classify reads.  use -h argument to see command line options.  
NOTE -- only use --protein 1 if you have already built a database for protein classification.

####Binner database construction

We have empirically found binning using k-mers of 21bp in length to be effective for our research applications.  This may or may not be the optimal size for your application, but the commands that follow demonstrate how to construct a binner database using 21bp k-mers.  

1.  Create k-mer count files using kanalze of fasta files of reference nucleotide sequences:
   
 `java -jar kanalyze.jar count -k21 -d 30 -l 30 -o \<output file\> -m hex -f fasta \<input fasta\>`

2.  Merge all k-mer count files using the script binner_merge_kanalyze.py.  The help shows the following: 
    
    positional arguments:
    
    kmer_counts_1   name of kmer count file in kanalyze hex output format
    kmer_counts_2   name of kmer count file to merge, must also be in kanalyze hex output format
    output          merged output file

    optional arguments:
    
    -h, --help     show this help message and exit
    -id_1 ID_1     integer value used to identify first database, must be a multiple of 2
    -id_2 ID_2     integer value used to identify second databse, must be a multiple of 2
  
  The first kmer count file, kmer_counts_1 will have kmer_counts_2 merged into it.  Only use the -id_1 or -id_2 tags if the   corresponding k-mer count file (kmer_counts_1 and/or kmer_counts_2) haven't already been merged using this script.  The reason is that the flags -id_1, -id_2 will overwrite the id field in the file.
  
  For example, suppose we have 3 kmer count files created by kanalyze to merge: k1.kc, k2.kc, k3.kc
  First, we determine the integer ids to assign to each kmer count file(must be a power of 2) -- n our example let k1.kc, k2.kc, k3.kc be 1, 2, 4 respectively.  The commands to achieve this are (only two files can be merged at once):
  
  `python binner_merge_kanalyze.py k1.kc k2.kc k1_k2.merged.kc -id_1 1 -id_2 2`
  
  `python binner_merge_kanalyze.py k1_k2.merged.kc k3.kc all_k.merged.kc -id_2 4`
  
  Now the merged file is called all_k.merged.kc
  
3. Build our example database using :
  
  `python build_db.py example_binner_db all_k.merged.kc fake.sti fake.fasta -kl 21 -binner`

Where example_binner_db is the prefix for the output files that enable running the binner.  fake.sti and fake.fasta are placeholders and are in the scripts/ folder.  They can be used in any binner database build.  



