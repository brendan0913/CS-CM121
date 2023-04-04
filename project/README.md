Usage: ./pa.py
Usage: ./pa.py reads_file transcriptome_file k    
    
    The psuedoalignment script, pa.py, can either run with no arugments or with 3 arguments.
The script creates a file called output.tsv that holds equivalence counts as specified in the
spec in tsv format.
    When run with no arguments, the k-mer length is 21, the reads file is "reads.fasta", and the
transcriptome file is "chr11_transcriptome.fasta". These two files must be in the same directory 
as the pa.py script.
    When run with 3 arguments, the first argument must be the reads file (in fasta format), 
the second argument must be the transcriptome file (in fasta fomrat), and the third argument 
must be the k-mer length, preferably between 21-36. 
    If more than 0 and less than 3 arguments are given or more than 3 arguments are given, 
an error message is printed and it exits.