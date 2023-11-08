# 16s_XTrakTA
Python executable tool for extracting 16S rRNA sequences following genome annotation. 

This simple tool has been designed around the output of a Prokka genome assembly, with the aim of easily and automatically iterating through the GenBank output files from Prokka and extracting the 16S rRNA sequence (if annonated).

The only required input is the directory which contains the Prokka output directories:

    python ./16s_XTrakTA.py /PATH/TO/PROKKA/DIRECTORIES

