#################################

Readme: FRI-Failure Analysis Tool 
Copyright Tyler Camp.
2015/7/16

#################################

Code for analyzing sequencing data from Freshman Research Initiative genetic failure mode experiments

Requirements:
-plasmids stored in GenBank format, preferably annotated according to the GenBank DDBJ/EMBL/GenBank Feature Table Definitions.
-Sanger reads stored in .ab1/.abi file format
-Python version: at least 2.7 (standard on most Linux distros)
-Biopython package installed: biopython.org/wiki/Main_Page for more info

A short explanation of the different Python source files:
-frifailurejob.py: Python script that controls analysis pipeline with (currently) hard-coded filepaths for input and output
-jobmanager.py: Wrapper class for running analysis jobs involving custom Python modules, input, and nested output directories
-Analysis.py: Modules containing all relevant analysis functions for the FRI-Failure tool

Usage:

john$ python frifailurejob.py

Tool outputs logfile.log to log errors.

