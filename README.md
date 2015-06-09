# FRI-failure
Code for analyzing sequencing data from Freshman Research Initiative genetic failure mode experiments

All python scripts located in the /python directory are written to work with sequencing files
in the /sequences directory, organized as follows:
   -folder with student initials containing:
   	   -Sanger sequencing reads, containing either "vf" or "vr"

seqproc.py creates files necessary for MAFFT to create alignments.

analysis.py does the actual alignment analysis and creates the GenomeDiff file.
