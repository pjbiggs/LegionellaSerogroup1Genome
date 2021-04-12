# LegionellaSerogroup1Genome
A repository to accompany our manuscript on different characteristics of the *Legionella longbeachae* serogroup 1 genome.

This repository contains the following files and scripts:
* `LegionellaCircos.pl`
* `file B`


## Requirements to run these scripts 

The following software is required to run these scripts.

* Perl -- a relatively recent version of Perl, 5.20 and above.  Packages required are `Getopt::Long`, `Bio::SeqIO`, `Bio::Seq` and `DBI` (for connections to MySQL).
* R version 3.5 and above.  This has to be directly accessible from the path. 
* MySQL
* Reputer -- available from https://bibiserv.cebitec.uni-bielefeld.de/reputer.  Within the suite, the program `repfind` is used.
* Circos -- a relatively recent version of Circos (available from http://circos.ca/software/), version 0.69.5 and above.


## Script description

### LegionellaCircos.pl

This a complex Perl script that uses data from a variety of sources stored in a MySQL database to generate all the required input files for the Circos plot of the *Legionella longbeachae* serogroup 1 genome F1157CHC.  

Once the script is run, and edits are made to the Circos configuration file, the images are generated using:

`cd ~/software/circos-0.69-5/`
`./bin/circos -conf circos_F1157CHC.conf -outputfile circosF1157CHC`


## file b


