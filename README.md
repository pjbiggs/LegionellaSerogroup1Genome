# LegionellaSerogroup1Genome

## Overview

A repository to accompany our manuscript on the *Legionella longbeachae* serogroup 1 genome, the most clinically significant *Legionella* species in New Zealand, causing around two-thirds of all notified cases of Legionnaires’ disease.

**Title**: "Extensive epigenetic modification with large-scale chromosomal and plasmid recombination characterise the *Legionella longbeachae* serogroup 1 genome"

**Authors**: Sandy Slow, Trevor Anderson, David R. Murdoch, Samuel Bloomfield, David Winter and Patrick J. Biggs

This repository contains the following files and scripts:
* `LegionellaCircos.pl`
* `circos_F1157CHC.conf`
* `methylome_analyses.Rmd`


## Requirements to run these scripts 

The following software is required to run these scripts.

* Perl -- a relatively recent version of Perl, 5.20 and above.  Packages required are `Getopt::Long`, `Bio::SeqIO`, `Bio::Seq` and `DBI` (for connections to MySQL).
* R -- version 3.5 and above.  This has to be directly accessible from the path. Packages required are `ggplot2`, `reshape`, `srtringr` and `ggpubr`. 
* MySQL
* Reputer -- available from https://bibiserv.cebitec.uni-bielefeld.de/reputer.  Within the suite, the program `repfind` is used.
* Circos -- a relatively recent version of Circos (available from http://circos.ca/software/), version 0.69.5 and above.


## Script description

### LegionellaCircos.pl

This a complex Perl script that uses data from a variety of sources stored in a MySQL database to generate all the required input files for the Circos plot of the *Legionella longbeachae* serogroup 1 genome F1157CHC.  

Once the script is run, and edits are made to the Circos configuration file (`circos_F1157CHC.conf`), the images are generated using:

`cd ~/software/circos-0.69-5/`

`./bin/circos -conf circos_F1157CHC.conf -outputfile circosF1157CHC`


### circos_F1157CHC.conf

This is the Circos configuration file to show the structure of the file needed to generate the Circos plot.  There are many tracks with each one requiring its own data file, as well as general files to set up the image, and the tickmarks amongst others.  If this level of detail is of interest, please contact the author (@pjbiggs) for more information. 


### methylome_analyses.Rmd

This document contains the code used to analysis the methylome data and to produce the individual elements of figure 6 in the main body. To run this code, you will need to install the R packages `ggplot2`, `reshape`, `srtringr` and `ggpubr`.  In order to produce the PDF document included as s supplementary file you will need the package `rmarkdown`.
