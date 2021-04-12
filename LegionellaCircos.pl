#!/usr/bin/env perl
#
#	the script to process data and make the data required for the Circos plot
#
#	created by pjb on: 		2017-05-25
#	last edited by pjb on:	2021-04-13
#
##################################################

use warnings;
use strict;
use Bio::SeqIO;
use Bio::Seq;
use DBI;
use Getopt::Long;

#######################

my ($curMod, $col1, $throw1, $modBase, $gStart, $gEnd, $score, $strand, $strand2, $phase, $rest, $curType, $curDist);
my ($statement, $joiner, $dbh, $sth, $datasource, $querystring, $rowcount, $count, $end);

my $root   	 	= ("/path/to/data/folder/");
my $cRoot		= ($root . "results/");
my $COGfolder	= ($root . "COGonCP020894/");
my $Droot   	= ($root . "referenceGenomeV3/");
my $modF		= ($root . "corrected_mod_data/");
my $cogRef		= ($root . "COG/");
my $log 		= ($root . "mappingLog.txt");

my @mods    = ('m4C', 'm6A');
my @recom	= ('N', 'R');
my @bases	= ('m4C', 'm6A');
my @strand	= ('plus', 'neg');
my @consist	= ('yes', 'no');
my @types 	= ('f', 'r', 'c', 'p');
my @dist	= ('0', '-1', '-2', '-3');
my @syn		= ('Synonymous','Nonsynonymous');

open (LOG, ">$log") or die ("couldn't open $log: $!\n");

print ("Process started at " . scalar(localtime) . ".\n");
print LOG ("Process started at " . scalar(localtime) . ".\n");


## connect to the db ##

&dbConnect();


## define the files for the input data ##

my $n1				= ("Legionella");
my $BEDcounts  		= ($n1 . "BED");
my $SNPtable    	= ($n1 . "SNP");
my $modBases		= ($n1 . "ModBases");
my $geneTable		= ($n1 . "Genes");
my $geneTable2		= ($n1 . "GenesOnly");
my $blockTable		= ($n1 . "Blocks");
my $repTable		= ($n1 . "Reputer");
my $refCOGName		= ("cognames2003_2014");
my $refFuncName		= ("fun2003_2014");
my $localCOG		= ($n1 . "COG");
my $COGcolT			= ($localCOG . "colours");
my $COGresults		= ("COGnitorResOn_Legionella");

my $SNPdata     	= ($root . "SNP_information_F1157CHC.txt");
my $baseBED     	= ($modF . "1kb_windowsMod.bed");
my $nucBED			= ($modF . "nuc_content.tsv");
my $geneGFF			= ($Droot . "CP020894_3.gff");
my $geneGFFMod		= ($Droot . "CP020894_3Mod.txt");

my $recBlocks		= ($root . "Legionella_longbeachae_recombination_areas.txt");
my $LegFa			= ($Droot . "CP020894_3.fasta");
my $refCOGtable		= ($cogRef . $refCOGName . ".txt");
my $refFuncTable	= ($cogRef . $refFuncName . ".txt");
my $COGcolours		= ($cogRef . "COGcoloursExtra.txt");
my $COGout			= ($COGfolder . "eggNOGdata.txt");


## create the tables ##

&tableCreate();


## load in the gff files ##

&gffTables();


## genome analysis for COG work ##

&COGparse();


## reputer work ##

&reputerCalc();


## work on the table ##

&tableMods1();


## start to categorise SNPs and modified bases ##

&tableMods2();


## prepare for some Circos output ##

&circosOut();


print ("All finished at " . scalar(localtime) . ".\n");
print LOG ("All finished at " . scalar(localtime) . ".\n");

close LOG;


#####################
#					#
#	subroutines 	#
#					#
#####################


sub COGparse {

	## genome GFF file to parse and load ##

	open (GFF, "<$geneGFF") or die ("couldn't open $geneGFF: $!\n");	
	open (GFF2, ">$geneGFFMod") or die ("couldn't open $geneGFFMod: $!\n");		

	while (<GFF>) {
		chomp;
		my @data	= split("\t");
		my $size	= @data;

		my ($contig, $method, $annoType, $genStart, $genEnd, $dot1, $strand, $dot2, $info)	= @data;
		
		if ($contig eq 'CP020894.3' && $size == 9) {

			my @longAnno	= split (";", $info);
			my $geneID		= $longAnno[0];

			pop(@data);
			push(@data, $geneID);

			print GFF2 (join("\t", @data), "\n");
		}
	}

	close GFF;
	close GFF2;	
	
	$sth = $dbh->prepare (qq{drop table if exists $geneTable});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $geneTable (altChrs varchar(30), method varchar(25), objectType varchar(25), gStart mediumint, gEnd mediumint, score mediumint, strand varchar(3), phase varchar(20), geneName varchar(100))});	$sth->execute();

	$sth = $dbh->prepare (qq{alter table $geneTable add index index1(gStart, gEnd, strand, geneName)});	$sth->execute();
	
	$sth = $dbh->prepare (qq{load data local infile '$geneGFFMod' into table $geneTable});	$sth->execute();		

	$sth = $dbh->prepare (qq{delete from $geneTable where objectType = 'exon'});	$sth->execute();
	$sth = $dbh->prepare (qq{update $geneTable set geneName = replace(geneName, 'ID=', '')});	$sth->execute();

	$sth = $dbh->prepare (qq{alter table $geneTable add column COG varchar(15) default 'nd'});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $geneTable add column COGcode varchar(10) default 'nd'});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $geneTable add column colChoice varchar(20) default '255,255,255'});	$sth->execute();

	$sth = $dbh->prepare (qq{alter table $geneTable add index index2(COGcode, COG)});	$sth->execute();
	
	$sth = $dbh->prepare (qq{update $geneTable set COG = 'ps' where objectType = 'pseudogene'});	$sth->execute();
	$sth = $dbh->prepare (qq{update $geneTable set COG = 'tR' where objectType = 'tRNA'});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $geneTable set COG = 'rR' where objectType = 'rRNA'});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $geneTable set COGcode = 'ps' where objectType = 'pseudogene'});	$sth->execute();
	$sth = $dbh->prepare (qq{update $geneTable set COGcode = 'tR' where objectType = 'tRNA'});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $geneTable set COGcode = 'rR' where objectType = 'rRNA'});	$sth->execute();
	
	
	## create data for outer rings ##

	$sth = $dbh->prepare (qq{update $COGresults set COGcode = substring_index(matchingOGs, ",", -1)});	$sth->execute();
	$sth = $dbh->prepare (qq{update $COGresults set COGcode = 'nd' where COGcode not like "COG%"});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $COGresults set COGcode = replace(COGcode, '\@1', '')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $COGresults set COGcode = replace(COGcode, '\@2', '')});	$sth->execute();

 	$sth = $dbh->prepare (qq{update $geneTable g, $COGresults l set g.COGcode = l.COGcode where concat("gene-", l.geneName) = concat(g.geneName, "_1") and COG not in ('ps', 'rR', 'tR')});	$sth->execute();
 	$sth = $dbh->prepare (qq{update $geneTable g, $COGresults l set g.COG = l.COGcat where concat("gene-", l.geneName) = concat(g.geneName, "_1")});	$sth->execute();

	$sth = $dbh->prepare (qq{update $geneTable g, (select b.* from (select * from Legion4Genes where COG not in ('nd', '', 'ps', 'tR', 'rR')) b where length(b.COG) != 1) a set g.COG = 'coCOG' where g.geneName = a.geneName});	$sth->execute();	
 	$sth = $dbh->prepare (qq{update $geneTable g, $COGcolT c set g.colChoice = c.colChoice where g.COG = c.COGtype});	$sth->execute();

    print ("\tNew COG analysis complete at " . scalar(localtime) . ".\n");
    print LOG ("New COG analysis complete at " . scalar(localtime) . ".\n");
}



sub reputerCalc {
	
	foreach $curType (@types) {
		my $data	= ($cRoot . "RepRes_" . $curType . ".txt");
		my $dataMod	= ($cRoot . "RepRes_" . $curType . "Mod.txt");	
		my $working	= ("-" . $curType);
		
		system "repfind $working -l 30 -h 3 -best 10000 $LegFa > $data";
		system "cat $data | perl -lpe 's/^\\s+//g;s/\\s+/\\t/g;' > $dataMod";
		
		$sth = $dbh->prepare (qq{load data local infile '$dataMod' into table $repTable ignore 2 lines});	$sth->execute();			
	}

	$sth = $dbh->prepare (qq{delete from $repTable where rep1length = 0});	$sth->execute();

	$sth = $dbh->prepare (qq{alter table $repTable add column rep1Start mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $repTable add column rep2Start mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $repTable add column rep1End mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $repTable add column rep2End mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $repTable add column ID mediumint auto_increment primary key first});	$sth->execute();

	$sth = $dbh->prepare (qq{update $repTable set rep1Start = rep1Loc where type in ('P', 'C')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $repTable set rep2Start = rep2Loc where type in ('P', 'C')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $repTable set rep1Start = (rep1Loc + 1) where type in ('R', 'F')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $repTable set rep2Start = (rep2Loc + 1) where type in ('R', 'F')});	$sth->execute();

	$sth = $dbh->prepare (qq{update $repTable set rep1End  = (rep1Loc + rep1Length - 1) where type in ('P', 'C')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $repTable set rep2End  = (rep2Loc + rep2Length - 1) where type in ('P', 'C')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $repTable set rep1End  = (rep1Loc + rep1Length) where type in ('R', 'F')});	$sth->execute();
	$sth = $dbh->prepare (qq{update $repTable set rep2End  = (rep2Loc + rep2Length) where type in ('R', 'F')});	$sth->execute();

	print ("Repeats processed at " . scalar(localtime) . ".\n");
	print LOG ("Repeats processed at " . scalar(localtime) . ".\n");	
}



sub circosOut {

	## file for checking ##
	
	my $BEDdata	= ($cRoot . "allBEDdata.txt");
	my @header	= ('chromosome', 'altChrs', 'binStart', 'binEnd', 'm4C', 'm6A', 'pctAT', 'pctGC', 'SNP_N', 'SNP_R', 'yes_m4C', 'no_m4C', 'yes_m6A', 'no_m6', 'delta_m4C', 'delta_m6A');
	
	open (OUT, ">$BEDdata") or die "$BEDdata not opened\n";
		
	print OUT (join("\t", @header), "\n");	
		
	$statement = ("select *, (m4C - (yes_m4C + no_m4C)) as delta_m4C, (m6A - (yes_m6A + no_m6A)) as delta_m6A from $BEDcounts order by binStart, chromosome");
	&statementPull ($statement, "\t", "\n");
	
	close OUT;	


	## recombination areas ##

	my $Crecom	= ($cRoot . "recombination2.txt");	
	
	open (OUT, ">$Crecom") or die "$Crecom not opened\n";

	# add in now columns due to reordering of chromosome for starting at dnaA #

	$statement = ("select 'CP020894.3', newBlockStart, newBlockEnd, 'fill_color=dgreen' from $blockTable order by blockStart");		
	&statementPull ($statement, " ", "\n");
	
	close OUT;

	
	## genes for outer ring ##
	
	foreach my $curStrand (@strand) {
		my $outfile	= ($cRoot . "F1157CHC_" . $curStrand . ".txt");
		
		if ($curStrand eq 'plus'){		$strand2 = "+";
		} elsif ($curStrand eq 'neg'){	$strand2 = "-";
		}
		
		open (OUT, ">$outfile") or die "$outfile not opened\n";
			
		$statement = ("select 'CP020894.3', gStart, gEnd, concat('fill_color=', colChoice) from $geneTable where geneName like \"\%B0B39\%\" and strand = '$strand2' order by gStart");	
		&statementPull ($statement, " ", "\n");	
		
		close OUT;		
	}
	
	
	## SNPs ##

	push (@recom, 'all');
	
	foreach my $curRecom (@recom) {
		my $outSNP	= ($cRoot . "SNPdata_" . $curRecom . ".txt");
		my $SNPR	= ("SNP_" . $curRecom);
		
		open (OUT, ">$outSNP") or die "$outSNP not opened\n";

		if ($curRecom eq 'all') {
			$statement = ("select 'CP020894.3', binStart, binEnd, round(log10(SNP_R + SNP_N), 4) from $BEDcounts where chromosome = 'CP020894.3' and (SNP_R + SNP_N) != 0 order by binStart");				
		} else {
			$statement = ("select 'CP020894.3', binStart, binEnd, round(log10($SNPR), 4) from $BEDcounts where chromosome = 'CP020894.3' and $SNPR != 0 order by binStart");							
		}

		&statementPull ($statement, " ", "\n");
			
		close OUT;
	}
	

	## synonymous SNPs ##

	foreach my $curSyn (@syn) {
		my $outSNP	= ($cRoot . "SNPdata_" . $curSyn . ".txt");
		my $SNPR	= ("SNP_" . $curSyn);
		
		open (OUT, ">$outSNP") or die "$outSNP not opened\n";

		$statement = ("select 'CP020894.3', binStart, binEnd, round(log10($SNPR), 4) from $BEDcounts where chromosome = 'CP020894.3' and $SNPR != 0 order by binStart");				
		
		&statementPull ($statement, " ", "\n");
			
		close OUT;
	}

	
	## GC content & methylation ##
	
	my @mod2    = ('m4C', 'm6A', 'pctGC');
	
	foreach my $mods2 (@mod2) {
		my $outFile	= ($cRoot . "data_" . $mods2 . ".txt");
		
		open (OUT, ">$outFile") or die "$outFile not opened\n";

		if ($mods2 eq 'pctGC') {
			$statement = ("select 'CP020894.3', binStart, binEnd, $mods2 from $BEDcounts where chromosome = 'CP020894.3' order by binStart");			
		} else {
			$statement = ("select 'CP020894.3', binStart, binEnd, round(log10($mods2), 4) from $BEDcounts where chromosome = 'CP020894.3' and $mods2 != 0 order by binStart");			
		}			

		&statementPull ($statement, " ", "\n");
			
		close OUT;		
	}	
	
	
	## generate methyl consistent data ##	
	
	foreach my $curBase (@bases) {
		foreach my $curConsist (@consist) {
			my $compact	= ($curConsist . "_" . $curBase);	
			my $methOut	= ($cRoot . $compact . "_out.txt");

			open (OUT, ">$methOut") or die "$methOut not opened\n";
			
			$statement = ("select 'CP020894.3', binStart, binEnd, round(log10($compact), 4) from $BEDcounts where chromosome = 'CP020894.3' and $compact != 0 order by binStart");
			&statementPull ($statement, " ", "\n");
				
			close OUT;			
		}
	}
	
	## links for reputer ##
	
	foreach my $curType (@types) {
		foreach my $curDist (@dist) {
			my $compact	= ($curType . "_" . $curDist);	
			my $methOut	= ($cRoot . "Rep_outFor_" . $compact . ".txt");			
					
			open (OUT, ">$methOut") or die "$methOut not opened\n";

			$statement = ("select ID, 'CP020894.3', rep1Start, rep1End from $repTable where type = upper('$curType') and eValue <= '1e-20' and distance = '$curDist'");
			&statementPull ($statement, " ", "\n");
			$statement = ("select ID, 'CP020894.3', rep2Start, rep2End from $repTable where type = upper('$curType') and eValue <= '1e-20' and distance = '$curDist'");
			&statementPull ($statement, " ", "\n");
				
			close OUT;			
		}
	}	
	
	print ("Circos output made at " . scalar(localtime) . ".\n");
	print LOG ("Circos output made at " . scalar(localtime) . ".\n");	
}



sub gffTables {

	# modified bases #
	
	foreach my $newBase (@bases) {
		my $in		= ($modF . $newBase . "_lifted.gff");
		my $inMod	= ($modF . $newBase . "_lifted_mod.gff");	
	
		open (GFF, "<$in") or die ("couldn't open $in: $!\n");	
		open (GFF2, ">$inMod") or die ("couldn't open $inMod: $!\n");		
		
		while (<GFF>) {
			chomp;
			my ($col1, $throw1, $modBase, $gStart, $gEnd, $strand, $phase, $rest)	= split;
			
			my @other	= split(';', $rest);
			my $size	= @other;
			my @head1;
			
			if ($size == 4) {
				@head1 = ($col1, $modBase, $gStart, $gEnd, $strand, $rest, @other);
			} elsif ($size == 7) {
				@head1 = ($col1, $modBase, $gStart, $gEnd, $strand, $rest, $other[0], $other[1], $other[2], $other[6]);
			}
			
			print GFF2 (join("\t", @head1), "\n");		
		}
		
		$sth = $dbh->prepare (qq{load data local infile '$inMod' into table $modBases});	$sth->execute();	
			
		close GFF;
		close GFF2;
		
		print ("\tModifed base $newBase loaded at " . scalar(localtime) . ".\n");
		print LOG ("\tModifed base $newBase loaded at " . scalar(localtime) . ".\n");		
	}	
}



sub tableMods2 {
	
	# alter $modBases #
	
    $sth = $dbh->prepare (qq{alter table $modBases add column objectType varchar(25) default 'n/a'});	$sth->execute();	
    $sth = $dbh->prepare (qq{alter table $modBases add column objectStrand varchar(5) default 'n/a'});	$sth->execute();	
    $sth = $dbh->prepare (qq{alter table $modBases add column objectConsistent varchar(5) default 'n/a'});	$sth->execute();
	
	# output the gene coordinates and update the table #
	
	my $geneCoord	= ($root . "geneCoords.txt");
	
	open (OUT, ">$geneCoord") or die "$geneCoord not opened\n";
		
	$statement = ("select objectType, gStart, gEnd, strand from $geneTable where objectType in ('CDS', 'tRNA', 'rRNA', 'pseudogene') order by gStart");	
	&statementPull ($statement, "\t", "\n");	
	
	close OUT;	

	# read in and make updates to tables # 

	open (IN, "<$geneCoord") or die "$geneCoord not opened\n";
	
	while (<IN>) {
		chomp;
		my ($obj1, $gStart, $gEnd, $str)	= split;
		
		$sth = $dbh->prepare (qq{update $modBases set objectType = '$obj1' where gStart between $gStart and $gEnd});	$sth->execute();
		$sth = $dbh->prepare (qq{update $modBases set objectStrand = '$str' where gStart between $gStart and $gEnd});	$sth->execute();		
	}	
	
	close IN;

	# set strand consistency #

	$sth = $dbh->prepare (qq{update $modBases set objectConsistent = 'yes' where (strand = objectStrand and objectStrand = '+') or (strand = objectStrand and objectStrand = '-')});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $modBases set objectConsistent = 'no' where (strand != objectStrand and objectStrand = '+') or (strand != objectStrand and objectStrand = '-')});	$sth->execute();		

	# set the consistency for the factors under consideration #

	foreach my $curBase (@bases) {
		foreach my $curConsist (@consist) {
			my $compact	= ($curConsist . "_" . $curBase);
			my $nanoT	= ($compact . "Table");

			$sth = $dbh->prepare (qq{drop table if exists $nanoT});	$sth->execute();			
							
			$sth = $dbh->prepare (qq{create table $nanoT select any_value(altChrs) as altChrs, modBin, count(*) as $compact from $modBases where altChrs != 'n/a' and objectConsistent = '$curConsist' and modBase = '$curBase' group by modBin});	$sth->execute();
			$sth = $dbh->prepare (qq{alter table $nanoT add column binStart mediumint});	$sth->execute();	
			$sth = $dbh->prepare (qq{alter table $nanoT add column binEnd mediumint});	$sth->execute();	
		
			$sth = $dbh->prepare (qq{update $nanoT set binStart = (modBin * 1000) + 1});	$sth->execute();
			$sth = $dbh->prepare (qq{update $nanoT set binEnd = (modBin + 1) * 1000});	$sth->execute();
		
			$sth = $dbh->prepare (qq{alter table $BEDcounts add column $compact mediumint default 0});	$sth->execute();					   
			$sth = $dbh->prepare (qq{update $BEDcounts b, $nanoT t set b.$compact = t.$compact where b.chromosome = t.altChrs and t.binStart = b.binStart and t.binEnd = b.binEnd});	$sth->execute();				
			
			print ("\t$compact data calculated at " . scalar(localtime) . ".\n");			
			print LOG ("\t$compact data calculated at " . scalar(localtime) . ".\n");				
		}
	}		

	# develop data for COG work #
	
	print ("Second round of table modification complete at " . scalar(localtime) . ".\n");	
	print LOG ("Second round of table modification complete at " . scalar(localtime) . ".\n");
}



sub tableMods1 {

    # SNP table first #
    
    $sth = $dbh->prepare (qq{alter table $SNPtable add column actualSNP char(1)});	$sth->execute();
    $sth = $dbh->prepare (qq{alter table $SNPtable add column leftFlank varchar(30)});	$sth->execute();
    $sth = $dbh->prepare (qq{alter table $SNPtable add column rightFlank varchar(30)});	$sth->execute();

    $sth = $dbh->prepare (qq{update $SNPtable set actualSNP = substring(kmerSNP, 21, 1)});	$sth->execute();
    $sth = $dbh->prepare (qq{update $SNPtable set leftFlank = left(kmerSNP, 20)});	$sth->execute();
    $sth = $dbh->prepare (qq{update $SNPtable set rightFlank = right(kmerSNP, 20)});	$sth->execute();
        
    foreach $curMod (@mods) {
		my $curBED	= ($modF . $curMod . "_hist.bed");
		my $tempBED	= ("BEDtemp");
					   
		$sth = $dbh->prepare (qq{drop table if exists $tempBED});	$sth->execute();
		$sth = $dbh->prepare (qq{create table $tempBED (chromosome varchar(40), binStart mediumint, binEnd mediumint, $curMod mediumint)});	$sth->execute();
		$sth = $dbh->prepare (qq{alter table $tempBED add index index1(binStart, binEnd)});	$sth->execute();
	
		$sth = $dbh->prepare (qq{load data local infile '$curBED' into table $tempBED ignore 1 lines});	$sth->execute();

		$sth = $dbh->prepare (qq{update $tempBED set binStart = binStart + 1});	$sth->execute();	
		$sth = $dbh->prepare (qq{alter table $BEDcounts add column $curMod mediumint default 0});	$sth->execute();
		$sth = $dbh->prepare (qq{update $BEDcounts b, $tempBED t set b.$curMod = t.$curMod where t.chromosome = b.chromosome and t.binStart = b.binStart and t.binEnd =  b.binEnd});	$sth->execute();
		
		$sth = $dbh->prepare (qq{drop table if exists $tempBED});	$sth->execute();
    }
    
    # add in the nucleotide data #
    
	$sth = $dbh->prepare (qq{drop table if exists nucCounts});	$sth->execute();
	$sth = $dbh->prepare (qq{create table nucCounts (chromosome varchar(40), binStart mediumint, binEnd mediumint, pctAT decimal(4,3), pctGC decimal(4,3), numA mediumint, numC mediumint, numG mediumint, numT mediumint, numN mediumint, numOth mediumint, seqLen mediumint)});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table nucCounts add index index1(binStart, binEnd)});	$sth->execute();

	$sth = $dbh->prepare (qq{load data local infile '$nucBED' into table nucCounts ignore 1 lines});	$sth->execute();

	$sth = $dbh->prepare (qq{update nucCounts set binStart = binStart + 1});	$sth->execute();
	
	$sth = $dbh->prepare (qq{alter table $BEDcounts add column pctAT decimal(4,3) default 0});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $BEDcounts add column pctGC decimal(4,3) default 0});	$sth->execute();		
	$sth = $dbh->prepare (qq{update $BEDcounts b, nucCounts t set b.pctAT = t.pctAT where b.chromosome = 'CP020894.3' and t.binStart = b.binStart and t.binEnd =  b.binEnd});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $BEDcounts b, nucCounts t set b.pctGC = t.pctGC where b.chromosome = 'CP020894.3' and t.binStart = b.binStart and t.binEnd =  b.binEnd});	$sth->execute();
	
	$sth = $dbh->prepare (qq{drop table if exists nucCounts});	$sth->execute();
    		
	# add in SNP data # 
	
	foreach my $curRecom (@recom) {
		my $localR	= ("SNPcount_" . $curRecom);
		my $SNPR	= ("SNP_" . $curRecom);
		
		$sth = $dbh->prepare (qq{drop table if exists $localR});	$sth->execute();
		$sth = $dbh->prepare (qq{create table $localR select SNPbin, count(*) as $SNPR from $SNPtable where recombination = "$curRecom" group by SNPbin});	$sth->execute();
		$sth = $dbh->prepare (qq{alter table $localR add column binStart mediumint});	$sth->execute();	
		$sth = $dbh->prepare (qq{alter table $localR add column binEnd mediumint});	$sth->execute();	
	
		$sth = $dbh->prepare (qq{update $localR set binStart = (SNPbin * 1000) + 1});	$sth->execute();
		$sth = $dbh->prepare (qq{update $localR set binEnd = (SNPbin + 1) * 1000});	$sth->execute();
	
		$sth = $dbh->prepare (qq{alter table $BEDcounts add column $SNPR mediumint default 0});	$sth->execute();					   
		$sth = $dbh->prepare (qq{update $BEDcounts b, $localR t set b.$SNPR = t.$SNPR where b.chromosome = 'CP020894.3' and t.binStart = b.binStart and t.binEnd =  b.binEnd});	$sth->execute();
		
		$sth = $dbh->prepare (qq{drop table if exists $localR});	$sth->execute();
	}

	# add in synonymous status data # 
	
	foreach my $curSyn (@syn) {
		my $localR	= ("SNPstuff_" . $curSyn);
		my $SNPR	= ("SNP_" . $curSyn);
		
		$sth = $dbh->prepare (qq{drop table if exists $localR});	$sth->execute();
		$sth = $dbh->prepare (qq{create table $localR select SNPbin, count(*) as $SNPR from $SNPtable where SNPtype = "$curSyn" group by SNPbin});	$sth->execute();
		$sth = $dbh->prepare (qq{alter table $localR add column binStart mediumint});	$sth->execute();	
		$sth = $dbh->prepare (qq{alter table $localR add column binEnd mediumint});	$sth->execute();	
	
		$sth = $dbh->prepare (qq{update $localR set binStart = (SNPbin * 1000) + 1});	$sth->execute();
		$sth = $dbh->prepare (qq{update $localR set binEnd = (SNPbin + 1) * 1000});	$sth->execute();
	
		$sth = $dbh->prepare (qq{alter table $BEDcounts add column $SNPR mediumint default 0});	$sth->execute();
		$sth = $dbh->prepare (qq{update $BEDcounts b, $localR t set b.$SNPR = t.$SNPR where b.chromosome = 'CP020894.3' and t.binStart = b.binStart and t.binEnd =  b.binEnd});	$sth->execute();
		
		$sth = $dbh->prepare (qq{drop table if exists $localR});	$sth->execute();
	}
		
	# add in name for Circos #
	
	$sth = $dbh->prepare (qq{alter table $BEDcounts add column altChrs varchar(30) default 'n/a' after chromosome});	$sth->execute();
	$sth = $dbh->prepare (qq{update $BEDcounts set altChrs = 'F1157CHC_pacbio' where chromosome = 'CP020894.3'});	$sth->execute();
	
	# tidy up the modified bases table #
	
	$sth = $dbh->prepare (qq{alter table $modBases drop column IPDratioN});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $modBases set coverage = replace(coverage, 'coverage=', '')});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $modBases set context = replace(context, 'context=', '')});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $modBases set idQv = replace(idQv, 'identificationQv=', '')});	$sth->execute();	

	$sth = $dbh->prepare (qq{alter table $modBases modify column idQv mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $modBases modify column coverage mediumint});	$sth->execute();

	# add in modbases bins #
	
	$sth = $dbh->prepare (qq{alter table $modBases add column modBin mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{update $modBases set modBin = floor(gStart/1000)});	$sth->execute();	
	
	print ("First round of table modification complete at " . scalar(localtime) . ".\n");	
	print LOG ("First round of table modification complete at " . scalar(localtime) . ".\n");	
}



sub tableCreate {

    print ("Base tables work started at " . scalar(localtime) . ".\n");
    print LOG ("Base tables work started at " . scalar(localtime) . ".\n");

	
	# gffModBases #

	$sth = $dbh->prepare (qq{drop table if exists $modBases});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $modBases (altChrs varchar(30), modBase varchar(5), gStart mediumint, gEnd mediumint, strand varchar(3), attributes varchar(160), coverage varchar(20), context varchar(60), IPDratioN varchar(20), idQv varchar(20))});	$sth->execute();
	
	$sth = $dbh->prepare (qq{alter table $modBases add index index1(gStart, gEnd, strand, context)});	$sth->execute();	

  
	# recombination blocks  #
  
 	$sth = $dbh->prepare (qq{drop table if exists $blockTable});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $blockTable (blockStart mediumint, blockEnd mediumint)});	$sth->execute();
 
 	$sth = $dbh->prepare (qq{load data local infile '$recBlocks' into table $blockTable ignore 1 lines});	$sth->execute();	
 
	$sth = $dbh->prepare (qq{alter table $blockTable add column newBlockStart mediumint});	$sth->execute(); 
 	$sth = $dbh->prepare (qq{alter table $blockTable add column newBlockEnd mediumint});	$sth->execute(); 
 
	$sth = $dbh->prepare (qq{update $blockTable set newBlockStart = mod((blockStart + 2566986), 4142881)});	$sth->execute(); 
	$sth = $dbh->prepare (qq{update $blockTable set newBlockEnd = mod((blockEnd + 2566986), 4142881)});	$sth->execute(); 
 
 
	# reputer table #
 
	$sth = $dbh->prepare (qq{drop table if exists $repTable});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $repTable (rep1Length smallint, rep1Loc mediumint, type enum ('F', 'R', 'C', 'P'), rep2Length smallint, rep2Loc mediumint, distance tinyint, eValue float)});	$sth->execute();

   
    # $SNPtable #

	$sth = $dbh->prepare (qq{drop table if exists $SNPtable});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $SNPtable (CP020894pos mediumint, SNPtype enum ('Synonymous','Nonsynonymous'), kmerSNP varchar(50))});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $SNPtable add index index1(CP020894pos, kmerSNP)});	$sth->execute();

	$sth = $dbh->prepare (qq{load data local infile '$SNPdata' into table $SNPtable ignore 1 lines});	$sth->execute();

	$sth = $dbh->prepare (qq{alter table $SNPtable add column newPos mediumint});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $SNPtable set newPos = mod((CP020894pos + 2566986), 4142881)});	$sth->execute();	
	
	$sth = $dbh->prepare (qq{alter table $SNPtable add column recombination char(2) default 'N'});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $SNPtable t, $blockTable b set t.recombination = 'R' where t.newPos = b.newBlockStart});	$sth->execute();	
	$sth = $dbh->prepare (qq{update $SNPtable t, $blockTable b set t.recombination = 'R' where t.newPos = b.newBlockEnd});	$sth->execute();
	$sth = $dbh->prepare (qq{update $SNPtable t, $blockTable b set t.recombination = 'R' where t.newPos between b.newBlockStart and b.newBlockEnd});	$sth->execute();
	
	$sth = $dbh->prepare (qq{alter table $SNPtable add column SNPbin mediumint});	$sth->execute();
	$sth = $dbh->prepare (qq{update $SNPtable set SNPbin = floor(newPos/1000)});	$sth->execute();
	
	
	# $BEDcounts #
	
	$sth = $dbh->prepare (qq{drop table if exists $BEDcounts});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $BEDcounts (chromosome varchar(40), binStart mediumint, binEnd mediumint)});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $BEDcounts add index index1(binStart, binEnd)});	$sth->execute();

	$sth = $dbh->prepare (qq{load data local infile '$baseBED' into table $BEDcounts ignore 1 lines});	$sth->execute();
	
	$sth = $dbh->prepare (qq{update $BEDcounts set binStart = (binStart + 1)});	$sth->execute();

	
	# refFuncName table #

	$sth = $dbh->prepare (qq{drop table if exists $refFuncName});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $refFuncName (code char(5),codeFunction varchar(100))});	$sth->execute();
	$sth = $dbh->prepare (qq{load data local infile '$refFuncTable' into table $refFuncName ignore 1 lines});	$sth->execute();


	# refCOGName table #

	$sth = $dbh->prepare (qq{drop table if exists $refCOGName});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $refCOGName (COG varchar(10), code char(5) references	$refFuncName.code, COGname varchar(100))});	$sth->execute();
	$sth = $dbh->prepare (qq{load data local infile '$refCOGtable' into table $refCOGName ignore 1 lines});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $refCOGName add index(COG)});	$sth->execute();


	# COG colours #
		
	$sth = $dbh->prepare (qq{drop table if exists $COGcolT});	$sth->execute();
	$sth = $dbh->prepare (qq{create table $COGcolT (COGtype varchar(10), colChoice varchar(20), hexChoice varchar(10), COGname varchar(100))});	$sth->execute();
	$sth = $dbh->prepare (qq{load data local infile '$COGcolours' into table $COGcolT ignore 1 lines});	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $COGcolT add index(COGtype)});	$sth->execute();

	
	# $COGresults #

	$sth = $dbh->prepare (qq{drop table if exists $COGresults});	$sth->execute();	
	$sth = $dbh->prepare (qq{create table $COGresults (geneName varchar(20), seedOrtholog varchar(50), eValue float, scoreVal float, bestTaxLvl varchar(30), preferredName varchar(20), GOterms text, ECnumber varchar(70), KEGGKO varchar(90), KEGGpathway text, KEGGmodule varchar(80), KEGGreaction varchar(250), KEGGrclass varchar(160), BRITE varchar(100), KEGGTC varchar(80), CAZy varchar(20), BiGGreaction text, annotLvl varchar(20), matchingOGs varchar(150), BestOG varchar(20), COGcat varchar(10), description text)});	$sth->execute();

	$sth = $dbh->prepare (qq{load data local infile '$COGout' into table $COGresults ignore 1 lines});	$sth->execute();

	$sth = $dbh->prepare (qq{alter table $COGresults add column COGcode varchar(50) default 'nd' after matchingOGs });	$sth->execute();
	$sth = $dbh->prepare (qq{alter table $COGresults add index(geneName, COGcat, COGcode)});	$sth->execute();
	
    print ("All base tables work complete at " . scalar(localtime) . ".\n");
    print LOG ("All base tables work complete at " . scalar(localtime) . ".\n");		
}



sub statementPull {
	($statement, $joiner, $end) = @_;

	$sth = $dbh->prepare (qq{$statement});	$sth->execute();
				
	$count++;
				
	while (my @row_items = $sth->fetchrow_array ()) {
		$rowcount++;
		print OUT (join ("$joiner", @row_items), "$end");
		} unless ($rowcount) {
		print OUT ("No data to display\n");
	}

	return ($statement, $joiner, $end);
}



sub dbConnect {
	$count 		= 0;
	$rowcount 	= 0;

	$datasource = "DBI:mysql:<<database>>;mysql_local_infile=1";
	$dbh = DBI->connect($datasource, '<<yourName>>', '<<yourPassword>>');
	$querystring = '';
}

