#! /usr/bin/perl -w 

## This script was writted in order to add layer of Drug Resistance-Succeptible information to variant position
## Input 
# 1. Annotated Variant in VCF format from snpEff 
# 2. Target Region of Drug Resistance Genes in BED format
# 3. Drug resistance Mutation List 

use strict;
use Cwd;
use File::Spec;
use File::Basename;
use Getopt::Long;

use POSIX qw(strftime);

my ($ARG_workingDir, $inputDir, $ARG_BED_DrugTarget, $ARG_DB);
## Get option Seesion
GetOptions(
	"i=s"	=> \$inputDir,

	"w=s"	=> \$ARG_workingDir,	
	"t=s"	=> \$ARG_BED_DrugTarget,
	"DB=s"	=> \$ARG_DB,
);

## Pre-Configure Variable for INSTALLATION PROCESS
my $INSTALL_DIR = "/home/apiluck/bin/perl/TBSNP";

my ($workingDir, $BED_DrugTarget, $DB) = setDefaultParam($ARG_workingDir, $INSTALL_DIR ,$ARG_BED_DrugTarget, $ARG_DB);
my $abs_input = checkInputDir($inputDir);


if (-e "ouput"){
	print "Found Previous Ouput directory\n";
	`rm -r output`;
	print "Remove Old Output Directory\n";
	`mkdir -p output`;
}else{
	`mkdir -p output`;
	print "Create New Output Directory\n"; sleep 1;
}

my $time = getTime();
my $abs_ResultDir = createResultDir($time);

## Main Loop for Each Snp
for my $VCF (`ls $abs_input/*.vcf`){
	chomp($VCF);
	
	open(VCF, $VCF);
	my @VCFFile = <VCF>;chomp(@VCFFile);
	my $VCFbase = basename($VCF, ".vcf");
	print "Searching on $VCFbase\n";
	
	mkdir "$abs_ResultDir/onBED"; mkdir "$abs_ResultDir/NotonBED";
	
	open (OutputBED, ">$abs_ResultDir/onBED/$VCFbase.txt");
	open (OutputNotBED, ">$abs_ResultDir/NotonBED/$VCFbase.txt");
		
	print OutputBED join("\t", 'Pos', 'Ref', 'Alt', 'Effect', 'Gene Name', 'LocusTag', 'Pos.On Gene', 'Ref', 'Alt', 'Condon', 'AA1', 'AA2', 'Drug Name', 'Report'),"\n";
	print OutputNotBED join("\t", 'Pos', 'Ref', 'Alt', 'Effect', 'Gene Name', 'LocusTag', 'Pos.On Gene', 'Ref', 'Alt', 'Condon', 'AA1', 'AA2', 'Drug Name', 'Report'),"\n";
	
	for my $variant (@VCFFile){
		if ($variant !~ /^#/){
			my $pos = getvarPos($variant);
			my $ref = getRefBase($variant);
			my $alt = getAltBase($variant);
			my @Ann = getVarAnno($variant); 

			my @Gene = getGeneinBED($pos, $BED_DrugTarget);
			my $TargetGene	= $Gene[1];
			my $NucCh	= $Ann[3];
			my $ProtCh	= $Ann[4];

			my @DrugDB = getDrugEff($pos, $Gene[1], @Ann, $DB);
			if($Gene[0] eq '-'){
				print OutputNotBED join("\t", $pos, $ref, $alt, @Ann, @DrugDB); #sleep 1;
				
			}else{
				print OutputBED join("\t", $pos, $ref, $alt, @Ann, @DrugDB); #sleep 1;
				
			}
			
		}else{		}
	}
}
print "CHECK THE OUTPUT FILE IN THE FOLLOWING DIRECTORY: [$abs_ResultDir]\n\n";

### sub-routein Session

sub checkInputDir{
	my $inputDir = $_[0];

	my $abs_inputDir = ();

	if (defined $inputDir){
		
		$abs_inputDir = $workingDir."/$inputDir";
		if(-e $abs_inputDir){

			my @vcfInput = `ls $abs_inputDir/*.vcf`;
			if(scalar(@vcfInput) == 0 ){	
				my $Err3 = "Cannot Find VCF in Input Directory ($abs_inputDir) \n";
				printHelp($Err3);
			}else{
			}
		}else{
			my $Err2 = "Recheck the working directory & input directory\n";
			printHelp($Err2);
		}
	}else{
		my $Err1 = "Input directory Not found: [Specify -i option]\n";
		printHelp($Err1);
	}
	return($abs_inputDir)
}

sub setDefaultParam {

	my $workingDir = $_[0];
	my $INSTALL_DIR = $_[1];
	my $BED_DrugTarget = $_[2];
	my $DB = $_[3];

	## set of default parameters
	if (!defined $workingDir){
	# use default current Directory for working directory
		$workingDir = cwd();	
		print "## use Default Current Working Directory:\t$workingDir\n";
	}else{ print "## use Default Current Working Directory:\t$workingDir\n";}

	if (!defined $BED_DrugTarget){
	# use default BEDTarget Drug Resistance file
		$BED_DrugTarget = $INSTALL_DIR.'/BED/DrugResistance_MTB.bed';
		print "## use Default BED Target file:\t$BED_DrugTarget\n";
	}else{print "## use Default BED Target file:\t$BED_DrugTarget\n";}

	if (!defined $DB){
	# use default DB Drug Resistance Report file
		$DB = $INSTALL_DIR.'/DB/DrugResistanceMutation.txt';
		print "## use Default Drug Resistance file:\t$DB\n";
	}else{print "## use Default Drug Resistance file:\t$DB\n";}

	return ($workingDir, $BED_DrugTarget, $DB);
}

sub createResultDir
{
	my $date = $_[0];
	 
	my $Rep = 1;
	my $Result_Dir = 'TBSNP_'.$date."_$Rep";        
    if(-e $Result_Dir){    
        $Rep++;
        $Result_Dir = 'TBSNP_'.$date."_$Rep";
              
	    while(-e $Result_Dir)
	    {	
            $Rep++;
            $Result_Dir = 'TBSNP_'.$date."_$Rep";
	    }
    }
    mkdir $Result_Dir;
    my $abs_resultDir = File::Spec->rel2abs($Result_Dir);
    return ($abs_resultDir);
}

sub getTime
{
	my $currentTime = strftime "%a_%b_%d", localtime;
	return ($currentTime);
}

sub getvarPos {
	my $variantInput = $_[0];
	my $pos = (split("\t", $variantInput))[1];
	return ($pos);
}

sub getRefBase {
	my $variantInput = $_[0];
	my $ref = (split("\t", $variantInput))[3];
	return ($ref);
}

sub getAltBase {
	my $variantInput = $_[0];
	my $alt = (split("\t", $variantInput))[4];
	return ($alt);
}

sub getVarAnno {
	my $variantInput = $_[0];
	my $ANN	= (split(';',(split("\t", $variantInput))[7]))[3];
	my $ann	= (split('=',(split(",", $ANN))[0]))[1];
	my $eff	= (split(/\|/, $ann))[1];
	my $gNa	= (split(/\|/, $ann))[3];
	my $GlT	= (split(/\|/, $ann))[4];
	my $nuc	= (split(/\|/, $ann))[9]; 
	my $pro = (split(/\|/, $ann))[10];
	my $codon = '-';	
	my $aa1	= '-';
	my $aa2 = '-';

	if ($pro !~ /^p/){
		$pro = 'AA-';
	}else{	
		$codon = substr $pro, 5, -3;		 
		$aa1 = substr $pro, 2, 3;
		$aa2 = substr $pro, -3;
	}

	my $PosOnGene = substr $nuc, 2, -3;
	my $RefB	= substr $nuc, -3, 1;
	my $MutB	= substr $nuc, -1; 

	my @VarAnno = ($eff, $gNa, $GlT, $PosOnGene, $RefB, $MutB, $codon, $aa1, $aa2);
	return (@VarAnno);
}

sub getGeneinBED{
	my $pos = $_[0];
	my $BED = $_[1];

	open (BED, $BED);
	my @BED = <BED>; chomp(@BED);


	my $geneTarget = '-';
	my $product = '-';

	for my $target (@BED){
		my @bedData = split("\t", $target);
		my $str = $bedData[1];
		my $stp = $bedData[2];
		my $gene = $bedData[3];
		my $pro	= $bedData[4];

		if($pos > $str && $pos < $stp){
			$geneTarget = $gene;
			$product = $pro;
		}else{

		}
	}
	
	my @BED_Search = ($geneTarget, $product);
	#print join("\t", @BED_Search),"\n";
	return (@BED_Search);
}

sub getDrugEff{

	my $TargetGene = $_[1];
	my $eff = $_[2];
	my $PosOnGene = $_[5];
	my $codon	= $_[8];
	my $RefB = $_[6];
	my $MutB = $_[7];	
	my $aa1 = $_[9];
	my $aa2 = $_[10];
	my $DrugMutationDB = $_[11];


	#Return Value
	my $FindGeneName 	= '-';
	my $FNucChange 		= '-';
	my $FAAChange		= '-';
	my $FChange			= '-';
	my $FDrugName		= 'NotinDB';
	my $FReport			= 'NotinDB';

	open (DrugMutation, $DrugMutationDB);
	my @MutationDB = <DrugMutation>; chomp(@MutationDB);

	for my $DB (@MutationDB){
		chomp($DB);
		my @drugData = split("\t", $DB);
		
		my $Target 	= $drugData[0];
		my $Change 	= $drugData[1];
		my $drugName = $drugData[2];
		my $report 	= $drugData[3];
		
		
		if ($Target =~ /promotor/){
			my $refB = substr $Change, 0, 1;
			my $mutB = substr $Change, -1;
			my $posG = substr $Change, 1, -1;
			
			if ($TargetGene eq $Target && $PosOnGene eq $posG && $RefB eq $refB && $MutB eq $mutB){
				$FDrugName = $drugName;
				$FReport = $report;	
			}else{}
		
		}else{
			my $codonDB = substr $Change, 3, -3;
			my $DBaa1 = substr $Change, 0, 3;
			my $DBaa2 = substr $Change, -3;
		 	
			if ($TargetGene eq $Target && $codon eq $codonDB && $aa1 eq $DBaa1 && $aa2 eq $DBaa2){
				$FDrugName = $drugName;
				$FReport = $report;	
			}else{}
		}
	}
	return ($FDrugName, $FReport),"\n";
}

sub printHelp
{
	my $ErrMasseage = $_[0];
	my @help = <DATA>;
	print STDERR @help;
	die("Plase Defined Valid Input: $ErrMasseage\n");
}

###################### USE OF TBVariant_RS_Annotation ############################

__DATA__

NAME
~~~~
	TBVariant_RS_Annotation integrate Annotated Variant Calling Data from SNPEff Software with Drug Resistance Database Return to the table of summary variant.

SYNOPSIS
~~~~~~~~
	TBVariant_RS_Annotation -i VCF_input_Directory [-t -w -DB: Optional]	

OPTION AGRUMENTS
~~~~~~~~~~~~~~~~

	-t  [BED] Sequence Target File for Drug Resistance Gene

		Format Specification: 5 columns BED format
		1st Column: Chromosome Name of Reference [NC_000962]
		2nd Column: Start Position on Reference Sequence
		3rd Column: Stop Position on Reference Sequence
		4th Column: LocusTag of Target Region
		5th Column: Gene Name/Name of Target Region
		
		In case of promotor Target Region the 5th Column should be code in the format of 
			"Gene Name"_"promotor" ex. katG_promotor

		!!CAUTION: Coding in Different Format will cause unsatify result

	-w [Absolute Path of Current Working Directory]
		Default is set by using CWD module of perl

	-DB [TAB] Drug Resistance Gene Target/ Changing Position which has been reported

		Format Specification: 4 Columns TAB file	
		1st Column: Target Name 
		2nd Column: Changing [Gly415Ala, C-15T]
		3rd Column: Drug Name
		4th Column: Report


DESCRIPTION
~~~~~~~~~~~
	This script was written in perl programming language. The script was written in order to 
	integrate Annotated Variant Calling Data from SNPEff Software with Drug Resistance Database 
	Return to the table of summary variant.

AUTHOR
~~~~~~
	ApilucK Musigkain, Field Application Specialist (FAS) Bioinfomatics, Geneplus.CO., Ltd
	apiluck@gene-plus.com, (+66)994-242-541 

VERSION
~~~~~~~
	Verion 1.0 [Production]

Copyright (C) 2016 by Gene-plus.CO.,Ltd. Ayothaya Tower, Ratchadaphisek Rd.,Huay Kwang, 
Bangkok 10310 Thailand. All rights reserved.

