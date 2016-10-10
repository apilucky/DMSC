#! /usr/bin/perl -w 

## This script was writted in order to add layer of Drug Resistance-Succeptible information to variant position
## Input 
# 1. Annotated Variant in VCF format from snpEff 
# 2. Target Region of Drug Resistance Genes in BED format
# 3. Drug resistance Mutation List 

use strict;
use Cwd;
use File::Basename;
use Getopt::Long;

my ($workingDir, $VCF_Dir, $BED_DrugTarget, $DB);
## Get option Seesion
GetOptions(
	"w=s"	=> \$workingDir,	
	"i=s"	=> \$inputDir,
	"t=s"	=> \$BED_DrugTarget,
	"DB=s"	=> \$DB,
);

## set of default parameters
if (!defined $workingDir){
# use default current Directory for working directory
	$workingDir = cwd();	
}else{

}

if (!defined $BED_DrugTarget){
# use default BEDTarget Drug Resistance file
	$BED_DrugTarget =  
	my $DrugBED = $VCF_Dir."/DrugResistance_MTB.bed";

}else{

}
my $VCF_Dir = "/home/apiluck/Shared/Apiluck_DMSc" ;

# Create Database for Target Region
my $DrugBED = $VCF_Dir."/DrugResistance_MTB.bed";
# Create Database for Drug Resistance Mutation List
my $DrugMutation = $VCF_Dir."/DrugResistanceMutation.txt";


## Main Loop for Each Snp
for my $VCF (`ls $VCF_Dir/*.vcf`){
	chomp($VCF);
	
	open(VCF, $VCF);
	my @VCFFile = <VCF>;chomp(@VCFFile);
	
	for my $variant (@VCFFile){
		if ($variant !~ /^#/){
			my $pos = getvarPos($variant);
			my $ref = getRefBase($variant);
			my $alt = getAltBase($variant);
			my @Ann = getVarAnno($variant); 

			my @Gene = getGeneinBED($pos, $DrugBED);
			my $TargetGene	= $Gene[1];
			my $NucCh	= $Ann[3];
			my $ProtCh	= $Ann[4];

			my @DrugDB = getDrugEff($pos, $Gene[1], @Ann, $DrugMutation);

#			print "Variant input is \n";
			print join("\t", $pos, $ref, $alt, @Ann, @Gene, @DrugDB); #sleep 1;

		}else{

		}
	}

	#sleep 1;
} 

### sub-routein Session
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

#	my @Drug = getDrugEff($pos, $Gene[1], @Ann, $DrugMutation);
#	my @VarAnno = ($eff, $gNa, $GlT, $PosOnGene, $RefB, $MutB, $codon, $aa1, $aa2);
	#Input Agrument
#	print "input in are\n";
#	print join ("\t", @_), "\n";

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

