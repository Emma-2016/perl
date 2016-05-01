#! /usr/bin/perl -w
use strict;

use File::Basename qw(dirname basename);
use Getopt::Long;

## prel $0 /path/to/sample.db /path/to/screen.configure /path/to/projectname /path/to/outputdir 3 CANONICAL
my ($sampleDB, $screenConfig, $inputPath, $outputPath, $tools, $trans, $help) = ('', '', '', '', 1, 'CANONICAL', '');
GetOptions(
	"sampleDB:s" => \$sampleDB,
	"screenConfigure:s" => \$screenConfig,
	"projectDir:s" => \$inputPath,
	"outputDir:s" => \$outputPath,
	"tools:i" => \$tools,
	"transcriptType:s" => \$trans,
	"help:s" => \$help
);

my $usage=<<USAGE;
usage:perl $0
  -samPleDB:		path to your sample.db.
  -screenConfigure:	Path to your screen.configure.
  -projectDir:		Path to your project directory, which contains all results of samples (list in the sample.db).
  -outputDir:		path of output directory to contain report results.
  -tools:		Number of variant calling tools that called concordantly on the same position. Choose among [1, 2, 3, 4]. Default is 1, which means variants that were called by not less than one tool would be taken into account.
  -transcriptType:	Transcript type on which you would like to know the consequence caused by mutation. Choose between ['CANONICAL', 'Most harmful']. Default is CANONICAL. 
  -help:		help
USAGE
die $usage if (!$sampleDB || !$screenConfig || !$inputPath || !$outputPath || $help);

`rm -r $outputPath` if -d "$outputPath";
`mkdir -p $outputPath`;

my @samplesID = readInSamples($sampleDB);
my @functionalElements = readInScreen($screenConfig);
foreach my $sampleID (@samplesID)
{
	`mkdir -p $outputPath/$sampleID` unless -d "$outputPath/$sampleID";
	getQC($sampleID, $inputPath, $outputPath);
	getAlignment($sampleID, $inputPath, $outputPath);
	getVCF($sampleID, $inputPath, $outputPath, $tools, \@functionalElements, $trans);
}

my $allSamples = `ls $outputPath`;
@samplesID = split /\n/, $allSamples;
@samplesID = grep (!/Items/, @samplesID);

open OUT, ">$outputPath/sampleInfo";
foreach $_(@samplesID)
{
	my @out = (split /\t/, `grep $_ $sampleDB`)[0, 1, 3, 4, 5, 8, 9, 11, 12];
	my $out = join("\t", @out);
	print OUT $out."\n";
}
close OUT;
combination("snp", $outputPath);
combination("indel", $outputPath);
combination("alignment", $outputPath);

#------------------------Make up html file
my $html=<<HTMLCODE;
<!DOCTYPE html>
<html>

<head>
	<meta charset='utf-8'>
	<title>Sample Report of Exome Capture sequencing in BGI</title>
	<meta name='viewport' content='width=device-width, initial-scale=1.0'>
	<meta name='description' content='The website of Project WES report '>
	<meta name='author' content='Mendelian Disease Group @ BGI'>

	<!-- Javascript -->
	<script type='text/javascript' src='./jquery-1.7.2.min.js'></script>
	<script type='text/javascript' src='./prettify.js'></script>
	<script type='text/javascript' src='./bootstrap.min.js'></script>

	<!-- Le styles -->
	<link type='text/css' rel='stylesheet' href='./bootstrap.min.css'>
	<link type='text/css' rel='stylesheet' href='./prettify.css'>

	<style type='text/css'>
	span.RED{color:Red;}
	span.BLUE{color:Blue;}
	span.CORNFLOWERBLUE{color:cornflowerblue;}

	.brand {
		padding-top: 0;
		padding-bottom: 0;
	}
	footer {
		background-color: #fafafa;
		color: #999;
		padding: 20px;
		box-shadow: 0 -1px 10px rgba(0,0,0,0.1);
	}
    .navbar-fixed-top .container {
		width: 820px;
    }
	#searchform {
		margin-top: 0 !important;
	}
	.container {
		width: 800px;
	}
	p {
		text-align: justify;
	}
	h1 {
		font-size: 26px;
	}
	.center {
		text-align: center;
	}
	.right {
		text-align: right;
	}
	.sidebar {
		width: 220px;
		position: fixed;
		top: 50%;
		margin-top: -200px;
		left: 50%;
		margin-left: -640px;
	}
	.well {
		padding: 12px 4px;
	}
	</style>
</head>

<body>

<header class='navbar navbar-fixed-top'>
	<div class='navbar-inner'>
		<div class='container'>
			<ul class='nav'>
				<li class="dropdown">
				<a href="javascript:void(0)" class="dropdown-toggle" data-toggle="dropdown">BGI<b class="caret"></b></a>
					<ul class="dropdown-menu">
						<li> <a href="http://xinxi.229andme.com/" target="_blank">229andme</a> </li>
						<li> <a href="http://www.genomics.cn/index" target="_blank">BGI HOME</a> </li>
						<li> <a href="http://www.bgitechsolutions.com/" target="_blank">BGI Tech</a> </li>
						<li> <a href="http://www.bgidx.cn/index" target="_blank">BGI Dx</a> </li>
					</ul>
				</li>
				<li class="dropdown">
				<a href="javascript:void(0)" class="dropdown-toggle" data-toggle="dropdown">Useful websites<b class="caret"></b></a>
					<ul class="dropdown-menu">
					<li> <a href="http://www.irdirc.org/" target="_blank">IRDIRC</a> </li>
					<li> <a href="http://www.biostars.org/" target="_blank">Biostar</a> </li>
					<li> <a href="http://seqanswers.com/" target="_blank">Seqanswer</a> </li>
					<li> <a href="https://www.broadinstitute.org/gatk/guide/topic?name=faqs" target="_blank">GATK FAQ</a> </li>
					</ul>
				</li>
			</ul>
		</div>  <!-- class='container' -->
	</div>  <!-- class='navbar-inner' -->
</header>
<!--    header ends     -->


<div class='container-fluid' id="top" style="margin-top: 60px;">
	<div class="row-fluid">
	<!--            side column             -->
		<div class="span3" style="position:fixed; margin: -25px 0 0 -20px;">
			<div class="well sidebar-nav">

				<p class="center"><b>Report for Panel</b></p>
				<ul class="nav nav-list">
                	<li style="margin-top:12px; list-style-type:none;"> <a href="#mark0">Sample Information</a> </li>
                	<li style="margin-top:12px; list-style-type:none;"> <a href="#mark1">Description of Workflow</a> </li>
						<ul>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#ExperimentalPipeline">1. Experimental pipeline</a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#BioinformaticsPipeline">2. Bioinformatics analysis pipeline</a></li>
						</ul>
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark2">Project Results</a> </li>
						<ul>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#AlignmentResults">1. Alignment Results</a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#SequencingDepthDistribution">2. Sequencing Depth Distribution</a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#snpResults">3. SNP Results</a></li>
                        	<li style=" margin-top:7px; list-style-type:none;"><a href="#indelResults">4. INDEL Results</a></li>
						</ul>
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark3">Notes on Analyzing Methods</a> </li>
						<ul>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#MethodsAndParameters">1. Brief description of bioinformatics methods and parameters<a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#Databases">2. Databases</a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#FileFormat">3. File format</a></li>
						</ul>
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark4">Glossary</a> </li>
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark5">References</a> </li>
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark6">FAQ</a> </li>
				</ul>

			</div> <!--/.well -->
		</div> <!--/span-->
<!--            side column             -->


<!--            main column             -->
<div class="span9" style="margin-left:25%">

<!-- div for sample information -->
<div>
	<div id="mark0"></div><hr /><h3> Sample Information </h3><hr />
	<h4></h4>
	<table border="1" cellpdding="5" cellspacing="0" bordercolordark="#DCDCDC" align="center">
HTMLCODE

#----------------------------------------------------------------------------
sub combination
{
	my ($searchFile, $outputPath) = @_;
	my $allFile = `find $outputPath/*/ -name $searchFile`;
	my @allFile = split /\n/, $allFile;
	$allFile = join(' ', @allFile);
	my $itemFile = $searchFile."Items";
	my $resultFile = $searchFile."Result";
	`paste $outputPath/$itemFile $allFile >	$outputPath/$resultFile`;
	#`rm $allFile`;
}

sub getQC
{
	my ($sampleID, $inputPath, $outputPath) = @_;
	my $sampleDir = "$inputPath/$sampleID";
	my $outputDir = "$outputPath/$sampleID";

	open OUT, ">$outputDir/QC";
	my $fqc1 = `find $sampleDir/*/lane/*/1.cleanFQ/SOAPnuke/pe_1.fq_fastqc/ -name fastqc_data.txt`;
	my $fqc1Result = `grep '>>' $fqc1`;
	my @fqc1Result = split /\n/, $fqc1Result;
	foreach $_(@fqc1Result)
	{
		next if ($_ =~ /END_MODULE/);
		$_ =~ s/>>//;
		print OUT $_."\n";
	}
	close OUT;


	open OUT, ">>$outputDir/QC";
	print OUT "\n\n";
	my $fqc2 = `find $sampleDir/*/lane/*/1.cleanFQ/SOAPnuke/pe_2.fq_fastqc/ -name fastqc_data.txt`;
	my $fqc2Result = `grep '>>' $fqc2`;
	my @fqc2Result = split /\n/, $fqc2Result;
	foreach $_(@fqc2Result)
	{
		next if ($_ =~ /END_MODULE/);
		$_ =~ s/>>//;
		print OUT $_."\n";
	}
	close OUT;
}

sub getAlignment
{
	my ($sampleID, $inputPath, $outputPath) = @_;
	my $sampleDir = "$inputPath/$sampleID";
	my $outputDir = "$outputPath/$sampleID";
	my $figures = `find $sampleDir/*/sample/1.alignment/BAM_Stat/ -name *.png`;
	my @figures = split /\n/, $figures;
	foreach my $figure (@figures){`cp $figure $outputDir/`;}

	my $alignment = `find $sampleDir/*/sample/1.alignment/BAM_Stat/ -name information.xls`;
	chomp $alignment;
	`awk -F "\t" '{print \$2}' $alignment > $outputDir/alignment`;
	unless (-e "$outputPath/alignmentItems")
	{`awk -F "\t" '{print \$1}' $alignment >$outputPath/alignmentItems`;}
}

sub getVCF
{
	my ($sampleID, $inputPath, $outputPath, $tools, $functionalElements, $trans) = @_;
	my $sampleDir = "$inputPath/$sampleID";
	my $outputDir = "$outputPath/$sampleID";
	my $vcf = `find $sampleDir/*/sample/3.annotation -name sample.unionVar.VEP.vcf`;
	chomp $vcf;
	my ($consequenceIndex , $rsNum, $Ti, $Tv, $TiTR, $TvTR)= (0, 0, 0, 0, 0, 0);
	my ($snp, $snpHomoTR, $snpHetTR, $snpHomo, $snpHet, $rsSNPnum) = (0, 0, 0, 0, 0, 0);
	my ($indel, $indelHomoTR, $indelHetTR, $indelHomo, $indelHet, $rsINDELnum) = (0, 0, 0, 0, 0, 0);
	my %SNPconsequence = ();
	my %INDELconsequence = ();
	my @vcf = ();
	my @SNP = ();
	my @INDEL = ();

	open IN, $vcf;
	while (<IN>)
	{
		chomp;
		next if $_ =~ /^#/;
		#chr6    144915196       .       T       C       0:0:1:0 TR      CSQ=C|
		my ($chr, $pos, $rs, $ref, $alt, $callingTool, $rg, $vep, $format, $gt) = split /\s+/, $_;
		my $toolNum = 0;
		$toolNum = 1 if $callingTool =~ /1:0:0:0/;
		$toolNum = 2 if $callingTool =~ /0:1:0:0/;
		$toolNum = 3 if $callingTool =~ /0:0:1:0/;
		$toolNum = 4 if $callingTool =~ /0:0:0:1/;
		next unless $toolNum >= $tools;
		if (length($ref) == 1 and length($alt) == 1)
		{
			$snp += 1;
			$rsSNPnum += 1 if $rs !~ /\./;

			my @base = ($ref, $alt); @base = sort @base; my $base = join('', @base);
			if ($base eq 'AG' | $base eq 'CT'){$Ti += 1; if ($rg eq 'TR'){$TiTR +=1;}}
			else {$Tv += 1; if ($rg eq 'TR'){$TvTR +=1;}}
			if ($gt =~ /^1\/1/){$snpHomo += 1;if ($rg eq 'TR'){$snpHomoTR += 1;}}
			else {$snpHet += 1;if ($rg eq 'TR'){$snpHetTR += 1;}}

			my $csq = getcsq($vep, $trans);
			if (defined $SNPconsequence{$csq}){$SNPconsequence{$csq} += 1;}
			else {$SNPconsequence{$csq} = 1;}
		}
		else
		{
			$indel += 1;
			$rsINDELnum += 1 if $rs !~ /\./;

			if ($gt =~ /^1\/1/){$indelHomo += 1;if ($rg eq 'TR'){$indelHomoTR += 1;}}
			else {$snpHet += 1;if ($rg eq 'TR'){$snpHetTR += 1;}}

			my $csq = getcsq($vep, $trans);
			if (defined $INDELconsequence{$csq}){$INDELconsequence{$csq} += 1;}
			else {$INDELconsequence{$csq} = 1;}
		}	

	}	
	close IN;

	foreach my $element(@{$functionalElements})
	{
		if (defined $SNPconsequence{$element}){push @SNP, $SNPconsequence{$element};}
		else {push @SNP, 0;}
	
		if (defined $INDELconsequence{$element}){push @INDEL, $INDELconsequence{$element};}
		else {push @INDEL, 0;}
	}

#---------------------------------------------
	my $snpTR = $snpHomoTR + $snpHetTR;
	my $TiTv = sprintf ("%.2f", $Ti/$Tv);
	my $TR_TiTv = sprintf ("%.2f", $TiTR/$TvTR);
	my $elements = join("\n", @SNP);	
	my $snpResult = <<snpResult;
$snp
$snpTR
$snpHomo
$snpHet
$rsSNPnum
$TiTv
$TR_TiTv
$elements
snpResult
	open OUT, ">$outputDir/snp";
	print OUT $snpResult;
	close OUT;

#---------------------------------------------
	my $indelTR = $indelHomoTR + $indelHetTR;
	$elements = join("\n", @INDEL);
	my $indelResult = <<indelResult;
$indel
$indelTR
$indelHomo
$indelHet
$rsINDELnum
$elements
indelResult
	open OUT, ">$outputDir/indel";
	print OUT $indelResult;
	close OUT;

#---------------------------------------------
	unless (-e "$outputPath/snpItems")
	{
		open OUT, ">$outputPath/snpItems";
		my $elements = join ("\n", @{$functionalElements});
		my $snpItems = <<snpItems;
Total number of SNPs
Number of SNPs in Target Regions
Number of Homozygous SNPs
Number of Heterozygous SNPs
Number of SNPs in dbSNP
Ti/Tv
Ti/Tv in Target Regions
$elements
snpItems
		print OUT $snpItems;
		close OUT;
	}

#---------------------------------------------
	unless (-e "$outputPath/indelItems")
	{
		open OUT, ">$outputPath/indelItems";
		my $elements = join ("\n", @{$functionalElements});
		my $indelItems = <<indelItems;
Total number of INDELs
Number of INDELs in Target Regions
Number of Homozygous INDELs
Number of Heterozygous INDELs
Number of INDELs in dbINDEL
$elements
indelItems
		print OUT $indelItems;
		close OUT;
	}
}

sub getcsq
{
	my ($vep, $trans) = @_;
	my @str = (split /,/, $vep);
	my %impact = ('transcript_ablation' => 4, 'splice_acceptor_variant' => 4, 'splice_donor_variant' => 4, 'stop_gained' => 4, 'frameshift_variant' => 4, 'stop_lost' => 4, 'start_lost' => 4, 'transcript_amplification' => 4, 'inframe_insertion' => 3, 'inframe_deletion' => 3, 'missense_variant' => 3, 'protein_altering_variant' => 3, 'splice_region_variant' => 2, 'incomplete_terminal_codon_variant' => 2, 'stop_retained_variant' => 2, 'synonymous_variant' => 2, 'coding_sequence_variant' => 1, 'mature_miRNA_variant' => 1, '5_prime_UTR_variant' => 1, '3_prime_UTR_variant' => 1, 'non_coding_transcript_exon_variant' => 1, 'intron_variant' => 1, 'NMD_transcript_variant' => 1, 'non_coding_transcript_variant' => 1, 'upstream_gene_variant' => 1, 'downstream_gene_variant' => 1, 'TFBS_ablation' => 1, 'TFBS_amplification' => 1, 'TF_binding_site_variant' => 1, 'regulatory_region_ablation' => 3, 'regulatory_region_amplification' => 1, 'feature_elongation' => 1, 'regulatory_region_variant' => 1, 'feature_truncation' => 1, 'intergenic_variant' => 1);

	if ($trans =~ /CANONICAL/i)
	{
		my $count = 0;
		foreach my $str(@str)
		{
			my @items = (split /\|/, $str);
			next if $#items <17;
			next unless $items[17] eq 'YES';
			my $csq = $items[4];
			return $csq;
		}
	}
	else
	{       
		my $maxHarmfulness = 0;
		my $csq = '';
		foreach my $str(@str)
		{
			my @items = split /\|/, $str;
			if (defined $impact{$items[4]} and $impact{$items[4]} > $maxHarmfulness)
			{
				$maxHarmfulness = $impact{$items[4]};
				$csq = $items[4];
			}
		}
		return $csq;
	}
}

sub readInScreen
{
	my ($screenConfig) = @_;
	my @elements = ();
	open IN, $screenConfig;
	while (<IN>)
	{
		chomp;
		next unless $_ =~ /^functional_element/;
		my $element = (split /=/)[1];
		$element =~ s/'//g;
		@elements = split /,/, $element;
		last;
	}
	close IN;
	return @elements;
}

sub readInSamples
{
	my ($sampleDB, $outputPath) = @_;
	my @sampleID = ();
	open IN, "$sampleDB";
	while (<IN>)
	{
		chomp;
		next if $_ =~/^#/;
		my $sampleID = (split /\t/)[1];
		push @sampleID, $sampleID;
	}
	close IN;
	return @sampleID;
}
