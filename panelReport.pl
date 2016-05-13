#! /usr/bin/perl -w
use strict;

use File::Basename qw(dirname basename);
use Getopt::Long;

my ($sampleDB, $screenConfig, $inputPath, $outputPath, $tools, $trans, $help) = ('', '', '', '', '', '', '');
GetOptions
(
	"sampleDB:s" => \$sampleDB,
	"screenConfigure:s" => \$screenConfig,
	"projectDir:s" => \$inputPath,
	"outputDir:s" => \$outputPath,
	"tools:i" => \$tools,
	"transcriptType:s" => \$trans,
	"help:s" => \$help
);

$tools ||= 1;
$trans ||= 'CANONICAL';
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

`mkdir -p $outputPath` unless -d "$outputPath";
`cp /ifshk1/BC_MD/PROJECT/jiangchongyi/htmlConfig/* $outputPath/`;

#------------------------Combine results
my @samplesID = readInSamples($sampleDB); @samplesID = sort @samplesID;
my @functionalElements = readInScreen($screenConfig);
foreach my $sampleID (@samplesID)
{
	`mkdir -p $outputPath/$sampleID` unless -d "$outputPath/$sampleID";
	getQC($sampleID, $inputPath, $outputPath);
	getAlignment($sampleID, $inputPath, $outputPath);
	getVCF($sampleID, $inputPath, $outputPath, $tools, \@functionalElements, $trans);
}
combination("alignment", $outputPath);
combination("snp", $outputPath);
combination("indel", $outputPath);
getSampleInfo($outputPath, \@samplesID);
cleanUpvepFeature("snp", $outputPath);	#discard features that are all zero;
cleanUpvepFeature("indel", $outputPath);
cleanUpalignment($outputPath);

#-------------------------draw distribution figures
my $cum = `cat $outputPath/cum.txt`;
my @cum = split /\n/, $cum;
@cum= sort{$a <=> $b}@cum;
my $cum_max = $cum[-1];
my $histogram = `cat $outputPath/histogram.txt`;
my @histogram = split /\n/, $histogram;
@histogram = sort{$a<=>$b}@histogram;
my $histogram_max = $histogram[-1];
foreach my $sampleID (@samplesID){getFigure($sampleID, $inputPath, $outputPath, $cum_max, $histogram_max);}

#------------------------Prepare html components
my $numInLine = 3;
my $tablesNum;
if((@samplesID) % ($numInLine) == 0){$tablesNum = @samplesID / $numInLine;}
else{$tablesNum = (int(@samplesID / $numInLine) + 1);}

my $sampleInfo = html_sub("sampleInfo", $outputPath, $tablesNum);
my $AlignmentResults = html_sub("alignmentResult", $outputPath, $tablesNum);
my $snpResults = html_sub("snpResult", $outputPath, $tablesNum);
my $indelResults = html_sub("indelResult", $outputPath, $tablesNum);
my $distribution = distributionHTML($outputPath, \@samplesID);

#------------------------Make up html file
if ($tools == 1){$tools = '1 tool'} else{$tools = "$tools tools";}
$trans = lc($trans);
my $html = <<HTMLCODE;
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
<!--
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark3">Notes on Analyzing Methods</a> </li>
						<ul>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#MethodsAndParameters">1. Brief description of bioinformatics methods and parameters</a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#Databases">2. Databases</a></li>
							<li style=" margin-top:7px; list-style-type:none;"><a href="#FileFormat">3. File format</a></li>
						</ul>
-->
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark4">Glossary</a> </li>

					<!--
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark5">References</a> </li>
					<li style="margin-top:12px; list-style-type:none;"> <a href="#mark6">FAQ</a> </li>
					-->

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
<!--	<table border="1" cellpdding="5" cellspacing="0" bordercolordark="#DCDCDC" align="center"> -->
		$sampleInfo
	<br />
</div>

<div>
	<div id="mark2"> </div><hr /><h3>Project Results</h3><hr />
		<div class='row-fluid accordion-heading'>
			<div id="AlignmentResults"></div>
			<div class='span11'> <h4>1. Alignment Results</h4> </div>
		</div>
			<div class='accordion-inner'>
			$AlignmentResults
			</div><br />

		<div class='row-fluid accordion-heading'>
			<div id="SequencingDepthDistribution"></div>
			<div class='span11'> <h4>2. Sequencing Depth Distribution</h4> </div>
		</div>
			<div class='accordion-inner'>
			$distribution
			</div><br />

		<div class='row-fluid accordion-heading'>
			<div id="snpResults"></div>
			<div class='span11'> <h4>3. SNP Results</h4> </div>
		</div>
			<div class='accordion-inner'>
			$snpResults
			</div><br />

		<div class='row-fluid accordion-heading'>
			<div id="indelResults"></div>
			<div class='span11'> <h4>4. INDEL Results</h4> </div>
		</div>
			<div class='accordion-inner'>
			$indelResults
		<p style="text-indent: 2em;">(1) Four variant calling tools were used in this project. SNPs were called by (<a href = 'https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php' target = '_blank'>GATK HaplotypeCaller</a>, <a href = 'http://soap.genomics.org.cn/soapsnp.html' target = '_blank'>SOAPsnp</a>, and <a href = 'http://www.well.ox.ac.uk/platypus' target = '_blank'>Platypus</a>). INEDLs were called by (<a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php" target="_blank">GATK HaplotypeCaller</a>, <a href="http://samtools.sourceforge.net/" target="_blank">SAMtools</a>,and <a href="http://www.well.ox.ac.uk/platypus" target="_blank">Platypus</a>). Tabulations above were based on variants that were called by at least $tools. 
		<p style="text-indent: 2em;">(2) Variants were annotated by <a href = 'http://www.ensembl.org/info/docs/tools/vep/index.html' target = '_blank'>VEP</a>. A single gene often has more than one transcript, and these transcripts can be of different classes (please refer to <a href = '#mark4'>Glossary</a> below). Tabulations above were based on variants that were located in $trans transcripts. 
		<p style="text-indent: 2em;">(3) VEP annotations are divided into <a href = 'http://asia.ensembl.org/info/genome/variation/predicted_data.html'>35 categories</a>. This panel pipeline focus on the following features: 'stop_gained', 'frameshift_variant', 'stop_lost', 'initiator_codon_variant', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant',  '5_prime_UTR_variant', '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'non_coding_exon_variant', 'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant',  'stop_retained_variant', 'synonymous_variant',  'mature_miRNA_variant', 'incomplete_terminal_codon_variant', 'coding_sequence_variant', and 'intron_variant'. Tabulations above only displays those features that at least one variant from any sample were classified into.</p>
			</div><br />
</div>
<div>
	<div id="mark4"> </div><hr /><h3>Glossary</h3><hr />
    <h4></h4>
			<p><b>1.&nbsp;<a href = 'http://www.ensembl.org/Help/Glossary?id=346' target = '_blank'>Canonical transcript</a></b></p>
			<p>For human, the canonical transcript for a gene is set according to the following hierarchy: 1. Longest CCDS translation with no stop codons. 2. If no (1), choose the longest Ensembl/Havana merged translation with no stop codons. 3. If no (2), choose the longest translation with no stop codons. 4. If no translation, choose the longest non-protein-coding transcript.</p><br />
</div>

<!-- the following </div> is for the main column -->
</div>
<!-- the following </div> is for both side and main column -->
</div>
</div>
HTMLCODE

open OUThtml, ">$outputPath/Panel_Report.html";
print OUThtml $html;
close OUThtml;

#------------------------sub methods
sub readInSamples
{
	my ($sampleDB, $outputPath) = @_;
	my @sampleID = ();
	open IN, "$sampleDB";
	while (<IN>)
	{
		chomp;
		next if $_ =~ /^\s*$/;
		next if $_ =~/^#/;
		my $sampleID = (split /\t/)[1];
		push @sampleID, $sampleID;
	}
	close IN;
	return @sampleID;
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

sub getQC
{
	my ($sampleID, $inputPath, $outputPath) = @_;
	my $sampleDir = "$inputPath/$sampleID";
	my $outputDir = "$outputPath/$sampleID";
	my @lanes = `ls $sampleDir/*/lane/`;
	
	open OUT, ">$outputDir/QC.html";
	my $html_sub = <<html_sub;
<!DOCTYPE html>
<html>
<head>
	<meta charset='utf-8'>
	<style type = 'text/css'>
	p {
		text-align: justify;
		font-family: verdana;
	}
	table 
	{ 
	border-collapse: collapse; 
	border: none; 
	} 
	td {
		font-family: verdana;
		border: solid #000 1px;
	}
		
	</style>
</head>
html_sub
	print OUT $html_sub."\n";
	if (@lanes == 1){print OUT "<h2>Reads for sample $sampleID are produced from 1 lane.</h2>";}
	else {print OUT "<h2>Reads for sample $sampleID are produced from ".@lanes." lanes.</h2>";}
	foreach my $lane (@lanes)
	{
		chomp $lane;
		my %hash = ();
		my $raw1 = `find $sampleDir/*/lane/$lane/1.cleanFQ/SOAPnuke/$lane\_1.fq_fastqc/ -name fastqc_data.txt`;
		my $raw1Result = `grep '>>' $raw1`;
		my @raw1Result = split /\n/, $raw1Result;
		foreach $_(@raw1Result)
		{
        	next if ($_ =~ /END_MODULE/ or $_ =~ /Basic Statistics/);
        	$_ =~ s/>>//;
        	my ($k, $v) = split /\t/;
        	$hash{$k} = $v;
		}

		my $raw2 = `find $sampleDir/*/lane/$lane/1.cleanFQ/SOAPnuke/$lane\_2.fq_fastqc/ -name fastqc_data.txt`;
		my $raw2Result = `grep '>>' $raw2`;
		my @raw2Result = split /\n/, $raw2Result;
		foreach $_(@raw2Result)
		{
			next if ($_ =~ /END_MODULE/ or $_ =~ /Basic Statistics/);
			$_ =~ s/>>//;
			my ($k, $v) = split /\t/;
			$hash{$k} = ucfirst($hash{$k})." | ".ucfirst($v) if defined $hash{$k};
		}

		my $clean1 = `find $sampleDir/*/lane/$lane/1.cleanFQ/SOAPnuke/pe_1.fq_fastqc/ -name fastqc_data.txt`;
		my $clean1Result = `grep '>>' $clean1`;
		my @clean1Result = split /\n/, $clean1Result;
		foreach $_(@clean1Result)
		{
			next if ($_ =~ /END_MODULE/ or $_ =~ /Basic Statistics/);
			$_ =~ s/>>//;
			my ($k, $v) = split /\t/;
			$hash{$k} = ucfirst($hash{$k})."\t".ucfirst($v) if defined $hash{$k};
		}

		my $clean2 = `find $sampleDir/*/lane/$lane/1.cleanFQ/SOAPnuke/pe_2.fq_fastqc/ -name fastqc_data.txt`;
		my $clean2Result = `grep '>>' $clean2`;
		my @clean2Result = split /\n/, $clean2Result;
		print OUT "Qualtiy of reads from lane <b>$lane</b>:<br /><br />";
		$html_sub = <<html_sub;
<table border = '0' cellpdding = '2' cellspacing = '0' bordercolordark = '#DCDCDC' align = 'center'>
<tr>
<td width = '295' align = 'left'>&nbsp;Features</td>
<td width = '145' align = 'center'>Raw reads</td>
<td width = '145' align = 'center'>clean reads</td>
</tr>
html_sub
		foreach $_(@clean2Result)
		{
			next if ($_ =~ /END_MODULE/ or $_ =~ /Basic Statistics/);
			$_ =~ s/>>//;
			my ($k, $v) = split /\t/;
			$hash{$k} = ucfirst($hash{$k})." | ".ucfirst($v) if defined $hash{$k};
			my ($rawQC, $cleanQC) = split /\t/, $hash{$k};
			$html_sub .= <<html_sub;
<tr>
<td width = '295' align = 'left'>&nbsp;$k</td>
<td width = '145' align = 'center'>$rawQC</td>
<td width = '145' align = 'center'>$cleanQC</td>
</tr>
html_sub
		}
		$html_sub .= "</table><br /><br />";
		print OUT $html_sub;
	}
my $note = <<note;
<br />
<p style="text-indent: 2em;">(1) All features and ranks (Pass, Warn, or Fail) for each category were acquird from FastQC which is widely used for reads quality control. Fastqc performs a quick evaluation of whether several analysis modules seem entirely normal (Pass), slightly abnormal (Warn) or very unusual (Fail). For detailed description about FastQC, please refer to <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">homepage of FastQC</a>.</p>
<p style="text-indent: 2em;">(2) For pair reads, ranks are seperated by "|". The first one refers to read1 of the pair, while the second one refers to read2 of the pair.</p>
note
	print OUT $note;
	close OUT;
}

sub getAlignment
{
	my ($sampleID, $inputPath, $outputPath) = @_;
	my $sampleDir = "$inputPath/$sampleID";
	my $outputDir = "$outputPath/$sampleID";

#Prepare for cumulative figure
	my $cumulativeFig = `find $sampleDir/*/sample/1.alignment/BAM_Stat/ -name depth_frequency.xls`;
	chomp $cumulativeFig;
	open IN, "$cumulativeFig";
	my $max = 0;
	while (<IN>)
	{
		chomp;
		my $ratio = (split /\s+/)[1];
		$ratio = $ratio * 100;
		if ($ratio > $max){$max = $ratio;}
	}
	close IN;
	if (-e "$outputPath/cum.txt"){open OUT, ">>$outputPath/cum.txt"; print OUT $max."\n";}
	else{open OUT, ">$outputPath/cum.txt"; print OUT $max."\n";}

	my $alignment = `find $sampleDir/*/sample/1.alignment/BAM_Stat/ -name information.xls`;
	chomp $alignment;

#Prepare for histogram figure
	if (-e "$outputPath/histogram.txt"){`awk -F "\t" '{if (\$1~/Average sequencing depth on target/)print \$2}' $alignment >>$outputPath/histogram.txt`;}
	else{`awk -F "\t" '{if (\$1~/Average sequencing depth on target/)print \$2}' $alignment >$outputPath/histogram.txt`;}

#alignment result
	open IN, ">$outputDir/alignment"; print IN "$sampleID\n"; close IN;
	`awk -F "\t" 'NR > 1{print \$2}' $alignment >> $outputDir/alignment`;
	unless (-e "$outputPath/alignmentItems")
	{
		open OUT, ">$outputPath/alignmentItems";
		print OUT "Features\n"; close OUT;
		`awk -F "\t" '{print \$1}' $alignment >>$outputPath/alignmentItems`;
	}
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

			if ($gt =~ /^1\/1/)
			{$indelHomo += 1;if ($rg eq 'TR'){$indelHomoTR += 1;}}
			else
			{$indelHet += 1;if ($rg eq 'TR'){$indelHetTR += 1;}}

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
$sampleID
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
$sampleID
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
Features
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
Features
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
	my %impact = (
	'transcript_ablation' => 4, 'splice_acceptor_variant' => 4, 'splice_donor_variant' => 4, 
	'stop_gained' => 4, 'frameshift_variant' => 4, 'stop_lost' => 4, 'start_lost' => 4, 
	'transcript_amplification' => 4, 'inframe_insertion' => 3, 'inframe_deletion' => 3, 
	'missense_variant' => 3, 'protein_altering_variant' => 3, 'splice_region_variant' => 2, 
	'incomplete_terminal_codon_variant' => 2, 'stop_retained_variant' => 2, 'synonymous_variant' => 2, 
	'coding_sequence_variant' => 1, 'mature_miRNA_variant' => 1, '5_prime_UTR_variant' => 1, 
	'3_prime_UTR_variant' => 1, 'non_coding_transcript_exon_variant' => 1, 'intron_variant' => 1, 
	'NMD_transcript_variant' => 1, 'non_coding_transcript_variant' => 1, 'upstream_gene_variant' => 1, 
	'downstream_gene_variant' => 1, 'TFBS_ablation' => 1, 'TFBS_amplification' => 1, 'TF_binding_site_variant' => 1, 
	'regulatory_region_ablation' => 3, 'regulatory_region_amplification' => 1, 'feature_elongation' => 1, 
	'regulatory_region_variant' => 1, 'feature_truncation' => 1, 'intergenic_variant' => 1
	);

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

sub combination
{
	my ($searchFile, $outputPath) = @_;
	my $allFile = `find $outputPath/*/ -name $searchFile`;
	my @allFile = split /\n/, $allFile;
	$allFile = join(' ', @allFile);
	my $itemFile = $searchFile."Items";
	my $resultFile = $searchFile."Result";
	`paste $outputPath/$itemFile $allFile >	$outputPath/$resultFile`;
	`rm $allFile $outputPath/$itemFile`;
}

sub getSampleInfo
{
	my ($outputPath, $samplesID) = @_;
	my @samplesID = (@{$samplesID});
	open OUT, ">$outputPath/sampleInfo";
	my $count = 0;
	my %out = ();
	my @items = ('SampleID', 'FamlilyID', 'Gender', 'Age', 'Population',
					'Disease', 'CaseOrControl', 'CaptureChip', 'Platform', 'Reads Quality');
	$out{$count} = \@items;
	foreach $_(@samplesID)
	{
		$count += 1;
    	my @out = (split /\t/, `grep $_ $sampleDB`)[1, 0, 3, 4, 5, 8, 9, 11, 12];
		push @out, "<a href = '$_/QC.html'>$_</a>";
		$out{$count} = \@out;
	}
	my @newItems = ();
	my $newItems = '';
	for $_(0 .. 9)
	{
		@newItems = ();
		for my $j (0 .. @samplesID)
		{
			push @newItems, ${$out{$j}}[$_];
		}
		$newItems = join ("\t", @newItems);
		print OUT $newItems."\n";
	}
	close OUT;
}

sub cleanUpvepFeature
{
	my ($name, $outputPath) = @_;
	$name = $name."Result";
	open IN, "$outputPath/$name";
	my $out = '';
	while (<IN>)
	{
		chomp;
		my @items = split /\t/;
		my $mark = 0;
		A: foreach my $i(1 .. $#items)
		{
			next if ($items[$i] eq '0');
			$mark = 1;
			last A;
		}
		if ($mark == 1){$out .= $_."\n";}
	}
	close IN;
	open OUT, ">$outputPath/$name";
	chomp $out;
	print OUT $out;
	close OUT;
}	

sub cleanUpalignment
{
	my ($outputPath) = @_;
	open IN, "$outputPath/alignmentResult";
	my $out = '';
	while (<IN>)
	{
		chomp;
		my @items = split /\t/;
		foreach my $i (1..$#items)
		{
			if (($items[$i] eq '0.00%') or ($items[$i] eq '0.00') 
			or ($items[$i] eq '0.0%') or ($items[$i] eq '0.0'))
			{
				$items[$i] = '0';
			}
		}
		my $items = join("\t", @items);
		$out .= $items."\n";
	}
	close IN;
	open OUT, ">$outputPath/alignmentResult";
	chomp $out;
	print OUT $out;
	close OUT;
}

sub getFigure
{
	my ($sampleID, $inputPath, $outputPath, $cum_max, $histogram_max) = @_;
	my $outputDir = $outputPath."/".$sampleID;	
	my $sampleDir = $inputPath."/".$sampleID;	
	my $depth_frequency = `find $sampleDir/*/sample/1.alignment/BAM_Stat/ -name depth_frequency.xls`;
	my $cumu = `find $sampleDir/*/sample/1.alignment/BAM_Stat/ -name cumu.xls`;
	chomp $depth_frequency;
	chomp $cumu;
	
	my $ylim = $cum_max;
	my ($xlim, $xbin, $ybin);
	$ylim = int($ylim) + 1;
	if($ylim <= 3){$ybin = 0.5;}else{$ybin=1;}

	if (int(($histogram_max * 2) / 100) % 2 == 0)
	{
		$xlim = int($histogram_max * 2 - (($histogram_max * 2) % 100));
	}
	else
	{
		$xlim = int($histogram_max * 2 - (($histogram_max * 2) % 100) + 100);
	}
	if ($xlim <= 600){$xbin = 100;}
	elsif ($xlim <= 1400){$xbin = 200;}
	elsif ($xlim % 4 == 0){$xbin = 400;}
	elsif ($xlim % 3 == 0){$xbin = 300;}
	else {die 'Average depth is more than 2000x';}

	histPlot($outputDir,"$depth_frequency",$ylim,$ybin,$xlim,$xbin);
	cumuPlot($outputDir,"$cumu",$xlim,$xbin);
}
sub cumuPlot
{
	my ($outdir, $dataFile, $xlim, $xbin) = @_;
	my $figFile = "$outdir/cumuPlot.pdf";
	my $Rline=<<Rline;
	png(filename="$outdir/cumuPlot.png",width = 480, height = 360)
	rt <- read.table("$dataFile",header=T)
	opar <- par()
	x <- rt\$Depth[1:($xlim+1)]
	y <- 100*rt\$TRPercent[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);
	system("/opt/blc/genome/biosoft/R-212/bin/R CMD BATCH $figFile.R");
	system("rm  $figFile.R  `pwd`/cumuPlot.pdf.Rout");
}

sub histPlot
{
	my ($outdir, $dataFile, $ylim, $ybin, $xlim, $xbin) = @_;
	my $figFile = "$outdir/histPlot.pdf";
	my $Rline=<<Rline; 
	png(filename="$outdir/histPlot.png",width = 480, height = 360)
	rt <- read.table("$dataFile")
	opar <- par()
	t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
	y=c(rt\$V2[1:$xlim],t)
	y <- y*100
	x <- rt\$V1[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("/opt/blc/genome/biosoft/R-212/bin/R CMD BATCH $figFile.R");
	system("rm  $figFile.R  `pwd`/histPlot.pdf.Rout");
}

sub html_sub
{
    my ($searchFile, $outputPath, $tablesNum) = @_;	
	my @tables = ();
	open IN, "$outputPath/$searchFile";
	while (<IN>)
	{
		chomp;
		my @items = split /\t/;
		for my $key(1..$tablesNum)
		{
			if ($searchFile eq 'sampleInfo')
			{$tables[$key] .= "<tr><td width = '160' align = 'left'>&nbsp;$items[0]</td>";}
			else{$tables[$key] .= "<tr><td width = '295' align = 'left'>&nbsp;$items[0]</td>";}
			if($numInLine * $key <= $#items)
			{
				if ($searchFile eq 'sampleInfo')
				{map {$tables[$key] .= "<td width = '160' align = 'center'>$items[$_]</td>";}
				($numInLine * ($key-1) + 1)..($numInLine * $key);}
				else
				{map {$tables[$key] .= "<td width = '145' align = 'center'>$items[$_]</td>";}
                ($numInLine * ($key-1) + 1)..($numInLine * $key);}
			}
			else
			{
				if ($searchFile eq 'sampleInfo')
				{map {$tables[$key] .= "<td width = '160' align = 'center'>$items[$_]</td>";}
				($numInLine * ($key-1) + 1)..$#items;}
				else
				{map {$tables[$key] .= "<td width = '145' align = 'center'>$items[$_]</td>";}
				($numInLine * ($key-1) + 1)..$#items;}
			}
			$tables[$key] .= "</tr>";
		}
	}
	close IN;
	my $html_sub = '';
	for (1 .. $#tables)
	{
		$html_sub .= "<h5></h5><table border = '1' cellpdding = '2' cellspacing = '0' bordercolordark = '#DCDCDC' align = 'center'>".$tables[$_]."</table><br />";
	}
	return $html_sub;
}

sub distributionHTML
{
	my ($outputPath, $samplesID) = @_;
	my $html_sub = "<table align = 'center' cellpadding = '5'>";
	my @samplesID2 = @{$samplesID};
	my $sampleNum = @samplesID2;
	my $totalLines = $sampleNum * 2;
	foreach my $line (0 .. $totalLines)
	{
		if ($line == 0)
        {
			$html_sub .= "<tr><th width = '360' align= 'center'>Histogram of depth distribution in Target Regions</th><th width = '360' align = 'center'>Evenness of exome capture sequencing</th></tr>";
        }
		if ($line % 2 > 0)
        {
			my $index = int ($line / 2);
			$html_sub .= <<html_sub;
<tr>
    <td><span style = 'width:350px;display:block;'><img src = './$samplesID2[$index]/histPlot.png' alt = '$samplesID2[$index] histPlot' /></span></td>
    <td><span style = 'width:350px;display:block;'><img src = './$samplesID2[$index]/cumuPlot.png' alt = '$samplesID2[$index] cumuPlot' /></span></td>
</tr>
html_sub
		}
        if (($line > 0) && ($line % 2 == 0))
        {
                my $index = ($line / 2) - 1;
				$html_sub .= <<html_sub;
<tr>
<td colspan = '3' align = 'center'> Sequencing Depth distribution for sample <b><i>$samplesID2[$index]</i></b></td>
</tr>
html_sub
		}
	}
	$html_sub .= '</table>';
	return $html_sub;
}

`rm $outputPath/*Result $outputPath/sampleInfo $outputPath/*.txt`;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$year=$year+1900;
$mon=$mon+1;
my $date = "$year-$mon-$mday";
