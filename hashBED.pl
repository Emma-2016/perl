#! /usr/bin/pelr -w
use strict;

use File::Basename qw(basename dirname);
use Getopt::Long;
use Data::Dumper;

my $samtools = "path/to/samtools";
my ($region, $bam, $RefCCDS, $outDir, $help);
my ($BestOff, $SecBestOff) = (15, 5); #alignment score assigned by BWA
my $flank = 200;  #length of flank region

GetOptions(
  "flank:i" => \$flank,
  "region:s" => \$region,
  "bam:s" => \$bam,
  "OutDir:s" => \$outDir,
  "BestOff:i" => \$BestOff,
  "SecBestOff:i" => \$SecBestOff,
  "help!" => \$help
);

die `pod2text $0` if ($help or !$region or !$bam or !$outDir);  #pod2text is to print out description
my $chrs = basename $region;
$chrs =~ s/\.bed// if $chr =~ /\.bed/;

#-----------------------------------------------------------
my $CCDSStatFile = "$outDir/$chr\_CCDS.stat";
my $CCDSGffFile = "$outDir/$chr\_CCDS.gff.gz";
my $TRDepthStatFile = "$outDir/$chr\_DepthByChr.stat";
my $TRDepthFile = "$outDir/$chr\_TRDepth.gz";
my $FRDepthFile = "$outDir/$chr\_FRDepth.gz";
my $TRDepthCntFile = "$outDir/$chr\_TRDepthCnt.stat";
my $FRDepthcntFile = "$outDir/$chr\_FRDepthCnt.stat";
my $CumDepthFile = "$outDir/$chr\_CumDepth.stat";


open RG, $region or die $!;
my (%TR, %FR) = ((),());
while (<RG>)
{
  chomp;
  my ($start, $end) = (split /[:-\t]/)[1,2];  #put the char into [] instead of "|";
  $TR{$_} = 0 foreach ($start..$end);
}
seek (RG, 0, 0);

while (<RG>)
{
  chomp;
  my ($start, $end) = (split /[:-\t]/)[1,2];
  foreach (($start-$flank)..$start)
  {
    next if exists $TR{$_};
    $FR{$_};
  }
  foreach (($end+1)..($end+$flank))
  {
    next if exists $TR{$_};
    $FR{$_};
  }
}
close RG;
print "This chromosome BED is hash over...\n";

#----------start hash the ccds file----------------#
#chr1    801943  802434  CCDS    CDS    -       name=NCRNA00115; ID=CCDS1.1; exonNumber=1; geneID=79854;
my (%TRCCDS, %ReadCCDS)=((),());
if ($RefCCDS) #chromosome ref gene file
{
  open CCDS, "$RefCCDS" or die $!;
  while (<CCDS>)
  {
    chomp;
    my ($start, $end, $lable) = (split /\t/)[1,2,4];
    next if $lable ne "CDS";
    foreach ($start..$end)
    {
      $TRCCDS[$_] = 0;
      $ReadCCDS[$_] = 0;  #CCDS covered by sequenced reads
      $TRCCDS[$_]++ if (exists $TRpos{$_}); #whether the CCDS is covered by bed; #TRCCDS[$_] = 1
    }
  }
  close CCDS;
}
print "CCDS file hash over...\n";

open IN,"$samtools view -F 4 $bam|" or die $!;
my ($ReadCnt, $TRReadCnt, $FRReadCnt, $UReadCnt, $UTRReadCnt, $UFRReadCnt) = (0, 0, 0, 0, 0, 0);
my ($Base, $UBase, $FRBase, $UFRBase, $TRBase, $UTRBase) = (0, 0, 0, 0, 0, 0);
my ($MissNum, $GCBase) = (0, 0);

##A819WHABXX:5:1102:20929:38038#GCCAATAT 65      chr21   9411193 11      2S48M40S        chr11   102573635       0TGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGAGTCACTCCACTAAAATTCACCAAGATTTCAAAGGGGATTG     ?B?>?BFECFCEFCDDDFBDDCCGBCGBDDDGFDGFHEFFDCGCFEBCGCGAGEFFGDEGFBEFFEEHEGDFCH@DEEHEFEGFEDA@BF      SA:Z:chr11,102573731,-,41M49S,11,0;     BD:Z:NNQNTSTPOOPPOPPNEOPOENNONPPQNFFOONOORQQRNNQQRQPPPQRRRQOPOPQOPPNFFOOOQOQQOPRRPGPRQHRUQRQQOP MD:Z:48 PG:Z:MarkDuplicates     RG:Z:RGID_110222_I288_FC819WHABXX_L5_EXONPEP00002020    BI:Z:RRSPVTUSPRRTQRTQISSRISPSRSTUQIIRRPRRRRPOPQTSTRTSTSUUSTSSRRTTSSPIJTPSTTUURTTTQJSURJSQOQUTQT NM:i:0  AS:i:48 XS:i:45
while (<IN>)
{
  chomp;
  ($start, $cigar, $ReadSeq) = (split /\t/)[3, 5, 9];
  my ($best, $SecBest) = (0, 0);
  $best = $1 if (/AS:i:(\d+)/);
  $SecBest = $1 if (/XS:i:(\d+)/;
  $MapTimes = 1 if ($best >= $BestOff and $SecBest < $SecBestOff);
  $miss = $1 =~tr/ATCG/ATCG/ if (/MD:Z:(\w+)/); #what if it is deletion or insertion?
  $MapLen = length $ReadSeq;
  next if $cigar=~/[IDS]/;
  
  $MissNum += $miss;
  $GCBase += ($ReadSeq =~ tr /GC/GC/);
  my ($TRMark, $FRMark) =(0, 0);
  if ($MapTimes == 1)
  {
    $UReadCnt++;
    $Ubase += $MapLen;
    for my $p ($start .. ($start+$MapLen-1))  #the start is 0-based; bed file is also 0-based, but gff may be not. !!
    {
      if (exists $TRpos{$p})
      {
        $TRpos{$p}++;
        $UTRBase{$p}++;
        $TRMark=1;
      }
      if (exists $FRpos{$p})
      {
        $FRpos{$p}++;
        $UFRBase++;
        $FRMark=1;
      }
      if ($RefCCDS and exists $ReadCCDS{$p}) #whether this base in CCDS cds
      {
        $ReadCCDS{$p}++;
      }
    }
    $UTRReadCnt++ if ($TRMark);
  }
  elsif ($MapTimes > 1)
  {
    $ReadCnt++;
    $Base += $MapLen;
    my $end = $start+$MapLen-1;
    for my $p ($start .. $end)
    {
      if (exists $TRpos{$p})
      {
        $TRpos{$p}++;
        $TRBase++;
        $TRMark = 1;
      }
      elsif (exists $FRpos{$p})
      {
        $FRpos{$p}++;
        $FRBase++;
        $FRMark=1;
      }
      if($RefCCDS and exists $ReadCCDS{$p})
      {
        $ReadCCDS{$p}++;
      }
    }
    $TRReadCnt++ if ($TRMark);
    $FRReadCnt++ if ($FRMark);
  }
  else
  {
    print "Something wrong with the bam file...\n";
  }
}
close IN;
$Readcnt += $UReadCnt;
$TRReadCnt += $UTRReadCnt;
$FRReadCnt += $UFRReadCnt;
$Base += $UBase;
$TRBase += $UTRBase;
$FRBase += $UFRBase;

my (%ID, %TRID, %ReadID, %Gene, %TRGene, %ReadGene);
my ($ExonCnt, $TRExonCnt, $ReadExonCnt) = (0, 0, 0);
my ($ccds4, $ccds10, $ccds20) = (0, 0, 0);
if ($RefCCDS)
{
  open CCDS, "$RefCCDS" or die $!;
  open CCDSGFF, "|gzip >$CCDSGffFile" or die $!;
  while (<CCDS>)
  {
    chomp;
    my ($start, $end, $lable, $express) = (split /\t/)[1,2,4,6];
    next if ($lable ne 'CDS');
    my ($ProbeCoved, $ReadCoved, $ReadBase) = (0, 0, 0);
    my $length = $end - $start + 1; #here may be wrong;
    foreach my $pos ($start .. $end)
    {
      $ProbeCoved++ if ($TRCCDS{$pos});
      if ($ReadCCDS{$pos})
      {
        $ReadCoved++;
        $ReadBase += $ReadCCDS{$pos};
      }
    }
    my $Probecovage = sprintf ("%u", 100*$ProbeCoved/$length);
    my $ReadCovage = sprintf ("%u", 100*$ReadCoved/$length);
    my $MeanDep = sprintf ("%.2f", $ReadBase/$length);
    print CCDSGFF "$_; probeCov=$Probecovage; readCov=$readcovage; $meanDepth=$MeanDep;\n"
  
    #chr1    801943  802434  CCDS    CDS    -       name=NCRNA00115; ID=CCDS1.1; exonNumber=1; geneID=79854;
    my ($gene, $id) = $express =~ /name=(\S+);\s+ID=(\S+);/;
    $Exoncnt++;
    $ID{$id} = 1;
    $Gene{$gene} = 1;
    if ($Probecovage > 0)
   {
      $TRExonCnt++;
      $TRID{$id} = 1;
      $TRGene{$gene} = 1;
    }
    if ($Readcovage > 0)
    {
      $ReadExonCnt++;
      $ReadID{$id} = 1;
      $ReadGene{$gene} = 1;
      $ccds4++ if ($MeanDep >= 4);
      $ccds10++ if ($MeanDep >= 10);
      $ccds20++ if ($MeanDep >= 20);
    }
  }
(%TRCCDS, %ReadCCDS) = ((), ());
close CCDS;
close CCDSGFF;
}
print "$chr CCDS 2 hash over ...\n\n";
########--------------End-----------------########

########--------------Start-----------------########
if($RefCCDS)
{       open CCDSSTAT,">$CCDSStatFile" or die $!;
        my $IDNum=keys %ID;
        my $IDNum1=keys %TRID;
        my $IDNum2=keys %ReadID;
        my $GeneNum=keys %Gene;
        my $GeneNum1=keys %TRGene;
        my $GeneNum2=keys %ReadGene;
	print CCDSSTAT "#Item\tTotalNum\tProbeCovNum\tReadCovNum\tDepth>=4x\tDepth>=10x\tDepth>=20x\n";
    	print CCDSSTAT "$chr\_Exon\t$ExonCnt\t$TRExonCnt\t$ReadExonCnt\t$ccds4\t$ccds10\t$ccds20\n";
	print CCDSSTAT "$chr\_Eon_ID\t$IDNum\t$IDNum1\t$IDNum2\t-\t-\t-\n";
	print CCDSSTAT "$chr\_Gene\t$GeneNum\t$GeneNum1\t$GeneNum2\t-\t-\t-\n";
	close CCDSSTAT;
}
########--------------End-----------------########

# -------------------------------------- output depth ----------------------------------

my ($TRCov,$FRCov,$TRBase)=(0,0,0);
my (%TRDepthCnt,%FRDepthCnt,%TRCumDepthCnt,%FRCumDepthCnt);
open TRD,"|gzip >$TRDepthFile" or die "Can not open $TRDepthFile\n";
foreach my $k( sort {$a<=>$b} keys %TRpos )
{	my $depth = $TRpos{$k};
	$TRCov++ if($depth);
	$TRDepthCnt{$depth}++;
	$TRBase += $depth;
	print TRD "$k\t$depth\n";
}
close TRD;
my  $TRSize = scalar keys %TRpos;
%TRpos=();
my $TRMeanDepth = sprintf("%.2f",$TRBase/$TRSize);
my $TRCoverage = sprintf("%.2f",100*$TRCov/$TRSize);
my  ($TRCount,$TRMedianDepth)=(0,0);
foreach (sort {$a<=>$b} keys %TRDepthCnt)
{       $TRCount +=$TRDepthCnt{$_};
	if($TRCount*2 >= $TRSize)
	{	$TRMedianDepth=$_;	
		last;
	}
}

open FRD,"|gzip >$FRDepthFile" or die "Can not open $FRDepthFile\n";
foreach my $k( sort {$a<=>$b} keys %FRpos)
{	my $depth = $FRpos{$k};
	$FRCov++ if($depth);
	$FRDepthCnt{$depth} ++;
	$FRBase += $depth;
	print FRD "$k\t$depth\n";
}
close FRD;
my  $FRSize = scalar keys %FRpos;
%FRpos=();
my $FRMeanDepth = sprintf("%.2f",$FRBase/$FRSize);
my $FRCoverage = sprintf("%.2f",100*$FRCov/$FRSize);
my  ($FRCount,$FRMedianDepth)=(0,0);
foreach (sort {$a<=>$b} keys %FRDepthCnt)
{       $FRCount +=$FRDepthCnt{$_};
	if($FRCount*2 >= $FRSize)
	{       $FRMedianDepth=$_;
		last;
	}
}

my $TRFRSize=$TRSize+$FRSize;
open CUM,">$CumDepthFile" or die $!;
print CUM "#Chr\tDepth\tTRCumDepthCnt\tFRCumDepthCnt\tTRFRCumDepthCnt\n";
$TRCumDepthCnt{0}+=$TRDepthCnt{$_} for (keys %TRDepthCnt);
$FRCumDepthCnt{0}+=$FRDepthCnt{$_} for (keys %FRDepthCnt);
my $TRFRCumDepthCnt=$TRCumDepthCnt{0}+$FRCumDepthCnt{0};
print CUM "$chr\t0\t$TRCumDepthCnt{0}\t$FRCumDepthCnt{0}\t$TRFRCumDepthCnt\n";
for(1..200)
{	my $j=$_-1;
	$TRDepthCnt{$j}=0 if (!exists $TRDepthCnt{$j}); #--------------------
	$FRDepthCnt{$j}=0 if (!exists $FRDepthCnt{$j});  #--------------------
	$TRCumDepthCnt{$_}=$TRCumDepthCnt{$j}-$TRDepthCnt{$j};	# cumulative depth;
	$FRCumDepthCnt{$_}=$FRCumDepthCnt{$j}-$FRDepthCnt{$j};
	$TRFRCumDepthCnt=$TRCumDepthCnt{$_}+$FRCumDepthCnt{$_};
	print CUM "$chr\t$_\t$TRCumDepthCnt{$_}\t$FRCumDepthCnt{$_}\t$TRFRCumDepthCnt\n";
}
close CUM;
open TRDC,">$TRDepthCntFile" or die $!;
print TRDC "$_\t$TRDepthCnt{$_}\n" for (sort {$a<=>$b} keys %TRDepthCnt);
close TRDC;
open FRDC,">$FRDepthCntFile" or die $!;
print FRDC "$_\t$FRDepthCnt{$_}\n" for (sort {$a<=>$b} keys %FRDepthCnt);
close FRDC;
(%TRDepthCnt,%FRDepthCnt,%TRCumDepthCnt,%FRCumDepthCnt)=((),(),(),());	 
open TRDSTAT,">$TRDepthStatFile" or die $!;
print TRDSTAT "#Item\tRegionSize\tReadNum\tUniqReadNum\tBaseNum\tUniqBaseNum\tCovPosNum\tCovrege\tMeanDepth\tMedianDepth\tMissNuM\tGCBaseNum\n";
print  TRDSTAT "$chr\_Region\t-\t$ReadCnt\t$UReadCnt\t$Base\t$UBase\t-\t-\t-\t-\t$MissNum\t$GCBase\n";
print  TRDSTAT "$chr\_TargetRegion\t$TRSize\t$TRReadCnt\t$UTRReadCnt\t$TRBase\t$UTRBase\t$TRCov\t$TRCoverage\t$TRMeanDepth\t$TRMedianDepth\t-\t-\n";
print  TRDSTAT "$chr\_FlankRegion\t$FRSize\t$FRReadCnt\t$UFRReadCnt\t$FRBase\t$UFRBase\t$FRCov\t$FRMeanDepth\t$FRCoverage\t$FRMedianDepth\t-\t-\n";
close  TRDSTAT;
`date`;
