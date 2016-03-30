#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);

my ($refseq, $out, $refSeqByChr) = @ARGV[0, 1, 2];

open OUT, ">$out" or die $!;
open IN, "$refSeq" or die $! if $refseq !~ /\.gz$/;
open IN, "gunzip -dc $refseq|" or die $! if $refseq =~ /\.gz$/;

while(<IN>)
{
  chomp;
  next if ($_ =~ /imcmpl/ or $_ =~ /none/ or $_ =~ /unk/);
  my ($bin, $id, $chr, $strand, $rnaStart, $rnaEnd, $cdsStart, $cdsEnd, $exonCnt, $exonStarts, $exonEnds, $score, $geneName, $cdsStartStat, $cdsEndStat, exonFrames) = (split /\t/);
  next unless ($chr =~ /^chr\d+$/ or $chr =~ /^chr[XY]$/);
  my @exonStarts = split /,/, $exonStarts;
  my @exonEnds = split /,/, $exonEnds;
  my ($exonNo, $intronNo) = (0, 0);
  my $utr = "5-UTR";
  $utr = "3-UTR" if $strand eq "-";
  
  for (my $index = 0; $index < $exonCount; $index++)
  {
    ($exonNo, $intronNo) = ($index+1, $index);
    ($exonNo, $intronNo) = ($exonCount-$index, $exonCount-$index) if $strand eq "-";
    if ($index > 0)
    {
      print OUT "$chr\t($exonEnds[$index-1]+1)\t$exonStarts[$index]\trefSeq\tintron\t$strand\tname=$geneName; ID=$id; intnNumber=$intronNum;\n";
    }
    if ($exonStarts[$index]


