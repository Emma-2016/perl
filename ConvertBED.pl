#! /usr/bin/perl -w
use strict;

my $usage=<<Usage;
This script is to convert original BED files into those that can be used directly for alignment or variant calling.
Usage

die ($usage) unless @ARGV>0;
my ($bed, $outDir, $chipName)=@ARGV;
`mkdir -p $outDir/$chipName\_0bp` unless (-e "$outDir/$chipName\_0bp"); #use back slash to start the constant
`mkdir -p $outDir/$chipName\_200bp` unless (-e "$outDir/$chipName\_200bp");

`awk -F "\t" '{if (\$1~/^chr/)print \$1"\t"\$2"\t"\$3 >"$outDir/"\$1".tmp"}' $bed`;#awk concact string immediately;put string inside "".
#awk -F "\t" '{if ($1~/^chr/)print $1"\t"$2"\t"$3 >$1".tmp"}' S07604514_Regions.bed 

foreach my $i (1..22,"X","Y")
{
    `awk -F "\t" '{print \$1,\$2-0,\$3+0}' "$outDir/chr$i.tmp" | awk 'BEGIN{b=0}{if(\$2>b){print \$0;b=\$3}}' | awk 'BEGIN{a=0;b=-1}{if(\$2>b+1){if(NR>1){print \$1"\t"a"\t"b >"$outDir/$chipName\_0bp/"\$1".bed"};a=\$2;b=\$3};b=\$3}END{print \$1"\t"a"\t"b >"$outDir/$chipName\_0bp/"\$1".bed"}'`;
    #awk do evaluation;
    `awk -F "\t" '{print \$1,\$2-200,\$3+200}' "$outDir/chr$i.tmp" | awk 'BEGIN{b=0}{if(\$2>b){print \$0;b=\$3}}' | awk 'BEGIN{a=0;b=-1}{if(\$2>b+1){if(NR>1){print \$1"\t"a"\t"b >"$outDir/$chipName\_200bp/"\$1".bed"};a=\$2;b=\$3};b=\$3}END{print \$1"\t"a"\t"b >"$outDir/$chipName\_200bp/"\$1".bed"}'`;
    `rm "$outDir/chr$i.tmp"`;
}
