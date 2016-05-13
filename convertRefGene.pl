#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

die "Usage:\nperl $0 <in|hg19_refGene.20110731.txt.gz> <out|refGene.1.gff> <outDir|refSeqByChr>\nFor Example:\nperl $0 $Bin/hg19_refGene.20110731.txt.gz $Bin/refGene.1.gff $Bin/refSeqByChr\n" unless @ARGV==3;

my ($refseq,$out,$refSeqByChr)=@ARGV[0,1,2];

$refseq||="/ifs2/BC_MD/GROUP/fujin/hg19_GeneData/CCDS_Refseq_Ensembl_snp/hg19_refGene.20110731.txt.gz";
$out||="$Bin/refGene.1.gff";

open OT,">$out" or die $!;
open IN,"<$refseq" or die $! if $refseq!~/\.gz$/;
open IN,"gunzip -cd <$refseq|" or die $! if $refseq=~/\.gz$/;
#671     NM_002761       chr16   -       11282193        11282693        11282349        11282596        2       11282193,11282484,      11282393,11282693,      0       PRM1    cmpl     cmpl    1,0,
#chr16   87451007        87451059        refSeq  5-UTR   +       name=TRAPPC2L; ID=NM_016209; exonNumber=1;
#chr16   87451060        87451092        refSeq  CDS     +       name=TRAPPC2L; ID=NM_016209; exonNumber=1;
#chr16   87451093        87452527        refSeq  intron  +       name=TRAPPC2L; ID=NM_016209; intnNumber=1;
#chr16   87452528        87452700        refSeq  CDS     +       name=TRAPPC2L; ID=NM_016209; exonNumber=2;
while(<IN>)
{
	chomp;
	###########shenyulan@genomics.cn##############

	next unless $_ !~/incmpl/;

	#############################################
	my ($ID,$chr,$strand,$mRNAsta,$mRNAend,$cds_st,$cds_ed,$exon_sts,$exon_eds,$geneName,$phases,$exonCount)=(split /\t/)[1,2,3,4,5,6,7,9,10,12,-1,8];
	next unless $chr=~/^chr\d+$/ || $chr=~/^chr[XY]$/;
	my @exon_st=split /,/,$exon_sts;
	my @exon_ed=split /,/,$exon_eds;
	my ($exonNum,$intronNum)=(0,0);
	my $utr="5-UTR";
	$utr="3-UTR" if $strand eq "-";
	for(my $index=0;$index<$exonCount;$index++)
	{
		($exonNum,$intronNum)=($index+1,$index);
		($exonNum,$intronNum)=($exonCount-$index,$exonCount-$index) if $strand eq "-";
		if($index>0)
		{
			print OT "$chr\t".($exon_ed[$index-1]+1)."\t$exon_st[$index]\trefSeq\tintron\t$strand\tname=$geneName; ID=$ID; intnNumber=$intronNum;\n";
		}
		if($exon_st[$index]<$cds_st)
		{
			if($exon_ed[$index]<$cds_st+1)
			{
				print OT "$chr\t".($exon_st[$index]+1)."\t$exon_ed[$index]\trefSeq\t$utr\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
			}elsif($exon_ed[$index]>=$cds_st+1 && $exon_ed[$index]<=$cds_ed)
			{
				print OT "$chr\t".($exon_st[$index]+1)."\t$cds_st\trefSeq\t$utr\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
				print OT "$chr\t".($cds_st+1)."\t$exon_ed[$index]\trefSeq\tCDS\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
			}else
			{
				print OT "$chr\t".($exon_st[$index]+1)."\t$cds_st\trefSeq\t$utr\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
				print OT "$chr\t".($cds_st+1)."\t$cds_ed\trefSeq\tCDS\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
				$utr="3-UTR";
				$utr="5-UTR" if $strand eq "-";
				print OT "$chr\t".($cds_ed+1)."\t$exon_ed[$index]\trefSeq\t$utr\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
			}
		}elsif($exon_st[$index]>=$cds_st && $exon_st[$index]<=$cds_ed)
		{
			if($exon_ed[$index]<=$cds_ed)
			{
				print OT "$chr\t".($exon_st[$index]+1)."\t$exon_ed[$index]\trefSeq\tCDS\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
			}else
			{
				$utr="3-UTR";
				$utr="5-UTR" if $strand eq "-";
				print OT "$chr\t".($exon_st[$index]+1)."\t$cds_ed\trefSeq\tCDS\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
				print OT "$chr\t".($cds_ed+1)."\t$exon_ed[$index]\trefSeq\t$utr\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
			}
		}else
		{
			$utr="3-UTR";
			$utr="5-UTR" if $strand eq "-";
			print OT "$chr\t".($exon_st[$index]+1)."\t$exon_ed[$index]\trefSeq\t$utr\t$strand\tname=$geneName; ID=$ID; exonNumber=$exonNum;\n";
		}
	}
}
close OT;
close IN;

print "$out\n";

$refSeqByChr||="/ifs2/BC_MD/GROUP/fujin/hg19_GeneData/CCDS_Refseq_Ensembl_snp/refSeqByChr";
`mkdir $refSeqByChr` unless -d "$refSeqByChr";
`awk '{print \$0 >"$refSeqByChr/"\$1}' $out`;
`if [[ -f $Bin/aaaa ]];then rm $Bin/aaaa;fi`;
`for i in {1..22} X Y;do cat $refSeqByChr/chr\$i >>$Bin/aaaa;done`;
`mv $Bin/aaaa $out`;
