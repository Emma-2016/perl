#! /usr/bin/perl -w
use strict;

$/="\n\n";
open IN,"$ARGV[0]";
my %relation=();
my %name=();
my @id=();
while (<IN>)
{
	next if $.==1;
	next if $_=~/is_obsolete:/;
	my @items=(split /\n/);
	my ($id,$parentID)=('','');
	foreach $_(@items)
	{
		if ($_=~s/^id:\s+//){$id=$_;}
		if ($_=~s/^name:\s+//){$name=$_;}
		if (($_=~s/^is_a:\s+//) and ($_=~s/\s+!.+//)){$parentID.=$_." ! ";}
	}
	$parentID=~s/\s+!\s+$//;
	$relation{$id}=$parentID;
	if($id=~/^HP:/){push @id,$id;}
}
close IN;
$/="\n";

A: foreach $_(@id)
{
	my $id=$_;
	my $out='';
	fac($id,$out,%relation);
}

sub fac
{
	my ($one,$out,%relation)=@_;
	$out.=$one."\t";
	if (defined $relation{$one})
	{
		if ($relation{$one}=~/\s+!\s+/)
		{
			my @p=split /\s+!\s+/,$relation{$one};
			foreach my $p(@p){fac($p,$out,%relation);}
		}
		else
		{
			$one=$relation{$one};
			fac($one,$out,%relation);
		}
	}   
	else
	{
	  $out=~s/\s+$//;
	  print $out."\n";
	}
}
