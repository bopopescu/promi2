#!/usr/bin/perl -w

use strict;

open(MFILE, $ARGV[0]) or die;
open(MLFILE, $ARGV[1]) or die;
open(CFILE, $ARGV[2]) or die;
open(FFILE, $ARGV[3]) or die;
open(LFILE, $ARGV[4]) or die;
my $datastart = $ARGV[5];

my %counts = ();
my %factors = ();
my $count;
my $factor;
while (<LFILE>)
{
	chomp;
	chomp($count = <CFILE>);
	chomp($factor = <FFILE>);
	$counts{$_} = $count;
	$factors{$_} = $factor;
}

## Precalculate tpm multiplicand
my @tpm = ();
my $tpm = 0;
while (<MLFILE>)
{
    chomp;
	$tpm = defined($counts{$_}) ? (1e6)/($counts{$_} * $factors{$_}) : 0;
    push(@tpm, $tpm);
}
close(MLFILE);
close(CFILE);
close(FFILE);
close(LFILE);

my @fields = ();
## Transform each tagcount per TC to tpm
my $line = "";
while (<MFILE>)
{
    chomp;
    @fields = split(/\t/, $_);
    for (my $i = $datastart; $i < $datastart + @tpm; $i++)
    {
		$fields[$i] = $tpm[$i-$datastart] != 0 ? $fields[$i] * $tpm[$i-$datastart] : "NA";
    }
	print $fields[$datastart];
	for (my $i=($datastart+1); $i<@fields; $i++)
	{
		print "\t".$fields[$i];
	}
	print "\n";
}
close(MFILE);
