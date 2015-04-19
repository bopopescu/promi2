#!/usr/bin/perl -w

use strict;

open(IN, $ARGV[0]) or die;

while (<IN>)
{
	if ($_ =~ /(chr[^\t]+):(\d+)\.\.(\d+),([+-]).*/)
	{
		print $1."\t".$2."\t".$3."\t".$1.":".$2."..".$3.",".$4."\t.\t".$4."\n";
	}
}
close(IN);
