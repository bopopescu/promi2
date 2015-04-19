#!/usr/bin/env perl

# This script takes a wiggle style binary file (i.e. phastcons track) and a GFF
# file and extracts the wiggle information for the regions in the GFF. To create
# the binary wiggle files needed by this program use the script wigFix2wib.pl or
# the Bio::Graphics::Wiggle::Loader module directly.
#
# !!Warning!! There are at least two different .wib formats! One is created by
# the BioPerl's Bio::Graphics::Wiggle::Loader module and the other is created by
# the Kent source tools that are used at UCSC. If you have directly downloaded a
# .wib file from UCSC it is likely UCSC/Kent-style and can not be used with
# Bio::Graphics::Wiggle (which is what this script uses). See for more details:
#
# http://search.cpan.org/~lds/Bio-Graphics-2.31/lib/Bio/Graphics/Wiggle.pm#___top
# http://genomewiki.ucsc.edu/index.php/Wiggle
#
# Example GFF file (MUST BE TAB-SEPARATED!!!):
# 
# chr1	TeleGene	enhancer	10918	11000	500	+	.	touch1
# chr1	TeleGene	promoter	11500	11700	900	+	.	touch1
#
# Note: I cut out a lot of unnecessary stuff from this script. In the future
# switch to using .bed files since they can contain only the information we need
# (chr, start and stop), and also allow one to look up a region by giving the
# coordinates directly on the command line. It should also be able to find the
# correct .wib file for the chromosome based on a simple naming convention,
# rather than having to specify the .wib file directly (just give a directory
# name instead).
#
# Original by ?
# Modified by Matt Huska (huska@molgen.mpg.de)
# Last updated 2012-10-18

use lib "/project/genetics/packages/x86_64/perl/Bio-Graphics-1.90/lib";
use Bio::Graphics::Wiggle;

$phastcons_file = $ARGV[0];
$region_file = $ARGV[1];

print STDERR "phastcons file: $phastcons_file\n";
print STDERR "region file: $region_file\n";

$wig = Bio::Graphics::Wiggle->new($phastcons_file,0);

# read through the region file (gff) and extract the values
open(REGIONS,"<".$region_file) || die "Can not open file $region_file : $!\n";;

$type = "";
$idlist = "";

while (<REGIONS>) {
    next if (/^\#/);

    @col=split /[\t\n]/;
    # $chr = $col[0]; # FIXME: Use this later to automatically open the correct .wib file
    $start = $col[3];
    $stop = $col[4];
    $id = $col[8];

    my @values = @{ $wig->values($start, $stop) };
    print "$id, ".join(", ", @values)."\n";
}
