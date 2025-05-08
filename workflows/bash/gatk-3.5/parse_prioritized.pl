#!/usr/bin/env perl
#
#   Script to aggregate values in columns from samples of MToolBox => VCF_file.vcf
#
#   Last Modified; March/05/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

#use JSON;

# Defining a few variables
my $filein = '';
my $debug  = 0;
my $man    = 0;
my $help   = 0;
my $delimiter = '|';

# Reading arguments
GetOptions(
    "input|i=s" => \$filein,    # numeric
    'help|?'    => \$help,      # flag
    "man"       => \$man,       # flag
    "debug"     => \$debug      # flag
) or pod2usage(2);
pod2usage(1) if !$filein;
pod2usage(1) if $help;
pod2usage( -exitval => 0, -verbose => 2 ) if $man;

# We will store the whole file in a bidimensional array
my @bidi = ();

# We print the header
print join( "\t", 'REF', 'ALT', 'GT', 'DP', 'HF' ) . "\n";

# Reading and parsing SG-Adviser file
open my $fh, '<', $filein;

# Ready-Go!
while ( defined( my $line = <$fh> ) ) {

# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MA0005602M-DNA_MIT	MA0005601P-DNA_MIT	MA0005603F-DNA_MIT
# chrMT	347	.	G	GC,C,T	.	PASS	AC=1,1,1;AN=6	GT:DP:HF:CILOW:CIUP	0	0/1/2/3:2069:0.004,0.006,0.002:0.002,0.003,0.001:0.008,0.01,0.006	0
    chomp $line;
    push( @bidi, [ split /\t/, $line ] );
}

# Get the dimensions of the array
my $nrow = scalar @bidi;
my $ncol = scalar @{ $bidi[0] };
print "ROWS:$nrow COLS:$ncol\n" if $debug;

# Start printing the elements
for ( my $i = 1 ; $i < $nrow ; $i++ ) {
    print "$bidi[$i][3]\t$bidi[$i][4]\t";
    my @tmp_GT = ();
    my @tmp_DP = ();
    my @tmp_HF = ();
    for ( my $j = 9 ; $j < $ncol ; $j++ ) {
        my $sample = $bidi[0][$j];
        $sample =~ s/-DNA_MIT//;
        $sample = substr( $sample, -3 );
        my @tmp_fields = ( split /:/, $bidi[$i][$j] );
        my $GT         = $tmp_fields[0];
        my $DP         = $tmp_fields[1] // 'N/A';
        my $HF         = $tmp_fields[2] //  'N/A';
        push @tmp_GT, $sample . ':' . $GT;
        push @tmp_DP, $sample . ':' . $DP;
        push @tmp_HF, $sample . ':' . $HF;
    }
    print join( $delimiter, sort @tmp_GT ) . "\t";    #sorted alphabetically
    print join( $delimiter, sort @tmp_DP ) . "\t";    #sorted alphabetically
    print join( $delimiter, sort @tmp_HF ) . "\n";    #sorted alphabetically

}

=head1 NAME

parse_prioritized: cript to aggregate values in columns from samples of MToolBox => VCF_file.vcf

=head1 SYNOPSIS


parse_prioritized.pl -i vcf_file [-options]

     Arguments:                       
       -i|input                       MToolBox file

     Options:
       -h|help                        Brief help message
       -man                           Full documentation
       -debug                         Print debugging (from 1 to 5, being 5 max)
