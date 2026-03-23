#!/usr/bin/env perl
#
# Script to aggregate values in columns from samples of MToolBox => VCF_file.vcf
# driven by a POS list, with NA filling for missing variants
#
# Last Modified: March/05/2025
#
# Copyright (C) 2025-2026 Manuel Rueda - CNAG

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Pod::Usage;

# -------------------------------
# Arguments
# -------------------------------
my $filein   = '';
my $varlist  = '';
my $debug    = 0;
my $man      = 0;
my $help     = 0;
my $delimiter = '|';

GetOptions(
    "input|i=s"    => \$filein,
    "variants|v=s" => \$varlist,
    "help|?"       => \$help,
    "man"          => \$man,
    "debug"        => \$debug
) or pod2usage(2);

pod2usage(
    -message => "ERROR: both --input and --variants are required\n",
    -exitval => 1,
    -verbose => 1
) if !$filein || !$varlist;

pod2usage( -exitval => 0, -verbose => 1 ) if $help;
pod2usage( -exitval => 0, -verbose => 2 ) if $man;

# -------------------------------
# Read variant list (POS only)
# -------------------------------
my @variants;
open my $vf, '<', $varlist;
while (<$vf>) {
    chomp;
    next if /^\s*$/;
    push @variants, $_;
}
close $vf;

# -------------------------------
# Read VCF into memory
# -------------------------------
my @bidi;
open my $fh, '<', $filein;
while ( my $line = <$fh> ) {
    next if $line =~ /^##/;          # skip meta headers
    chomp $line;
    push @bidi, [ split /\t/, $line ];
}
close $fh;

my $nrow = scalar @bidi;
my $ncol = scalar @{ $bidi[0] };

print "ROWS:$nrow COLS:$ncol\n" if $debug;

# -------------------------------
# Build NA template (once)
# -------------------------------
my @na_samples;
for ( my $j = 9 ; $j < $ncol ; $j++ ) {
    my $sample = $bidi[0][$j];
    $sample =~ s/-DNA_MIT//;
    $sample = substr( $sample, -3 );
    push @na_samples, "$sample:NA";
}

my $na_block = join( $delimiter, sort @na_samples );

# -------------------------------
# Index VCF by POS (duplicate check)
# -------------------------------
my %vcf_data;

for ( my $i = 1 ; $i < $nrow ; $i++ ) {

    my $pos = $bidi[$i][1];   # POS column

    die "ERROR: duplicated POS $pos found in VCF\n"
        if exists $vcf_data{$pos};

    my @tmp_GT;
    my @tmp_DP;
    my @tmp_HF;

    for ( my $j = 9 ; $j < $ncol ; $j++ ) {
        my $sample = $bidi[0][$j];
        $sample =~ s/-DNA_MIT//;
        $sample = substr( $sample, -3 );

        my @fields = split /:/, $bidi[$i][$j];
        my $GT = $fields[0] // 'NA';
        my $DP = $fields[1] // 'NA';
        my $HF = $fields[2] // 'NA';

        push @tmp_GT, "$sample:$GT";
        push @tmp_DP, "$sample:$DP";
        push @tmp_HF, "$sample:$HF";
    }

    $vcf_data{$pos} = {
        REF => $bidi[$i][3],
        ALT => $bidi[$i][4],
        GT  => join( $delimiter, sort @tmp_GT ),
        DP  => join( $delimiter, sort @tmp_DP ),
        HF  => join( $delimiter, sort @tmp_HF ),
    };
}

# -------------------------------
# Output driven by POS list
# -------------------------------
print join( "\t", 'REF', 'ALT', 'GT', 'DP', 'HF' ) . "\n";

for my $pos (@variants) {

    if ( exists $vcf_data{$pos} ) {
        print join "\t",
            $vcf_data{$pos}{REF},
            $vcf_data{$pos}{ALT},
            $vcf_data{$pos}{GT},
            $vcf_data{$pos}{DP},
            $vcf_data{$pos}{HF};
        print "\n";
    }
    else {
        print join "\t",
            'NA',
            'NA',
            $na_block,
            $na_block,
            $na_block;
        print "\n";
    }
}

exit 0;

__END__

=head1 NAME

parse_prioritized.pl - Aggregate per-sample VCF fields driven by a POS list

=head1 SYNOPSIS

parse_prioritized.pl -i input.vcf -v variants.txt [options]

=head1 DESCRIPTION

This script reads a MToolBox VCF file and aggregates per-sample
GT, DP, and HF fields. Output order is driven by a one-column
list of variant positions (POS).

If a POS from the list is not present in the VCF, the script
outputs NA values for all samples.

If a POS occurs more than once in the VCF, the script exits
with an error.

=head1 ARGUMENTS

=over 4

=item B<-i, --input>

Input VCF file (MToolBox format)

=item B<-v, --variants>

One-column file with POS values in the desired output order

=back

=head1 OPTIONS

=over 4

=item B<-h, --help>

Print brief help message

=item B<--man>

Print full documentation

=item B<--debug>

Print debugging information

=back

=head1 AUTHOR

Manuel Rueda <manuel.rueda@cnag.eu>

=cut

