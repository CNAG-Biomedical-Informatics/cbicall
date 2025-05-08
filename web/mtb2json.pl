#!/usr/bin/env perl
#
#   Script to read and filter MToolBox prioritized variants
#   The output can be HASH | JSON | JSON4HTML | TSV
#
#   Last Modified; Dec/08/2016
#
#   Version: 1.0.5
#
#   Copyright (C) 2016 Manuel Rueda (mrueda@scripps.edu)
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, see <https://www.gnu.org/licenses/>. 
#
#   If this program helps you in your research, please cite

use strict;
use warnings;
use autodie;
use FindBin qw($Bin);
use lib $Bin;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use JSON::XS;
use Text::CSV::Slurp;

# Defining a few variables
my $filein     = '';
my $debug      = 0;
my $denovo     = 0;
my $hf_cutoff  = 0;
my $maf_cutoff = 0;
my $man        = 0;
my $help       = 0;
my $format     = 'hash';
my $sga_dir    = $Bin;

# Reading arguments
GetOptions(
    "input|i=s"  => \$filein,        # numeric
    "denovo"     => \$denovo,        # flag
    "HF=f"       => \$hf_cutoff,     # real
    "MAF=f"      => \$maf_cutoff,    # real
    "format|f=s" => \$format,        # string
    'help|h'     => \$help,          # flag
    "man"        => \$man,           # flag
    "debug"      => \$debug          # flag
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitval => 0, -verbose => 2 ) if $man;
pod2usage(1) if !$filein;

# Defining the position of the fields in the header
# head -1 prioritized_variants.txt   |sed 's/ /_/g' |t2n |nl -v 0 | awk '{print $2,$1}' | sed "s/^/'/" | sed "s/ /' => /" | sed 's/$/,/'
my %mtb_data_loc = (
    'Variant_Allele'                => 0,
    'Sample'                        => 1,
    'Locus'                         => 2,
    'Nt_Variability'                => 3,
    'Codon_Position'                => 4,
    'Aa_Change'                     => 5,
    'Aa_Variability'                => 6,
    'tRNA_Annotation'               => 7,
    'Disease_Score'                 => 8,
    'RNA_predictions'               => 9,
    'Mitomap_Associated_Disease(s)' => 10,
    'Mitomap_Homoplasmy'            => 11,
    'Mitomap_Heteroplasmy'          => 12,
    'Somatic_Mutations'             => 13,
    'SM_Homoplasmy'                 => 14,
    'SM_Heteroplasmy'               => 15,
    'ClinVar'                       => 16,
    'OMIM_link'                     => 17,
    'dbSNP_ID'                      => 18,
    'Mamit-tRNA_link'               => 19,
    'PhastCons20Way'                => 20,
    'PhyloP20Way'                   => 21,
    'AC/AN_1000_Genomes'            => 22,
    '1000_Genomes_Homoplasmy'       => 23,
    '1000_Genomes_Heteroplasmy'     => 24,
    'REF'                           => 25,
    'ALT'                           => 26,
    'GT'                            => 27,
    'DP'                            => 28,
    'HF'                            => 29
);

# Defining the keys to be reported in an array
my @keys2report =
  qw ( Variant_Allele   Sample  Locus   Nt_Variability  Codon_Position  Aa_Change       Aa_Variability  REF ALT GT DP HF tRNA_Annotation Disease_Score   RNA_predictions Mitomap_Associated_Disease(s)   Mitomap_Homoplasmy      Mitomap_Heteroplasmy    Somatic_Mutations       SM_Homoplasmy   SM_Heteroplasmy ClinVar OMIM_link       dbSNP_ID        Mamit-tRNA_link PhastCons20Way  PhyloP20Way     AC/AN_1000_Genomes      1000_Genomes_Homoplasmy 1000_Genomes_Heteroplasmy );

# The information about the variants will be stored in a Hash named %hash_out
my %hash_out = ();

###############################################
#  START READING AND PARSING MTOOLBOXR FILE   #
###############################################
open my $fh, '<', $filein;

# Ready-Go!
while ( defined( my $line = <$fh> ) ) {
    next if $. == 1;    # Skip header
    chomp $line;
    my @mtb_fields = split /\t/, $line;

    $mtb_fields[1] = sort_sample_alphabetically( $mtb_fields[1] );
    my $index = $mtb_fields[1] . '_' . $mtb_fields[2] . '_' . $mtb_fields[0];

    ###########################################
    #  The section below includes our filters #
    ###########################################

    # Discarding synonymous variants
    next if $mtb_fields[ $mtb_data_loc{'Aa_Change'} ] eq 'syn';

    # Discarding variants with AC/AN_1000_Genomes > $maf_cutoff (Note that we're keeping void values 'AC/AN_1000_Genomes' => '')
    next
      if ( $maf_cutoff
        && $mtb_fields[ $mtb_data_loc{'AC/AN_1000_Genomes'} ]
        && $mtb_fields[ $mtb_data_loc{'AC/AN_1000_Genomes'} ] >= $maf_cutoff );

    # Discarding variants that don't met hf_cutoff (if arg exists > 0)
    if ( $hf_cutoff > 0 ) {
        my $max_HF = 0;
        my $tmp_DP = $mtb_fields[ $mtb_data_loc{'HF'} ];

        # single mode => 0.02,0.08
        # cohort mode => 01P:0.02,0.8|02M:N/A|03F:1.0
        my @tmp_fields_sample = split /\|/, $tmp_DP;

        # Since we only check Probands (01P) we will apply the same operation to single/cohort modes
        $tmp_fields_sample[0] =~ s/01P\://;    # Get rid of 01P
        $tmp_fields_sample[0] =~ s#N/A#0#g;
        my @tmp_fields_bis = split /,/, $tmp_fields_sample[0];
        my @sorted_tmp_fields = sort { $b cmp $a } @tmp_fields_bis;

        # Load the highesy value
        $max_HF = $sorted_tmp_fields[0];

        # discard non-passing entries
        next unless $max_HF > $hf_cutoff;
    }

    # Those variants who passed all filters are loaded
    foreach my $key (@keys2report) {
        $hash_out{$index}{$key} = $mtb_fields[ $mtb_data_loc{$key} ];
    }

}

#############################################
#  END READING AND PARSING MTOLBOOX FILE  #
##############################################

#############################################
#  PRINTING ACCORDING TO USER PARAMETERS    #
#############################################

if ( $format eq 'hash' ) {
    $Data::Dumper::Sortkeys = 1;
    print Dumper \%hash_out;
}

if ( $format eq 'json' ) {
    my $json = encode_json \%hash_out;
    print "$json\n";
}

if ( $format eq 'json4html' ) {
    my $array_ref = hash2array( \%hash_out );
    my $json      = encode_json $array_ref;
    $json =~ s/^/{"data":/;
    $json =~ s/$/\}/;
    print "$json\n";
}

if ( $format eq 'tsv' ) {

#  Text::CSV::Slurp only works with @AoH (not HoH)/ Note that we miss the 1D of the HoH.
    my @array4tsv = map { $hash_out{$_} } keys %hash_out; #NB: keys (columns) will end up sorted naturally but rows will not
    my $tsv = 'None';
    if (scalar @array4tsv) {  # AoH needs to have elemnts otherwise error
        $tsv = Text::CSV::Slurp->create( input => \@array4tsv, sep_char => "\t" );
    }
    print "$tsv\n";
}

##################################################
#  END OF PRINTING ACCORDING TO USER PARAMETERS  #
##################################################

sub hash2array {

    # dataTables it's incredibly nitpicky
    my $hash_ref = shift;
    my @array    = ();      # Bidimensional array
    my $row      = 0;
    for my $key1 ( keys %{$hash_ref} ) {

        my @keys4html =
          qw ( Sample Locus Variant_Allele    REF ALT   Aa_Change   GT DP HF  tRNA_Annotation   Disease_Score   RNA_predictions  Mitomap_Associated_Disease(s)   Mitomap_Homoplasmy      Mitomap_Heteroplasmy    ClinVar OMIM_link       dbSNP_ID        Mamit-tRNA_link   AC/AN_1000_Genomes   1000_Genomes_Homoplasmy 1000_Genomes_Heteroplasmy );

        my $tmp_var = '';
        for my $key2 (@keys4html) {

            $tmp_var = $hash_ref->{$key1}{$key2};

            # Replacing values with URLs
            $tmp_var =~ s#$tmp_var#<a target="_blank" href="https://ghr.nlm.nih.gov/gene/$tmp_var\#conditions">$tmp_var</a>#
              if $key2 eq 'Locus';
            $tmp_var =~ s#$tmp_var#<a target="_blank" href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=$tmp_var">$tmp_var</a>#
              if $key2 eq 'dbSNP_ID';    # Won't work with multiple rs
            $tmp_var =~ s#$tmp_var#<a target="_blank" href="$tmp_var">$tmp_var</a>#
              if $key2 eq 'OMIM_link';    # Won't work with multiple rs

            # Trimming name of samples
            $tmp_var =~ s#-DNA_MIT##g if $key2 eq 'Sample';

            # For non-empty strings
            if ($tmp_var) {

                # Adding color to suspicius cells providing they have contents
                if ( $key2 eq 'Mitomap_Associated_Disease(s)' ) {
                    $tmp_var =~
                      s/\+/_plus_/g;  #Na+/K+/Ca++  =~ s/// won't work otherwise
                    $tmp_var =~
                      s#$tmp_var#<p class="alert-danger">$tmp_var</p>#;
                }
                $tmp_var =~ s#$tmp_var#<p class="alert-warning">$tmp_var</p>#
                  if ( $key2 eq 'Disease_Score(s)' && $tmp_var >= 0.8 );

                # Changing sample delimiter to ",<br>"
                $tmp_var =~ s#\|#,<br />#g if $key2 eq 'GT';
                $tmp_var =~ s#\|#,<br />#g if $key2 eq 'DP';
                $tmp_var =~ s#\|#,<br />#g if $key2 eq 'HF';
            }

            # Fill
            push @{ $array[$row] }, $tmp_var;
        }
        $row++;
    }
    return \@array;
}

sub sort_sample_alphabetically {
    my $str_in = shift;
    my @samples = split /,/, $str_in;
    @samples = sort (@samples);
    my $str_out = join( ',', @samples );
    return $str_out;
}

=head1 NAME

mtb2json: A script for parsing MToolBox prioritized variants output

=head1 SYNOPSIS

mtb2json.pl -i mit_prioritized_variants.txt [-options]

     Arguments:                       
       -i|input                       MToolBox file
       -f|format                      Output format [>hash|json|json4html|tsv]
       -denovo                        Denovo Mutations
       -hf                            Heteroplasmic Fraction cutoff for 01P [0.0]
       -maf                           maf filter [0.01] (we keep the variant if 1000G MAF < maf value)

     Options:
       -h|help                        Brief help message
       -man                           Full documentation
       -debug                         Print debugging (from 1 to 5, being 5 max)
