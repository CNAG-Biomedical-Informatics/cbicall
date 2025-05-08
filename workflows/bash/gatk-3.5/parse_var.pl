#!/usr/bin/env perl
#
#   Script for parsing MToolbox prioritized_variants.txt to match VCF_file.tmp
#
#   Last Modified; March/05/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

use strict;
use warnings;

while (<>) {
    chomp;
    s/[A-Z]//g;
    s/\.//g;

    # Possibilites I've seen
    #    310.C => 310
    #    286-287d => 285
    #    8270-8278d => 8269
    #    Insertions keep the same numbering
    my $del = $_ =~ /d/ ? '1' : '0';
    s/[a-z]//g;    # getting rid of d
    my @tmp_fields = split /-/;
    my $var        = $tmp_fields[0] - $del;
    print "$var\n";
}
