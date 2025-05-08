#!/usr/bin/env perl
#
#   Script to calculate mean and standard deviation
#
#   Last Modified; Feb/08/2017
#
#   Copyright (C) 2017 Manuel Rueda (mrueda@scripps.edu)
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#   If this program helps you in your research, please cite.

use strict;
use warnings;

# https://www.thoughtco.com/population-vs-sample-standard-deviations-3126372
my ( $sum, $sumsq, $mean, $counter, $desv );

# Per fer _Population_ Standard Desviation (not Sample SD)
while (<>) {
    next unless /[0-9]/;
    die "You have commmas\n" if /\,/;    # perl does not like commas, takes integer
    chomp;
    $sum   += $_;
    $sumsq += $_ * $_;
    $counter++;
}

$mean = $sum / $counter;
$desv = sqrt( ( $sumsq - $sum**2 / $counter ) / $counter );
#printf "%12.5f %8.4f \n ", $mean,$desv;
#printf "%12.5f  ", $mean;
printf "%8.2f %8.2f \n ", $mean, $desv;
#printf "%15.10f %15.10f \n ", $mean,$desv;
