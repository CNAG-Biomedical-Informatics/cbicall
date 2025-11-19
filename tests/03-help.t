#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use File::Temp qw(tempfile);
use FindBin;
use lib "$FindBin::Bin/../lib";
use CBICALL::Help;  # Module under test


### Test GoodBye message

# Since the file defines package GoodBye (not Help::GoodBye), call it directly.
my $goodbye_obj = GoodBye->new();
my $bye = $goodbye_obj->say_goodbye();
ok( defined $bye && $bye ne '', 'GoodBye returns a non-empty string' );

### Test usage with valid arguments

# Create a dummy YAML file to simulate a parameters file.
my ($fh, $dummy_yaml) = tempfile( SUFFIX => '.yaml' );
print $fh "dummy: dummy\n";
close $fh;

# Simulate command-line arguments.
local @ARGV = ('-t', '4', '-p', $dummy_yaml);

# usage should return a hash reference with the parsed options.
my $args;
eval { $args = Help::usage('TestVersion') };
ok( $args->{threads} == 4, 'Usage returns correct threads' );
like( $args->{paramfile}, qr/\Q$dummy_yaml\E/, 'Usage returns correct paramfile' );
