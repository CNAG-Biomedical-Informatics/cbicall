#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use File::Spec::Functions qw(catdir catfile);
use File::Temp qw(tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use CBICALL::DNAseq;  # Module under test

# Create a temporary directory for the dummy project.
my $tempdir = tempdir( CLEANUP => 1 );

### Test variant_calling for bash workflow_engine

# Create a dummy bash script.
my $dummy_script = catfile($tempdir, 'dummy.sh');
open my $fh, '>', $dummy_script or die $!;
print $fh "#!/bin/bash\nexit 0\n";
close $fh;
chmod 0755, $dummy_script;

# Build a dummy object for the bash mode.
my $obj_bash = {
    pipeline      => 'wes',
    workflow_engine => 'bash',
    mode          => 'single',
    projectdir    => $tempdir,
    threads       => 4,
    gatk_version  => 'gatk-4.6',
    id            => 'testid',
    debug         => 0,
    'bash_wes_single' => $dummy_script,
};

# Override submit_cmd to capture the command
my $cmd_bash;
{
    no warnings 'redefine';
    *DNAseq::submit_cmd = sub {
        my ($cmd, $dir, $log, $id, $debug) = @_;
        $cmd_bash = $cmd;
        return 1;
    };
}

my $dnaseq_bash = DNAseq->new($obj_bash);
$dnaseq_bash->variant_calling();
like( $cmd_bash, qr/cd .* && .*bash_wes_single/, 'Bash command is constructed properly' );

### Test variant_calling for snakemake workflow_engine

# Create a dummy snakemake file.
my $dummy_smk = catfile($tempdir, 'smk_wes_single.smk');
open $fh, '>', $dummy_smk or die $!;
print $fh "# dummy snakemake file\n";
close $fh;
chmod 0755, $dummy_smk;

my $cmd_smk;
my $obj_smk = {
    pipeline      => 'wes',
    workflow_engine => 'snakemake',
    mode          => 'single',
    projectdir    => $tempdir,
    threads       => 2,
    gatk_version  => 'gatk-4.6',
    id            => 'testid2',
    debug         => 0,
    smk_wes_single => $dummy_smk,
};

{
    no warnings 'redefine';
    *DNAseq::submit_cmd = sub {
        my ($cmd, $dir, $log, $id, $debug) = @_;
        $cmd_smk = $cmd;
        return 1;
    };
}

my $dnaseq_smk = DNAseq->new($obj_smk);
$dnaseq_smk->variant_calling();
like( $cmd_smk, qr/sn.* -s .*smk_wes_single/, 'Snakemake command is constructed properly' );

### Test that an invalid workflow_engine throws an error

my $obj_invalid = {
    pipeline      => 'wes',
    workflow_engine => 'invalid_mode',
    mode          => 'single',
    projectdir    => $tempdir,
    threads       => 2,
    gatk_version  => 'gatk-4.6',
    id            => 'testid3',
    debug         => 0,
};
my $dnaseq_invalid = DNAseq->new($obj_invalid);
eval {
    $dnaseq_invalid->variant_calling();
};
like( $@, qr/Invalid workflow_engine/, 'Invalid workflow_engine dies as expected' );

done_testing();
