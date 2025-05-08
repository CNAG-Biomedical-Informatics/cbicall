#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 10;
use File::Temp qw(tempfile tempdir);
use File::Path qw(make_path);
use File::Spec::Functions qw(catfile catdir);
use Cwd qw(abs_path);
use FindBin;
use lib "$FindBin::Bin/../lib";
use CBICALL::Config;  # Module under test

### Test read_param_file with valid YAML content

# Create a temporary YAML file with valid parameters.
my ($fh, $yaml_file) = tempfile( SUFFIX => '.yaml' );
print $fh <<"EOF";
mode: cohort
pipeline: wes
organism: "Mus musculus"
technology: NovaSeq
workflow_engine: snakemake
gatk_version: gatk-3.5
sample: /tmp/some_sample_dir
EOF
close $fh;

# Read the file.
my $params = Config::read_param_file($yaml_file);
is( $params->{mode},          'cohort',         'Mode should be cohort' );
is( $params->{pipeline},      'wes',            'Pipeline should be wgs' );
is( $params->{organism},      'Mus musculus',   'Organism should be Mus musculus' );
is( $params->{technology},    'NovaSeq',        'Technology should be NovaSeq' );
is( $params->{workflow_engine}, 'snakemake',      'Workflow_mode should be snakemake' );
is( $params->{gatk_version}, 'gatk-3.5',      'Workflow_mode should be gatk4.2' );

### Test that an invalid parameter key causes an error

($fh, my $yaml_invalid) = tempfile( SUFFIX => '.yaml' );
print $fh <<"EOF";
invalid_key: some_value
EOF
close $fh;
{
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    eval { Config::read_param_file($yaml_invalid) };
    like( $@, qr/Parameter 'invalid_key' does not exist/, 'Invalid parameter key should die' );
}

### Test set_config_values with dummy workflows directories

# Create a temporary directory structure.
my $tempdir = tempdir( CLEANUP => 1 );
# Create a "bin" directory and a parallel "workflows" directory.
my $bin_dir = catdir($tempdir, 'bin');
mkdir $bin_dir;
my $workflows_dir  = catdir($tempdir, 'workflows');
my $bash_dir       = catdir($workflows_dir, 'bash',  $params->{gatk_version});
my $snakemake_dir  = catdir($workflows_dir, 'snakemake',  $params->{gatk_version});
mkdir $workflows_dir;
make_path($bash_dir);
make_path($snakemake_dir);

# Create dummy bash files in workflows/bash and mark them executable.
my @bash_files = qw(parameters.sh wes_single.sh wes_cohort.sh mit_single.sh mit_cohort.sh coverage.sh jaccard.sh vcf2sex.sh);
foreach my $file (@bash_files) {
    my $path = catfile($bash_dir, $file);
    open my $fh, '>', $path or die "Cannot create $path: $!";
    print $fh "#!/bin/bash\n";
    close $fh;
    chmod 0755, $path;
}

# Create dummy snakemake files in workflows/snakemake.
my @smk_files = qw(wes_single.smk wes_cohort.smk mit_single.smk mit_cohort.smk config.yaml);
foreach my $file (@smk_files) {
    my $path = catfile($snakemake_dir, $file);
    open my $fh, '>', $path or die "Cannot create $path: $!";
    print $fh "# dummy file\n";
    close $fh;
}

# Set $main::Bin to the bin directory so that set_config_values finds the workflows.
{
    no strict 'refs';
    ${"main::Bin"} = $bin_dir;
}

# Call set_config_values using our $params (which already contains a sample).
my $config = Config::set_config_values($params);
ok( defined $config->{projectdir}, 'projectdir is defined in config' );
ok( -x $config->{bash_parameters}, 'bash_parameters file is executable' );
ok( exists $config->{smk_config},    'smk_config key exists in config' );

done_testing();
