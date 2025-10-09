package Config;

use strict;
use warnings;
use autodie;
use Cwd qw(abs_path);
use Sys::Hostname;
use File::Spec::Functions qw(catdir catfile);
use YAML::XS              qw(LoadFile);
use List::Util            qw(all);
use Types::Standard       qw(Str Int);
use Type::Utils           qw(enum);

#use Data::Dumper;

# Define custom enum types for allowed parameter values
my $mode_type            = enum [qw(single cohort)];
my $pipeline_type        = enum [qw(wes wgs mit)];
my $organism_type        = enum [ 'Homo Sapiens',   'Mus musculus' ];
my $technology_type      = enum [ 'Illumina HiSeq', 'NovaSeq' ];
my $workflow_engine_type = enum [qw(bash nextflow snakemake)];
my $gatk_version         = enum [qw(gatk-3.5 gatk-4.6)];

# Map each parameter key to its type constraint
my %param_types = (
    mode            => $mode_type,
    pipeline        => $pipeline_type,
    organism        => $organism_type,
    technology      => $technology_type,
    workflow_engine => $workflow_engine_type,
    gatk_version    => $gatk_version,
);

# Define default values
my %default = (
    mode            => 'single',
    sample          => undef,
    sample_map      => undef,
    output_basename => undef,
    pipeline        => 'wes',
    organism        => 'Homo Sapiens',
    technology      => 'Illumina HiSeq',
    workflow_engine => 'bash',
    gatk_version    => 'gatk-3.5',
    projectdir      => 'cbicall',
    cleanup_bam     => 0
);

# Define allowed pipeline-mode combinations per GATK version
my %allowed_combos = (
    'gatk-3.5' => {
        wes => [qw(single cohort)],
        mit => [qw(single cohort)],
    },
    'gatk-4.6' => {
        wes => [qw(single cohort)],
        wgs => [qw(single cohort)],
    },
);

sub read_param_file {
    my $yaml_file = shift;

    # Keeping booleans as 'true' or 'false'. Perl still handles 0 and 1 internally.
    $YAML::XS::Boolean = 'JSON::PP';

    my $param = LoadFile($yaml_file);

    # Merge provided parameters with defaults, and validate allowed values
    foreach my $key ( keys %$param ) {
        if ( exists $default{$key} ) {
            my $value = $param->{$key};
            if ( exists $param_types{$key} && defined $value ) {
                $param_types{$key}->assert_valid($value);
            }
            $default{$key} = $value;
        }
        else {
            die "Parameter '$key' does not exist (typo?)\n";
        }
    }

    # Add full path to 'sample_map'
    if ( defined $param->{sample_map} ) {
        $default{sample_map} = abs_path( $param->{sample_map} );
    }

    # Validate pipeline-mode combination for the selected GATK version
    {
        my $version  = $default{gatk_version};
        my $pipeline = $default{pipeline};
        my $mode     = $default{mode};
        unless ( exists $allowed_combos{$version}{$pipeline}
            && grep { $_ eq $mode } @{ $allowed_combos{$version}{$pipeline} } )
        {
            die
"Pipeline-mode '${pipeline}_${mode}' is not supported for GATK version ${version}\n";
        }
    }

    return wantarray ? %default : \%default;
}

sub set_config_values {
    my $param = shift;
    my $user  = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);

    # Base directories for workflows
    my $workflows_bash_dir =
      abs_path( catdir( $main::Bin, '..', 'workflows', 'bash' ) );
    my $workflows_snakemake_dir =
      abs_path( catdir( $main::Bin, '..', 'workflows', 'snakemake' ) );

    # Select the versioned subdirectories
    my $bash_version_dir =
      catdir( $workflows_bash_dir, $param->{gatk_version} );
    my $snakemake_version_dir =
      catdir( $workflows_snakemake_dir, $param->{gatk_version} );

    # Populate config with paths to scripts in the correct versioned directories
    my %config = ( user => $user, );

    # Common bash scripts
    $config{bash_parameters} = catfile( $bash_version_dir, 'parameters.sh' );
    $config{bash_coverage}   = catfile( $bash_version_dir, 'coverage.sh' );
    $config{bash_jaccard}    = catfile( $bash_version_dir, 'jaccard.sh' );
    $config{bash_vcf2sex}    = catfile( $bash_version_dir, 'vcf2sex.sh' );

    # Selected bash workflow (only one per pipeline-mode)
    my $bash_script = $param->{pipeline} . '_' . $param->{mode} . '.sh';
    my $bash_key    = 'bash_' . $param->{pipeline} . '_' . $param->{mode};
    $config{$bash_key} = catfile( $bash_version_dir, $bash_script );

    # Snakemake workflows (only if chosen engine)
    if ( $param->{workflow_engine} eq 'snakemake' ) {
        my $smk_script = $param->{pipeline} . '_' . $param->{mode} . '.smk';
        my $smk_key    = 'smk_' . $param->{pipeline} . '_' . $param->{mode};
        $config{$smk_key} = catfile( $snakemake_version_dir, $smk_script );
        $config{smk_config} = catfile( $snakemake_version_dir, 'config.yaml' );
    }

    # Internal settings
    $config{id}   = time . substr( "00000$$", -5 );
    $config{date} = localtime();
    my $tmp_str = join '_',
      (
        $param->{projectdir}, $param->{workflow_engine},
        $param->{pipeline},   $param->{mode}, $param->{gatk_version},
        $config{id},
      );

    # The directory is ALWAYS created below $param->{sample}
    if ( $param->{sample} ) {
        $config{projectdir} =
          catdir( abs_path( $param->{sample} ), $tmp_str );
        my @parts = split m{/}, $param->{sample};
        $config{output_basename} = $parts[-1];
    }
    else {
        $config{projectdir} = catdir( abs_path($tmp_str) );
    }
    $config{hostname} = hostname;
    chomp( my $threadshost = qx{/usr/bin/nproc} // 1 );
    $config{threadshost} = 0 + $threadshost;
    $config{threadsless} = $threadshost > 1 ? $threadshost - 1 : 1;
    my $th_less = $config{threadsless};
    $config{zip} =
      -x '/usr/bin/pigz'
      ? "/usr/bin/pigz -p $th_less"
      : '/bin/gunzip';

    # Set the genome and capture
    $config{genome} = $param->{gatk_version} eq 'gatk-4.6' ? 'b37' : 'hg19';
    $config{capture} =
      $param->{gatk_version} eq 'gatk-4.6'
      ? 'GATK_bundle_b37'
      : 'Agilent SureSelect';

    # Architecture detection
    chomp( my $uname = `uname -m` );
    my $arch =
        $uname eq 'x86_64'  ? 'x86_64'
      : $uname eq 'aarch64' ? 'arm64'
      :                       $uname;
    $config{arch} = $arch;

    # Validate executable permissions
    my @exe_keys = (
        qw(bash_parameters bash_coverage bash_jaccard bash_vcf2sex),
        'bash_' . $param->{pipeline} . '_' . $param->{mode},
    );

    die "Missing +x on one or more workflow scripts"
      unless all { -x $config{$_} } @exe_keys;

    return wantarray ? %config : \%config;
}

1;
