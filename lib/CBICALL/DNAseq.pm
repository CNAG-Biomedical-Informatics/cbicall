package DNAseq;

use strict;
use warnings;
use Carp;
use File::Spec::Functions qw(catfile);
use feature               qw(say);

#use Data::Dumper;

sub new {
    my ( $class, $self ) = @_;
    bless $self, $class;
    return $self;
}

sub variant_calling {
    my $self = shift;
    my ( $pipeline, $mode, $dir, $threads, $engine, $gatk_version, $id, $debug )
      = @{$self}{
        qw/pipeline mode projectdir threads workflow_engine gatk_version id debug cleanup_bam/
      };

    # Build the script name and base command
    my $suffix = "${pipeline}_$mode";
    my ( $script, $cmd_base );
    if ( $engine eq 'bash' ) {
        $script = $self->{"bash_$suffix"}
          // die "Missing bash script for mode ‘$mode’";
        $cmd_base = "$script -t $threads";
    }
    elsif ( $engine eq 'snakemake' ) {
        $script = $self->{"smk_$suffix"}
          // die "Missing Snakefile for mode ‘$mode’";
        $cmd_base = "snakemake --forceall all -s $script --cores $threads";
    }
    else {
        die "Invalid workflow_engine: $engine";
    }

    # Add the pipeline/config flag only for non-gatk-3.5
    if ( $gatk_version ne 'gatk-3.5' ) {
        if ( $engine eq 'bash' ) {
            $cmd_base .= " --pipeline $pipeline";

            # Add cleanup_bam if set to true
            $cmd_base .= " --cleanup-bam" if $self->{cleanup_bam};
        }
        else {
            $cmd_base .= " --config pipeline=$pipeline";
        }
    }
    say "$cmd_base" if $self->{debug};

    my $log = "${engine}_${suffix}.log";
    my $cmd = "cd $dir && $cmd_base > $log 2>&1";
    submit_cmd( $cmd, $dir, $log, $id, $debug );
    return 1;
}

# Helper to run the command and handle errors
sub submit_cmd {
    my ( $cmd, $dir, $log, $id, $debug ) = @_;
    my $path_log = catfile( $dir, $log );
    my $msg = "Failed to execute: $id\nPlease check this file:\n$path_log\n";
    {
        local $SIG{__DIE__} = 'DEFAULT';
        system($cmd) == 0 or ( $debug ? confess($msg) : die($msg) );
    }
    return 1;
}

1;

