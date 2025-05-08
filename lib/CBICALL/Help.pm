package Help;

use strict;
use warnings;
use feature qw(say);
use Pod::Usage;
use Getopt::Long;

sub usage {

    # http://www.gsp.com/cgi-bin/man.cgi?section=3&topic=Getopt::Long
    my $version = shift;
    my %arg     = ();
    GetOptions(
        'v'           => sub { print "$version\n"; exit },
        'debug=i'     => \$arg{debug},                       # numeric (integer)
        'verbose'     => \$arg{verbose},                     # flag
        'h|help'      => \$arg{help},                        # flag
        'man'         => \$arg{man},                         # flag
        't|threads=i' => \$arg{threads},                     # numeric (integer)
        'p|param=s'   => \$arg{paramfile},                   # string
        'no-color|nc' => \$arg{nocolor},
        'no-emoji|ne' => \$arg{noemoji},
    ) or pod2usage( -exitval => 0, -verbose => 1 );

    # Control check
    pod2usage( -exitval => 0, -verbose => 2 ) if $arg{man};
    pod2usage( -exitval => 0, -verbose => 1 ) if $arg{help};
    pod2usage( -exitval => 1, -verbose => 1 )
      if ( !$arg{threads} || !$arg{paramfile} );
    pod2usage(
        -exitval => 1,
        -verbose => 1,
        -message => 'Option --i requires a parameters file'
    ) if ( !-s $arg{paramfile} );
    pod2usage(
        -exitval => 1,
        -verbose => 1,
        -message => 'Option --n requires a positive integer'
    ) if ( $arg{threads} <= 0 );    # Must be positive integer

    # Initialize undefs
    $arg{debug} = 0 if !$arg{debug};

    # Apply global settings (if necessary)
    $ENV{'ANSI_COLORS_DISABLED'} = 1 if $arg{nocolor};

    return wantarray ? %arg : \%arg;
}

package GoodBye;

sub new {
    my $class = shift;
    return bless {}, $class;
}

=head2 goodbye

    About   : Well, the name says it all :-)
    Usage   :         
    Args    : 

=cut

sub say_goodbye {
    my @words = ( <<"EOF" =~ m/^\s*(.+)/gm );
      Aavjo
      Abar Dekha-Hobe
      Adeus
      Adios
      Aloha
      Alvida
      Ambera
      Annyong hi Kashipshio
      Arrivederci
      Auf Wiedersehen
      Au Revoir
      Ba'adan Mibinamet
      Dasvidania
      Donadagohvi
      Do Pobatchenya
      Do Widzenia
      Eyvallah
      Farvel
      Ha Det
      Hamba Kahle
      Hooroo
      Hwyl
      Kan Ga Waanaa
      Khuda Hafiz
      Kwa Heri
      La Revedere
      Le Hitra Ot
      Ma'as Salaam
      Mikonan
      Na-Shledanou
      Ni Sa Moce
      Paalam
      Rhonanai
      Sawatdi
      Sayonara
      Selavu
      Shalom
      Totsiens
      Tot Ziens
      Ukudigada
      Vale
      Zai Geen
      Zai Jian
      Zay Gesunt
EOF
    my $random_word = $words[ rand @words ];
    return $random_word;
}
1;
