use strict;
use warnings;
use utf8;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;
use File::Find;

install_share dist => 'share';

my $GITHUB_URL = 'https://github.com/CNAG-Biomedical-Informatics/cbicall';

WriteMakefile(
    NAME          => 'CBICall',
    ABSTRACT_FROM => 'lib/CBICALL/Help.pm',
    AUTHOR        => q{Manuel Rueda <mrueda@cpan.org>},
    LICENSE       => 'perl',
    VERSION_FROM  => 'bin/cbicall',
    META_MERGE    => {
        'dynamic_config' => 0,
        'meta-spec' => { version   => 0 },
        'no_index'  => { directory => [qw(docs t nb)] },
        'prereqs'   => { runtime   => { requires => { perl => '5.036000' } } },
        'resources' => {
            bugtracker => { web => "$GITHUB_URL/issues" },
            homepage   => $GITHUB_URL,
            license    =>
              ['http://www.opensource.org/licenses/artistic-license-2.0'],
            repository =>
              { type => 'git', url => "$GITHUB_URL.git", web => $GITHUB_URL },
        },
        'x_contributors' => [],
    },
    PREREQ_PM => {
        'JSON::XS'        => '4.03',
        'YAML::XS'        => '0.88',
        'Path::Tiny'      => '0.144',
        'Term::ANSIColor' => '5.01',
        'Types::Standard' => '2.004000',
        'File::ShareDir::ProjectDistDir' => '1.000009'

    },
    EXE_FILES => [
        'bin/cbicall'

    ],
    TEST_REQUIRES => {
        'Test::More'            => '1.30',
        'File::Path'            => '2.18',
        'File::Spec::Functions' => '3.75'
    },
    test               => { TESTS => 't/*.t' },
    CONFIGURE_REQUIRES => {
        'File::ShareDir::Install' => '0.14',
        'ExtUtils::MakeMaker'     => '7.70'
    }
);

package MY;
use File::ShareDir::Install 'postamble';
