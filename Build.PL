use 5.006;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'cworld::dekker',
    license             => 'perl',
    dist_author         => q{Bryan R. Lajoie <bryan.lajoie@gmail.com>},
    dist_version_from   => 'lib/cworld/dekker.pm',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => 0,
    },
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'GD'              => 0,
    },
    add_to_cleanup     => [ 'cworld-dekker-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
