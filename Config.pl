#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component       = 'PC';
our $Code            = 'IPIC3D';
our $MakefileDefOrig = 'Makefile.def.ipic3d';
our @Arguments       = @ARGV;

my $config     = "share/Scripts/Config.pl";

my $GITCLONE = "git clone"; my $GITDIR = "herot:/GIT/FRAMEWORK/";
if (not -f $config and not -f "../../$config"){
    `$GITCLONE $GITDIR/share; $GITCLONE $GITDIR/util`;
}

if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;
our $NewGridSize;
our $ShowGridSize;
our $Hdf5;

my $compiler;
my $debug;


# my $makefileamrex = "share/amrex/GNUmakefile";

# foreach (@Arguments){
#     if(/^-s$/)                 {$Show=1;  next};
#     if(/^-compiler=(.*)/i)     {$compiler="$1";  next};
#     if(/^-debug=(.*)/i)        {$debug=uc("$1");  next};
#     warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
# }


&set_options;

&print_help if $Help;
&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{
     # open(FILE, $makefileamrex) or die "$ERROR could not open $makefileamrex\n";
     #  while(<FILE>){
     #  	next if /^\s*!/; # skip commented out lines
     #  	$compiler = $1 if/^COMP = (\S*)/i;
     # 	$debug = $1 if/^DEBUG = (\S*)/i; 
     #  }
     #  close $makefileamrex;
}

################################################################################
sub set_options{

    # die "$ERROR File $makefileamrex does not exist!\n" unless -f $makefileamrex;

    #  @ARGV = ($makefileamrex);
    # while(<>){
    # 	if($compiler  ne "")     {s/^COMP =.*/COMP = $compiler/;};
    # 	if($debug     ne "")     {s/^DEBUG =.*/DEBUG = $debug/;};
    # 	print;
    # }    
}

################################################################################
sub show_settings{
    &get_settings;
    
    # print "Compiler = $compiler\n";
    # print "Debug    = $debug\n";
}
################################################################################

sub print_help{

    print "
N/A

";
    exit -0;
}
