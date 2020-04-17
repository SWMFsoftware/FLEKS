#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component       = 'PC';
our $Code            = 'FLEKS';
our $MakefileDefOrig = 'Makefile.def.FLEKS';
our $MakefileConf    = 'Makefile.conf';
our @Arguments       = @ARGV;

my $config     = "share/Scripts/Config.pl";

my $GITCLONE = "git clone"; my $GITDIR = "git\@gitlab.umich.edu:swmf_software";
if (not -f $config and not -f "../../$config"){
    `$GITCLONE $GITDIR/share; $GITCLONE $GITDIR/util`;    
}

if(-f $config){
    require $config;
}else{
    require "../../$config";
}

my $AmrexDir = "util/AMREX";
if(not -d $AmrexDir){
    $AmrexDir = "../../util/AMREX"; 
    if(not -d $AmrexDir) {
        die "Error: AMReX has not been installed!";
        }
}

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;

 foreach (@Arguments){
     if(/^-s$/)                 {$Show=1;  next};
     if(/^-h$/)                 {$Help=1;  next};
#     if(/^-compiler=(.*)/i)     {$compiler="$1";  next};
#     if(/^-debug=(.*)/i)        {$debug=uc("$1");  next};
#     warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
 }

my $AmrexComp;
my $AmrexDebug;
my $AmrexTinyProfile; 

&set_options;

&print_help if $Help;
&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{
    print "get_settings\n";
    my $AmrexMakefile = "${AmrexDir}/GNUmakefile";
    print $AmrexMakefile;
     open(FILE, $AmrexMakefile) or die "$ERROR could not open $AmrexMakefile \n";
       while(<FILE>){
       	next if /^\s*!/; # skip commented out lines
       	$AmrexComp = $1 if/^COMP = (\S*)/i;
      	$AmrexDebug = $1 if/^DEBUG = (\S*)/i; 
        $AmrexTinyProfile = $1 if/^TINY_PROFILE = (\S*)/i; 
       }
    close $AmrexMakefile;
}

################################################################################
sub set_options{
    # print "set_options\n"; 
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

    print "AMReX compiler     = $AmrexComp \n";
    print "AMReX debug        = $AmrexDebug \n";
    print "AMReX tiny profile = $AmrexTinyProfile \n";    
}
################################################################################

sub print_help{

    print "
There is not any configuration setting needed for FLEKS so far. 
Show the configuration for AMReX:
    ./Config.pl 
";
    exit -0;
}
