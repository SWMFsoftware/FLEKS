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

my $mypseudorand;
my $myng;
my $myBATSRUS;
my $ng;
my $pseudrand;
my $BATSRUS;

my $makefilelocal = "src/Makefile.local";

#print "IPIC3D argv = @Arguments \n";

foreach (@Arguments){
    if(/^-s$/)                 {$Show=1;  next};
    if(/^-ng=(.*)/i)           {$ng="$1";  next};
    if(/^-pseudrand=(.*)/i)    {$pseudrand="$1";  next};
    if(/^-BATSRUS=(.*)/i)      {$BATSRUS="$1";  next};
    if(/^-sethdf5/)            {next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&get_settings;

# set local Makefile options
&set_cpp_options;

&print_help if $Help;
&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{
    $mypseudorand = "false";

    # Read size of the grid from topology file
    open(FILE, $makefilelocal) or die "$ERROR could not open $makefilelocal\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
	$mypseudorand="true" if/^.*-DPSEUDRAND\b/i;
	$myng=$1 if/^.*-DNG_P=(\d+)/i;
	$myBATSRUS="true" if/^.*-DBATSRUS\b/i;
    }
    close $makefilelocal;
}

################################################################################
sub set_cpp_options{

    die "$ERROR File $makefilelocal does not exist!\n" unless -f $makefilelocal;

    @ARGV = ($makefilelocal);
    while(<>){
	if($ng        ne "")     {s/DNG_P=[0-9*]/DNG_P=$ng/;};
	if($pseudrand ne "")     {s/\-DPSEUDRAND($|\s)//g;}; # remove if -pseudrand=* is set
	if($pseudrand eq "true") {s/FLAGC_LOCAL=/FLAGC_LOCAL=-DPSEUDRAND /;};
	if($BATSRUS   ne "")     {s/-DBATSRUS($|\s)//g;}; # remove if -BATSRUS=* is set
	if($BATSRUS   eq "true") {s/FLAGC_LOCAL=/FLAGC_LOCAL=-DBATSRUS /;};
	# $Hdf5 is inherited from share/Scripts/Config.pl and it is always set to "yes" or "no"
	s/-DHDF5($|\s)//g;
	if($Hdf5       eq "yes")  {s/FLAGC_LOCAL=/FLAGC_LOCAL=-DHDF5 /};
	print;
    }
}

################################################################################
sub show_settings{
    print "Number of ghost cell layers for the particles: $myng\n";
    print "Use pseudo-random numbers instead of rand()  : $mypseudorand\n";
    print "Couple with BATSRUS                          : $myBATSRUS\n";
}
################################################################################

sub print_help{

    print "
Additional options for IPIC3D/Config.pl:

-ng=NG               Use NG ghost cell layers for particles. Default is 3.

-pseudrand=BOOLEAN   Use pseudo-random numbers or the built-in random numbers.
                     BOOLEAN value is true or false. Default is true.

-BATSRUS=BOOLEAN     Couple iPIC3D with BATSRUS or not. 
                     BOOLEAN value is true or false. Default is true.


-s                   Show IPIC3D settings

Examples:

Set number of ghostcells to 1 and switch on pseudo-random numbers:

   Config.pl -ng=1 -pseudorand=true

Show settings:

   Config.pl -s

";
    exit -0;
}
