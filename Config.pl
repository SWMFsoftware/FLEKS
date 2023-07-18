#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

# Add local directory to search
push @INC, ".";

use strict;

our $Component       = 'PC';
our $Code            = 'FLEKS';
our $MakefileDefOrig = 'Makefile.def.FLEKS';
our $MakefileConf    = 'Makefile.conf';
our @Arguments       = @ARGV;

my $config     = "share/Scripts/Config.pl";

my $GITCLONE = "git clone"; my $GITDIR = "git\@github.com:SWMFsoftware";

if (not -f $config and not -f "../../$config"){
    # Stand-alone
    `$GITCLONE $GITDIR/share; $GITCLONE $GITDIR/util`;    
}

my $AmrexDir = "util/AMREX";
if(-f $config){
    #Stand-alone FLEKS. Turn on amrex automatically. 
    push @Arguments, "-amrex";
    require $config;
    die "Error: AMReX doest not exist!\n" unless -d $AmrexDir;
}else{
    require "../../$config";
    $AmrexDir = "../../util/AMREX"; 
}


# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;

my $TPSave = "P";
my $NewTPSave;
my %nTPString=(7=>'P', 10=>'PB', 13=>'PBE');
my %nTPSave=('P' => 7, 'PB' => 10, 'PBE' => 13);
my %TPInfo=('P' => "Particle", 'PB' => "Particle+B", 'PBE' => "Particle+B+E");

my $AmrexDim;

my $nLevMax=1;
my $NewLevMax;

foreach (@Arguments){
     if(/^-s$/)                 {$Show=1;       next};
     if(/^-h$/)                 {$Help=1;       next};     
     if(/^-tp=(.*)$/)           {$NewTPSave=$1; next};
     if(/^-lev=(.*)$/)          {$NewLevMax=$1; next};           
     warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
 }

die "$ERROR -tp input should be 'P', 'PB', or 'PBE'.\n" 
    if $NewTPSave and not exists $nTPSave{$NewTPSave};

my $AmrexComp;
my $AmrexDebug;
my $AmrexTinyProfile; 

&set_test_particle if $NewTPSave;

&set_grid if $NewLevMax;

&print_help if $Help;

&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{
    my $AmrexMakefile = "${AmrexDir}/GNUmakefile";
     open(FILE, $AmrexMakefile) or die "$ERROR could not open $AmrexMakefile \n";
       while(<FILE>){
       	next if /^\s*!/; # skip commented out lines
       	$AmrexComp = $1 if/^COMP = (\S*)/i;
      	$AmrexDebug = $1 if/^DEBUG = (\S*)/i; 
        $AmrexTinyProfile = $1 if/^TINY_PROFILE = (\S*)/i; 
        $AmrexDim = $1 if/^DIM = (\S*)/i; 
       }
    close $AmrexMakefile;
    
    my $NameConstFile = "include/Constants.h";
    `make $NameConstFile`;
    my $nSize;
    open(FILE, $NameConstFile) or die "$ERROR could not open $NameConstFile\n";   
    while(<FILE>){                                                                      
        $nSize=$2           if /\b(ptRecordSize\s*=\s*)(\d+)/i;
	    $nLevMax=$2        if /\b(nLevMax\s*=\s*)(\d+)/i;
    }            
    $TPSave = $nTPString{$nSize};
    close FILE;

    
}

################################################################################
sub set_grid{
    $nLevMax = $NewLevMax if $NewLevMax;

    my $NameConstFile = "include/Constants.h";
    @ARGV = ($NameConstFile);
    while(<>){
        s/\b(nLevMax\s*=[^0-9]*)(\d+)/$1$nLevMax/i;
        print;
    } 
}

################################################################################
sub set_test_particle{
    $TPSave = $NewTPSave if $NewTPSave;

    my $NameConstFile = "include/Constants.h";
    my $nSize = $nTPSave{$TPSave};     
    @ARGV = ($NameConstFile);
    while(<>){
        s/\b(ptRecordSize\s*=[^0-9]*)(\d+)/$1$nSize/i;
        print;
    } 

}

################################################################################
sub show_settings{
    &get_settings;    

    print "AMReX compiler     = $AmrexComp \n";
    print "AMReX debug        = $AmrexDebug \n";
    print "AMReX nDim         = $AmrexDim \n";
    print "AMReX nLevMax      = $nLevMax \n";    
    print "AMReX tiny profile = $AmrexTinyProfile \n";    
    print "Test Particle info = $TPInfo{$TPSave} \n";
}
################################################################################

sub print_help{

    print "
-s            Show the configuration for AMReX.

-lev          Number of maximum grid levels. It is a uniform grid without AMR if lev=1.	

-tp=P,PB,PBE  Test particle output information. 
              P: only save particle velocity + location
              PB: particle + magnetic field. 
              PBE: particle + magnetic field + electric field

\n";
    exit -0;
}
