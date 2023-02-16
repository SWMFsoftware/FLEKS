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

my $GITCLONE = "git clone"; my $GITDIR = "git\@gitlab.umich.edu:swmf_software";

if (not -f $config and not -f "../../$config"){
    # Stand-alone
    `$GITCLONE $GITDIR/share; $GITCLONE $GITDIR/util`;    
}

if(-f "../../$config"){
    # Not stand-alone
    `$GITCLONE $GITDIR/srcBATL srcBATL_orig` unless -d "srcBATL_orig";
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
our $NewGridSize;
our $Help;
our $ERROR;
our $WARNING;

my $TPSave = "P";
my $NewTPSave;
my %nTPString=(7=>'P', 10=>'PB', 13=>'PBE');
my %nTPSave=('P' => 7, 'PB' => 10, 'PBE' => 13);
my %TPInfo=('P' => "Particle", 'PB' => "Particle+B", 'PBE' => "Particle+B+E");

# BATL Grid size variables
my $NameGridFile = "srcBATL/BATL_size.f90";
my $GridSize;
my ($nI, $nJ, $nK, $iRatio, $jRatio, $kRatio);
my $GhostCell;
my $NewGhostCell;

foreach (@Arguments){
     if(/^-s$/)                 {$Show=1;       next};
     if(/^-h$/)                 {$Help=1;       next};
     if(/^-ng=(.*)$/)           {$NewGhostCell=$1; next};
     if(/^-tp=(.*)$/)           {$NewTPSave=$1;    next};    
     warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
 }

die "$ERROR -tp input should be 'P', 'PB', or 'PBE'.\n" 
    if $NewTPSave and not exists $nTPSave{$NewTPSave};

my $AmrexComp;
my $AmrexDebug;
my $AmrexTinyProfile; 

&set_options;

&set_test_particle if $NewTPSave;

&print_help if $Help;

&get_batl_settings;

&set_batl_grid_size;

&show_settings if $Show;

exit 0;
#############################################################################

sub get_batl_settings{

    # Read size of the grid from $NameGridFile
    open(MODSIZE,$NameGridFile) or die "$ERROR could not open $NameGridFile\n";
    while(<MODSIZE>){
        next if /^\s*!/; # skip commented out lines
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
        $nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
        $nK=$1           if /\bnK\s*=\s*(\d+)/i;
	$iRatio=$1       if /\biRatio\s*=\s*min\(\s*(\d)/;
	$jRatio=$1       if /\bjRatio\s*=\s*min\(\s*(\d)/;
	$kRatio=$1       if /\bkRatio\s*=\s*min\(\s*(\d)/;
	$GhostCell=$1    if /\bnG\s*=\s*(\d)/;
    }
    close MODSIZE;

    die "$ERROR could not read nI from $NameGridFile\n" unless length($nI);
    die "$ERROR could not read nJ from $NameGridFile\n" unless length($nJ);
    die "$ERROR could not read nK from $NameGridFile\n" unless length($nK);

    die "$ERROR could not read iRatio from $NameGridFile\n" 
	unless length($iRatio);
    die "$ERROR could not read jRatio from $NameGridFile\n" 
	unless length($jRatio);
    die "$ERROR could not read kRatio from $NameGridFile\n" 
	unless length($kRatio);

    die "$ERROR could not read nG from $NameGridFile\n" 
	unless length($GhostCell);

    $GridSize  = "$nI,$nJ,$nK";
}

#############################################################################

sub set_batl_grid_size{

    $GridSize = $NewGridSize if $NewGridSize;

    $GhostCell = $NewGhostCell if $NewGhostCell;
    
    if($GridSize =~ /^[1-9]\d*,[1-9]\d*,[1-9]\d*$/){
	($nI,$nJ,$nK) = split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR ".
	    "-g=$GridSize should be 3 positive integers separated by commas\n";
    }

    @ARGV = ($NameGridFile);

    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nI\s*=[^0-9]*)\d+/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)\d+/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)\d+/$1$nK/i;
	s/\b(nG\s*=[^0-9]*)\d+/$1$GhostCell/i;
	print;
    }

}

sub get_settings{
    my $AmrexMakefile = "${AmrexDir}/GNUmakefile";
     open(FILE, $AmrexMakefile) or die "$ERROR could not open $AmrexMakefile \n";
       while(<FILE>){
       	next if /^\s*!/; # skip commented out lines
       	$AmrexComp = $1 if/^COMP = (\S*)/i;
      	$AmrexDebug = $1 if/^DEBUG = (\S*)/i; 
        $AmrexTinyProfile = $1 if/^TINY_PROFILE = (\S*)/i; 
       }
    close $AmrexMakefile;
    
    my $NameConstFile = "include/Constants.h";
    `make $NameConstFile`;
    my $nSize;
    open(FILE, $NameConstFile) or die "$ERROR could not open $NameConstFile\n";   
    while(<FILE>){                                                                      
        $nSize=$2           if /\b(ptRecordSize\s*=\s*)(\d+)/i;                                
    }            
    $TPSave = $nTPString{$nSize};
    close FILE;

    
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
    print "Test Particle info = $TPInfo{$TPSave} \n";
    print "FLEKS BATL block size       : nI=$nI, nJ=$nJ, nK=$nK \n";
    print "FLEKS BATL ghost cell layers: nG=$GhostCell \n";
}
################################################################################

sub print_help{

    print "
-s            Show the configuration for AMReX.

-tp=P,PB,PBE  Test particle output information. 
              P: only save particle velocity + location
              PB: particle + magnetic field. 
              PBE: particle + magnetic field + electric field

\n";
    exit -0;
}
