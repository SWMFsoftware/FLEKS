#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# Allow in-place editing
$^I = "";

# Add local directory to search
push @INC, ".";

use strict;
use File::Copy qw(copy);

our $Component       = 'PC';
our $Code            = 'FLEKS';
our $MakefileDefOrig = 'Makefile.def.FLEKS';
our $MakefileConf    = 'Makefile.conf';
our @Arguments       = @ARGV;

my $config     = "share/Scripts/Config.pl";
my $ConstantsFile = "include/Constants.h";
my $UserSourceFile = "include/UserSource.h";
my $UserSourceDir = "userfiles";

if (@Arguments == 1 and $Arguments[0] eq "-u") {
    &print_user_source_options;
    exit 0;
}

my $GITDIR = "git\@github.com:SWMFsoftware";

if (not -f $config and not -f "../../$config"){
    # Stand-alone
    if (not -f $config) {
        print "--- Cloning SWMFsoftware/share ---\n";
        if (-d "share") {
            system("git clone $GITDIR/share share_tmp && cp -rp share_tmp/. share/ && rm -rf share_tmp") == 0 or die "Error: could not clone share\n";
        } else {
            system("git clone $GITDIR/share") == 0 or die "Error: could not clone share\n";
        }
    }
    if (not -d "util") {
        print "--- Cloning SWMFsoftware/util ---\n";
        system("git clone $GITDIR/util") == 0 or die "Error: could not clone util\n";
    }
    if (not -d "util/AMREX") {
        print "--- Cloning SWMFsoftware/AMREX ---\n";
        system("git clone $GITDIR/AMREX util/AMREX") == 0 or die "Error: could not clone AMREX\n";
    }
}

my $AmrexDir = (-f $config) ? "util/AMREX" : "../../util/AMREX";
if (-f "Makefile.def" and not grep { /^-install(=s)?$/i } @Arguments) {
    if (open(my $mf, "<", "Makefile.def")) {
        while (<$mf>) {
            if (/^\s*include\s+(.*)\/Makefile\.def\s*$/) {
                $AmrexDir = "$1/util/AMREX" if -d "$1/util/AMREX";
                last;
            }
        }
        close $mf;
    }
}

# GPU (CUDA) and CPU OpenMP builds. -cuda[=ARCH] links FLEKS with a CUDA build
# of AMReX (InstallDir<N>D_cuda, built on demand), -openmp with an OpenMP build
# (InstallDir<N>D_omp, built on demand), and -nocuda / -noopenmp with the
# serial CPU build (InstallDir<N>D). Without these flags the current choice is
# kept: it is the library InstallDir points to.
my $NewCuda;
my $CudaArch;
my $UseCuda = &amrex_uses_cuda("$AmrexDir/InstallDir");
my $NewOmp;
my $UseOmp = &amrex_uses_omp("$AmrexDir/InstallDir");
foreach (@Arguments){
    if(/^-cuda(?:=([\d,]+))?$/i){ $NewCuda = 1; $CudaArch = $1 if defined $1; }
    if(/^-nocuda$/i)            { $NewCuda = 0; }
    if(/^-openmp$/i)            { $NewOmp  = 1; }
    if(/^-noopenmp$/i)          { $NewOmp  = 0; }
}
die "Error: -cuda and -openmp cannot be enabled at the same time.\n"
    if $NewCuda and $NewOmp;
@Arguments = grep { !/^-(no)?(cuda(=\S+)?|openmp)$/i } @Arguments;
if (defined $NewCuda) {
    $UseCuda = $NewCuda;
    $UseOmp = 0 if $UseCuda;
}
if (defined $NewOmp) {
    $UseOmp = $NewOmp;
    $UseCuda = 0 if $UseOmp;
}

my $HasDimArg = scalar grep { /^-amrex(2d|3d)?$/i } @Arguments;
my $AmrexDimArg = "3d";
foreach (@Arguments){
    if (/^-amrex(2d|3d)?$/i){
        $AmrexDimArg = lc($1 || "3d");
        last;
    }
}
unless ($HasDimArg){
    my $AmrexConfig = "${AmrexDir}/InstallDir/include/AMReX_Config.H";
    if (-f $AmrexConfig) {
        open(FILE, $AmrexConfig) or warn "Config.pl WARNING: could not open $AmrexConfig \n";
        while(<FILE>){
            $AmrexDimArg = "$1d" if /^#define\s+AMREX_SPACEDIM\s+(\d+)/i;
        }
        close FILE;
    }
}

if(-f $config){
    # Local FLEKS dependency tree. Turn on AMReX automatically on install, or
    # if the base CPU AMReX build is missing. When only switching dimension on
    # an existing install, link_amrex below updates the InstallDir symlink.
    my $needAmrex = (!-f $MakefileConf or
                     !-e "$AmrexDir/InstallDir\U$AmrexDimArg\E/lib/libamrex.a" or
                     grep { /^-install(=.*)?$/i } @Arguments);
    if ($needAmrex) {
        &set_c_compiler_conf("gcc_mpicc", undef, $UseOmp) if -f $MakefileConf;
        push @Arguments, "-amrex$AmrexDimArg"
            unless grep { /^-amrex(2d|3d)?$/i } @Arguments;
    } else {
        @Arguments = grep { !/^-amrex(2d|3d)?$/i } @Arguments;
        if (defined $NewCuda or defined $NewOmp) {
            &set_c_compiler_conf($UseCuda ? "nvcc" : "gcc_mpicc",
                                 "$AmrexDir/InstallDir\U$AmrexDimArg\E_cuda",
                                 $UseOmp);
        }
    }
    push @Arguments, "-show" unless @ARGV;
    push @Arguments, "-nodebug", "-install" unless -f $MakefileConf;
    require $config;

    # Compile share lib if it's not there
    if (not -f "lib/libSHARE.a") {
        print "Building libSHARE.a...\n";
        mkdir "lib" unless -d "lib";
        system("cd share/Library/src; make LIB") == 0 or die "Error: could not build libSHARE.a\n";
    }

    die "Error: AMReX does not exist!\n" unless -d $AmrexDir;

    if (($needAmrex and ($UseCuda or $UseOmp)) or $HasDimArg or
        defined $NewCuda or defined $NewOmp) {
        &link_amrex($AmrexDimArg, $UseCuda, $UseOmp);
    }
}else{

    my $IsInstall = 0;
    foreach (@Arguments) {
        if (/^-install$/) {
            $_ = "-install=c";
        }
        $IsInstall = 1 if /^-install(=.*)?$/;
    }
    push @Arguments, "-install=c" unless -f $MakefileConf or $IsInstall;
    @Arguments = grep { !/^-amrex(2d|3d)?$/i } @Arguments
        if -e "$AmrexDir/InstallDir\U$AmrexDimArg\E/lib/libamrex.a";
    require "../../$config";
    if ($UseCuda or $UseOmp or $HasDimArg or
        defined $NewCuda or defined $NewOmp) {
        &link_amrex($AmrexDimArg, $UseCuda, $UseOmp);
    }
}

# These are inherited from $config
our %Remaining;   # Arguments not handled by share/Scripts/Config.pl
our $Show;
our $Help;
our $ERROR;
our $WARNING;

my $TPSave = "PBED";
my $NewTPSave;
my %nTPString=(7=>'P', 10=>'PB', 13=>'PBE', 19=>'PBED', 22=>'PBEG');
my %nTPSave=('P' => 7, 'PB' => 10, 'PBE' => 13, 'PBED' => 19, 'PBEG' => 22);
my %TPInfo=('P' => "Particle", 'PB' => "Particle+B", 'PBE' => "Particle+B+E", 'PBED' => "Particle+B+E+Drifts", 'PBEG' => "Particle+B+E+GradB");

my $AmrexDim;

my $nLevMax=1;
my $NewLevMax;
my $NewUserSource;
my $ListUserSources;
my $VerboseGmres;
my $SimdFlags;

foreach (@Arguments){
     if(/^-s$/)                 {$Show=1;       next};
     if(/^-h$/)                 {$Help=1;       next};     
     if(/^-tp=(.*)$/)           {$NewTPSave=$1; next};
     if(/^-lev$/)               {$NewLevMax=1;  next};
     if(/^-lev=(.*)$/)          {$NewLevMax=$1; next};
     if(/^-u$/)                 {$ListUserSources=1; next};
     if(/^-u=(.*)$/)            {$NewUserSource=$1; next};
     if(/^-debug$/i)            {$VerboseGmres="yes"; next};
     if(/^-nodebug$/i)          {$VerboseGmres="no";  next};
     if(/^-simd(?:=(\S+))?$/i)  {$SimdFlags = defined($1) ? $1 : "yes"; next};
     if(/^-nosimd$/i)           {$SimdFlags="no";  next};
     warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

die "$ERROR -tp input should be 'P', 'PB', 'PBE', or 'PBEG'.\n"
    if $NewTPSave and not exists $nTPSave{$NewTPSave};

die "$ERROR -u input should use letters, numbers, or underscores.\n"
    . &get_user_source_options_text
    if defined $NewUserSource and $NewUserSource !~ /^[A-Za-z0-9_]+$/;

my $AmrexComp;
my $AmrexDebug;
my $AmrexTinyProfile; 

if ($ListUserSources) {
    &print_user_source_options;
    exit 0;
}

&get_settings;

# During install, seed a local UserSource.h from the default template only if
# the user did not select one and no local UserSource.h exists yet.
if (!defined $NewUserSource && !-f $UserSourceFile &&
    scalar grep { /^-install(=.*)?$/ } @Arguments) {
    $NewUserSource = "Default";
}

&set_user_source if defined $NewUserSource;

&set_test_particle if $NewTPSave and $NewTPSave ne $TPSave;

&set_grid if $NewLevMax and $NewLevMax ne $nLevMax;

&set_verbose_gmres if defined $VerboseGmres;

&set_simd_flags if defined $SimdFlags;

&print_help if $Help;

&show_settings if $Show;

exit 0;
#############################################################################

sub get_settings{
    my $AmrexMakefile = "${AmrexDir}/GNUmakefile";
    open(FILE, $AmrexMakefile) or warn "$WARNING could not open $AmrexMakefile \n";
    while(<FILE>){
       	next if /^\s*!/; # skip commented out lines
       	$AmrexComp = $1 if/^COMP = (\S*)/i;
    }
    close FILE;

    $AmrexDebug = "False";
    $AmrexTinyProfile = "False";
    my $AmrexConfig = "${AmrexDir}/InstallDir/include/AMReX_Config.H";
    open(FILE, $AmrexConfig) or warn "$WARNING could not open $AmrexConfig \n";
    while(<FILE>){
      	$AmrexDebug = "True" if/^#define AMREX_DEBUG/i; 
        $AmrexTinyProfile = "True" if/^#define AMREX_TINY_PROFILING/i; 
        $AmrexDim = $1 if/^#define AMREX_SPACEDIM (\S*)/i; 
    }
    close FILE;

    system("make", "-s", $ConstantsFile) == 0
        or die "$ERROR could not create $ConstantsFile\n";
    my $nSize;
    open(FILE, $ConstantsFile) or die "$ERROR could not open $ConstantsFile\n";   
    while(<FILE>){                                                                      
        $nSize = $2   if /\b(ptRecordSize\s*=\s*)(\d+)/i;
	$nLevMax = $2 if /\b(nLevMax\s*=\s*)(\d+)/i;
    }
    $TPSave = $nTPString{$nSize};
    close FILE;
}
################################################################################
sub set_grid{
    $nLevMax = $NewLevMax if $NewLevMax;

    @ARGV = ($ConstantsFile);
    while(<>){
        s/\b(nLevMax\s*=[^0-9]*)(\d+)/$1$nLevMax/i;
        print;
    }
}
################################################################################
sub set_verbose_gmres{
    # Set the verbose GMRES residual logging flag in src/Makefile:
    #   -debug   -> VERBOSE_GMRES = -DGMRES_VERBOSE (log residual of every matvec)
    #   -nodebug -> VERBOSE_GMRES is empty (only iteration count + final residual)
    my $VERBOSE_GMRES = "";
    $VERBOSE_GMRES = '-DGMRES_VERBOSE' if $VerboseGmres eq "yes";

    print "Setting verbose GMRES logging to '$VerboseGmres' in src/Makefile\n";

    @ARGV = ("src/Makefile");
    while(<>){
        s/^(\s*VERBOSE_GMRES\s*=).*/$1 $VERBOSE_GMRES/;
        print;
    }
}
################################################################################
sub set_simd_flags{
    # Set SIMD compiler flags in src/Makefile:
    #   -simd        -> SIMD_FLAGS = -march=native -fopenmp-simd (default)
    #   -simd=ARCH   -> SIMD_FLAGS = -march=ARCH -fopenmp-simd
    #   -nosimd      -> SIMD_FLAGS is empty
    my $SIMD_FLAGS = "";
    if ($SimdFlags eq "yes") {
        $SIMD_FLAGS = '-march=native -fopenmp-simd';
    } elsif ($SimdFlags ne "no") {
        $SIMD_FLAGS = "-march=$SimdFlags -fopenmp-simd";
    }

    print "Setting SIMD flags to '$SIMD_FLAGS' in src/Makefile\n";

    @ARGV = ("src/Makefile");
    while(<>){
        s/^(\s*SIMD_FLAGS\s*=).*/$1 $SIMD_FLAGS/;
        print;
    }
}
################################################################################
sub set_test_particle{
    $TPSave = $NewTPSave if $NewTPSave;

    my $nSize = $nTPSave{$TPSave};
    @ARGV = ($ConstantsFile);
    while(<>){
        s/\b(ptRecordSize\s*=[^0-9]*)(\d+)/$1$nSize/i;
        print;
    } 
}
################################################################################
sub set_user_source{
    my $sourceFile = "$UserSourceDir/${NewUserSource}Source.h";

    die "$ERROR could not find $sourceFile\n" . &get_user_source_options_text
        unless -f $sourceFile;

    copy($sourceFile, $UserSourceFile)
        or die "$ERROR could not copy $sourceFile to $UserSourceFile: $!\n";
}
################################################################################
sub print_user_source_options{
    print &get_user_source_options_text;
}
################################################################################
sub get_user_source_options_text{
    my $text = "Available user source options:\n";

    my @sourceFiles = sort glob("$UserSourceDir/*Source.h");
    unless (@sourceFiles) {
        return $text . "  None found in $UserSourceDir\n";
    }

    foreach my $sourceFile (@sourceFiles) {
        my $name = $sourceFile;
        $name =~ s#^\Q$UserSourceDir\E/##;
        $name =~ s/Source\.h$//;
        $text .= "  $name\n";
    }

    return $text;
}
################################################################################
sub show_settings{
    &get_settings;    

    print "AMReX: compiler     = $AmrexComp \n";
    print "AMReX: debug        = $AmrexDebug \n";
    print "AMReX: nDim         = $AmrexDim \n";
    print "AMReX: nLevMax      = $nLevMax \n";
    print "AMReX: tiny profile = $AmrexTinyProfile \n";
    print "AMReX: OpenMP       = " .
        (&amrex_uses_omp("$AmrexDir/InstallDir") ? "True" : "False") . " \n";
    print "AMReX: GPU          = " .
        (&amrex_uses_cuda("$AmrexDir/InstallDir") ? "CUDA" : "none") . " \n";
    print "Test Particle info  = $TPInfo{$TPSave} \n";
    print "User Source         = " . &get_user_source_info . " \n";
}
################################################################################
sub get_user_source_info{
    open(my $file, "<", $UserSourceFile)
        or return "Unavailable ($UserSourceFile not found)";

    local $/;
    my $content = <$file>;
    close $file;

    if($content =~ /\binfo\s*=\s*"((?:\\.|[^"\\])*)"\s*;/s){
        my $info = $1;
        $info =~ s/\\(["\\])/$1/g;
        return $info;
    }

    return "Unavailable (info not set in $UserSourceFile)";
}
################################################################################
sub amrex_uses_cuda{
    # True if the AMReX installation in directory $_[0] was built with CUDA.
    my $config = "$_[0]/include/AMReX_Config.H";
    open(my $file, "<", $config) or return 0;
    my $useCuda = grep { /^#define\s+AMREX_USE_CUDA\s+1/ } <$file>;
    close $file;
    return $useCuda;
}
################################################################################
sub amrex_uses_omp{
    # True if the AMReX installation in directory $_[0] was built with OpenMP.
    my $config = "$_[0]/include/AMReX_Config.H";
    open(my $file, "<", $config) or return 0;
    my $useOmp = grep { /^#define\s+AMREX_USE_OMP\s+1/ } <$file>;
    close $file;
    return $useOmp;
}
################################################################################
sub target_makefile_conf{
    my $conf = $MakefileConf;
    if (-f $conf and open(my $fh, "<", $conf)) {
        while (<$fh>) {
            if (/^\s*include\s+(.*Makefile\.conf)\s*$/) {
                $conf = $1 if -f $1;
                last;
            }
        }
        close $fh;
    }
    return $conf;
}
################################################################################
sub set_c_compiler_conf{
    my ($comp, $amrexInstallDir, $useOmp) = @_;
    my $cMake = (-f "share/build/Makefile.$comp")
        ? "share/build/Makefile.$comp" : "../../share/build/Makefile.$comp";
    my $targetConf = &target_makefile_conf();
    return unless -f $targetConf and -f $cMake;
    my $mpicxxShow = `mpicxx -show 2>/dev/null`;
    my ($cudaArch, @genCode);
    if ($comp eq "nvcc" and $amrexInstallDir and
        open(my $pc, "<", "$amrexInstallDir/lib/pkgconfig/amrex.pc")) {
        while (<$pc>) {
            push @genCode, /(--generate-code\s+\S+)/g;
        }
        close $pc;
        $cudaArch = join(",", map { /arch=compute_(\d+)/ ? $1 : () } @genCode);
    }
    my ($skip, $isDebug) = (0, 0);
    @ARGV = ($targetConf);
    while (<>) {
        $isDebug = 1 if !$skip and /^\s*DEBUG\s*=\s*\$\{DEBUGFLAG\}/;
        if (!$skip and defined $useOmp) {
            s/^#(OPENMPFLAG\s*=)/$1/ if $useOmp;
            s/^(OPENMPFLAG\s*=)/#$1/ if !$useOmp;
        }
        if (/^# C language related part of Makefile\.conf:/) {
            $skip = 1;
            open(my $f, "<", $cMake) or die "$ERROR could not open $cMake\n";
            while (my $line = <$f>) {
                if ($line =~ /^CPPLIB/) {
                    $line =~ s/ -lmpi_cxx// if $mpicxxShow !~ /-lmpi_cxx/;
                    $line =~ s/ -lmpi\+\+// if $mpicxxShow !~ /-lmpi\+\+/;
                }
                $line =~ s/^#(DEBUGC\s*=)/$1/ if $isDebug;
                $line =~ s/^CUDA_ARCH\s*:=.*/CUDA_ARCH := $cudaArch/ if $cudaArch;
                $line =~ s/--generate-code\s+\S+/join(" ", @genCode)/e if @genCode;
                print $line;
            }
            close $f;
        }
        print unless $skip;
    }
}
################################################################################
sub link_amrex{
    # Point $AmrexDir/InstallDir to the serial CPU, OpenMP CPU, or CUDA build
    # of AMReX, building the requested library on demand, and select
    # Makefile.nvcc or Makefile.gcc_mpicc in Makefile.conf.
    my ($dimArg, $useCuda, $useOmp) = @_;
    my $nDim = substr($dimArg, 0, 1);
    my $InstallDir = $useCuda ? "InstallDir${nDim}D_cuda"
                   : $useOmp  ? "InstallDir${nDim}D_omp"
                              : "InstallDir${nDim}D";
    if ($useCuda) {
        my $arch = $CudaArch ||
            `nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -dc 0-9` ||
            "80";
        (my $archList = $arch) =~ s/,/ /g;
        my @libArch = (`cat $AmrexDir/$InstallDir/lib/pkgconfig/amrex.pc 2>/dev/null`
                       =~ /arch=compute_(\d+)/g);
        if (not -e "$AmrexDir/$InstallDir/lib/libamrex.a" or not @libArch or
            (defined $NewCuda and join(" ", @libArch) ne $archList)) {
            my $flags = "USE_CUDA=TRUE USE_OMP=FALSE CUDA_ARCH=\"$archList\" USE_CUDA_FAST_MATH=FALSE " .
                "NVCC_CCBIN=mpicxx DIM=$nDim AMREX_INSTALL_DIR=$InstallDir";
            print "Building the CUDA AMReX library in $AmrexDir/$InstallDir (sm_$arch)\n";
            system("cd $AmrexDir && rm -rf $InstallDir && make clean && make -j 8 $flags && make install $flags") == 0
                or die "$ERROR could not build $AmrexDir/$InstallDir\n";
        }
        &set_c_compiler_conf("nvcc", "$AmrexDir/$InstallDir", 0);
    } elsif ($useOmp) {
        if (not -e "$AmrexDir/$InstallDir/lib/libamrex.a" or
            not &amrex_uses_omp("$AmrexDir/$InstallDir")) {
            my $flags = "USE_OMP=TRUE USE_CUDA=FALSE DIM=$nDim AMREX_INSTALL_DIR=$InstallDir";
            print "Building the OpenMP AMReX library in $AmrexDir/$InstallDir\n";
            system("cd $AmrexDir && rm -rf $InstallDir && make clean && make -j 8 $flags && make install $flags") == 0
                or die "$ERROR could not build $AmrexDir/$InstallDir\n";
        }
        &set_c_compiler_conf("gcc_mpicc", undef, 1);
    } else {
        if (not -e "$AmrexDir/$InstallDir/lib/libamrex.a") {
            my $flags = "USE_OMP=FALSE USE_CUDA=FALSE DIM=$nDim AMREX_INSTALL_DIR=$InstallDir";
            print "Building the CPU AMReX library in $AmrexDir/$InstallDir\n";
            system("cd $AmrexDir && rm -rf $InstallDir && make clean && make -j 8 $flags && make install $flags") == 0
                or die "$ERROR could not build $AmrexDir/$InstallDir\n";
        }
        &set_c_compiler_conf("gcc_mpicc", undef, 0);
    }
    if ((defined $NewCuda or defined $NewOmp) and -f $MakefileConf) {
        system("make -s clean");
    }
    system("cd $AmrexDir && rm -f InstallDir && ln -s $InstallDir InstallDir && touch $InstallDir/include/AMReX_Config.H") == 0
        or die "$ERROR could not link $AmrexDir/InstallDir to $InstallDir\n";
    my $label = $useCuda ? "CUDA" : ($useOmp ? "OpenMP CPU" : "CPU");
    print "Linking FLEKS with the $label AMReX library $InstallDir\n";
}
################################################################################
sub print_help{

    print "
-s            Show the configuration for AMReX.

-debug        Compile with debug flags. Verbose GMRES residual logging is
              switched on: VERBOSE_GMRES = -DGMRES_VERBOSE in src/Makefile.
-nodebug      Do not use debug flags, and switch off the verbose GMRES
              residual logging (only the iteration count and the final
              residual are reported by the solver).

-simd[=ARCH]  Enable SIMD compiler flags (-march=native -fopenmp-simd by default,
              or -march=ARCH -fopenmp-simd if ARCH is given, e.g. x86-64-v3) in
              src/Makefile for C++ compilation. Enabled by default.
-nosimd       Disable extra SIMD compiler flags in src/Makefile.

-lev          Number of maximum grid levels. It is a uniform grid without AMR if lev=1.

-amrex2d      Link with the 2D AMReX library (AMReX_SPACEDIM = 2).
-amrex3d      Link with the 3D AMReX library (AMReX_SPACEDIM = 3). This is the default.
              AMReX is enabled automatically in the standalone FLEKS setup; by default the
              dimension of the currently-linked library is kept unless -amrex2d/-amrex3d is given.

-openmp       Build with OpenMP multithreading: enable OPENMPFLAG in Makefile.conf
              and link with an OpenMP build of AMReX in util/AMREX/InstallDir<N>D_omp
              (compiled on first use with USE_OMP=TRUE).
-noopenmp     Disable OpenMP flags and link with the serial CPU build of AMReX.

-cuda[=ARCH]  Build for NVIDIA GPUs: link with a CUDA build of AMReX in
              util/AMREX/InstallDir<N>D_cuda, compiled on first use for compute
              capability ARCH (default: auto-detected via nvidia-smi, or 80 if
              no GPU is present). The choice is kept by later Config.pl calls
              until -nocuda is given.
-nocuda       Link with the CPU build of AMReX.

-u            Show available user source options.

-u=NAME       Select user source file userfiles/NAMESource.h.

-tp=P,PB,PBE,PBED,PBEG  Test particle output information.
              P: only save particle velocity + location
              PB: particle + magnetic field. 
              PBE: particle + magnetic field + electric field
              PBED: particle + B-field + E-field + gradient/curvature drift velocities (default)
              PBEG: particle + B-field + E-field + raw gradient of B-field

\n";
    exit -0;
}
