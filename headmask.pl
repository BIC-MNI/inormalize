#!xPERLx
#---------------------------------------------------------------------------
#@COPYRIGHT :
#             Copyright 1996, Alex P. Zijdenbos
#             McConnell Brain Imaging Centre,
#             Montreal Neurological Institute, McGill University.
#             Permission to use, copy, modify, and distribute this
#             software and its documentation for any purpose and without
#             fee is hereby granted, provided that the above copyright
#             notice appear in all copies.  The author and McGill University
#             make no representations about the suitability of this
#             software for any purpose.  It is provided "as is" without
#             express or implied warranty.
#---------------------------------------------------------------------------- 
#$RCSfile: headmask.pl,v $
#$Revision: 1.1 $
#$Author: jason $
#$Date: 2002-04-09 17:20:07 $
#$State: Exp $
#---------------------------------------------------------------------------

require "ctime.pl";
require "file_utilities.pl";
require "minc_utilities.pl";
require ParseArgs;

use Startup;
use JobControl;

&Startup();

&Initialize();

&HeadMask($InFile, $OutFile);

print "Updating history\n" if $Verbose;
my(@history) = &get_history($InFile);
push(@history, $HistoryLine);
&put_history($OutFile, @history);

&Cleanup(1);

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &SetHelp
#@INPUT      : none
#@OUTPUT     : none
#@RETURNS    : nothing
#@DESCRIPTION: Sets the $Help and $Usage globals, and registers them
#              with ParseArgs so that user gets useful error and help
#              messages.
#@METHOD     : 
#@GLOBALS    : $Help, $Usage
#@CALLS      : 
#@CREATED    : 95/08/25, Greg Ward (from code formerly in &ParseArgs)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub SetHelp
{
   $Usage = <<USAGE;
Usage: $ProgramName [options] <in.mnc> <out.mnc>

USAGE

   $Help = <<HELP;

$ProgramName 
   creates a rough binary head mask by thresholding <in.mnc> at a
   value derived from the histogram under the assumption that the
   histogram is bimodal.

HELP

   &ParseArgs::SetHelpText ($Help, $Usage);
}

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &Initialize
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Sets global variables, parses command line, finds required 
#              programs, and sets their options.  Dies on any error.
#@METHOD     : 
#@GLOBALS    : general: $Verbose, $Execute, $Clobber, $KeepTmp
#              mask:    $Mask
#@CALLS      : &SetupArgTables
#              &ParseArgs::Parse
#              
#@CREATED    : 96/04/27, Alex Zijdenbos
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub Initialize
{
    chop ($ctime = &ctime(time));
    $HistoryLine = "$ctime>>> $0 @ARGV";

    $Clobber  = 0;
    $Execute  = 1;
    $Verbose  = 1;

    $Mask = ();

    &SetHelp;

    ($HeadMaskArgsTbl) = &SetupArgTables;

    &ParseArgs::Parse ([@DefaultArgs, @$HeadMaskArgsTbl,], \@ARGV) || exit 1;

    if ($#ARGV != 1) {
	die "Please supply two file arguments!\n";
    }
    $InFile  = shift(@ARGV);
    $OutFile = shift(@ARGV);

    &CheckFiles($InFile) || &Fatal();

    $OutFile =~ s/\.(Z|gz|z)$//;
    if (!$Clobber && -e $OutFile) {
	die "Output file $OutFile exists; use -clobber to overwrite\n";
    }

    &SetSpawnOptions ("Verbose", $Verbose, "Execute", $Execute,
		      "ErrorAction", "&cleanup_and_die()");
    
    $SearchPath = defined($ENV{'MS_LESION'}) ? 
	"$ENV{'MS_LESION'}/bin:$ENV{'MS_LESION'}/extern/bin" :
	    "$ENV{'PATH'}:/usr/people/alex/Release/bin";	    
    
    ($MincMath       = &FindProgram ("mincmath",     $SearchPath)) || &Fatal();
    ($VolumeStats    = &FindProgram ("volume_stats", $SearchPath)) || &Fatal();
}

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &SetupArgTables
#@INPUT      : none
#@OUTPUT     : none
#@RETURNS    : References to the option tables:
#                @HeadMaskArgs
#@DESCRIPTION: Defines the tables of command line (and config file) 
#              options that we pass to ParseArgs.  There are four
#              separate groups of options, because not all of them
#              are valid in all places.  See comments in the routine
#              for details.
#@METHOD     : 
#@GLOBALS    : makes references to many globals (almost all of 'em in fact)
#              even though most of them won't have been defined when
#              this is called
#@CALLS      : 
#@CREATED    : 96/08/08, Alex Zijdenbos
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub SetupArgTables
{
    my (@HeadMaskArgs);

    # Preferences -- these may be given in the configuration file
    # or the command line

    @HeadMaskArgs = 
	(["Mask options", "section"],
	 ["-mask", "string", 1, \$Mask,
	  "calculate the threshold value over the voxels included in the specified mask only."]);

    (\@HeadMaskArgs);
}

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &HeadMask
#@INPUT      : $file
#@OUTPUT     : none
#@RETURNS    : $mask
#@DESCRIPTION: Creates a rough non-BG mask
#@METHOD     : 
#@GLOBALS    : Standard ($Execute, ...)
#@CALLS      : 
#@CREATED    : 96/08/08, Alex Zijdenbos
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub HeadMask {
    my($file, $mask) = @_;

    # Obtain BG threshold
    $VolumeStatsCmd = "$VolumeStats -quiet -biModalT $file";
    if (defined($Mask)) {
	$VolumeStatsCmd .= " -mask $Mask";
    }
	
    my($result, $T) = &Spawn($VolumeStatsCmd);
    chop($T);

    &Spawn(&AddOptions("$MincMath -byte -const $T -ge $file $mask", $Verbose, 1));
}

sub AddOptions {
    my($string, $verbose, $clobber) = @_;

    if ($clobber) {
	$string =~ s/^([^\s]+)(.*)$/$1 -clobber$2/;
    }

    if ($verbose) {
	$string =~ s/^([^\s]+)(.*)$/$1 -verbose$2/;
    }
    else {
	$string =~ s/^([^\s]+)(.*)$/$1 -quiet$2/;
    }

    $string;
}