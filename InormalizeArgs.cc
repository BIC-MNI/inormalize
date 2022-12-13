/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: InormalizeArgs.cc,v $
$Revision: 1.3 $
$Author: claude $
$Date: 2006-06-01 21:14:55 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include "InormalizeArgs.h"
#include <iostream>
using namespace std;

int    InormalizeArgs::clobber = FALSE;
int    InormalizeArgs::verbose = TRUE;
int    InormalizeArgs::printValues = FALSE;
int    InormalizeArgs::cache = TRUE;
int    InormalizeArgs::compress = FALSE;
int    InormalizeArgs::normalize = TRUE;
int    InormalizeArgs::useWorldCoord = FALSE;

char  *InormalizeArgs::modelString = 0;
int    InormalizeArgs::zNormalize = FALSE;
int    InormalizeArgs::yNormalize = FALSE;
int    InormalizeArgs::xNormalize = FALSE;

int    InormalizeArgs::method = MEDIAN_OF_RATIOS;
int    InormalizeArgs::minVoxels = 1000;
int    InormalizeArgs::constrained = FALSE;

char  *InormalizeArgs::maskString = 0;
#ifdef HAVE_MATLAB
char  *InormalizeArgs::matlabOutputString = 0;
#endif
double InormalizeArgs::thresholds[2] = { -MAXDOUBLE, MAXDOUBLE };
double InormalizeArgs::constants[2]  = { -MAXDOUBLE, -MAXDOUBLE };
double InormalizeArgs::rangePct      = -MAXDOUBLE;

// All this method stuff should be streamlined somehow
const int InormalizeArgs::RATIO_OF_MEANS     = 0;
const int InormalizeArgs::RATIO_OF_MEDIANS   = 1;
const int InormalizeArgs::MEAN_OF_RATIOS     = 2;
const int InormalizeArgs::MEAN_OF_LOG_RATIOS = 3;
const int InormalizeArgs::MEDIAN_OF_RATIOS   = 4;
const int InormalizeArgs::RMS                = 5;
const int InormalizeArgs::VR                 = 6;
const int InormalizeArgs::RANGE              = 7;
const char *InormalizeArgs::methods[] = { "ratioOfMeans", "ratioOfMedians",
					  "meanOfRatios", "meanOfLogRatios", 
					  "medianOfRatios", "rms", "vr" , "range" };
const int InormalizeArgs::nMethods = 8;


ArgvInfo InormalizeArgs::argTable[] = {
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "General options:"},
  {"-clobber", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::clobber, 
   "Overwrite existing file."},
  {"-noclobber", ARGV_CONSTANT, (char *) (int) FALSE, (char *) &InormalizeArgs::clobber, 
   "Do not overwrite existing file (default)."},
  {"-verbose", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::verbose, 
   "Print out log messages as processing is being done (default)."},
  {"-quiet", ARGV_CONSTANT, (char *) (int) FALSE, (char *) &InormalizeArgs::verbose, 
   "Do not print out any log messages."},
  {"-print", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::printValues, 
   "Calculate and print the values obtained with every solution method."},
  {"-compress", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::compress, 
   "Compress <outfile>."},
  {"-nocompress", ARGV_CONSTANT, (char *) (int) FALSE, (char *) &InormalizeArgs::compress, 
   "Do not compress <outfile> (default)."},
  {"-nocache", ARGV_CONSTANT, (char *) (int) FALSE, (char *) &InormalizeArgs::cache, 
   "Do not use volume caching."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nNormalization options:"},
  {"-const", ARGV_FLOAT, (char *) 1, (char *) InormalizeArgs::constants,
   "Specify a constant value to normalize to."},
  {"-const2", ARGV_FLOAT, (char *) 2, (char *) InormalizeArgs::constants,
   "Specify two constant values (for -range)."},
  {"-range", ARGV_FLOAT, (char *) 1, (char *) &InormalizeArgs::rangePct, 
   "Normalize the range of <infile> to const values or model.\n\t\t   Requires a float argument specifying the top- and bottom % to exclude, e.g., \"-range 5\""},
  {"-model", ARGV_STRING, (char *) 1, (char *) &InormalizeArgs::modelString,
   "Normalize <infile> to <model>."},
  {"-znormalize", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::zNormalize,
   "Normalize <infile> in z."},
  {"-ynormalize", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::yNormalize,
   "Normalize <infile> in y."},
  {"-xnormalize", ARGV_CONSTANT, (char *) (int) TRUE, (char *) &InormalizeArgs::xNormalize,
   "Normalize <infile> in x."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nSolution options:"},
  {"-minvoxels", ARGV_INT, (char *) 1, (char *) &InormalizeArgs::minVoxels,
   "Minimum number of voxels needed for estimation."},
  {"-constrained", ARGV_CONSTANT, (char *) (int) TRUE, (char *)&InormalizeArgs::constrained,
   "Constrain the ratio array to [10e-0.1, 10e+0.1]"},
  {"-rms", ARGV_CONSTANT, (char *)(int)InormalizeArgs::RMS, (char *) &InormalizeArgs::method,
   "Minimize RMS voxel difference."},
  {"-vr", ARGV_CONSTANT, (char *)(int)InormalizeArgs::VR, (char *) &InormalizeArgs::method,
   "Minimize variance of ratios."},
  {"-ratioOfMeans", ARGV_CONSTANT, (char *)(int)InormalizeArgs::RATIO_OF_MEANS, 
   (char *) &InormalizeArgs::method, "Use the ratio of voxel means."},
  {"-ratioOfMedians", ARGV_CONSTANT, (char *)(int)InormalizeArgs::RATIO_OF_MEDIANS, 
   (char *) &InormalizeArgs::method, "Use the ratio of voxel medians."},
  {"-meanOfRatios", ARGV_CONSTANT, (char *)(int)InormalizeArgs::MEAN_OF_RATIOS, 
   (char *) &InormalizeArgs::method, "Use the mean of voxel ratios."},
  {"-meanOfLogRatios", ARGV_CONSTANT, (char *)(int)InormalizeArgs::MEAN_OF_LOG_RATIOS, 
   (char *) &InormalizeArgs::method, "Use the mean of voxel ratios in the log domain."},
  {"-medianOfRatios", ARGV_CONSTANT, (char *)(int)InormalizeArgs::MEDIAN_OF_RATIOS, 
   (char *) &InormalizeArgs::method, "Use the median of voxel ratios (default)."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nMask options:"},
  {"-mask", ARGV_STRING, (char *) 1, (char *) &InormalizeArgs::maskString,
   "Ignore voxels for which <mask> is zero."},
  {"-useVoxelCoord", ARGV_CONSTANT, (char *) (int) FALSE,(char *)&InormalizeArgs::useWorldCoord,
   "Use voxel coordinates when matching mask (default)."},
  {"-useWorldCoord", ARGV_CONSTANT, (char *) (int) TRUE,(char *)&InormalizeArgs::useWorldCoord,
   "Use world coordinates when matching mask."},
  {"-threshold", ARGV_FLOAT, (char *) 2, 
   (char *) InormalizeArgs::thresholds,
   "Ignore voxels outside [<min>, <max>]."},

  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
   "\nOutput options:"},
#ifdef HAVE_MATLAB
  {"-matlab", ARGV_STRING, (char *) 1, (char *)&InormalizeArgs::matlabOutputString,
   "Save a sample (max 10000 points) of the voxels used in a matlab (.mat) file (-model only)."},
#endif
  {"-nonormalize", ARGV_CONSTANT, (char *) (int) FALSE, (char *) &InormalizeArgs::normalize,
   "Do not actually normalize the input volume (calculate mapping(s) only).\n"},

  {NULL, ARGV_END, NULL, NULL, NULL}
};
   
/* ----------------------------- MNI Header -----------------------------------
@NAME       : Constructor of InormalizeArgs
@INPUT      : argc - number of command-line arguments
              argv - command-line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Builds a InormalizeArgs object from commandline arguments
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 29, 1995 (Alex Zijdenbos)
@MODIFIED   : 
---------------------------------------------------------------------------- */

InormalizeArgs::InormalizeArgs(int argc, char **argv)
: command(argv[0])
{
  for (unsigned i = 1; i < argc; i++) {
    command += " ";
    command += argv[i];
  }

  // Call ParseArgv
  if (ParseArgv(&argc, argv, argTable, 0) || 
      (((argc != 3) && normalize) || ((argc != 2) && !normalize))) {
    cerr << endl << "Usage: " << argv[0] 
	 << " [<options>] <infile> [<outfile>]" << endl;
    cerr <<         "       " << argv[0] << " [-help]" << endl << endl;
    exit(EXIT_FAILURE);
  }

  // Check and expand paths
  MString extension;
  Boolean useModel = FALSE;
  if (modelString) {
    modelPath = modelString;
    if (!modelPath.expanded().existsCompressed(&extension)) {
      cerr << "Can't find " << modelPath << endl;
      exit(EXIT_FAILURE);
    }
    modelPath = modelPath.expanded() + extension;
    useModel = TRUE;
  }

  Boolean useMask = FALSE;
  if (maskString) {
    maskPath = maskString;
    if (!maskPath.expanded().existsCompressed(&extension)) {
      cerr << "Can't find " << maskPath << endl;
      exit(EXIT_FAILURE);
    }
    maskPath = maskPath.expanded() + extension;
    useMask = TRUE;
  }

  inputPath = argv[1];
  if (!inputPath.expanded().existsCompressed(&extension)) {
    cerr << "Can't find " << inputPath << endl;
    exit(EXIT_FAILURE);
  }
  inputPath = inputPath.expanded() + extension;

  if (normalize) {
    outputPath = argv[2];
    outputPath = outputPath.expanded().removeCompressedExtension();
    if (!clobber) {
      MString extension;
      if (outputPath.existsCompressed(&extension)) {
	cerr << outputPath + extension << " exists!" << endl;
	exit(EXIT_FAILURE);
      }
    }
  }

  // Verify other args
  Boolean sliceNorm = (xNormalize || yNormalize || zNormalize);
  oneConst          = (constants[0] != -MAXDOUBLE) && !(constants[1] != -MAXDOUBLE);
  twoConst          = (constants[0] != -MAXDOUBLE) && (constants[1] != -MAXDOUBLE);
  useRange          = (rangePct != -MAXDOUBLE);

  if (!sliceNorm && !useModel && !oneConst && !twoConst) {
    cerr << "Please specify one of -{x,y,z}normalize, -model, -const, or -const2" << endl;
    exit(EXIT_FAILURE);
  }

  if (sliceNorm) {
    if (useModel) {
      cerr << "You cannot use -model with -{x,y,z}normalize" << endl;
      exit(EXIT_FAILURE);
    }
    if (oneConst || twoConst) {
      cerr << "You cannot use -const or -const2 with -{x,y,z}normalize" << endl;
      exit(EXIT_FAILURE);
    }
    if (useRange) {
      cerr << "-range with -{x,y,z}normalize is unreliable!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  if (useModel) {
    if (oneConst || twoConst) {
      cerr << "You cannot use -const or -const2 with -model" << endl;
      exit(EXIT_FAILURE);
    }
  }    

  if (useRange) {
    if ((rangePct < 0) || (rangePct >= 50)) {
      cerr << "Range % must be in [0, 50>" << endl;
      exit(EXIT_FAILURE);
    }

    method = RANGE;
  }

  if (twoConst && !useRange) {
    cerr << "You must specify -range <pct> with -const2" << endl;
    exit(EXIT_FAILURE);
  }

  if (oneConst &&
      !(method == RATIO_OF_MEANS) &&
      !(method == RATIO_OF_MEDIANS)) {
    cerr << "You must specify either -ratioOfMeans or -ratioOfMedians with -const" << endl;
    exit(EXIT_FAILURE);
  }

  // Do not cache main volume, since range modification is not yet implemented 
  // with volume caching.
  volume = loadVolume(inputPath, verbose);
  if (!volume) {
    cerr << "Couldn't load " << inputPath << endl;
    exit(EXIT_FAILURE);
  }

  if (cache) {
    set_n_bytes_cache_threshold(0);           // Always cache volume
    set_default_max_bytes_in_cache(524288);   // 0.5M cache size
    //set_default_max_bytes_in_cache(0);      // Cache 1 block at a time
    set_cache_block_sizes_hint(SLICE_ACCESS); // Cache volume slices
  }

  model = 0;
  if (useModel) {
    model = loadVolume(modelPath, verbose);
    if (!model) {
      cerr << "Couldn't load " << modelPath << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  mask = 0;
  if (useMask) {
    mask = loadVolume(maskPath, verbose);
    if (!mask) {
      cerr << "Couldn't load " << maskPath << endl;
      exit(EXIT_FAILURE);
    }
    //    int blockSizes[3];
    //    blockSizes[0] = 180;
    //    blockSizes[1] = 216;
    //    blockSizes[2] = 1;
    //    set_volume_cache_block_sizes(mask, blockSizes);
  }
}

//
//
//
VIO_Volume
loadVolume(const Path& path, int verbose)
{
  VIO_Real amountDone;
  volume_input_struct inputInfo;
  VIO_Volume volume = 0;

  if (verbose)
    cout << "Reading volume " << path << flush;
  if (start_volume_input((char *) path.string(), 3, NULL, NC_UNSPECIFIED, 
			 TRUE, 0.0, 0.0, TRUE, &volume, (minc_input_options *) NULL, 
			 &inputInfo) != VIO_OK)
    return 0;

  while (input_more_of_volume(volume, &inputInfo, &amountDone))
    if (verbose)
      cout << "." << flush;
  if (verbose)
    cout << "Done" << endl;
  delete_volume_input(&inputInfo);

  return volume;
}

