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
$RCSfile: InormalizeArgs.h,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2005-08-17 23:14:26 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _I_NORMALIZE_ARGS_H
#define _I_NORMALIZE_ARGS_H

extern "C" {			/* (bert)-use C linkage */
#include <ParseArgv.h>
}
#include "EBTKS/Minc.h"
#include "EBTKS/Path.h"

class InormalizeArgs {
public:
  MString command;
  Path    inputPath;
  Path    modelPath;
  Path    maskPath;
  Path    outputPath;

  VIO_Volume  volume, model, mask;
  Boolean oneConst, twoConst, useRange;

  static int    clobber;
  static int    verbose;
  static int    printValues;
  static int    compress;
  static int    cache;
  static int    normalize;

  static int    useWorldCoord;

  static char  *modelString;
  static int    zNormalize, yNormalize, xNormalize;

  static int    method;
  static int    minVoxels;
  static int    constrained;

  static char  *maskString;
  static double thresholds[2];
  static double constants[2];
  static double rangePct;

#ifdef HAVE_MATLAB
  static char  *matlabOutputString;
#endif

  static ArgvInfo argTable[];

  // All this methods stuff should be streamlined somehow
  static const int RMS, VR;
  static const int RATIO_OF_MEANS, RATIO_OF_MEDIANS;
  static const int MEAN_OF_RATIOS, MEAN_OF_LOG_RATIOS, MEDIAN_OF_RATIOS;
  static const int RANGE;
  static const char *methods[];
  static const int nMethods;

  // Constructors/destructor
  InormalizeArgs(int argc, char **argv);
  ~InormalizeArgs() {}
};

VIO_Volume loadVolume(const Path& path, int verbose = 0);

#endif

