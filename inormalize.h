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
$RCSfile: inormalize.h,v $
$Revision: 1.2 $
$Author: jason $
$Date: 2002-03-28 23:53:20 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _I_NORMALIZE_H
#define _I_NORMALIZE_H

#include "EBTKS/SimpleArray.h"
#include "InormalizeArgs.h"
#include "EBTKS/Minc.h"
#include "EBTKS/ValueMap.h"

typedef SimpleArray<float>   FloatArray;
typedef SimpleArray<Boolean> BoolArray;
//typedef SimpleArray<float>   FloatArray;

void      scanVoxelRange(const VIO_Volume volume, double *voxelMin, double *voxelMax);
unsigned  getValidVoxels(const InormalizeArgs& args, BoolArray& VV);
void      floatArrayFromVolume(FloatArray& array, const VIO_Volume volume, 
			       const BoolArray& validVoxels, unsigned *N = 0, 
			       int verbose = 0);
void      floatArraysFromSlices(const VIO_Volume volume, const BoolArray& validVoxels,
				unsigned axis, unsigned slice1, unsigned slice2, 
				FloatArray& array1, FloatArray& array2);
void      selfNormalizeMain(char *dimension, const BoolArray& validVoxels,
			    const InormalizeArgs& args, MString& history);
Array<LinearMap> 
          selfNormalize(VIO_Volume volume, const BoolArray& validVoxels, 
			char *dimension, const InormalizeArgs& args, int method = -1);
void      reMapVolume(VIO_Volume volume, const LinearMap& iMap, int verbose = 0);
void      reMapVolume(VIO_Volume volume, int axis, const Array<LinearMap>& iMaps, 
		      int verbose = 0);
LinearMap determineMap(const FloatArray& modelArray, 
		       const FloatArray& dataArray, 
		       const InormalizeArgs& args, int method = -1);
Boolean   saveVolume(const VIO_Volume volume, const Path& path, const Path& mincModel,
		     const MString& history, int compress = 1, int verbose = 0);
MString&  appendToString(MString& string, const LinearMap& iMap);
MString&  appendToString(MString& string, const Array<LinearMap>& iMaps);
double    evaluateRMS(void *data, float *f);
double    evaluateVR(void *data, float *f);
FloatArray pctExtrema(const FloatArray& array, double pct, int verbose = 0);

#endif
