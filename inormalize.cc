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
$RCSfile: inormalize.cc,v $
$Revision: 1.4 $
$Author: claude $
$Date: 2007-10-09 21:56:35 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : inormalize
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : error status
@DESCRIPTION: Program for the intensity normalization of a minc volume,
                following a (registered) model volume.
@METHOD     : Minimization of the MSE voxel intensity difference by exhaustive
                search
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 11, 1995 (Alex Zijdenbos)
@MODIFIED   : $Log: inormalize.cc,v $
@MODIFIED   : Revision 1.4  2007-10-09 21:56:35  claude
@MODIFIED   : fixed inormalize for float data type
@MODIFIED   :
@MODIFIED   : Revision 1.3  2005/08/17 23:14:26  bert
@MODIFIED   : Minor changes for some C++ warnings and issues
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/28 23:53:19  jason
@MODIFIED   : switched to SimpleArray from CachedArray as there appear to be IO problems in the cached array implementation
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  2002/03/27 18:36:59  jason
@MODIFIED   : First import of inormalize sources. Does not work perfectly yet,
@MODIFIED   : especially the median estimation appears to fail.
@MODIFIED   :
@MODIFIED   :
@MODIFIED   : Revision 1.6  1998/04/01 18:29:49  alex
@MODIFIED   : Made OpTimers quiet if -quiet was specified
@MODIFIED   :
@MODIFIED   : Revision 1.5  1998/02/28 20:04:34  alex
@MODIFIED   : Added -range and -const options
@MODIFIED   :
@MODIFIED   : Revision 1.4  1997/09/11 14:26:31  alex
@MODIFIED   : Various changes directed at g++ compatibility
@MODIFIED   :
@MODIFIED   : Revision 1.3  1996/12/19 20:47:16  alex
@MODIFIED   : Various changes directed at g++ compatibility
@MODIFIED   :
 * Revision 1.2  1996/09/04  20:56:35  alex
 * Checked Inormalize into main branch
 *
 * Revision 1.1.1.1  1996/08/29  19:05:04  alex
 * Source for inormalize
 *
@COPYRIGHT  :
---------------------------------------------------------------------------- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "inormalize.h"
#include <iostream>
using namespace std;
#include <time.h>
//  #include "EBTKS/CachedArray.h"
#include "EBTKS/Dictionary.h"
#include "EBTKS/FileIO.h"
#include "EBTKS/OpTimer.h"
#include "EBTKS/amoeba.h"
#include "EBTKS/trivials.h"

class OptData {
public:
  const FloatArray& modelArray;
  const FloatArray& dataArray;
  unsigned n;
  Boolean verbose;

  OptData(const FloatArray& ma, const FloatArray& da, unsigned nEl,
	  Boolean vb = FALSE) 
    : modelArray(ma),
      dataArray(da)
  { n = nEl; verbose = vb; }
};

//
// Main program
//

int 
main(int argc, char *argv[])
{
  // Set timers
  OpTimer usr(OpTimer::USR, "inormalize");
  OpTimer sys(OpTimer::SYS, "inormalize");
  OpTimer cpu(OpTimer::CPU, "inormalize");

  // Get argument information
  InormalizeArgs args(argc, argv);

  if (!args.verbose) {
    usr.verbose(FALSE);
    sys.verbose(FALSE);
    cpu.verbose(FALSE);
  }

  // Determine valid voxels from the mask and threshold settings
  SimpleArray<Boolean> validVoxels;
  unsigned nVoxels = getValidVoxels(args, validVoxels);
  if (args.verbose)
    cout << "Considering " << nVoxels << " voxels (" 
	 << 100.0*nVoxels/validVoxels.size() << "% of total)" << endl;

  if (!nVoxels) {
    cerr << "Error! voxel selection criterai removed all voxels. Not normalizing volume" << endl;
    exit(EXIT_FAILURE);
  }

  if (args.verbose)
    cout << "Selected normalization: " << args.methods[args.method] << endl;

  // Can delete mask volume here
  if (args.mask) {
    delete_volume(args.mask);
    args.mask = 0;
  }

  // Create history to append to MINC file
  time_t  timeVal = time(0);
  MString history(ctime(&timeVal));
  history.chop();
  history += ">>> ";
  history += args.command;
  history += "; ";

  // Self-normalize the input volume
  if (args.zNormalize)
    selfNormalizeMain(MIzspace, validVoxels, args, history);

  if (args.yNormalize)
    selfNormalizeMain(MIyspace, validVoxels, args, history);

  if (args.xNormalize)
    selfNormalizeMain(MIxspace, validVoxels, args, history);

  if (args.oneConst) {

    SimpleArray<float> volumeArray;

    floatArrayFromVolume(volumeArray, args.volume, validVoxels, &nVoxels, args.verbose);

    LinearMap iMap;
    if (args.method == InormalizeArgs::RATIO_OF_MEANS)
      iMap.factor(args.constants[0] / mean(volumeArray));
    else if (args.method == InormalizeArgs::RATIO_OF_MEDIANS) {
      iMap.factor(args.constants[0] / median(volumeArray));
      cout << median(volumeArray) << endl;
      cout << mean(volumeArray) << endl;
    }
    else {
      cerr << "Invalid method " << args.methods[args.method] << endl;
      exit(EXIT_FAILURE);
    }
    
    MString mapString;
    appendToString(mapString, iMap);

    if (args.normalize) {
      history += mapString;
      reMapVolume(args.volume, iMap, args.verbose);
    }
  }

  else if (args.twoConst) {
    SimpleArray<float> volumeArray;
    floatArrayFromVolume(volumeArray, args.volume, validVoxels, &nVoxels, args.verbose);

    FloatArray extrema(pctExtrema(volumeArray, args.rangePct, args.verbose));

    LinearMap iMap(extrema[(unsigned int)0], 
                   extrema[(unsigned int)1], 
                   args.constants[(unsigned int)0], 
                   args.constants[(unsigned int)1]);
    
    MString mapString;
    appendToString(mapString, iMap);
    
    if (args.normalize) {
      history += mapString;
      reMapVolume(args.volume, iMap, args.verbose);
    }
  }

  else if (args.model) {

    // Make sure volume and model have compatible sizes.
    int volume_sizes[3], model_sizes[3];
    get_volume_sizes(args.volume, volume_sizes);
    get_volume_sizes(args.model, model_sizes);
    if( ( model_sizes[0] != volume_sizes[0] ) || 
        ( model_sizes[1] != volume_sizes[1] ) ||
        ( model_sizes[2] != volume_sizes[2] ) ) {
      cerr << "VIO_Volume and model dimensions do not match!" << endl;
      exit(EXIT_FAILURE);
    }

    // Convert entire volume to float arrays (cached to conserve memory)
    SimpleArray<float> volumeArray;
    SimpleArray<float> modelArray;

    floatArrayFromVolume(volumeArray, args.volume, validVoxels, &nVoxels, args.verbose);
    floatArrayFromVolume(modelArray, args.model, validVoxels, &nVoxels, args.verbose);

    // Can delete model volume here
    delete_volume(args.model);
    args.model = 0;

#ifdef HAVE_MATLAB
    if (args.matlabOutputString) {
      SimpleArray<float> samples = SimpleArray<float>(volumeArray.sample(10000));
      samples.saveMatlab(args.matlabOutputString, "volume");
      samples = SimpleArray<float>(modelArray.sample(10000));
      samples.saveMatlab(args.matlabOutputString, "model");
    }
#endif

    history += "global:";

    for (unsigned method = 0; method < args.nMethods; method++)
      if ((method == args.method) || args.printValues) {
	if (args.printValues)
	  cout << "  " << args.methods[method] << ":" << flush;

	LinearMap iMap = determineMap(modelArray, volumeArray, args, method);
	MString mapString;
	appendToString(mapString, iMap);
	if (args.printValues)
	  cout << mapString << endl;

#ifdef HAVE_MATLAB
	if (args.matlabOutputString) {
	  SimpleArray<float> mapArray(2);
	  LinearMap tempMap(iMap);
	  tempMap.inv();
	  mapArray[0] = tempMap.factor();
	  mapArray[1] = tempMap.offset();
	  mapArray.saveMatlab(args.matlabOutputString, args.methods[method]);
	}
#endif
	if (args.normalize && (args.method == method)) {
	  history += mapString;
	  reMapVolume(args.volume, iMap, args.verbose);
	}
      }
  }

  if (args.normalize) {
    history += "\n";
    if (!saveVolume(args.volume, args.outputPath, args.inputPath, history, 
		    args.compress, args.verbose))
      exit(EXIT_FAILURE);
  }
  
  return(EXIT_SUCCESS);
}

//
//
//
void
scanVoxelRange(const VIO_Volume volume, double *voxelMin, double *voxelMax)
{
  int sizes[3];
  get_volume_sizes(volume, sizes);
  unsigned D1 = sizes[0];
  unsigned D2 = sizes[1];
  unsigned D3 = sizes[2];

  VIO_Real minVal = get_volume_voxel_value(volume, 0, 0, 0, 0, 0);
  VIO_Real maxVal = minVal;

  for (unsigned d1 = 0; d1 < D1; d1++)
    for (unsigned d2 = 0; d2 < D2; d2++)
      for (unsigned d3 = 0; d3 < D3; d3++) {
	VIO_Real value = get_volume_voxel_value(volume, d1, d2, d3, 0, 0);
	if (value < minVal)
	  minVal = value;
	else if (value > maxVal)
	  maxVal = value;
      }
  
  *voxelMin = minVal;
  *voxelMax = maxVal;
}

//
//
//
unsigned
getValidVoxels(const InormalizeArgs& args, BoolArray& validVoxels)
{
  int sizes[3];
  get_volume_sizes(args.volume, sizes);
  unsigned D1 = sizes[0];
  unsigned D2 = sizes[1];
  unsigned D3 = sizes[2];
  unsigned nVoxels = D1*D2*D3;
  validVoxels.newSize(nVoxels);

  double minT = args.thresholds[0];
  double maxT = args.thresholds[1];

  Boolean useThreshold = 
    ((maxT >= minT) && ((minT > -MAXDOUBLE) || (maxT < MAXDOUBLE)));

  if (args.mask || useThreshold) {
    validVoxels.clear(FALSE);
    unsigned maskD1, maskD2, maskD3;
    if (args.mask) {
      if (args.verbose)
	cout << "Scanning mask" << flush;
      int maskSizes[3];
      get_volume_sizes(args.mask, maskSizes);
      maskD1 = maskSizes[0];
      maskD2 = maskSizes[1];
      maskD3 = maskSizes[2];

      if (!args.useWorldCoord && ((maskD1 != D1) || (maskD2 != D2) || (maskD3 != D3))) {
	cerr << "VIO_Volume and mask dimensions do not match!" << endl;
	exit(EXIT_FAILURE);
      }
    }
    
    nVoxels    = 0;
    unsigned i = 0;

    for (unsigned d1 = 0; d1 < D1; d1++) {
      for (unsigned d2 = 0; d2 < D2; d2++) {
	for (unsigned d3 = 0; d3 < D3; d3++) {
	  Boolean valid = TRUE;
	  if (useThreshold) {
	    VIO_Real value = get_volume_real_value(args.volume, d1, d2, d3, 0, 0);
	    if ((value < minT) || (value > maxT))
	      valid = FALSE;
	  }
	  
	  if (valid && args.mask) {
	    VIO_Real maskValue = 0;
	    if (args.useWorldCoord) {
	      VIO_Real xWorld, yWorld, zWorld;
	      convert_3D_voxel_to_world(args.volume, d1, d2, d3, 
					&xWorld, &yWorld, &zWorld);
	      VIO_Real voxel1, voxel2, voxel3;
	      convert_3D_world_to_voxel(args.mask, xWorld, yWorld, zWorld,
					&voxel1, &voxel2, &voxel3);
	      int v1 = VIO_ROUND(voxel1);
	      int v2 = VIO_ROUND(voxel2);
	      int v3 = VIO_ROUND(voxel3);
	      if ((v1 >= 0) && (v2 >= 0) && (v3 >= 0) &&
		  (v1 < maskD1) && (v2 < maskD2) && (v3 < maskD3))
		maskValue = get_volume_real_value(args.mask, v1, v2, v3, 0, 0);
	    }
	    else
	      maskValue = get_volume_real_value(args.mask, d1, d2, d3, 0, 0);

	    if (!maskValue)
	      valid = FALSE;
	  }
	  
	  if (validVoxels[i++] = valid)
	    nVoxels++;
	}
      }
      if (args.verbose)
	cout << "." << flush;
    }
    if (args.verbose)
      cout << "Done" << endl;
  }
  else
    validVoxels.clear(TRUE);
  
  return nVoxels;
}

//
//
//
void
floatArrayFromVolume(FloatArray& array, const VIO_Volume volume, 
		     const BoolArray& validVoxels, unsigned *N, 
		     int verbose)
{
  unsigned nVoxels = (N && *N) ? *N : validVoxels.occurrencesOf(TRUE);
  if (N && !*N)
    *N = nVoxels;

  if (verbose)
    cout << "Converting volume" << flush;

  array.newSize(nVoxels);

  int sizes[3];
  get_volume_sizes(volume, sizes);
  unsigned D1 = sizes[0];
  unsigned D2 = sizes[1];
  unsigned D3 = sizes[2];

  unsigned voxelCtr = 0;
  unsigned valueCtr = 0;

  for (unsigned d1 = 0; d1 < D1; d1++) {
    for (unsigned d2 = 0; d2 < D2; d2++) {
      for (unsigned d3 = 0; d3 < D3; d3++)
	if (validVoxels[voxelCtr++]) 
	  array[valueCtr++] = (float) get_volume_real_value(volume, d1, d2, d3, 0, 0);
    }
    
    if (verbose)
      cout << "." << flush;
  }

  if (verbose)
    cout << "Done" << endl;
}

//
//
//
void
floatArraysFromSlices(const VIO_Volume volume, const BoolArray& validVoxels,
		      unsigned axis, unsigned slice1, unsigned slice2, 
		      FloatArray& array1, FloatArray& array2)
{
  int sizes[3];
  get_volume_sizes(volume, sizes);
  unsigned D1 = sizes[0];
  unsigned D2 = sizes[1];
  unsigned D3 = sizes[2];

  switch(axis) {
  case 0: {
    unsigned N = D2*D3;
    array1.newSize(N);
    array1.newSize(0); // Allocated but empty
    array2.newSize(N);
    array2.newSize(0); // Allocated but empty
    unsigned valid1 = slice1*N;
    unsigned valid2 = slice2*N;
    //    const Boolean *valid1 = validVoxels.contents() + slice1*N;
    //    const Boolean *valid2 = validVoxels.contents() + slice2*N;
    for (unsigned d2 = 0; d2 < D2; d2++)
      for (unsigned d3 = 0; d3 < D3; d3++, valid1++, valid2++)
	if (validVoxels[valid1] && validVoxels[valid2]) {
	  VIO_Real value1 = get_volume_real_value(volume, slice1, d2, d3, 0, 0);
	  VIO_Real value2 = get_volume_real_value(volume, slice2, d2, d3, 0, 0);
	  if (value1 && value2) {
	    array1.append(value1);
	    array2.append(value2);
	  }
	}
    break;
  }
  case 1: {
    unsigned N = D1*D3;
    array1.newSize(N);
    array1.newSize(0); // Allocated but empty
    array2.newSize(N);
    array2.newSize(0); // Allocated but empty
    //    const Boolean *validBase = validVoxels.contents();
    for (unsigned d1 = 0; d1 < D1; d1++) {
      unsigned valid1 = (d1*D2 + slice1)*D3;
      unsigned valid2 = (d1*D2 + slice2)*D3;
      //      const Boolean *valid1 = validBase + (d1*D2 + slice1)*D3;
      //      const Boolean *valid2 = validBase + (d1*D2 + slice2)*D3;
      for (unsigned d3 = 0; d3 < D3; d3++, valid1++, valid2++)
	if (validVoxels[valid1] && validVoxels[valid2]) {
	  VIO_Real value1 = get_volume_real_value(volume, d1, slice1, d3, 0, 0);
	  VIO_Real value2 = get_volume_real_value(volume, d1, slice2, d3, 0, 0);
	  if (value1 && value2) {
	    array1.append(value1);
	    array2.append(value2);
	  }
	}
    }
    break;
  }
  case 2: {
    unsigned N = D1*D2;
    array1.newSize(N);
    array1.newSize(0); // Allocated but empty
    array2.newSize(N);
    array2.newSize(0); // Allocated but empty
    //    const Boolean *validBase = validVoxels.contents();
    for (unsigned d1 = 0; d1 < D1; d1++) {
      //      const Boolean *valid1 = validBase + d1*D2*D3 + slice1;
      //      const Boolean *valid2 = validBase + d1*D2*D3 + slice2;
      unsigned valid1 = d1*D2*D3 + slice1;
      unsigned valid2 = d1*D2*D3 + slice2;
      for (unsigned d2 = 0; d2 < D2; d2++) {
	if (validVoxels[valid1] && validVoxels[valid2]) {
	  VIO_Real value1 = get_volume_real_value(volume, d1, d2, slice1, 0, 0);
	  VIO_Real value2 = get_volume_real_value(volume, d1, d2, slice2, 0, 0);
	  if (value1 && value2) {
	    array1.append(value1);
	    array2.append(value2);
	  }
	}
	valid1 += D3;
	valid2 += D3;
      }
    }
    break;
  }
  }
}

//
//
//
void
selfNormalizeMain(char *dimension, const BoolArray& validVoxels, 
		  const InormalizeArgs& args, MString& history)
{
  for (unsigned method = 0; method < args.nMethods; method++)
    if ((method == args.method) || args.printValues) {
      Array<LinearMap> iMaps(selfNormalize(args.volume, validVoxels, dimension, 
					   args, method));
      MString mapString(dimension);
      mapString += ":";
      appendToString(mapString, iMaps);
      
      if (args.method == method)
	history += mapString;
      if (args.printValues)
	cout << "  " << args.methods[method] << ":" << mapString << endl;
    }
}

//
//
//
Array<LinearMap>
selfNormalize(VIO_Volume volume, const BoolArray& validVoxels, char *dimension, 
	      const InormalizeArgs& args, int method)
{
  if (method < -1)
    method = args.method;

  int sizes[3];
  get_volume_sizes(volume, sizes);
  int axis;

  if (!convert_dim_name_to_spatial_axis(dimension, &axis)) {
    cerr << "inormalize: couldn't convert dim_name " << dimension 
	 << " to axis" << endl;
    exit(EXIT_FAILURE);
  }

  axis = 2 - axis; // ??? Why is this necessary?

  unsigned sliceDim = sizes[axis];

  Array<LinearMap> iMaps(sliceDim);

  if (sliceDim < 2)
    return iMaps;

  unsigned centerSlice = sliceDim/2;

  unsigned int slice;
  FloatArray refArray, array;
  for (slice = centerSlice; slice < sliceDim - 1; slice++) {
    floatArraysFromSlices(volume, validVoxels, axis, slice, slice+1, refArray, array);
    LinearMap sliceMap(determineMap(refArray, array, args, method));

    iMaps[slice + 1] = iMaps[slice];
    iMaps[slice + 1].concat(sliceMap);
  }
  
  for (slice = centerSlice; slice > 0; slice--) {
    floatArraysFromSlices(volume, validVoxels, axis, slice, slice-1, refArray, array);
    LinearMap sliceMap(determineMap(refArray, array, args, method));

    iMaps[slice - 1] = iMaps[slice];
    iMaps[slice - 1].concat(sliceMap);
  }

  for (slice = 0; slice < sliceDim; slice++)
    cout << slice << ": " << iMaps[slice].factor() << endl;

  reMapVolume(volume, axis, iMaps, args.verbose);

  return iMaps;
}

//
//
//
void
reMapVolume(VIO_Volume volume, const LinearMap& iMap, int verbose)
{
  // Scan voxel range, just in case the volume attributes are messed up
  double voxelMin, voxelMax;
  scanVoxelRange(volume, &voxelMin, &voxelMax); 

  // Get real extrema
  double realMin = CONVERT_VOXEL_TO_VALUE(volume, voxelMin);
  double realMax = CONVERT_VOXEL_TO_VALUE(volume, voxelMax);

  // cout << "(" << voxelMin << ", " << voxelMax << ") -> (" 
  //      << realMin << ", " << realMax << ")" << endl;

  if (verbose)
    cout << "Normalizing (" << iMap.factor() << ", " << iMap.offset() << ")" << endl;

  // Concatenate original and derived mapping
  realMin = iMap(realMin);
  realMax = iMap(realMax);

  //LinearMap fullMap(voxelMin, voxelMax, realMin, realMax);
  
  //cout << "FullMap: " << fullMap << endl;

  // cout << "(" << voxelMin << ", " << voxelMax << ") -> (" 
  // << realMin << ", " << realMax << ")" << endl;

  // Scale the voxel values if not working with ranges (for FLOAT and DOUBLE).
  if( !volume->real_range_set ) {
    cout << "Fixing voxel values for real data type..." << endl;
    int sizes[3];
    get_volume_sizes(volume, sizes);
    for (unsigned d1 = 0; d1 < sizes[0]; d1++) {
      for (unsigned d2 = 0; d2 < sizes[1]; d2++) {
        for (unsigned d3 = 0; d3 < sizes[2]; d3++) {
	  VIO_Real value = get_volume_voxel_value(volume, d1, d2, d3, 0, 0);
	  set_volume_voxel_value(volume, d1, d2, d3, 0, 0, 
			         ::clamp(iMap(value), realMin, realMax));
        }
      }
    }
  } else {
    // Reset the voxel range and set the desired real range
    set_volume_voxel_range(volume, voxelMin, voxelMax);
  }
  set_volume_real_range(volume, realMin, realMax);

}

//
//
//
void
reMapVolume(VIO_Volume volume, int axis, const Array<LinearMap>& iMaps, int verbose)
{
  int sizes[3];
  get_volume_sizes(volume, sizes);
  unsigned D1 = sizes[0];
  unsigned D2 = sizes[1];
  unsigned D3 = sizes[2];

  // Scan voxel range, just in case the volume attributes are messed up
  double voxelMin = MAXDOUBLE;
  double voxelMax = -MAXDOUBLE;

  double mappedRealMin = MAXDOUBLE;
  double mappedRealMax = -MAXDOUBLE;

  // Traverse volume to find new global min and max
  LinearMap iMap;
  unsigned d1;
  for (d1 = 0; d1 < D1; d1++) {
    if (!axis)
      iMap = iMaps[d1];
    for (unsigned d2 = 0; d2 < D2; d2++) {
      if (axis == 1)
	iMap = iMaps[d2];
      for (unsigned d3 = 0; d3 < D3; d3++) {
	if (axis == 2)
	  iMap = iMaps[d3];
	VIO_Real value = get_volume_voxel_value(volume, d1, d2, d3, 0, 0);

	if (value < voxelMin)
	  voxelMin = value;
	if (value > voxelMax)
	  voxelMax = value;

	value = iMap(CONVERT_VOXEL_TO_VALUE(volume, value));

	if (value < mappedRealMin)
	  mappedRealMin = value;
	if (value > mappedRealMax)
	  mappedRealMax = value;
      }
    }
  }

  // Get real extrema
  double realMin = CONVERT_VOXEL_TO_VALUE(volume, voxelMin);
  double realMax = CONVERT_VOXEL_TO_VALUE(volume, voxelMax);

  // Retrieve original voxel->real mapping
  LinearMap voxelToRealMap(voxelMin, voxelMax, CONVERT_VOXEL_TO_VALUE(volume, voxelMin),
			   CONVERT_VOXEL_TO_VALUE(volume, voxelMax));

  // Reset the voxel range to the full range of the data type
  // And set the desired real range
  set_volume_voxel_range(volume, voxelMin, voxelMax);
  set_volume_real_range(volume, mappedRealMin, mappedRealMax);

  // Retrieve range of data type
  //  get_volume_voxel_range(volume, &voxelMin, &voxelMax);

  // Determine final real->voxel map
  LinearMap realToVoxelMap(mappedRealMin, mappedRealMax, voxelMin, voxelMax);

  /*
  cout << "Axis: " << axis << endl
       << "voxelMin: " << voxelMin << " voxelMax: " << voxelMax << endl
       << "Voxel->real map: " << voxelToRealMap << endl
       << "VIO_Real->voxel map: " << realToVoxelMap << endl;
       */

  // Remap volume
  if (verbose)
    cout << "Remapping slices" << flush;
  for (d1 = 0; d1 < D1; d1++) {
    if (!axis) {
      iMap = voxelToRealMap;
      iMap.concat(iMaps[d1]);
      iMap.concat(realToVoxelMap);
    //      iMap = realToVoxelMap(iMaps[d1](voxelToRealMap));
    }
    for (unsigned d2 = 0; d2 < D2; d2++) {
      if (axis == 1) {
	iMap = realToVoxelMap;
	iMap.concat(iMaps[d2]);
	iMap.concat(voxelToRealMap);
	//	iMap = realToVoxelMap(iMaps[d2](voxelToRealMap));
      }
      for (unsigned d3 = 0; d3 < D3; d3++) {
	if (axis == 2) {
	  iMap = realToVoxelMap;
	  iMap.concat(iMaps[d3]);
	  iMap.concat(voxelToRealMap);
	    //	  iMap = realToVoxelMap(iMaps[d3](voxelToRealMap));
	}
	VIO_Real value = get_volume_voxel_value(volume, d1, d2, d3, 0, 0);
	set_volume_voxel_value(volume, d1, d2, d3, 0, 0, 
			       ::clamp(iMap(value), voxelMin, voxelMax));
      }
    }
    if (verbose)
      cout << "." << flush;
  }
  if (verbose)
    cout << "Done" << endl;
}

//
//
//
LinearMap 
determineMap(const FloatArray& modelArray, const FloatArray& dataArray, 
	     const InormalizeArgs& args, int method)
{
  unsigned nVoxels = size(modelArray);
  if (size(dataArray) != nVoxels) {
    cerr << "determineMap: sizes of data- and modelarray are not equal!" << endl;
    exit(EXIT_FAILURE);
  }

  if (method < 0)
    method = args.method;

  LinearMap map;

  if (nVoxels < args.minVoxels)
    return map;

  if (method == InormalizeArgs::RMS) {
    OptData optData(modelArray, dataArray, nVoxels, args.verbose);

    FloatArray temp(modelArray);

    double f = temp.medianVolatile();
    temp = dataArray;
    f /= temp.medianVolatile();

    double fStep = 0.1;

    amoeba_struct amoeba;

    initialize_amoeba(&amoeba, 1, &f, &fStep, evaluateRMS, (void *) &optData, 1e-6);
    
    unsigned iter = 0;
    while ((iter < 100) && perform_amoeba(&amoeba))
      iter++;

    double minRMS = get_amoeba_parameters(&amoeba, &f);

    terminate_amoeba(&amoeba);

    map.factor() = f;
  }

  else if (method == InormalizeArgs::VR) {
    OptData optData(modelArray, dataArray, nVoxels, args.verbose);

    FloatArray temp(modelArray);

    double f = temp.medianVolatile();
    temp = dataArray;
    f /= temp.medianVolatile();

    double fStep = 0.1;

    amoeba_struct amoeba;

    initialize_amoeba(&amoeba, 1, &f, &fStep, evaluateVR, (void *)&optData, 1e-6);
    
    unsigned iter = 0;
    while ((iter < 100) && perform_amoeba(&amoeba))
      iter++;

    double minVR = get_amoeba_parameters(&amoeba, &f);

    terminate_amoeba(&amoeba);

    map.factor() = f;
  }

  else if (method == InormalizeArgs::RATIO_OF_MEANS)
    map.factor() = mean(modelArray)/mean(dataArray);
  
  else if (method == InormalizeArgs::RATIO_OF_MEDIANS) {
    FloatArray temp(modelArray);
    map.factor() = temp.medianVolatile();
    temp = dataArray;
    map.factor() /= temp.medianVolatile();
  }

  else if (method == InormalizeArgs::MEAN_OF_RATIOS) {
    FloatArray temp(modelArray);
    temp /= dataArray;
    prune(temp);

    if (args.constrained) {
      float ceil  = pow(10, 0.1);
      float floor = pow(10, -0.1);
      temp.removeAllNotIn(floor, ceil);
    }
    if (size(temp) >= args.minVoxels) {
      map.factor() = mean(temp);
      if (args.constrained && args.verbose)
	cout << "+";
    }
  }

  else if (method == InormalizeArgs::MEAN_OF_LOG_RATIOS) {
    FloatArray temp(modelArray);
    temp /= dataArray;
    temp = log(temp);
    prune(temp);

    if (args.constrained) {
      float ceil  = 0.1;
      float floor = -0.1;
      temp.removeAllNotIn(floor, ceil);
    }
    if (size(temp) >= args.minVoxels) {
      map.factor() = pow(10, mean(temp));
      if (args.constrained && args.verbose)
	cout << "+";
    }
  }

  else if (method == InormalizeArgs::MEDIAN_OF_RATIOS) {
    FloatArray temp(modelArray);
    temp /= dataArray;
    prune(temp);

    if (args.constrained) {
      float ceil  = pow(10, 0.1);
      float floor = pow(10, -0.1);
      temp.removeAllNotIn(floor, ceil);
    }
    if (size(temp) >= args.minVoxels) {
      map.factor() = temp.medianVolatile();
      if (args.constrained && args.verbose)
	cout << "+";
    }
  }

  else if (method == InormalizeArgs::RANGE) {
    FloatArray dataExtrema(pctExtrema(dataArray, args.rangePct, args.verbose));
    FloatArray modelExtrema(pctExtrema(modelArray, args.rangePct, args.verbose));

    map(dataExtrema[(unsigned int)0], 
        dataExtrema[(unsigned int)1], 
        modelExtrema[(unsigned int)0], 
        modelExtrema[(unsigned int)1]);
  }

  else {
    cerr << "Unknown method " << method << " encountered" << endl;
    exit(EXIT_FAILURE);
  }

  return map;
}

//
//
//
Boolean
saveVolume(const VIO_Volume volume, const Path& path, const Path& mincModel,
	   const MString& history, int compress, int verbose)
{
  char *mincModelString = 0;
  if (!mincModel.isEmpty())
    mincModelString = (char *) mincModel.string();

  char *historyString = 0;
  if (!history.isEmpty())
    historyString = (char *) history.string();

  Path outputPath(path.expanded());
  outputPath.removeCompressedExtension();
  if (output_modified_volume(outputPath, NC_UNSPECIFIED, FALSE, 0.0, 0.0, 
			     volume, mincModelString, historyString,
			     (minc_output_options *) NULL) != VIO_OK) {
    cerr << "Couldn't save " << outputPath << endl;
    return FALSE;
  }

  if (compress) {
    if (verbose) 
      cout << "Compressing " << outputPath << "..." << flush;
    MString zipCommand("gzip -f " + outputPath); 
    system(zipCommand);
    if (verbose)
      cout << "Done" << endl;
  }

  return TRUE;
}

//
//
//
MString&
appendToString(MString& string, const LinearMap& iMap)
{
  string += " (";
  string += iMap.factor();
  string += ", ";
  string += iMap.offset();
  string += ")";
  
  return string;
}

//
//
//
MString&
appendToString(MString& string, const Array<LinearMap>& iMaps)
{
  unsigned n = iMaps.size();
  if (n)
    for (unsigned i = 0; i < n; i++) {
      string += " (";
      string += iMaps[i].factor();
      string += ", ";
      string += iMaps[i].offset();
      string += ")";
    }
  
  return string;
}

double
evaluateRMS(void *data, float *f)
{
  OptData& optData = *(OptData *) data;

  double RMS = 0;
  optData.modelArray.resetIterator();
  optData.dataArray.resetIterator();
  for (unsigned i = optData.n; i; i--) {
    double value = optData.modelArray++ - (optData.dataArray++) * (*f);
    RMS += value*value;
  }

  RMS = sqrt(RMS/optData.n);

  if (optData.verbose)
    cout << *f << ": " << RMS << endl;

  return RMS;
}

double
evaluateVR(void *data, float *f)
{
  OptData& optData = *(OptData *) data;

  double sum  = 0.0;
  double sum2 = 0.0;
  
  optData.modelArray.resetIterator();
  optData.dataArray.resetIterator();
  for (unsigned i = optData.n; i; i--) {
    double value = optData.modelArray++ / ((optData.dataArray++) * (*f));
    sum  += value;
    sum2 += value*value;
  }

  double VR = sum2/optData.n - SQR(sum/optData.n);

  if (optData.verbose)
    cout << *f << ": " << VR << endl;

  return VR;
}

FloatArray
pctExtrema(const FloatArray& array, double pct, int verbose)
{
  FloatArray extrema(2);

  unsigned n = size(array);

  if (!n) {
    cerr << "Warning! attempt to calculate pctExtrema on empty array" << endl;
    return extrema;
  }

  unsigned nKill = unsigned(ceil(pct / 100 * n));

  if (verbose)
    cout << "Calculating " << pct << "," << 100-pct << "% extrema (excluding " 
	 << 2*nKill << " voxels): " << flush;

  if (!nKill)
    array.extrema(&extrema[(unsigned int)0], &extrema[(unsigned int)1]);
  else {
    FloatArray temp(array);
    temp.qsort();
    extrema[(unsigned int)0] = (temp[(unsigned int)nKill] + 
                                temp[(unsigned int)nKill + 1])/2;
    extrema[(unsigned int)1] = (temp[(unsigned int)n - nKill] + 
                                temp[(unsigned int)n - nKill - 1])/2;
  }

  if (verbose)
    cout << extrema << endl;

  return extrema;  
}
