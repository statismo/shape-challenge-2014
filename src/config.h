/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/


#ifndef __CONFIG_H
#define __CONFIG_H


struct ConfigParameters {
 static const unsigned numSamplesForSpecificityComputations;
 static const unsigned numSamplingPointsSpecificity;
 static const unsigned numSamplingPointsGeneralization;
};



class FittingConfigParameters {
public:
// Some configuration parameters
static const short maxNumberOfIterations; // the maximum number of iterations to use in the optimization
static const double maxStepLength; // the maximum number of iterations to use in the optimization
static const double translationScale; // dynamic range of translations
static const double rotationScale; // dynamic range of rotations
static const double scalingScale; // dynamic range of rotations
static const double smScale; // dynamic range of statistical model parameters
};


#endif
