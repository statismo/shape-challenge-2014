/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/


#include "config.h"

const unsigned ConfigParameters::numSamplesForSpecificityComputations = 200;
const unsigned ConfigParameters::numSamplingPointsSpecificity = 1000;
const unsigned ConfigParameters::numSamplingPointsGeneralization = 50000;

const short FittingConfigParameters::maxNumberOfIterations = 5000;
const double FittingConfigParameters::defaultStepLength = 1;
const double FittingConfigParameters::gradientConvergenceTolerance = 1e-18;
const double FittingConfigParameters::lineSearchAccuracy = 0.1;
double const FittingConfigParameters::translationScale = 1;
double const FittingConfigParameters::rotationScale = 0.01;
double const FittingConfigParameters::scalingScale = 0.001;
double const FittingConfigParameters::smScale = 1;
