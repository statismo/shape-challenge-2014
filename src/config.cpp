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
const double FittingConfigParameters::maxStepLength = 1;
double const FittingConfigParameters::translationScale = 1;
double const FittingConfigParameters::rotationScale = 0.1;
double const FittingConfigParameters::smScale = 3;
