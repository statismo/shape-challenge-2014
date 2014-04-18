/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/

#include "shapemodelvalidation.h"

#include <string>
#include <vector>
#include <algorithm>

// the simplest compactness measure is simply the number of parameters of the model .

float compactness(StatisticalModelType::Pointer model) {
     return model->GetNumberOfPrincipalComponents();
}
