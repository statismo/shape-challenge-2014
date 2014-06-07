/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/


#include "shapemodelvalidation.h"
#include "config.h"

#include <string>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImage.h>
#include "itkMeshFileWriter.h"

typedef itk::LinearInterpolateImageFunction<DistanceImageType> InterpolatorType;

double l2DistToTrainingImage(DistanceImageType::Pointer trainingDistImage, MeshType::Pointer sample);

float specificity(Logger& logger, StatisticalModelType::Pointer model, const MeshDataList& testMeshes, unsigned numberOfShapes) {


	// draw a number of samples and compute its distance to the closest training dataset
    double accumulatedDistToClosestTrainingShape = 0;
    for (unsigned i = 0; i < numberOfShapes; i++) {
        MeshType::Pointer sample = model->DrawSample();

        double minDist = std::numeric_limits<double>::max();
        for (MeshDataList::const_iterator it = testMeshes.begin(); it != testMeshes.end(); ++it) {
            MeshType::Pointer testMesh = it->first;
            double dist = computeAverageDistance(testMesh, sample, ConfigParameters::numSamplingPointsSpecificity);
            logger.Get(logINFO) << "distance " << dist << std::endl;
            if (dist < minDist) {
                minDist = dist;
            }
        }
        logger.Get(logINFO) << "closest distance for sample " << i << ": " << minDist << std::endl;

        accumulatedDistToClosestTrainingShape += minDist;
    }
    double avgDist = accumulatedDistToClosestTrainingShape / numberOfShapes;
    logger.Get(logINFO) << "average distance " << avgDist << std::endl;
    return avgDist;
}



