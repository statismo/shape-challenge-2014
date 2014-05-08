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

float specificity(Logger& logger, StatisticalModelType::Pointer model, unsigned numberOfShapes) {

    // get the training meshes (that need to be provided in the scores   
    statismo::MatrixType S = model->GetModelInfo().GetScoresMatrix(); 
	vnl_matrix<float> scores = vnl_matrix<float>(S.data(), S.rows(), S.cols());

    if (scores.cols() == 0 || scores.rows() != model->GetNumberOfPrincipalComponents()) {
        throw std::runtime_error("An invalid scores matrix was provided");
    }

    typedef std::vector<MeshType::Pointer> TrainingSamplesType;
    TrainingSamplesType trainingSamples;

	// create for each training sample a distance image (to speed up distance computations)
    for (unsigned i  = 0; i < scores.cols(); i++) {
        vnl_vector<float> coeffsForIthSample = scores.get_column(i);
        MeshType::Pointer trainingSample = model->DrawSample(coeffsForIthSample);
        trainingSamples.push_back(trainingSample);
    }

	// draw a number of samples and compute its distance to the closest training dataset
    double accumulatedDistToClosestTrainingShape = 0;
    for (unsigned i = 0; i < numberOfShapes; i++) {
        MeshType::Pointer sample = model->DrawSample();

        double minDist = std::numeric_limits<double>::max();
        for (TrainingSamplesType::const_iterator it = trainingSamples.begin(); it != trainingSamples.end(); ++it) {
            MeshType::Pointer x = sample;
            double dist = computeAverageDistance(*it, sample, ConfigParameters::numSamplingPointsSpecificity);
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



