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


typedef itk::LinearInterpolateImageFunction<DistanceImageType> InterpolatorType;

double l2DistToTrainingImage(DistanceImageType::Pointer trainingDistImage, MeshType::Pointer sample);

float specificity(StatisticalModelType::Pointer model, unsigned numberOfShapes) {

    // get the training meshes (that need to be provided in the scores   
    statismo::MatrixType S = model->GetModelInfo().GetScoresMatrix(); 
	vnl_matrix<float> scores = vnl_matrix<float>(S.data(), S.rows(), S.cols());

    typedef std::vector<DistanceImageType::Pointer> TrainingSamplesType;
    TrainingSamplesType trainingSamples;

	// create for each training sample a distance image (to speed up distance computations)
    for (unsigned i  = 0; i < scores.cols(); i++) {
        vnl_vector<float> coeffsForIthSample = scores.get_column(i);
        MeshType::Pointer trainingSample = model->DrawSample(coeffsForIthSample);
		BinaryImageType::Pointer trainingSampleAsImage = meshToBinaryImage(trainingSample);
		DistanceImageType::Pointer trainingSampleAsDistImage = binaryImageToDistanceImage(trainingSampleAsImage);
        trainingSamples.push_back(trainingSampleAsDistImage);

		std::cout << "created distance image " << i << std::endl;
    }

	// draw a number of samples and compute its distance to the closest training dataset
    double accumulatedDistToClosestTrainingShape = 0;
    for (unsigned i = 0; i < numberOfShapes; i++) {
        MeshType::Pointer sample = model->DrawSample();
        double minDist = std::numeric_limits<double>::max();
        for (TrainingSamplesType::const_iterator it = trainingSamples.begin(); it != trainingSamples.end(); ++it) {
            double dist = l2DistToTrainingImage(*it, sample);
            if (dist < minDist) {
                minDist = dist;
            }
        }
		std::cout << "closest distance for sample " << i << ": " << minDist << std::endl;

        accumulatedDistToClosestTrainingShape += minDist;
    }
    double avgDist = accumulatedDistToClosestTrainingShape / numberOfShapes;
	std::cout << "average distance " << avgDist << std::endl;
    return avgDist;
}


// computes the average euclidean distance of the sample (mesh), to the training dataset, using 
// the distanceMap representation of the training dataset to speed up computations.  
double  l2DistToTrainingImage(DistanceImageType::Pointer trainingDistImage, MeshType::Pointer sample) {

    InterpolatorType::Pointer trainingImageInterpolated = InterpolatorType::New();
    trainingImageInterpolated->SetInputImage(trainingDistImage);

    float accumulatedDist = 0;
    for (unsigned i = 0; i < sample->GetNumberOfPoints();  i++) {
		PointType pt = sample->GetPoint(i);
		float distAtPoint = 0;
		if (trainingImageInterpolated->IsInsideBuffer(pt)) {
			distAtPoint = trainingImageInterpolated->Evaluate(pt);
		} else {
			distAtPoint = DMConfigParameters::imageMargin;
		}
        accumulatedDist += distAtPoint * distAtPoint;
    }
    return accumulatedDist / sample->GetNumberOfPoints();
}


