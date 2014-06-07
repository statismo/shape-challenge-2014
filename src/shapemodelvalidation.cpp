/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/



#include "shapemodelvalidation.h"
#include "config.h"
#include <iostream>
#include <string>
#include <list>
#include "logger.h"




int main(int argc, char* argv[]) {

    if (argc < 5) {
        std::cout << "usage: shapeModelValidation shapemodel trainingdatadir testdatadir logfile" << std::endl;
		return -1;
	}

    char* modelFn = argv[1];
    char* trainingdatadir = argv[2];
    char* testdatadir = argv[3];
    char* logfile = argv[4];

	try { 

		Logger logger(logfile, logDEBUG);


        logger.Get(logINFO) << "Reading statistical model " << modelFn << std::endl;

        RepresenterType::Pointer representer = RepresenterType::New();
        StatisticalModelType::Pointer model = StatisticalModelType::New();
        model->Load(representer, modelFn);

        ImageDataList testImages = getImagesInDir(logger, testdatadir);
        MeshDataList testMeshes = establishCorrespondenceAndAlignImages(logger, model, testImages);
        ImageDataList trainingImages = getImagesInDir(logger, trainingdatadir);
        MeshDataList trainingMeshes = establishCorrespondenceAndAlignImages(logger, model, trainingImages);

        GeneralizationResult generalizationScore = generalization(logger, model, testMeshes);
		logger.Get(logINFO) << "generalizationScore: avg = " << generalizationScore.averageDistance << " hd = " << generalizationScore.hausdorffDistance << std::endl;

        float specificityValue = specificity(logger, model, trainingMeshes, ConfigParameters::numSamplesForSpecificityComputations);
		logger.Get(logINFO) << "specificity value: " << specificityValue << std::endl;
    
        float compactnessScore = compactness(logger, model);
		logger.Get(logINFO) << "compactness score: " << compactnessScore  << std::endl;


	} catch (std::exception& e) { 
		std::cerr << "An error occurred : " << e.what() << std::endl;
		return -1;
	}

    return 0;
}

