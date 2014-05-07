/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/



#include "shapemodelvalidation.h"
#include <iostream>
#include <string>
#include <list>
#include "logger.h"


int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cout << "usage: shapeModelValidation shapemodel testdatadir logfile" << std::endl;
		return -1;
	}

    char* modelFn = argv[1];
    char* testdatadir = argv[2];
    char* logfile = argv[3];

	try { 

		Logger logger(logfile, logDEBUG);
	
		FileList testfiles = getTestImagesInDir(testdatadir);
		TestImageList testImages;
		
		for (FileList::const_iterator it = testfiles.begin(); it != testfiles.end(); ++it) { 
            std::string filename = std::string(testdatadir) +"/" + *it;
            logger.Get(logINFO) << "reading image " << filename << std::endl;
            BinaryImageType::Pointer testImage = readBinaryImage(filename);
			testImages.push_back(testImage);
		}

		logger.Get(logINFO) << "Reading statistical model " << modelFn << std::endl;

		RepresenterType::Pointer representer = RepresenterType::New();
		StatisticalModelType::Pointer model = StatisticalModelType::New();
		model->Load(representer, modelFn);

        GeneralizationResult generalizationScore = generalization(logger, model, testImages);
		logger.Get(logINFO) << "generalizationScore: avg = " << generalizationScore.averageDistance << " hd = " << generalizationScore.hausdorffDistance << std::endl;

        float specificityValue = specificity(logger, model, 1000);
		logger.Get(logINFO) << "specificity value: " << specificityValue << std::endl;
    
        float compactnessScore = compactness(logger, model);
		logger.Get(logINFO) << "compactness score: " << compactnessScore  << std::endl;


	} catch (std::exception& e) { 
		std::cerr << "An error occurred : " << e.what() << std::endl;
		return -1;
	}

    return 0;
}

