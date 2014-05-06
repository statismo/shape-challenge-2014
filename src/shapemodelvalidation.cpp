/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/



#include "shapemodelvalidation.h"

#include "itkDirectory.h"
#include "itkImageFileReader.h"
#include <iostream>
#include <string>
#include <list>

typedef std::list<std::string> FileList;
typedef itk::ImageFileReader<BinaryImageType> ImageReaderType;

FileList getTestImagesInDir (std::string dir)
{
	FileList files;

	itk::Directory::Pointer directory = itk::Directory::New();
	directory->Load(dir.c_str());

	for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
		std::string filename(directory->GetFile(i));
		if (filename.find(".vtk") != std::string::npos || filename.find(".mhd") != std::string::npos)
            files.push_back(filename);
	}

    return files;
}


void writeResult(const char* resultfile, unsigned objId, const char* modelFn, const GeneralizationResult& gndummy, float specificity, float compactness) { 
	std::ofstream resFile(resultfile);

	// TODO error handling
	if (!resFile) { 
		std::ostringstream msgos;
		msgos << "cannot create output file " << resultfile;
		throw std::exception(msgos.str().c_str());	
	}

	resFile << "{" << std::endl;
	resFile << "\"object-id\" : " << objId << "," << std::endl;
	resFile << "\"shape-model\" : " << "\"" << modelFn << "\"" << "," << std::endl;
	resFile << "\"error-code\" : " << 0 << "," << std::endl;
	resFile << "\"error-message\" : " << "\"\"" << "," << std::endl;
	resFile << "\"generalization\" : {" << std::endl;
	resFile << "\t\"average-distance\" : " << gndummy.averageDistance << "," << std::endl;
	resFile << "\t\"hausdorff-distance\" : " << gndummy.hausdorffDistance <<  std::endl;
	resFile << "}," << std::endl;
	resFile << "\"specificity\" : " << specificity << "," << std::endl;
	resFile << "\"compactness\" : " << compactness << std::endl;
	resFile << "}" << std::endl;
	resFile.close();
}



int main(int argc, char* argv[]) {

	if (argc < 6) { 
		std::cout << "usage: shapeModelValidation obj-id shapemodel testdatadir resultfile logfile" << std::endl;
		return -1;
	}

	unsigned int objId = static_cast<unsigned>(atoi(argv[1]));
	char* modelFn = argv[2];
	char* testdatadir = argv[3];
	char* resultfile = argv[4]; 
	char* logfile = argv[5];


	GeneralizationResult gdummy(0.1, 1.2);
	writeResult(resultfile, objId, modelFn, gdummy, 0.1, 0.1);
	// 
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// TODO remove me
	/////////////////////////////////
	exit(0);


	FileList testfiles = getTestImagesInDir(testdatadir);
	TestImageList testImages;

	for (FileList::const_iterator it = testfiles.begin(); it != testfiles.end(); ++it) { 
		std::cout << "reading image "<< *it << std::endl;

		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName(std::string(testdatadir) +"/" + *it);
		reader->Update();
		BinaryImageType::Pointer testImage = reader->GetOutput();
		testImages.push_back(testImage);
	}


    RepresenterType::Pointer representer = RepresenterType::New();
    StatisticalModelType::Pointer model = StatisticalModelType::New();
    model->Load(representer, modelFn);

	GeneralizationResult generalizationScore = generalization(model, testImages);
	std::cout << "generalizationScore: avg = " << generalizationScore.averageDistance << " hd = " << generalizationScore.hausdorffDistance << std::endl;

    float specificityValue = specificity(model, 1000);
    std::cout << "specificity value: " << specificityValue << std::endl;
    
	float compactnessScore = compactness(model);
    std::cout << "compactness score: " << compactnessScore  << std::endl;



    return 0;
}

