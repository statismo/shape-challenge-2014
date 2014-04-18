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


int main(int argc, char* argv[]) {

	if (argc < 3) { 
		std::cout << "usage: shapeModelValidation shapemodel testdatadir" << std::endl;
		return -1;
	}

	char* modelFn = argv[1];
	char* testdatadir = argv[2];

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

