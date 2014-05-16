/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/

#ifndef __SHAPE_MODEL_VALIDATION_H
#define __SHAPE_MODEL_VALIDATION_H

#include "statismo_ITK/itkStatisticalModel.h"
#include "Representers/ITK/itkStandardMeshRepresenter.h"
#include "itkImage.h"
#include "itkMesh.h"
#include "logger.h"


// Common definitions
const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions > MeshType;
typedef MeshType::PointType PointType;	
typedef itk::StandardMeshRepresenter<float, Dimensions> RepresenterType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::Image<float, 3> DistanceImageType;
typedef std::list<BinaryImageType::Pointer> TestImageList;
typedef std::list<std::string> FileList;

// holds the result for the two generalization metrics
struct GeneralizationResult { 
	GeneralizationResult(double a, double h) : averageDistance(a), hausdorffDistance(h) {}

	double averageDistance;
	double hausdorffDistance;
};

// the main validation functions
GeneralizationResult generalization(Logger& logger, StatisticalModelType::Pointer model, const TestImageList& testImages);
float specificity(Logger& logger, StatisticalModelType::Pointer model, unsigned numberOfShapes);
float compactness(Logger& logger, StatisticalModelType::Pointer model);

// some commonly used helper functions
double computeAverageDistance(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints);
double computeSymmetricAverageDistance(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints);
double computeHausdorffDistance(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints);
BinaryImageType::Pointer meshToBinaryImage(MeshType::Pointer mesh, unsigned imageResolution, double imageMargin);
MeshType::Pointer binaryImageToMesh(BinaryImageType* testImage);
DistanceImageType::Pointer binaryImageToDistanceImage(BinaryImageType::Pointer binaryImage);
FileList getTestImagesInDir (std::string dir);
BinaryImageType::Pointer readBinaryImage(const std::string& filename);

void writeBinaryImage(BinaryImageType* image, const char* filename);
void writeMesh(MeshType* mesh, const char* filename);

#endif 
