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


// Common definitions
const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions > MeshType;
typedef MeshType::PointType PointType;	
typedef itk::StandardMeshRepresenter<float, Dimensions> RepresenterType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::Image<float, 3> DistanceImageType;
typedef std::list<BinaryImageType::Pointer> TestImageList;


// holds the result for the two generalization metrics
struct GeneralizationResult { 
	GeneralizationResult(double a, double h) : averageDistance(a), hausdorffDistance(h) {}

	double averageDistance;
	double hausdorffDistance;
};

// the main validation functions
GeneralizationResult generalization(StatisticalModelType::Pointer model, const TestImageList& testImages);
float specificity(StatisticalModelType::Pointer model, unsigned numberOfShapes);
float compactness(StatisticalModelType::Pointer model);

// some commonly used helper functions
double computeAverageDistance(BinaryImageType::Pointer image1, BinaryImageType::Pointer image2);
double computeHausdorffDistance(BinaryImageType::Pointer image1, BinaryImageType::Pointer image2);
BinaryImageType::Pointer meshToBinaryImage(MeshType::Pointer mesh, BinaryImageType::Pointer infoImage = 0);
DistanceImageType::Pointer binaryImageToDistanceImage(BinaryImageType::Pointer binaryImage);


#endif 