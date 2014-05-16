/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/


#include "config.h"
#include "shapemodelvalidation.h"

#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMeshFileWriter.h>
#include <itkDirectory.h>
#include <itkBinaryMask3DMeshSource.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkPointsLocator.h>
#include <itkVersion.h>
#include <string>
#include <iostream>
#include <vnl/vnl_random.h>
#include <algorithm>
#include <numeric>

typedef std::list<std::string> FileList;
typedef itk::ImageFileReader<BinaryImageType> ImageReaderType;

FileList getTestImagesInDir (std::string dir)
{
    FileList files;

    itk::Directory::Pointer directory = itk::Directory::New();
    directory->Load(dir.c_str());

    for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
        std::string filename(directory->GetFile(i));
        if (filename.find(".vtk") != std::string::npos || filename.find(".mha") != std::string::npos)
            files.push_back(filename);
    }

    return files;
}


void writeBinaryImage(BinaryImageType* image, const char* filename) {
    typedef itk::ImageFileWriter<BinaryImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(image);
    writer->SetFileName(filename);
    writer->Update();
}

void writeMesh(MeshType* mesh, const char* filename) {
    typedef itk::MeshFileWriter<MeshType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(mesh);
    writer->SetFileName(filename);
    writer->Update();
}




BinaryImageType::Pointer readBinaryImage(const std::string& filename) {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // we make sure that the image really has value 0 and 1
    typedef itk::BinaryThresholdImageFilter<BinaryImageType, BinaryImageType> ThresholdFilterType;
    ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetLowerThreshold(1);
    thresholdFilter->SetUpperThreshold(255);
    thresholdFilter->SetInput(reader->GetOutput());
    thresholdFilter->Update();
    BinaryImageType::Pointer img = thresholdFilter->GetOutput();
    img->DisconnectPipeline();
    return img;
}


BinaryImageType::Pointer meshToBinaryImage(MeshType::Pointer mesh, unsigned imageResolution, double imageMargin) {
typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> TriangleMeshToBinaryImageFilterType;

    // We transform the mesh to a binary image
    TriangleMeshToBinaryImageFilterType::Pointer meshToBinImageFilter = TriangleMeshToBinaryImageFilterType::New();


    // we choose the image slightly larger than the bounding box of the mesh
    const MeshType::BoundingBoxType* boundingBox = mesh->GetBoundingBox();
    PointType minPt = boundingBox->GetMinimum();
    for (unsigned d = 0; d < 3; d++) {minPt.SetElement(d, minPt.GetElement(d) - imageMargin); }
    PointType maxPt = boundingBox->GetMaximum();
    for (unsigned d = 0; d < 3; d++) {maxPt.SetElement(d, maxPt.GetElement(d) + imageMargin); }

    meshToBinImageFilter->SetOrigin(minPt);
    TriangleMeshToBinaryImageFilterType::SizeType size;
    for (unsigned d = 0; d < 3; d++) {
    size[d] = imageResolution;
    }

    meshToBinImageFilter->SetSize(size);
    TriangleMeshToBinaryImageFilterType::SpacingType spacing;
    for (unsigned d = 0; d < 3; d++) { spacing[d] = (maxPt[d] - minPt[d]) / size[d]; }
   meshToBinImageFilter->SetSpacing(spacing);

    meshToBinImageFilter->SetInput(mesh);
    meshToBinImageFilter->Update();
    return meshToBinImageFilter->GetOutput();
}



DistanceImageType::Pointer binaryImageToDistanceImage(BinaryImageType::Pointer binaryImage) { 

	typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, DistanceImageType> DistanceImageFilterType;
	
    // The binary image is now converted to a distance map
    DistanceImageFilterType::Pointer dm = DistanceImageFilterType::New();
    dm->SetInput(binaryImage);
    dm->Update();
    DistanceImageType::Pointer distanceImage = dm->GetOutput();
	return distanceImage;
}


MeshType::Pointer binaryImageToMesh(BinaryImageType* testImage) {

    typedef itk::BinaryMask3DMeshSource< BinaryImageType, MeshType > MeshSourceType;
    MeshSourceType::Pointer meshSource = MeshSourceType::New();


    meshSource->SetObjectValue( 1 );
    meshSource->SetInput( testImage );
    try
    {
        meshSource->Update();
    }
    catch( itk::ExceptionObject & exp )
    {
        throw std::runtime_error("Could not extract surface from binary image.");
    }
    return meshSource->GetOutput();
}

double computeEuclideanPointDist(MeshType::PointType pt1, MeshType::PointType pt2) {
    double sumSquares = 0;
    for (unsigned i = 0; i < Dimensions; i++)  {
        double dist = pt1[i] - pt2[i];
        sumSquares += dist * dist;
    }
    return std::sqrt(sumSquares);
}


std::vector<double> computeDistances(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints) {

    #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
    typedef itk::PointsLocator< MeshType::PointsContainer > PointsLocatorType;
    #else
    typedef itk::PointsLocator<int, 3, double, MeshType::PointsContainer > PointsLocatorType;
    #endif

    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(mesh2->GetPoints());
    ptLocator->Initialize();


    // we assume that the mesh points are approximately uniformely distributed.

    std::vector<double> distanceValues;

    double totalDist = 0;
    vnl_random randGen;

    for (unsigned i = 0; i < numberOfSamplingPoints; i++) {
        unsigned ptId = randGen.lrand32(mesh1->GetNumberOfPoints() - 1);

        MeshType::PointType sourcePt = mesh1->GetPoint(ptId);
        int closestPointId = ptLocator->FindClosestPoint(sourcePt);
        MeshType::PointType targetPt = mesh2->GetPoint(closestPointId);
        distanceValues.push_back(computeEuclideanPointDist(sourcePt, targetPt));
    }
    return distanceValues;
}


double computeAverageDistance(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints) {
    std::vector<double> distanceValues = computeDistances(mesh1, mesh2, numberOfSamplingPoints);
    return std::accumulate(distanceValues.begin(), distanceValues.end(), 0.0) / numberOfSamplingPoints;
}

double computeSymmetricAverageDistance(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints) {
    double dist1 = computeAverageDistance(mesh1, mesh2, numberOfSamplingPoints);
    double dist2 = computeAverageDistance(mesh2, mesh1, numberOfSamplingPoints);
    return (dist1 + dist2) / 2;
}


double computeHausdorffDistance(MeshType::Pointer mesh1, MeshType::Pointer mesh2, unsigned numberOfSamplingPoints) {

    std::vector<double> distanceValues1 = computeDistances(mesh1, mesh2, numberOfSamplingPoints);
    std::vector<double> distanceValues2 = computeDistances(mesh2, mesh1, numberOfSamplingPoints);

    std::vector<double>::iterator maxElemIt1 = std::max_element(distanceValues1.begin(), distanceValues1.end());
    std::vector<double>::iterator maxElemIt2 = std::max_element(distanceValues2.begin(), distanceValues2.end());

    return  std::max(*maxElemIt1, *maxElemIt2);
}


