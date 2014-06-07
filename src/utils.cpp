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
#include <itkTransformMeshFilter.h>
#include <itkIdentityTransform.h>

typedef std::list<std::string> FileList;
typedef itk::ImageFileReader<BinaryImageType> ImageReaderType;



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


ImageDataList getImagesInDir (Logger& logger, std::string dirname)
{
    ImageDataList images;

    itk::Directory::Pointer directory = itk::Directory::New();
    directory->Load(dirname.c_str());

    for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
        std::string filename(directory->GetFile(i));
        if (filename.find(".vtk") != std::string::npos || filename.find(".mha") != std::string::npos) {
            std::string fullpath = std::string(dirname) +"/" + filename;
            logger.Get(logINFO) << "reading image " << fullpath << std::endl;
            BinaryImageType::Pointer image = readBinaryImage(fullpath);
            images.push_back(std::make_pair(image, fullpath));
        }
    }

    return images;
}





MeshType::Pointer cloneMesh(const MeshType* mesh) {

    // cloning is cumbersome - therefore we let itk do the job for, and use perform a
    // Mesh transform using the identity transform. This should result in a perfect clone.

    typedef itk::IdentityTransform<MeshType::PixelType, 3> IdentityTransformType;
    typedef itk::TransformMeshFilter<MeshType, MeshType, IdentityTransformType> TransformMeshFilterType;

    typename TransformMeshFilterType::Pointer tf = TransformMeshFilterType::New();
    tf->SetInput(mesh);
    typename IdentityTransformType::Pointer idTrans = IdentityTransformType::New();
    tf->SetTransform(idTrans);
    tf->Update();

    typename MeshType::Pointer clone = tf->GetOutput();
    clone->DisconnectPipeline();
    return clone;
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

