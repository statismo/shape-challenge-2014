/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/


#include "config.h"
#include "shapemodelvalidation.h"

#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkImageFileReader.h>
#include <itkDirectory.h>

#include <string>
#include <iostream>

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


BinaryImageType::Pointer readBinaryImage(const std::string& filename) {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    BinaryImageType::Pointer img = reader->GetOutput();
    img->DisconnectPipeline();
    return img;
}


BinaryImageType::Pointer meshToBinaryImage(MeshType::Pointer mesh, BinaryImageType::Pointer infoImage) { 
	typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> TriangleMeshToBinaryImageFilterType;

	// We transform the mesh to a binary image
    TriangleMeshToBinaryImageFilterType::Pointer meshToBinImageFilter = TriangleMeshToBinaryImageFilterType::New();
	if (infoImage.IsNull()) { 

		// we choose the image slightly larger than the bounding box of the mesh
		const MeshType::BoundingBoxType* boundingBox = mesh->GetBoundingBox();
		PointType minPt = boundingBox->GetMinimum();
		for (unsigned d = 0; d < 3; d++) {minPt.SetElement(d, minPt.GetElement(d) - DMConfigParameters::imageMargin); }
		PointType maxPt = boundingBox->GetMaximum();
		for (unsigned d = 0; d < 3; d++) {maxPt.SetElement(d, maxPt.GetElement(d) + DMConfigParameters::imageMargin); }

		meshToBinImageFilter->SetOrigin(minPt);
		TriangleMeshToBinaryImageFilterType::SizeType size;
		for (unsigned d = 0; d < 3; d++) {
			size[d] = DMConfigParameters::DistanceImageResolution; 
		}
	   
		meshToBinImageFilter->SetSize(size);
		TriangleMeshToBinaryImageFilterType::SpacingType spacing;
		for (unsigned d = 0; d < 3; d++) { spacing[d] = (maxPt[d] - minPt[d]) / size[d]; }
		meshToBinImageFilter->SetSpacing(spacing);
	}
	else {
		meshToBinImageFilter->SetInfoImage(infoImage);
	}

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
