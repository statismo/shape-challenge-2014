/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/



#include "shapemodelvalidation.h"
#include "config.h"


#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkMeanSquaresPointSetToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkLBFGSOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkCommand.h"
#include "itkMesh.h"
#include "itkPointSet.h"
#include "itkTransformMeshFilter.h"
#include "itkCompositeTransform.h"
#include "itkRigid3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCompositeTransform.h"
#include "itkCenteredVersorTransformInitializer.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"




GeneralizationResult generalization(Logger& logger, StatisticalModelType::Pointer model, const MeshDataList& testMeshes) {

	double totalAvgDistance = 0;
	double totalHdDistance = 0;
    unsigned numTestImages = testMeshes.size();

    for (MeshDataList::const_iterator it = testMeshes.begin(); it != testMeshes.end(); ++it) {
        std::string filename = it->second;
        logger.Get(logINFO) << "computing best reconstruction " << filename << std::endl;
        MeshType::Pointer testMesh = it->first;

        // project test mesh into model
        MeshType::Pointer bestReconstruction = model->DrawSample(model->ComputeCoefficientsForSample(testMesh));
        std::ostringstream pos;
        unsigned id = std::distance(it, testMeshes.begin());

#ifdef DEBUG
        pos << "/tmp/projected/" << "projected-" << id << ".vtk";
        writeMesh(testMesh, pos.str().c_str());
        std::ostringstream fos;
        fos << "/tmp/fitted/" << "fitted-" << id << ".vtk";
        writeMesh(fittedMesh, fos.str().c_str());
#endif
        double avgDistance = computeSymmetricAverageDistance(testMesh,  bestReconstruction, ConfigParameters::numSamplingPointsGeneralization);
        totalAvgDistance += avgDistance;
        double hdDistance = computeHausdorffDistance(testMesh, bestReconstruction, ConfigParameters::numSamplingPointsGeneralization);
        totalHdDistance += hdDistance;

        logger.Get(logINFO) << "avg distance " << avgDistance << std::endl;
        logger.Get(logINFO) << "hd distance " << hdDistance << std::endl;
	}
	return GeneralizationResult(totalAvgDistance / numTestImages, totalHdDistance / numTestImages);
}



