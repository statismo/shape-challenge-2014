#include "config.h"
#include "shapemodelvalidation.h"

#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkPointsLocator.h>
#include <itkVersion.h>
#include <string>
#include <iostream>
#include <vnl/vnl_random.h>
#include <algorithm>
#include <numeric>

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



