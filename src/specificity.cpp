/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/


#include "shapemodelvalidation.h"
#include "config.h"

#include <string>
#include <itkImage.h>
#include "itkMeshFileWriter.h"
#include <itkScaleTransform.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_rank.h>
#include <vnl/vnl_trace.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>



// typedefs and forward declarations of functions that are only used here.
//
typedef vnl_matrix<float> vnlMatrixType;
typedef vnl_vector<float> vnlVectorType;
typedef float ElementType;
typedef std::vector<MeshType::PointType> PointListType;

MeshType::Pointer normalizeScale(MeshType* mesh, MeshType* target);
double computeScalingFactor(MeshType* mesh1, MeshType* mesh2);
vnlVectorType calcMean(const vnlMatrixType& M);
ElementType calcVar(const vnlMatrixType& M, const vnlVectorType& mu) ;
vnlMatrixType calcCov(const vnlMatrixType& X,
                                            const vnlMatrixType& Y,
                                            const vnlVectorType& mu_x,
                                            const vnlVectorType& mu_y);
vnlMatrixType createMatrixFromPoints(const PointListType& points) ;



// Computes specificity of the model by comparing random samples of the model with the test mashes.
float specificity(Logger& logger, StatisticalModelType::Pointer model, const MeshDataList& testMeshes, unsigned numberOfSamples) {


	// draw a number of samples and compute its distance to the closest training dataset
    double accumulatedDistToClosestTrainingShape = 0;
    for (unsigned i = 0; i < numberOfSamples; i++) {
        MeshType::Pointer sample = model->DrawSample();

        double minDist = std::numeric_limits<double>::max();
        for (MeshDataList::const_iterator it = testMeshes.begin(); it != testMeshes.end(); ++it) {
            MeshType::Pointer testMesh = it->first;

            // before we compute the distances between the meshes, we normalize the scale by scaling them
            // to optimally match the mean. This makes sure that models that include scale and those that have them normalized
            // ar etreated the same.
             MeshType::Pointer sampledScaledToMean = normalizeScale(sample, model->DrawMean());
             MeshType::Pointer testScaledToMean = normalizeScale(testMesh, model->DrawMean());

             double dist = computeAverageDistance(testScaledToMean, sampledScaledToMean, ConfigParameters::numSamplingPointsSpecificity);
            logger.Get(logINFO) << "distance " << dist << std::endl;
            if (dist < minDist) {
                minDist = dist;
            }
        }
        logger.Get(logINFO) << "closest distance for sample " << i << ": " << minDist << std::endl;

        accumulatedDistToClosestTrainingShape += minDist;
    }
    double avgDist = accumulatedDistToClosestTrainingShape / numberOfSamples;
    logger.Get(logINFO) << "average distance " << avgDist << std::endl;
    return avgDist;
}


// returns a scaled version of the mesh, that is optimally scaled with respect to the target mesh.
MeshType::Pointer normalizeScale(MeshType* mesh, MeshType* target) {
    double scalingFactor = computeScalingFactor(mesh, target);
    itk::ScaleTransform<double, 3>::Pointer scaleTransform = itk::ScaleTransform<double,3>::New();
    scaleTransform->SetScale(scalingFactor);

    typedef itk::TransformMeshFilter<MeshType, MeshType, itk::Transform<double, 3> > TransformMeshFilterType;
    TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
    transformMeshFilter->SetInput(mesh);
    transformMeshFilter->SetTransform(scaleTransform);
    transformMeshFilter->Update();
    MeshType::Pointer scaledMesh = transformMeshFilter->GetOutput();
    return scaledMesh;

}


// compute optimal pose and scale using a procedure described by Umeyame,
// Least-squares estimation of transformation parameters between two point patterns (IEEE Pami 1991)
double computeScalingFactor(MeshType* mesh1, MeshType* mesh2)
{
    if (mesh1->GetNumberOfPoints() != mesh2->GetNumberOfPoints()) {
        throw std::runtime_error("meshes should be in correspondence when computing specificity");
    }

    PointListType points1;
    PointListType points2;
    for (unsigned i = 0; i < mesh1->GetNumberOfPoints(); i++) {

        MeshType::PointType pt1;
        mesh1->GetPoint(i, &pt1);
        points1.push_back(pt1);

        MeshType::PointType pt2;
        mesh2->GetPoint(i, &pt2);
        points2.push_back(pt2);
    }

    vnlMatrixType X  = createMatrixFromPoints(points1);
    vnlMatrixType Y =  createMatrixFromPoints(points2);

    vnlVectorType mu_x = calcMean(X);
    vnlVectorType mu_y = calcMean(Y);

    ElementType sigma2_x = calcVar(X, mu_x);
    ElementType sigma2_y = calcVar(Y, mu_y);

    vnlMatrixType Sigma_xy = calcCov(X, Y, mu_x, mu_y);

    vnl_svd<ElementType> SVD(Sigma_xy);

    unsigned m = X.rows();

    vnlMatrixType S(m,m);
    S.set_identity();
    ElementType detU = vnl_qr<ElementType>(SVD.U()).determinant();
    ElementType detV = vnl_qr<ElementType>(SVD.V()).determinant();
    if ( detU * detV == -1) {
        S[m-1][m-1] = -1;
    }


    // the procedure actually computes the optimal rotation, translation and scale. We only
    // use the scaling factor
    vnlMatrixType R = SVD.U() * S * SVD.V().transpose();

    ElementType c = 1/sigma2_x * vnl_trace(SVD.W()*S);
    vnlVectorType t = mu_y - c * R * mu_x;
    ElementType epsilon2 = sigma2_y -  pow(vnl_trace(SVD.W()*S), 2) / sigma2_x;

    return c;
}


vnlVectorType calcMean(const vnlMatrixType& M)  {
    vnlVectorType mu(M.rows(), 0.0);

    for (unsigned j = 0; j < M.cols(); j++) {
        mu += M.get_column(j);
    }
    mu /= M.cols();

    return mu;
}

ElementType calcVar(const vnlMatrixType& M, const vnlVectorType& mu)  {
    ElementType sigma2 = 0;
    for (unsigned j = 0; j < M.cols(); j++) {
        sigma2 += pow((M.get_column(j) - mu).two_norm(),2);
    }
    return sigma2/M.cols();
}

vnlMatrixType calcCov(const vnlMatrixType& X,
                                const vnlMatrixType& Y,
                                const vnlVectorType& mu_x,
                                const vnlVectorType& mu_y)
{

    vnlMatrixType Sigma_xy(X.rows(), X.rows(), 0.0);
    for (unsigned j = 0; j < X.cols(); j++) {
        Sigma_xy += outer_product(Y.get_column(j) - mu_y, X.get_column(j) - mu_x);
    }
    Sigma_xy /= X.cols();
    return Sigma_xy;
}


vnlMatrixType createMatrixFromPoints(const PointListType& points) {
    unsigned numberOfPoints = points.size();
    vnlMatrixType M(3, numberOfPoints);

    for (unsigned i = 0; i < numberOfPoints; i++) {
        vnlVectorType p = points[i].GetVnlVector();
        M.set_column(i, p);
    }

    return M;
}


