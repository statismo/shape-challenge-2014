/*
* Created by Marcel Luethi
* Copyright (c) 2011 University of Basel
*
* Licensed under the BSD, 3-Clause license
*
**/

#include "shapemodelvalidation.h"
#include "config.h"


#include "itkStatisticalShapeModelTransform.h"
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
#include "itkSimilarity3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCompositeTransform.h"
#include "itkCenteredVersorTransformInitializer.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <itkPointsLocator.h>


typedef itk::Similarity3DTransform<double> SimilarityTransformType;
typedef std::pair<MeshType::Pointer, SimilarityTransformType::Pointer> FittingResultType;


FittingResultType fitModelToTestImage(Logger& logger, StatisticalModelType::Pointer model, BinaryImageType::Pointer testImage);


/*
 * Returns a new mesh with each vertex of the mesh projected to its closet point on the target
 */
MeshType::Pointer projectMesh(MeshType::Pointer mesh, MeshType::Pointer target) {
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
    typedef itk::PointsLocator< MeshType::PointsContainer > PointsLocatorType;
#else
    typedef itk::PointsLocator<int, 3, double, MeshType::PointsContainer > PointsLocatorType;
#endif

    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(target->GetPoints());
    ptLocator->Initialize();

    MeshType::Pointer projectedMesh = cloneMesh(mesh);
    for (unsigned ptId = 0; ptId < mesh->GetNumberOfPoints(); ptId++) {

        MeshType::PointType sourcePt = mesh->GetPoint(ptId);
        int closestPointId = ptLocator->FindClosestPoint(sourcePt);
        MeshType::PointType targetPt = target->GetPoint(closestPointId);
        projectedMesh->SetPoint(ptId, targetPt);
    }
    return projectedMesh;
}

/*
 * Returns a mesh representation of the shape encoded in the test image, which is in correspondence with the model. Furthermore, the
 * pose and scale of the mesh are adjusted to the model
 */
MeshType::Pointer establishCorrespondenceAndAlignData(Logger& logger, StatisticalModelType::Pointer model, BinaryImageType::Pointer testImage) {

    FittingResultType fittedMeshAndTransform = fitModelToTestImage(logger, model, testImage);

    MeshType::Pointer testMesh = binaryImageToMesh(testImage);

    // we align the mesh to the model (i.e. correct pose and scale)
    typedef itk::TransformMeshFilter<MeshType, MeshType, itk::Transform<double, 3> > TransformMeshFilterType;
    TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
    transformMeshFilter->SetInput(testMesh);   
    itk::Transform<double, 3>::Pointer inverseTransform = fittedMeshAndTransform.second->GetInverseTransform();
    transformMeshFilter->SetTransform(inverseTransform);
    transformMeshFilter->Update();
    MeshType::Pointer alignedTestMesh = transformMeshFilter->GetOutput();

    MeshType::Pointer fittedMesh = fittedMeshAndTransform.first;
    MeshType::Pointer testMeshInCorrespondence = projectMesh(fittedMesh, alignedTestMesh);
    return testMeshInCorrespondence;

}



/*
 * Establishes correspondence for each image in a list (see establishCorrespondenceAndAlignImage for details)
 */
MeshDataList establishCorrespondenceAndAlignData(Logger& logger, StatisticalModelType::Pointer model, const ImageDataList& images) {
    MeshDataList alignedMeshes;

    for (ImageDataList::const_iterator it = images.begin(); it != images.end(); ++it) {
        BinaryImageType::Pointer image  = it->first;
        MeshType::Pointer alignedMesh = establishCorrespondenceAndAlignData(logger, model, image);
        std::string filename  = it->second;
        alignedMeshes.push_back(std::make_pair(alignedMesh, filename));
    }
    return alignedMeshes;
}


// Helper class needed by ITK to be able to log the status of the fittign in each iteration
template<class OptimizerType>
class IterationStatusObserver : public itk::Command
{
public:
  typedef IterationStatusObserver Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  itkNewMacro( Self );


   void SetLogger(Logger& logger) { m_logger = &logger; }

   void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const OptimizerType* optimizer =  dynamic_cast< const OptimizerType* >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
    {
      return;
    }

    if (m_logger != 0) {
       m_logger->Get(logINFO) << "Iteration: " << ++m_iter_no << "  value " <<  optimizer->GetValue() << std::endl; //<< " smodel parameters " << optimizer->GetCurrentPosition() << std::endl;
    }
  }


protected:
  IterationStatusObserver():
     m_iter_no(0), m_logger(0) {};

  virtual ~IterationStatusObserver(){};

private:
  int m_iter_no;
  Logger* m_logger;

};


// Fits a model to the given test image.
FittingResultType fitModelToTestImage(Logger& logger, StatisticalModelType::Pointer model, BinaryImageType::Pointer testImage)
{
    const unsigned Dimensions = 3;
    typedef itk::PointSet<float, Dimensions > PointSetType;
    typedef itk::Mesh<float, Dimensions > MeshType;
    typedef itk::Point<double, 3> PointType;
    typedef itk::MeanSquaresPointSetToImageMetric<PointSetType, DistanceImageType> MetricType;
    typedef itk::PointSetToImageRegistrationMethod<PointSetType, DistanceImageType> RegistrationFilterType;
    typedef itk::LBFGSOptimizer OptimizerType;
    //typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::LinearInterpolateImageFunction<DistanceImageType, double> InterpolatorType;
    typedef itk::CompositeTransform<double, 3> CompositeTransformType;
    typedef itk::CenteredVersorTransformInitializer<BinaryImageType, BinaryImageType> TransformInitializerType;
    typedef itk::StatisticalShapeModelTransform<MeshType, double, Dimensions> StatisticalModelTransformType;


    double meanImageMargin = 40.0;
    unsigned meanImageResolution = 256;

    // we compute a binary image of the model mean to use it for initialization of the shape
    BinaryImageType::Pointer refAsBinImage = meshToBinaryImage(model->GetRepresenter()->GetReference(), meanImageResolution, meanImageMargin);
    SimilarityTransformType::Pointer similarityTransform = SimilarityTransformType::New();
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetFixedImage(refAsBinImage);
    initializer->SetMovingImage(testImage);
    initializer->SetTransform(similarityTransform);
    initializer->ComputeRotationOff();
    //initializer->MomentsOn();
    initializer->InitializeTransform();
    similarityTransform->SetScale(1);


    // optimize only similarity transform in a first step

    CompositeTransformType::Pointer transform = CompositeTransformType::New();
    transform->AddTransform(similarityTransform);
    transform->SetAllTransformsToOptimizeOn(); // optimize shape and pose parameters

    // Setting up the fitting
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetDefaultStepLength(FittingConfigParameters::defaultStepLength);
    //optimizer->SetMaximumStepLength(FittingConfigParameters::maxStepLength);
    optimizer->MinimizeOn();
    //optimizer->SetNumberOfIterations(FittingConfigParameters::maxNumberOfIterations);
    optimizer->SetMaximumNumberOfFunctionEvaluations(FittingConfigParameters::maxNumberOfIterations);
    optimizer->SetGradientConvergenceTolerance(FittingConfigParameters::gradientConvergenceTolerance);
    optimizer->SetLineSearchAccuracy(FittingConfigParameters::lineSearchAccuracy);
    // set the scales of the optimizer, to compensate for potentially different scales of translation, rotation and shape parameters
    OptimizerType::ScalesType scales( similarityTransform->GetNumberOfParameters() );
    for (unsigned i = 0; i < 3; i++) {
        scales[i] = 1.0 / (FittingConfigParameters::rotationScale);
    }

    for (unsigned i =  3; i <  6; i++) {
        scales[i] = 1.0 / (FittingConfigParameters::translationScale);
    }
    scales[6] =  1.0 / (FittingConfigParameters::scalingScale);
    optimizer->SetScales(scales);


    // set up the observer to keep track of the progress
    typedef IterationStatusObserver<OptimizerType> ObserverType;
    ObserverType::Pointer observer = ObserverType::New();
    observer->SetLogger(logger);
    optimizer->AddObserver( itk::IterationEvent(), observer );

    // set up the metric and interpolators
    MetricType::Pointer metric = MetricType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    // connect all the components
    RegistrationFilterType::Pointer registration = RegistrationFilterType::New();
    registration->SetInitialTransformParameters(transform->GetParameters());
    registration->SetInterpolator(interpolator);
    registration->SetMetric(metric);
    registration->SetOptimizer( optimizer);
    registration->SetTransform( transform );





    // we create the fixedPointSet of the registration from the reference mesh of our model.
    // As we are fitting to the 0 level set of a distance image, we set the value of the pointdata to 0.
    PointSetType::Pointer fixedPointSet = PointSetType::New();
    fixedPointSet->SetPoints(model->GetRepresenter()->GetReference()->GetPoints());
    PointSetType::PointDataContainer::Pointer points = PointSetType::PointDataContainer::New();
    points->Reserve(model->GetRepresenter()->GetReference()->GetNumberOfPoints());
    for (PointSetType::PointDataContainer::Iterator it = points->Begin(); it != points->End(); ++it) {
        it->Value() = 0;
    }
    fixedPointSet->SetPointData(points);


    DistanceImageType::Pointer targetDistanceImage = binaryImageToDistanceImage(testImage);
    registration->SetFixedPointSet( fixedPointSet);
    registration->SetMovingImage(targetDistanceImage);

    try {
        logger.Get(logINFO) << "starting model fitting of similarity part only" << std::endl;
        registration->Update();
    } catch ( itk::ExceptionObject& o ) {
        logger.Get(logINFO) << "caught exception " << o << std::endl;
    }



    // now also optimize statistical model
    StatisticalModelTransformType::Pointer statModelTransform = StatisticalModelTransformType::New();
    statModelTransform->SetStatisticalModel(model);
    statModelTransform->SetIdentity();

    // compose the two transformation
    transform->AddTransform(statModelTransform);

    transform->SetAllTransformsToOptimizeOn(); // optimize shape and pose parameters

    unsigned numStatmodelParameters = statModelTransform->GetNumberOfParameters();
    unsigned totalNumParameters = similarityTransform->GetNumberOfParameters() + numStatmodelParameters;

    // set the scales of the optimizer, to compensate for potentially different scales of translation, rotation and shape parameters
    OptimizerType::ScalesType allScales( totalNumParameters );
    for (unsigned i = 0; i < numStatmodelParameters; i++) {
        allScales[i] = 1.0 / (FittingConfigParameters::smScale);
    }
    for (unsigned i = numStatmodelParameters; i < numStatmodelParameters + 3; i++) {
        allScales[i] = 1.0 / (FittingConfigParameters::rotationScale);
    }

    for (unsigned i = numStatmodelParameters + 3; i < statModelTransform->GetNumberOfParameters() + 6; i++) {
        allScales[i] = 1.0 / (FittingConfigParameters::translationScale);
    }
    allScales[numStatmodelParameters + 6] =  1.0 / (FittingConfigParameters::scalingScale);
    optimizer->SetScales(allScales);

    registration->SetInitialTransformParameters(transform->GetParameters());


    try {
        logger.Get(logINFO) << "starting model full fitting" << std::endl;
        registration->Update();
    } catch ( itk::ExceptionObject& o ) {
        logger.Get(logINFO) << "caught exception " << o << std::endl;
    }


    MeshType::Pointer fittedMesh = model->DrawSample(statModelTransform->GetCoefficients());

    return std::make_pair(fittedMesh, similarityTransform);
}



