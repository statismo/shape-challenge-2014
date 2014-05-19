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

MeshType::Pointer fitModelToTestImage(Logger& logger, StatisticalModelType::Pointer model, BinaryImageType::Pointer testImage);



GeneralizationResult generalization(Logger& logger, StatisticalModelType::Pointer model, const TestImageList& testImages) {

	double totalAvgDistance = 0;
	double totalHdDistance = 0;
	unsigned numTestImages = testImages.size();

	for (TestImageList::const_iterator it = testImages.begin(); it != testImages.end(); ++it) { 
        std::string filename = it->first;
        logger.Get(logINFO) << "fitting image " << filename << std::endl;
        BinaryImageType::Pointer testImage = it->second;
        MeshType::Pointer fittedMesh = fitModelToTestImage(logger, model, testImage);
		
        MeshType::Pointer testMesh = binaryImageToMesh(testImage);
        double avgDistance = computeSymmetricAverageDistance(testMesh,  fittedMesh, ConfigParameters::numSamplingPointsGeneralization);
        totalAvgDistance += avgDistance;
        double hdDistance = computeHausdorffDistance(testMesh, fittedMesh, ConfigParameters::numSamplingPointsGeneralization);
        totalHdDistance += hdDistance;

        logger.Get(logINFO) << "avg distance " << avgDistance << std::endl;
        logger.Get(logINFO) << "hd distance " << hdDistance << std::endl;
	}
	return GeneralizationResult(totalAvgDistance / numTestImages, totalHdDistance / numTestImages);
}


// resize image2 such that it has the same size as refImage
BinaryImageType::Pointer resizeToImage(BinaryImageType::Pointer refImage, BinaryImageType::Pointer image2) { 
	typedef itk::ResampleImageFilter<BinaryImageType, BinaryImageType> ResampleImageFilterType;
	typedef itk::NearestNeighborInterpolateImageFunction<BinaryImageType> NNInterpolatorType;
	itk::IdentityTransform<double, 3>::Pointer identTransform = itk::IdentityTransform<double, 3>::New();

	ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
	resampler->SetOutputParametersFromImage(refImage);
	resampler->SetInput(image2);
	resampler->SetTransform(identTransform);
	NNInterpolatorType::Pointer nnInterpolator = NNInterpolatorType::New();
	resampler->SetInterpolator(nnInterpolator);
	resampler->SetDefaultPixelValue(0);
	resampler->Update();
	return resampler->GetOutput();
}




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
       m_logger->Get(logINFO) << "Iteration: " << ++m_iter_no << "  value " <<  optimizer->GetValue() << "model parameters " << optimizer->GetCurrentPosition() << std::endl;
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


MeshType::Pointer fitModelToTestImage(Logger& logger, StatisticalModelType::Pointer model, BinaryImageType::Pointer testImage)
{
	const unsigned Dimensions = 3;
	typedef itk::PointSet<float, Dimensions > PointSetType;
	typedef itk::Mesh<float, Dimensions > MeshType;
	typedef itk::Point<double, 3> PointType;
	typedef itk::MeanSquaresPointSetToImageMetric<PointSetType, DistanceImageType> MetricType;
	typedef itk::PointSetToImageRegistrationMethod<PointSetType, DistanceImageType> RegistrationFilterType;
    //typedef itk::LBFGSOptimizer OptimizerType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	typedef itk::LinearInterpolateImageFunction<DistanceImageType, double> InterpolatorType;
    typedef itk::VersorRigid3DTransform<double> RigidTransformType;
	typedef itk::CompositeTransform<double, 3> CompositeTransformType;
    typedef itk::TransformMeshFilter<MeshType, MeshType, CompositeTransformType > TransformMeshFilterType;
    typedef itk::CenteredVersorTransformInitializer<BinaryImageType, BinaryImageType> TransformInitializerType;
	typedef itk::StatisticalShapeModelTransform<MeshType, double, Dimensions> StatisticalModelTransformType;


    double meanImageMargin = 40.0;
    unsigned meanImageResolution = 256;

	// we compute a binary image of the model mean to use it for initialization of the shape 
    BinaryImageType::Pointer refAsBinImage = meshToBinaryImage(model->GetRepresenter()->GetReference(), meanImageResolution, meanImageMargin);
	RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
	TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetFixedImage(refAsBinImage);
	initializer->SetMovingImage(testImage);
	initializer->SetTransform(rigidTransform);
    initializer->ComputeRotationOn();
    initializer->MomentsOn();
	initializer->InitializeTransform();

	StatisticalModelTransformType::Pointer statModelTransform = StatisticalModelTransformType::New();
	statModelTransform->SetStatisticalModel(model);
	statModelTransform->SetIdentity();

	// compose the two transformation
	CompositeTransformType::Pointer transform = CompositeTransformType::New();
	transform->AddTransform(rigidTransform);
	transform->AddTransform(statModelTransform);
	transform->SetAllTransformsToOptimizeOn(); // optimize shape and pose parameters

	// Setting up the fitting
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetMaximumStepLength(FittingConfigParameters::maxStepLength);
    optimizer->MinimizeOn();
    optimizer->SetNumberOfIterations(FittingConfigParameters::maxNumberOfIterations);


	unsigned numStatmodelParameters = statModelTransform->GetNumberOfParameters();
	unsigned totalNumParameters = rigidTransform->GetNumberOfParameters() + numStatmodelParameters;

	// set the scales of the optimizer, to compensate for potentially different scales of translation, rotation and shape parameters
	OptimizerType::ScalesType scales( totalNumParameters );
	for (unsigned i = 0; i < numStatmodelParameters; i++) {
		scales[i] = 1.0 / (FittingConfigParameters::smScale);
	}
	for (unsigned i = numStatmodelParameters; i < numStatmodelParameters + 3; i++) {
		scales[i] = 1.0 / (FittingConfigParameters::rotationScale);
	}

	for (unsigned i = numStatmodelParameters + 3; i < statModelTransform->GetNumberOfParameters() + 6; i++) {
		scales[i] = 1.0 / (FittingConfigParameters::translationScale);
	}
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
        logger.Get(logINFO) << "starting model fitting" << std::endl;
		registration->Update();
	} catch ( itk::ExceptionObject& o ) {
        logger.Get(logINFO) << "caught exception " << o << std::endl;
	}


	TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
	transformMeshFilter->SetInput(model->GetRepresenter()->GetReference());
    transformMeshFilter->SetTransform(transform);
	transformMeshFilter->Update();

	return transformMeshFilter->GetOutput();
}




