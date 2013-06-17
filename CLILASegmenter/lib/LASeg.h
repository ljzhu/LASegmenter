#ifndef LASEG_H
#define LASEG_H

// ITK
#include "itkImageMaskSpatialObject.h"
#include "itkAndImageFilter.h"
#include "itkOrImageFilter.h"

// VTK
#include "vtkImageGaussianSmooth.h"
#include "vtkMarchingCubes.h"
#include "vtkDecimatePro.h"
#include "vtkTriangleFilter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPLYReader.h"
#include "vtkPLYWriter.h"

// local
#include "./rss/labelMapPreprocessor.h"
#include "./rss/SFLSRobustStatSegmentor3DLabelMap_single.h"
//#include "./PoissonSurface/vtkPoissonReconstruction.h"

#include "utilities.hxx"
#include "ZernikeDescriptor.h"
#include "cmakeConfig.h"


namespace LASeg {

void LASegWithRSSandZM(const char* fnSrc, const char* fnMask);

template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(const char* fnSrc, const char* fnLab, const char* fnZM, double wIntensityIni, double wCurvature, \
                double wRSS, double exVolume, int nIterAll, int nIterRSS = 1);


template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(const char* fnSrc, const char* fnLab, const char* fnZM, double wIntensityIni, double wCurvature, \
                double wRSS, double exVolume, int nIterAll, int nIterRSS = 1);

template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(typename itk::Image<PixelType, 3>::Pointer imSrc, typename ImMaskType::Pointer imLab, const char* fnZM, double seedRad, double wIntensityIni, double wCurvature, \
                double wRSS, double exVolume, int nIterAll, int nIterRSS, std::string fnInter = "");

// RSS with Shape Constraint
void RSSWithShapeConstraint(const char* fnSrc, const char* fnSeed, const char* fnLA);

// Compute Zernike Moments for binary mask(s)
void ComputeZernikeMoments();


///*
//    Possion Surface Reconstruction
//*/
//template<typename ImMaskType>
//vtkPolyData* PossionSurfaceRecon(typename ImMaskType::Pointer imMask, float depth = 10, float deciReduction = 0.2, \
//                                 float deciError = 0.5, float scale = 1.25, float solDivide = 8, \
//                                 float isoDivide = 8, float samNode = 1.0);

/*
            Convert VTK to PLY
*/
vtkPolyData* VTK2PolyWithLargestConnectedComponent(vtkPolyData* surface, std::string fnPolySave = "");


///*
//            Extract surface reconstructed from Zernike Moments
//*/
//template< typename ImMaskType>
//vtkPolyData* ExtractZernikeSurface(typename ImMaskType::Pointer imMask, float depth = 8,  std::string fnSave = "");

///*
//            Extract surfaces reconstructed from Zernike Moments
//*/
//void ExtractZernikeSurfaces();

/*
            Convert VTK to PLY with Largest Connected Component
*/
void VTK2PLYWithLargestLargestConnectedComponent();

/*
            Resample Images
*/
void ReSampleImages();

/*
        Compute Volume
*/
template<typename LabImType>
double CalVolume(typename LabImType::Pointer maskImg);

/*
        Compute Volume Distribution
*/
void ComputeVolumeDistribution();
}


#include "LASeg.txx"

#endif // LASEG_H
