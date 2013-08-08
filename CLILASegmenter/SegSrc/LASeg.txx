#include "LASeg.h"

// ITK
#include "itkImageMaskSpatialObject.h"
#include "itkAndImageFilter.h"
#include "itkOrImageFilter.h"

// local
#include "./rss/labelMapPreprocessor.h"
#include "./rss/SFLSRobustStatSegmentor3DLabelMap_single.h"

namespace LASeg {

/************************************************************


************************************************************/
template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(typename itk::Image<PixelType, 3>::Pointer imSrc, typename ImMaskType::Pointer imLab, double seedRad, \
                double wIntensityIni, double wCurvature, double wRSS, double exVolume, int nIterAll, int nIterRSS,std::string fnInter){

    typedef CSFLSRobustStatSegmentor3DLabelMap< PixelType > SFLSRobustStatSegmentor3DLabelMap_c;

    // Do segmentation
    SFLSRobustStatSegmentor3DLabelMap_c seg;
    seg.setImage(imSrc);
    seg.setInputLabelImage(imLab);
    seg.setSeedRad(seedRad);
    seg.SetSaveInter(fnInter);

    seg.setNumIter(nIterAll); // a large enough number, s.t. will not be stoped by this creteria.
    seg.SetRSSNumIter(nIterRSS);

    double expectedVolume = exVolume; // the mean volume of the training pre-endo is 108
    seg.setMaxVolume(expectedVolume);

    double maxRunningTime = 30; // max running time 30 min
    seg.setMaxRunningTime(maxRunningTime);

    seg.setIntensityHomogeneity(wIntensityIni/2.0);
    seg.setCurvatureWeight(wCurvature/1.5);
    seg.SetRSSWeight(wRSS);

    seg.doSegmenation();

    // Get final mask
    typedef itk::BinaryThresholdImageFilter<itk::Image<PixelType, 3>, ImMaskType> binaryThresholdImageFilter_t;

    typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
    thlder->SetInput(seg.mp_phi);
    thlder->SetInsideValue(1);
    thlder->SetOutsideValue(0);
    thlder->SetLowerThreshold(-5);
    thlder->SetUpperThreshold(vnl_math::eps);
    thlder->Update();

    return thlder->GetOutput();
}

/************************************************************


************************************************************/
void RSSWithShapeConstraint(const char* fnSrc, const char* fnSeed, const char* fnLA, double wIntensity, double wCurvature) {

    const double wRSS = 0.3;
    const double exVolRSS = 60;
    const double exVolRSSZM = 300;
    const int iterRSS = 100;
    const int iter = 300;

    typedef float PixelType;
    typedef itk::Image< PixelType, 3 > ImSrcType;
    typedef itk::Image< char, 3 > ImMaskType;
    ImSrcType::Pointer imSrc, imSrcSmooth, imSrcIso;
    ImMaskType::Pointer imLab, imLabIso, imLA;

    imSrc = LASeg::readImage<ImSrcType>(fnSrc);
    imLab = LASeg::readImage<ImMaskType>(fnSeed);
    double spacingOri[3];
    spacingOri[0] =  imSrc->GetSpacing()[0];
    spacingOri[1] =  imSrc->GetSpacing()[1];
    spacingOri[2] =  imSrc->GetSpacing()[2];

    imSrcSmooth = LASeg::AnisotropicSmoothing<ImSrcType>(imSrc);
    float laSeedRad;
    laSeedRad = EstimateSeedRadius<ImMaskType>(imLab);

    // Isotropic resampling required by Zernike Moments comuptation
    const double spacingX = std::max(imSrc->GetSpacing()[0]*2, 1.0);
    double spacing[] = {spacingX, spacingX, spacingX};
    imSrcIso = LASeg::ResampleImage<ImSrcType, ImSrcType>(imSrcSmooth, spacing, 1);
    imLabIso = LASeg::ResampleImage<ImMaskType, ImMaskType>(imLab, spacing, 0);

    imLA = LASeg::SegWithRSSandZM<PixelType, ImMaskType>(imSrcIso, imLabIso, laSeedRad, wIntensity, wCurvature, wRSS, exVolRSS, iterRSS, iterRSS);
    imLA = LASeg::SegWithRSSandZM<PixelType, ImMaskType>(imSrcIso, imLA,  laSeedRad, wIntensity, wCurvature, wRSS, exVolRSSZM, iter, 0);

    imLA = LASeg::ResampleImage<ImMaskType, ImMaskType>(imLA, spacingOri, 0);

    LASeg::writeImage<ImMaskType>(imLA, fnLA);
}

/***************************************************************

***************************************************************/
template<typename ImSrcType>
double EstimateSeedRadius(typename ImSrcType::Pointer imSrc) {

    const double *spacing = imSrc->GetSpacing().GetDataPointer();
    typedef itk::ImageRegionConstIterator<ImSrcType> ImageRegionConstIterator;
    ImageRegionConstIterator it(imSrc, imSrc->GetLargestPossibleRegion() );

   long count = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
        count += it.Get() > 0 ? 1 : 0;
    }

    return std::sqrt(double(count)/PI)*spacing[0]; // assume the seed region locates at the Z/Axial direction only
}

}


