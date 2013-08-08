#ifndef LASEG_H
#define LASEG_H

#include "utilities.hxx"
#include "ZernikeDescriptor.h"
#include "cmakeConfig.h"


namespace LASeg {

// RSS with Shape Constraint
void RSSWithShapeConstraint(const char* fnSrc, const char* fnSeed, const char* fnLA);

template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(typename itk::Image<PixelType, 3>::Pointer imSrc, typename ImMaskType::Pointer imLab, double seedRad, double wIntensityIni, double wCurvature, \
                double wRSS, double exVolume, int nIterAll, int nIterRSS, std::string fnInter = "");

template<typename ImSrcType>
double EstimateSeedRadius(typename ImSrcType::Pointer imSrc);

}


#include "LASeg.txx"

#endif // LASEG_H
