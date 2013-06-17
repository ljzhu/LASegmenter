#ifndef SFLSRobustStatSegmentor3DLabelMap_single_txx_
#define SFLSRobustStatSegmentor3DLabelMap_single_txx_

#include "SFLSRobustStatSegmentor3DLabelMap_single.h"

#include <algorithm>
#include <ctime>

#include <limits>

// //debug//
// #include "cArrayOp.h"
#include <fstream>
// //DEBUG//

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"

/* ============================================================   */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::basicInit()
{
  SuperClassType::basicInit();

  m_statNeighborX = 1;
  m_statNeighborY = 1;
  m_statNeighborZ = 1;

  m_kernelWidthFactor = 10.0;

  m_inputImageIntensityMin = 0;
  m_inputImageIntensityMax = 0;

  m_laSeedRad = 0.0;

  m_numIterRSS = 1;
  m_fnInter = "";

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::setInputLabelImage(TLabelImagePointer l)
{
    // only works for isotropic images
    const double* spacing = l->GetSpacing().GetDataPointer();
    if( !(spacing[0] == spacing[1] && spacing[1] == spacing[2] && spacing[0] == spacing[2]) ) {
        std::cerr << "Error: do isotropic resampling first\n";
        raise(SIGABRT);
    }

  m_inputLabelImage = l;

  TSize size = m_inputLabelImage->GetLargestPossibleRegion().GetSize();

  TIndex start = m_inputLabelImage->GetLargestPossibleRegion().GetIndex();
  TIndex origin = {{0, 0, 0}};
  if( start != origin )
    {
    std::cout << "Warrning: Force mask start to be (0, 0, 0)\n";

    TRegion region = m_inputLabelImage->GetLargestPossibleRegion();
    region.SetIndex(origin);

    m_inputLabelImage->SetRegions(region);
    }

  if( this->m_nx + this->m_ny + this->m_nz == 0 )
    {
    this->m_nx = size[0];
    this->m_ny = size[1];
    this->m_nz = size[2];
    }
  else if( this->m_nx != (long)size[0] || this->m_ny != (long)size[1] || this->m_nz != (long)size[2] )
    {
    std::cerr << "Error: image sizes do not match with label image size.\n";
    raise(SIGABRT);
    }

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::computeForce()
{
      double kappaMax = std::numeric_limits<double>::min();
      double rssMax = std::numeric_limits<double>::min();

      long n = this->m_lz.size();
      double* kappaOnZeroLS = new double[ n ];
      double* zmForce = new double[ n ];
      double* rssForce = new double[ n ];


      std::vector<typename CSFLSLayer::iterator> m_lzIterVct( n );
      {
        long iiizzz = 0;
        for (typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz)
          m_lzIterVct[iiizzz++] = itz;
      }



    // #ifndef NDEBUG
    //     std::ofstream ff("/tmp/force.txt");
    // #endif


      for (long i = 0; i < n; ++i)
        {
          typename CSFLSLayer::iterator itz = m_lzIterVct[i];

          long ix = (*itz)[0];
          long iy = (*itz)[1];
          long iz = (*itz)[2];

          long ixs = (*itz)[0] - m_ZMOrigin[0];
          long iys = (*itz)[1] - m_ZMOrigin[1];
          long izs = (*itz)[2] - m_ZMOrigin[2];

          zmForce[i] = -m_ZM.GetMomentVariation(ixs, iys, izs);

           TIndex idx = {{ix, iy, iz}};

          kappaOnZeroLS[i] = this->computeAreaVolumeVariation(ix, iy, iz);

          std::vector<double> f(m_numberOfFeature);

          computeFeatureAt(idx, f);

          double a = -kernelEvaluationUsingPDF(f);
          rssForce[i] = a;

          rssMax = rssMax>fabs(a)?rssMax:fabs(a);
          kappaMax = kappaMax>fabs(kappaOnZeroLS[i])?kappaMax:fabs(kappaOnZeroLS[i]);
        }

        double ZMNorm, KappaNorm;
        ZMNorm = 0.0; KappaNorm = 0.0;
        for(long i = 0; i < n; i++) {
            ZMNorm += zmForce[i]*zmForce[i];
        }
          ZMNorm = std::sqrt(ZMNorm);

      this->m_force.resize(n);
      for (long i = 0; i < n; ++i)
        {
          this->m_force[i] = m_rssWeight*(rssForce[i]/(rssMax + vnl_math::eps)) + \
                  (1 - m_rssWeight)*(((1 - (this->m_curvatureWeight))*(zmForce[i]/(ZMNorm + 1e-10) + \
                             (this->m_curvatureWeight)*kappaOnZeroLS[i])/(kappaMax + 1e-10)));
        }

        delete[] kappaOnZeroLS;
        delete[] zmForce;
        delete[] rssForce;

}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::ComputeRSSForce()
{
  double fmax = std::numeric_limits<double>::min();
  double kappaMax = std::numeric_limits<double>::min();

  long    n = this->m_lz.size();
  double* kappaOnZeroLS = new double[n];
  double* cvForce = new double[n];

  std::vector<typename CSFLSLayer::iterator> m_lzIterVct( n );
    {
    long iiizzz = 0;
    for( typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz )
      {
      m_lzIterVct[iiizzz++] = itz;
      }
    }

// #ifndef NDEBUG
//     std::ofstream ff("/tmp/force.txt");
// #endif
  for( long i = 0; i < n; ++i )
    {
    typename CSFLSLayer::iterator itz = m_lzIterVct[i];

    long ix = (*itz)[0];
    long iy = (*itz)[1];
    long iz = (*itz)[2];

    TIndex idx = {{ix, iy, iz}};

    kappaOnZeroLS[i] = this->computeKappa(ix, iy, iz);

    std::vector<double> f(m_numberOfFeature);

    computeFeatureAt(idx, f);

    double a = -kernelEvaluationUsingPDF(f);

    fmax = fmax > fabs(a) ? fmax : fabs(a);
    kappaMax = kappaMax > fabs(kappaOnZeroLS[i]) ? kappaMax : fabs(kappaOnZeroLS[i]);

    cvForce[i] = a;
    }

  this->m_force.resize(n);
  for( long i = 0; i < n; ++i )
    {
    this->m_force[i] = (1 - (this->m_curvatureWeight) ) * cvForce[i] / (fmax + 1e-10) \
      +  (this->m_curvatureWeight) * kappaOnZeroLS[i] / (kappaMax + 1e-10);
    }

  delete[] kappaOnZeroLS;
  delete[] cvForce;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::ComputeRSSZMForce()
{
      double kappaMax = std::numeric_limits<double>::min();
      double rssMax = std::numeric_limits<double>::min();

      long n = this->m_lz.size();
      double* kappaOnZeroLS = new double[ n ];
      double* zmForce = new double[ n ];
      double* rssForce = new double[ n ];


      std::vector<typename CSFLSLayer::iterator> m_lzIterVct( n );
      {
        long iiizzz = 0;
        for (typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz)
          m_lzIterVct[iiizzz++] = itz;
      }



    // #ifndef NDEBUG
    //     std::ofstream ff("/tmp/force.txt");
    // #endif


      for (long i = 0; i < n; ++i)
        {
          typename CSFLSLayer::iterator itz = m_lzIterVct[i];

          long ix = (*itz)[0];
          long iy = (*itz)[1];
          long iz = (*itz)[2];

          long ixs = (*itz)[0] - m_ZMOrigin[0];
          long iys = (*itz)[1] - m_ZMOrigin[1];
          long izs = (*itz)[2] - m_ZMOrigin[2];

          zmForce[i] = -m_ZM.GetMomentVariation(ixs, iys, izs);

          TIndex idx = {{ix, iy, iz}};

          kappaOnZeroLS[i] = this->computeKappa(ix, iy, iz);

          std::vector<double> f(m_numberOfFeature);

          computeFeatureAt(idx, f);

          double a = -kernelEvaluationUsingPDF(f);
          rssForce[i] = a;

          rssMax = rssMax>fabs(a)?rssMax:fabs(a);
          kappaMax = kappaMax>fabs(kappaOnZeroLS[i])?kappaMax:fabs(kappaOnZeroLS[i]);
        }

        double ZMNorm, KappaNorm;
        ZMNorm = 0.0; KappaNorm = 0.0;
        for(long i = 0; i < n; i++) {
            ZMNorm += zmForce[i]*zmForce[i];
        }
        ZMNorm = std::sqrt(ZMNorm);

      // update coefficients
      // NOTE: the terms containing the derivative of RSS_Weight were dropped because of the
      // individual normalization of the RSS, ZM, and Kappa terms !!
      EvaluateRSSCof();

      this->m_force.resize(n);
      for (long i = 0; i < n; ++i) {
          this->m_force[i] = m_rssWeight*(rssForce[i]/(rssMax + vnl_math::eps)) + \
                             this->m_curvatureWeight*kappaOnZeroLS[i]/(kappaMax + 1e-10) +
                  (1.0 - m_rssWeight - this->m_curvatureWeight)*zmForce[i]/(ZMNorm + 1e-10);
        }

        delete[] kappaOnZeroLS;
        delete[] zmForce;
        delete[] rssForce;
}




/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::inputLableImageToSeeds()
{
  typedef itk::ImageRegionConstIteratorWithIndex<TLabelImage> ImageRegionConstIteratorWithIndex_t;
  ImageRegionConstIteratorWithIndex_t it(m_inputLabelImage, m_inputLabelImage->GetLargestPossibleRegion() );
  it.GoToBegin();

  std::ofstream sf("_seeds.txt");

    {
    std::vector<long> thisSeed(3);
    for( ; !it.IsAtEnd(); ++it )
      {
      if( it.Get() != 0 )  // 0 for bkgd, every label else is obj
        {
        TIndex idx = it.GetIndex();
        thisSeed[0] = idx[0];
        thisSeed[1] = idx[1];
        thisSeed[2] = idx[2];

        m_seeds.push_back(thisSeed);

        sf << thisSeed[0] << ", " << thisSeed[1] << ", " << thisSeed[2] << std::endl;
        }
      }
    }

  sf.close();

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::getThingsReady()
{
  /*
   1. Generate mp_mask from seeds
   2. Compute feature at each point
   3. Extract feature at/around the seeds
  */

  inputLableImageToSeeds();

  seedToMask();

  // dialteSeeds();

  initFeatureComputedImage();
  initFeatureImage();

  // computeFeature();
  getFeatureAroundSeeds();
  estimateFeatureStdDevs();

  estimatePDFs();

  return;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::initFeatureImage()
{
  if( !(this->mp_img) )
    {
    std::cerr << "Error: set input image first.\n";
    raise(SIGABRT);
    }
  for( long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature )
    {
    // TDoubleImagePointer fimg = TDoubleImage::New();

    TFloatImagePointer fimg = TFloatImage::New();
    fimg->SetRegions(this->mp_img->GetLargestPossibleRegion() );
    fimg->Allocate();
    fimg->CopyInformation(this->mp_img);

    m_featureImageList.push_back(fimg);
    }

  return;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::initFeatureComputedImage()
{
  if( !(this->mp_img) )
    {
    std::cerr << "Error: set input image first.\n";
    raise(SIGABRT);
    }

  m_featureComputed = TLabelImage::New();
  m_featureComputed->SetRegions(this->mp_img->GetLargestPossibleRegion() );
  m_featureComputed->Allocate();
  m_featureComputed->CopyInformation(this->mp_img);
  m_featureComputed->FillBuffer(0);

  return;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::computeFeatureAt(TIndex idx, std::vector<double>& f)
{
  f.resize(m_numberOfFeature);

  if( m_featureComputed->GetPixel(idx) )
    {
    // the feature at this pixel is computed, just retrive
    for( long i = 0; i < m_numberOfFeature; ++i )
      {
      f[i] = (m_featureImageList[i])->GetPixel(idx);
      }
    }
  else
    {
    // compute the feature
    std::vector<double> neighborIntensities;

    long ix = idx[0];
    long iy = idx[1];
    long iz = idx[2];
    for( long iiz = iz - m_statNeighborZ; iiz <= iz + m_statNeighborZ; ++iiz )
      {
      for( long iiy = iy - m_statNeighborY; iiy <= iy + m_statNeighborY; ++iiy )
        {
        for( long iix = ix - m_statNeighborX; iix <= ix + m_statNeighborX; ++iix )
          {
          if( 0 <= iix && iix < this->m_nx    \
              && 0 <= iiy && iiy < this->m_ny    \
              && 0 <= iiz && iiz < this->m_nz )
            {
            TIndex idxa = {{iix, iiy, iiz}};
            neighborIntensities.push_back(this->mp_img->GetPixel(idxa) );
            }
          }
        }
      }

    getRobustStatistics(neighborIntensities, f);
    for( long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature )
      {
      m_featureImageList[ifeature]->SetPixel(idx, f[ifeature]);
      }

    m_featureComputed->SetPixel(idx, 1);   // mark as computed
    }

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::EvaluateERSS() {

    long n = this->m_lz.size();

    std::vector<typename CSFLSLayer::iterator> m_lzIterVct( n );
      {
      long iiizzz = 0;
      for( typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz )
        {
        m_lzIterVct[iiizzz++] = itz;
        }
      }

  // #ifndef NDEBUG
  //     std::ofstream ff("/tmp/force.txt");
  // #endif
    m_ERSSVal = 0.0;
    for( long i = 0; i < n; ++i ) {

      typename CSFLSLayer::iterator itz = m_lzIterVct[i];

      long ix = (*itz)[0];
      long iy = (*itz)[1];
      long iz = (*itz)[2];

      TIndex idx = {{ix, iy, iz}};

      std::vector<double> f(m_numberOfFeature);
      computeFeatureAt(idx, f);

      // double a = -kernelEvaluation(f);
      m_ERSSVal  += kernelEvaluationUsingPDF(f);
    }
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::UpdateRSSZMShiftFactor() {

    if(m_laSeedRad < vnl_math::eps)
        return;

    double shift;

    shift = (m_laSeedRad - m_LARMean)/m_LARSTD;

    if(shift >= -1.0 && shift <= 1.0) {
        shift = 0;
    }
    else if(shift > 1.0 && shift <= 2.0) {
        shift = 1.0;
    }
    else if(shift > 2.0 && shift <= 3.0) {
        shift = 2.0;
    }
    else if(shift > 3.0) {
        shift = 3.0;
    }
    else if(shift < -1.0 && shift >= -2.0) {
        shift = -1.0;
    }
    else if(shift < -2.0 && shift >= -3.0) {
         shift = -2.0;
    }
    else if(shift < -3.0) {
        shift = -2.5;
    }


    m_RSSZMShiftFactor = m_RSSZMShiftFactor + shift;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::EvaluateRSSCof() {

    m_rssWeight = (1.0 - this->m_curvatureWeight)/(1.0 \
              + exp(m_RSSZMFactor*(this->m_insideVolume/1000.0 - (m_LAVolMean + m_RSSZMShiftFactor*m_LAVolSTD))));
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::EvaluateDevRSSCof() {

    double expVal;

    expVal = exp(m_RSSZMFactor*(this->m_insideVolume -m_LAVolMean));
    m_devRssWeight = -(1.0 - this->m_curvatureWeight)*m_RSSZMFactor*expVal/(1.0 + expVal)/(1.0 + expVal);
}


/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::doSegmenation()
{
    double startingTime = clock();

    getThingsReady();

    UpdateRSSZMShiftFactor();

    /*============================================================
     * From the initial mask, generate: 1. SFLS, 2. mp_label and
     * 3. mp_phi.
     */
    this->initializeSFLS();

  // #ifndef NDEBUG
  //   std::ofstream dbgf("/tmp/dbgo.txt", std::ios_base::app);
  // #endif
    const double STOPSCALE = 1.1; // default 1.1
    std::list<double> eValList;
    const unsigned int eWinSize = 10;
    double eVal, eValMed, eValNxt, eValFirst, eVal0, eValMax = 0.0;

    if(m_numIterRSS > 0) {
      this->m_maxVolume = std::max((m_LAVolMean + m_RSSZMShiftFactor*m_LAVolSTD)*1000, 6000.0);
    }

    // TODO: Detection spliting and set boundary constraint for it !!!
    for (unsigned int it = 0; it < this->m_numIter; ++it) {
        m_zmInter = it;

        //dbg//
//        std::cout<<"In iteration "<<it<<std::endl<<std::flush;
        //DBG//


        // keep current zero contour as history is required
        if (this->m_keepZeroLayerHistory) {
            (this->m_zeroLayerHistory).push_back(this->m_lz);
          }

        double oldVoxelCount = this->m_insideVoxelCount;

        // RSS
        if(it < m_numIterRSS) {            
            ComputeRSSForce();
        }
        // RSS with ZM Constraint
        else {
            // Compute Zernike Moments
            ComputeZM();
            eVal = m_ZM.GetZernikeDistance();
            if(it == m_numIterRSS) {
                eVal0 = eVal;
            }

            // Stop condition
            if(it < eWinSize + m_numIterRSS) {
              eValList.push_front(eVal);
              if(eVal > eValMax) {
                  eValMax = eVal;
              }
            }
            else {
                // find the Median between the three latest eVals
                std::list<double>::iterator it;
                it = eValList.begin();
                eValFirst = *it++;
                eValNxt = *it;

                if(eValFirst > std::min(eVal, eValNxt) && eValFirst < std::max(eVal, eValNxt)) {
                    eValMed = eValFirst;
                }
                else {
                    eValMed = std::min(eVal, eValNxt);
                }

                if(eValMed > eValMax*STOPSCALE && eValMed < eVal0) {
//                    std::cout << "stop condition reached!\n";
                    break;    // stop if it is greater than the maximum in period of previous iterations
                }
                else {
                  eValList.pop_back();
                  eValList.push_front(eVal);
                  eValMax = 0.0;

                  it = eValList.begin();
                  for(std::advance(it,2); it != eValList.end(); ++it) {
                      eValMax = eValMax < *it ? *it : eValMax;
                  }
                }
            }

            ComputeRSSZMForce();
        }

        this->normalizeForce();

        this->oneStepLevelSetEvolution();

        /*----------------------------------------------------------------------
          If the level set stops growing, stop */
        this->updateInsideVoxelCount();
        if( it > m_numIterRSS + 2 && oldVoxelCount >= this->m_insideVoxelCount )
          {
          std::ofstream f("/tmp/o.txt");
          f << "In the " << it << "-th iteration\n";

          f << "stop grow\n";
          f << "oldVoxelCount = " << oldVoxelCount << std::endl;
          f << "m_insideVoxelCount = " << this->m_insideVoxelCount << std::endl;

          f << "m_kernelWidthFactor = " << m_kernelWidthFactor << std::endl;
          f << "m_maxRunningTime = " << this->m_maxRunningTime << std::endl;
          f.close();

          break;
          }

        /* If the level set stops growing, stop
           ----------------------------------------------------------------------*/

        /*----------------------------------------------------------------------
          If the inside physical volume exceed expected volume, stop */
        double volumeIn = (this->m_insideVoxelCount) * (this->m_dx) * (this->m_dy) * (this->m_dz);
        if( volumeIn > (this->m_maxVolume) )
          {
          //          std::fstream f("/tmp/o.txt", std::ios_base::app);
          std::ofstream f("/tmp/o.txt");
          f << "In the " << it << "-th iteration\n";
          f << "reach max volume\n";

          f << "m_maxVolume = " << this->m_maxVolume << std::endl;
          f << "volumeIn = " << volumeIn << std::endl;

          f.close();

          break;
          }
        /*If the inside physical volume exceed expected volume, stop
          ----------------------------------------------------------------------*/

        double ellapsedTime = (clock() - startingTime) / static_cast<double>(CLOCKS_PER_SEC);
        if( ellapsedTime > (this->m_maxRunningTime) )
          {
          std::ofstream f("/tmp/o.txt");
          f << "running time = " << ellapsedTime << std::endl;
          f << "m_maxRunningTime = " << this->m_maxRunningTime << std::endl;
          f.close();

          break;
          }

    }

//    double ellapsedTime = (clock() - startingTime)/static_cast<double>(CLOCKS_PER_SEC);
//    std::cout << "avg. time: " << ellapsedTime/this->m_numIter << std::endl;

// #ifndef NDEBUG
//   dbgf.close();
// #endif

  return;
}

/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::LoadZMReference() {

    if(m_fnZMRef != "") {

        std::ifstream infile(m_fnZMRef.c_str());
        if (infile.good()) {

            m_ZMRef.clear();

            std::stringstream ss;
            std::string line;
            float real, imag;
            std::getline(infile, line);

            ss.str(line);
            ss >> m_ZMOrder;

            while(std::getline(infile, line))  {
                if(!line.empty()) {
                    ss.clear();
                    ss.str(line);
                    ss >> real >> imag;
                    ComplexT ZM(real, imag);
                    m_ZMRef.push_back(ZM);
                }
                else break;
              }

            m_ZM.SetZMReference(m_ZMRef);
          }
        else
          {
            std::cerr<<"Error: can not open file:"<<m_fnZMRef<<std::endl;
            raise(SIGABRT);
          }

        infile.close();

    }
    else {

        std::cout << "the reference ZM file doesn't exist:" << m_fnZMRef << std::endl;
    }
}

/* ============================================================ */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::LoadRadVolPrior() {

    m_LAVolMean =  106.376;
    m_LAVolSTD = 45.0909;
    m_LARMean =  11.58;
    m_LARSTD =  3.13;
    m_RSSZMFactor = 0.03;
    m_RSSZMShiftFactor =  -1.0;

//    std::cout << "prior loaded: " << m_fnRadVolPrior << std::endl;
//    if(m_fnRadVolPrior != "") {
//        std::string strID;

//        std::ifstream configFile(m_fnRadVolPrior.c_str());
//        if(configFile.is_open()) {
//           configFile >> strID >> m_LAVolMean;
//           configFile >> strID >> m_LAVolSTD;
//           configFile >> strID >> m_LARMean;
//           configFile >> strID >> m_LARSTD;
//           configFile >> strID >> m_RSSZMFactor;
//           configFile >> strID >> m_RSSZMShiftFactor;

//           std::cout << m_LAVolMean << "," << m_LAVolSTD << std::endl;

//           configFile.close();

//       }
//        else {
//            std::cout << "can not open " << m_fnRadVolPrior << std::endl;
//            return;
//        }

//    }
//    else {

//        std::cout << "the RadVolPrior file doesn't exist:" << m_fnRadVolPrior << std::endl;
//    }
}


/* ============================================================ */
template< typename TPixel >
void CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::ComputeZM() {

    // Extract the region within the contour
    long n = this->m_lz.size();
    long bound[6] = {100000, 0, 100000, 0, 100000, 0};

    std::vector<typename CSFLSLayer::iterator> m_lzIterVct( n );
    {
      long iiizzz = 0;
      for (typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz)
        m_lzIterVct[iiizzz++] = itz;
    }

    for (long i = 0; i < n; ++i) {
        typename CSFLSLayer::iterator itz = m_lzIterVct[i];

        long ix = (*itz)[0];
        long iy = (*itz)[1];
        long iz = (*itz)[2];

        if(ix < bound[0]) {
            bound[0] = ix;
        }
        if(ix > bound[1]) {
            bound[1] = ix;
        }
        if(iy < bound[2]) {
            bound[2] = iy;
        }
        if(iy > bound[3]) {
            bound[3] = iy;
        }
        if(iz < bound[4]) {
            bound[4] = iz;
        }
        if(iz > bound[5]) {
            bound[5] = iz;
        }
    }

    m_ZMOrigin[0] = bound[0];
    m_ZMOrigin[1] = bound[2];
    m_ZMOrigin[2] = bound[4];

    TSize size = {{bound[1] - bound[0] + 1, bound[3] - bound[2] + 1, bound[5] - bound[4]+ 1}};

    long dim = std::max(std::max(size[0], size[1]), size[2]);
    long total;//, widLine, areaPlane;
    total = dim*dim*dim;
    float* pData = new float[total];
    std::fill(pData, pData+total, 0.0);

    typedef itk::Image<unsigned char, 3> ImMaskType;
    typename ImMaskType::Pointer imLab = ImMaskType::New();

    if(m_fnInter != "") {
        typedef ImMaskType::IndexType TIndex;
        typename ImMaskType::IndexType start;
        typename ImMaskType::RegionType region;
         start[0] = 0;
         start[1] = 0;
         start[2] = 0;

         region.SetSize(this->mp_img->GetLargestPossibleRegion().GetSize());
         region.SetIndex(start);

         imLab->SetRegions(region);
         imLab->SetSpacing(this->mp_img->GetSpacing());
         imLab->SetOrigin(this->mp_img->GetOrigin());
         imLab->Allocate();
         imLab->FillBuffer(0);
    }

    TIndex idx;
    unsigned int i,j,k;
    for(i = 0; i < size[0]; i++) {
        for(j = 0; j < size[1]; j++) {
            for(k = 0; k < size[2]; k++) {
                idx[0] = i + m_ZMOrigin[0];
                idx[1] = j + m_ZMOrigin[1];
                idx[2] = k + m_ZMOrigin[2];
                pData[i + dim*(j + k*dim)] = this->mp_phi->GetPixel(idx) < vnl_math::eps ? 1.0 : 0.0;

                if(m_fnInter != "")
                    imLab->SetPixel(idx, pData[i + dim*(j + k*dim)]);
            }
        }
    }

   // Compute Zernike Moments
   m_ZM.ComputeZernikeMoments(pData, dim, m_ZMOrder);

    delete []pData;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::getRobustStatistics(std::vector<double>& samples, std::vector<double>& robustStat)
{
  /* note, sample is sorted, so the order is changed */
  robustStat.resize(m_numberOfFeature);

  std::sort(samples.begin(), samples.end() );

  double n = samples.size();

  double q1 = n / 4.0;
  double q1_floor;
  double l1 = modf(q1, &q1_floor);

  double q2 = n / 2.0;
  double q2_floor;
  double l2 = modf(q2, &q2_floor);

  double q3 = 3.0 * n / 4.0;
  double q3_floor;
  double l3 = modf(q3, &q3_floor);

  double median = (1 - l2) * samples[static_cast<long>(q2_floor)] + l2 * samples[static_cast<long>(q2_floor) + 1];

  double iqr = ( (1 - l3) * samples[static_cast<long>(q3_floor)] + l3 * samples[static_cast<long>(q3_floor) + 1] ) \
    - ( (1 - l1) * samples[static_cast<long>(q1_floor)] + l1 * samples[static_cast<long>(q1_floor) + 1] );

  robustStat[0] = median;
  robustStat[1] = iqr;

  /* next compute MAD */
  long                nn = samples.size();
  std::vector<double> samplesDeMedian(nn);
  for( long i = 0; i < nn; ++i )
    {
    samplesDeMedian[i] = fabs(samples[i] - median);
    }

  std::sort(samplesDeMedian.begin(), samplesDeMedian.end() );

  double mad =
    (1 - l2) * samplesDeMedian[static_cast<long>(q2_floor)] + l2 * samplesDeMedian[static_cast<long>(q2_floor) + 1];
  robustStat[2] = mad;

  return;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::seedToMask()
{
  if( !(this->mp_img) )
    {
    std::cerr << "Error: set input image first.\n";
    raise(SIGABRT);
    }

  if( this->mp_mask )
    {
    /* Sometimes, the mask is not corresponding to the seed, like
       using the mean shape as the mask and just some other seed as
       the samples for feature. In such cases, do not touch mask
       and just return. */

    return;
    }

  long n = m_seeds.size();
  if( n == 0 )
    {
    std::cerr << "Error: No seeds specified." << std::endl;
    raise(SIGABRT);
    }

  this->mp_mask = TMaskImage::New();
  this->mp_mask->SetRegions(this->mp_img->GetLargestPossibleRegion() );
  this->mp_mask->Allocate();
  this->mp_mask->CopyInformation(this->mp_img);
  this->mp_mask->FillBuffer(0);
  for( long i = 0; i < n; ++i )
    {
    if( 3 != m_seeds[i].size() )
      {
      std::cerr << "Error: 3 != m_seeds[i].size()\n";
      raise(SIGABRT);
      }

    long ix = m_seeds[i][0];
    long iy = m_seeds[i][1];
    long iz = m_seeds[i][2];
    for( long iiz = iz - 1; iiz <= iz + 1; ++iiz )
      {
      for( long iiy = iy - 1; iiy <= iy + 1; ++iiy )
        {
        for( long iix = ix - 1; iix <= ix + 1; ++iix )
          {
          if( 0 <= iix && iix < this->m_nx    \
              && 0 <= iiy && iiy < this->m_ny    \
              && 0 <= iiz && iiz < this->m_nz )
            {
            TIndex idx = {{iix, iiy, iiz}};
            this->mp_mask->SetPixel(idx, 1);
            }
          }
        }
      }
    }

  return;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::dialteSeeds()
{
  /* For each seed, add its 26 neighbors into the seed list. */

  if( !(this->mp_img) )
    {
    std::cerr << "Error: set input image first.\n";
    raise(SIGABRT);
    }

  long                            n = m_seeds.size();
  std::vector<std::vector<long> > newSeeds;

  if( n == 0 )
    {
    std::cerr << "Error: No seeds specified." << std::endl;
    raise(SIGABRT);
    }
  for( long i = 0; i < n; ++i )
    {
    if( 3 != m_seeds[i].size() )
      {
      std::cerr << "Error: 3 != m_seeds[i].size()\n";
      raise(SIGABRT);
      }

    long ix = m_seeds[i][0];
    long iy = m_seeds[i][1];
    long iz = m_seeds[i][2];
    for( long iiz = iz - 1; iiz <= iz + 1; ++iiz )
      {
      for( long iiy = iy - 1; iiy <= iy + 1; ++iiy )
        {
        for( long iix = ix - 1; iix <= ix + 1; ++iix )
          {
          if( 0 <= iix && iix < this->m_nx    \
              && 0 <= iiy && iiy < this->m_ny    \
              && 0 <= iiz && iiz < this->m_nz )
            {
            /* Some locations may be added multiple times,
               if the original seeds are close. But I think
               this is fine */

            std::vector<long> s(3);
            s[0] = iix;
            s[1] = iiy;
            s[2] = iiz;

            newSeeds.push_back(s);
            }
          }
        }
      }
    }

  m_seeds.assign(newSeeds.begin(), newSeeds.end() );

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::getFeatureAroundSeeds()
{
  if( !m_featureImageList[m_numberOfFeature - 1] )
    {
    // last feature image is not constructed
    std::cerr << "Error: construct feature images first.\n";
    raise(SIGABRT);
    }

  long n = m_seeds.size();
  if( n == 0 )
    {
    std::cerr << "Error: No seeds specified." << std::endl;
    raise(SIGABRT);
    }

//   short ax = 0;
//   short ay = 0;
//   short az = 0;

// #ifndef NDEBUG
//   std::ofstream intensityAtSeeds("/tmp/intenSeeds.txt");
// #endif
  for( long i = 0; i < n; ++i )
    {
    if( 3 != m_seeds[i].size() )
      {
      std::cerr << "Error: 3 != m_seeds[i].size()\n";
      raise(SIGABRT);
      }

    long ix = m_seeds[i][0];
    long iy = m_seeds[i][1];
    long iz = m_seeds[i][2];

    TIndex idx = {{ix, iy, iz}};

    std::vector<double> featureHere(m_numberOfFeature);
    computeFeatureAt(idx, featureHere);

    m_featureAtTheSeeds.push_back(featureHere);

// #ifndef NDEBUG
//       //intensityAtSeeds<<"inensity at ( "<<ix<<", "<<iy<<", "<<iz<<") is "<<this->mp_img->GetPixel(idx)<<std::endl;
//       intensityAtSeeds<<this->mp_img->GetPixel(idx)<<std::endl;
// #endif
    }

// #ifndef NDEBUG
//   intensityAtSeeds.close();
// #endif

  return;
}

/* ============================================================ */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::estimateFeatureStdDevs()
{
  m_kernelStddev.assign(m_numberOfFeature, 0.0);

  long n = m_seeds.size(); // == m_featureAtTheSeeds.size()
  for( long i = 0; i < m_numberOfFeature; ++i )
    {
    double m = 0;
    for( long ii = 0; ii < n; ++ii )
      {
      m += m_featureAtTheSeeds[ii][i];
      }
    m /= n;
    for( long ii = 0; ii < n; ++ii )
      {
      m_kernelStddev[i] += (m_featureAtTheSeeds[ii][i] - m) * (m_featureAtTheSeeds[ii][i] - m);
      }

    m_kernelStddev[i] /= (n - 1);
    m_kernelStddev[i] = sqrt(m_kernelStddev[i]);
    }

// // #ifndef NDEBUG
// //   std::ofstream dbgf("/tmp/dbgo.txt", std::ios_base::app);
//   for (long i = 0; i < m_numberOfFeature; ++i)
//     {
//       //dbgf<<"Feature "<<i<<" has var = "<<m_kernelStddev[i]<<std::endl;
//       std::cout<<"Feature "<<i<<" has var = "<<m_kernelStddev[i]<<std::endl;
//     }
// // #endif

  return;
}

template <typename TPixel>
double
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::kernelEvaluationUsingPDF(const std::vector<double>& newFeature)
{
  double p = 1;

  for( long i = 0; i < m_numberOfFeature; ++i )
    {
    long idx = static_cast<long>(newFeature[i] - m_inputImageIntensityMin);

    double probOfThisFeature = m_PDFlearnedFromSeeds[i][idx];

    p *= probOfThisFeature;
    }

  return p;
}

/* ============================================================  */
template <typename TPixel>
double
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::kernelEvaluation(const std::vector<double>& newFeature)
{
  long n = m_seeds.size(); // == m_featureAtTheSeeds.size()

  double p = 1;

  // double p = 0;
  for( long i = 0; i < m_numberOfFeature; ++i )
    {
    double pp = 0.0;

    double stdDev = m_kernelStddev[i] / m_kernelWidthFactor; // /10 as in Eric's appendix

    double var2 = -1.0 / (2 * stdDev * stdDev);
    double c = 1.0 / sqrt(2 * (vnl_math::pi) ) / stdDev;
    for( long ii = 0; ii < n; ++ii )
      {
      pp += exp(var2 * (newFeature[i] - m_featureAtTheSeeds[ii][i]) * (newFeature[i] - m_featureAtTheSeeds[ii][i]) );
      }

    pp *= c;
    pp /= n;

    p *= pp;
    // p = p>pp?p:pp;
    }

  return p;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::setKernelWidthFactor(double f)
{
  if( f < 0.1 )
    {
    m_kernelWidthFactor = 0.1;
    }

  if( f > 20.0 )
    {
    m_kernelWidthFactor = 20.0;
    }

  m_kernelWidthFactor = f;

//   std::ofstream fil("/tmp/d.txt", std::ios_base::app);
//   fil<<"m_kernelWidthFactor = "<<m_kernelWidthFactor<<std::endl;
//   fil.close();

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::setIntensityHomogeneity(double h)
{
//   std::ofstream fil("/tmp/d.txt", std::ios_base::app);
//   fil<<"intensity homogeneity = "<<h<<std::endl;
//   fil.close();

  double f = h * (20.0 - 0.1) + 0.1;

  setKernelWidthFactor(f);

  return;
}

/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::SetZMReference(const char* fnZMRef){

  m_fnZMRef = fnZMRef;

  LoadZMReference();

  return;
}

/* ============================================================  */
template< typename TPixel >
void
CSFLSRobustStatSegmentor3DLabelMap< TPixel >
::SetRadVolPrior(const char* fn){

  m_fnRadVolPrior = fn;

  LoadRadVolPrior();

  return;
}


/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::estimatePDFs()
{
  m_PDFlearnedFromSeeds.clear();

  computeMinMax(); // so we have the range of all pdfs

// #ifndef NDEBUG
//   std::cout<<"m_inputImageIntensityMin = "<<m_inputImageIntensityMin<<std::endl;
//   std::cout<<"m_inputImageIntensityMax = "<<m_inputImageIntensityMax<<std::endl;
// #endif

  long n = m_seeds.size();
  for( long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature )
    {
    std::vector<double> thisPDF(m_inputImageIntensityMax - m_inputImageIntensityMin + 1);
    // assumption: TPixel are of integer types.

    double stdDev = m_kernelStddev[ifeature] / m_kernelWidthFactor; // /10 as in Eric's appendix

    // std::cout<<"kernel sigma of "<<ifeature<<"-th feature is "<<stdDev<<std::endl;

// #ifndef NDEBUG
//       std::cout<<"stdDev of "<<ifeature<<"-th feature = "<<stdDev<<std::endl;
// #endif

//       //#ifndef NDEBUG
//       std::ofstream df("/tmp/detail.txt");
//       //#endif

    double var2 = -1.0 / (2 * stdDev * stdDev);
    double c = 1.0 / sqrt(2 * (vnl_math::pi) ) / stdDev;
    for( TPixel a = m_inputImageIntensityMin; a <= m_inputImageIntensityMax; ++a )
      {
      long ia = static_cast<long>(a - m_inputImageIntensityMin);

      double pp = 0.0;
      for( long ii = 0; ii < n; ++ii )
        {
        pp += exp(var2 * (a - m_featureAtTheSeeds[ii][ifeature]) * (a - m_featureAtTheSeeds[ii][ifeature]) );

        }

      pp *= c;
      pp /= n;

      thisPDF[ia] = pp;
      }

    m_PDFlearnedFromSeeds.push_back(thisPDF);
    }

  return;
}

/* ============================================================  */
template <typename TPixel>
void
CSFLSRobustStatSegmentor3DLabelMap<TPixel>
::computeMinMax()
{
  if( !(this->mp_img) )
    {
    std::cerr << "Error: set input image first.\n";
    raise(SIGABRT);
    }

  typedef itk::Image<TPixel, 3> itkImage_t;

  typedef itk::ImageRegionConstIterator<itkImage_t> itkImageRegionConstIterator_t;

  itkImageRegionConstIterator_t it( (this->mp_img), (this->mp_img)->GetLargestPossibleRegion() );
  it.GoToBegin();

  m_inputImageIntensityMin = std::numeric_limits<unsigned>::max(); // yes, it's twisted so easity to compute.
  m_inputImageIntensityMax = std::numeric_limits<unsigned>::min();
  for( ; !it.IsAtEnd(); ++it )
    {
    TPixel v = it.Get();

    m_inputImageIntensityMin = m_inputImageIntensityMin < v ? m_inputImageIntensityMin : v;
    m_inputImageIntensityMax = m_inputImageIntensityMax > v ? m_inputImageIntensityMax : v;
    }

  return;
}

#endif
