#include "LASeg.h"

namespace LASeg {

/************************************************************


************************************************************/
template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(const char* fnSrc, const char* fnLab, const char* fnZM, double wIntensityIni, double wCurvature, \
                double wRSS, double exVolume, int nIterAll, int nIterRSS,std::string fnInter){

    typedef CSFLSRobustStatSegmentor3DLabelMap< PixelType > SFLSRobustStatSegmentor3DLabelMap_c;

    // Read raw image
    typedef typename SFLSRobustStatSegmentor3DLabelMap_c::TImage Image_t;
    typename Image_t::Pointer imSrc = LASeg::readImage<Image_t>(fnSrc);
//    // Create an image
//    typename Image_t::Pointer imSrc = Image_t::New();
//    typename Image_t::RegionType region;
//    typename Image_t::IndexType start;
//     start[0] = 0;
//     start[1] = 0;
//     start[2] = 0;

//     typename Image_t::SizeType size;
//     unsigned int NumRows = 120;
//     unsigned int NumCols = 120;
//     unsigned int NumDeps = 120;
//     size[0] = NumRows;
//     size[1] = NumCols;
//     size[2] = NumDeps;

//     region.SetSize(size);
//     region.SetIndex(start);

//     imSrc->SetRegions(region);
//     imSrc->Allocate();
//     imSrc->FillBuffer(1);

    // Read input label image
    typedef typename SFLSRobustStatSegmentor3DLabelMap_c::TLabelImage LabelImage_t;
    typename LabelImage_t::Pointer imLab = LASeg::readImage<LabelImage_t>(fnLab);

    // Do segmentation
    SFLSRobustStatSegmentor3DLabelMap_c seg;
    seg.setImage(imSrc);
    seg.setInputLabelImage(imLab);
    seg.SetZMReference(fnZM);
    seg.SetSaveInter(fnInter);

    seg.setNumIter(nIterAll); // a large enough number, s.t. will not be stoped by this creteria.
    seg.SetRSSNumIter(nIterRSS);

    double expectedVolume = exVolume; // the mean volume of the training pre-endo is 108
    seg.setMaxVolume(expectedVolume);

    double maxRunningTime = 30; // max running time 10 min
    seg.setMaxRunningTime(maxRunningTime);

//    double curvatureWeight = 0.2;
    seg.setIntensityHomogeneity(wIntensityIni);
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

    // Get final mask
//    typename itk::Image<PixelType, 3>::Pointer imSeg = seg.mp_phi;
//    typename ImMaskType::Pointer imMask = ImMaskType::New();
//    int labOutValue = 1;

//    typename ImMaskType::SizeType size = imSeg->GetLargestPossibleRegion().GetSize();;
//    typename ImMaskType::IndexType start = {{0, 0, 0}};
//    long nx = size[0];
//    long ny = size[1];
//    long nz = size[2];

//    typename ImMaskType::RegionType region;
//    region.SetSize( size );
//    region.SetIndex( start );

//    imMask->SetRegions( region );
//    imMask->SetSpacing(imSeg->GetSpacing());
//    imMask->SetOrigin(imSeg->GetOrigin());
//    imMask->Allocate();
//    imMask->FillBuffer(0);

//    for (long ix = 0; ix < nx; ++ix)
//      {
//        for (long iy = 0; iy < ny; ++iy)
//          {
//            for (long iz = 0; iz < nz; ++iz)
//              {
//                typename ImMaskType::IndexType idx = {{ix, iy, iz}};
//                PixelType v = imSeg->GetPixel(idx);

//                imMask->SetPixel(idx, v<= threshold?labOutValue:0);
//              }
//          }
//      }

//    return imMask;

    return NULL;
}

/************************************************************


************************************************************/
template<typename PixelType, typename ImMaskType>
typename ImMaskType::Pointer
SegWithRSSandZM(typename itk::Image<PixelType, 3>::Pointer imSrc, typename ImMaskType::Pointer imLab, const char* fnZM, double seedRad, \
                double wIntensityIni, double wCurvature, double wRSS, double exVolume, int nIterAll, int nIterRSS,std::string fnInter){

    typedef CSFLSRobustStatSegmentor3DLabelMap< PixelType > SFLSRobustStatSegmentor3DLabelMap_c;

 //     LASeg::writeImage<Image_t>(imSrc, "../../Results/Test/UniformSrc.nrrd");
    std::string fnRSSZMCofs =  std::string(Prior3D_DIR) + "RSSZMCofs.txt"; //"../LASeg/config/fnRSSZMCofs.txt";
    // Do segmentation
    SFLSRobustStatSegmentor3DLabelMap_c seg;
    seg.setImage(imSrc);
    seg.setInputLabelImage(imLab);
    seg.SetZMReference(fnZM);
    seg.SetRSSZMCofs(fnRSSZMCofs.c_str());
    seg.setSeedRad(seedRad);
    seg.SetSaveInter(fnInter);

    seg.setNumIter(nIterAll); // a large enough number, s.t. will not be stoped by this creteria.
    seg.SetRSSNumIter(nIterRSS);

    double expectedVolume = exVolume; // the mean volume of the training pre-endo is 108
    seg.setMaxVolume(expectedVolume);

    double maxRunningTime = 30; // max running time 30 min
    seg.setMaxRunningTime(maxRunningTime);

    seg.setIntensityHomogeneity(wIntensityIni);
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

    // Get final mask
//    typename itk::Image<PixelType, 3>::Pointer imSeg = seg.mp_phi;
//    typename ImMaskType::Pointer imMask = ImMaskType::New();
//    int labOutValue = 1;

//    typename ImMaskType::SizeType size = imSeg->GetLargestPossibleRegion().GetSize();;
//    typename ImMaskType::IndexType start = {{0, 0, 0}};
//    long nx = size[0];
//    long ny = size[1];
//    long nz = size[2];

//    typename ImMaskType::RegionType region;
//    region.SetSize( size );
//    region.SetIndex( start );

//    imMask->SetRegions( region );
//    imMask->SetSpacing(imSeg->GetSpacing());
//    imMask->SetOrigin(imSeg->GetOrigin());
//    imMask->Allocate();
//    imMask->FillBuffer(0);

//    for (long ix = 0; ix < nx; ++ix)
//      {
//        for (long iy = 0; iy < ny; ++iy)
//          {
//            for (long iz = 0; iz < nz; ++iz)
//              {
//                typename ImMaskType::IndexType idx = {{ix, iy, iz}};
//                PixelType v = imSeg->GetPixel(idx);

//                imMask->SetPixel(idx, v<= threshold?labOutValue:0);
//              }
//          }
//      }

//    return imMask;

    return NULL;
}

/************************************************************


************************************************************/
void RSSWithShapeConstraint(const char* fnSrc, const char* fnSeed, const char* fnLA, float laSeedRad) {

    std::string fnRSSZMCof, fnLAZMPrior, fnRadVolPrior;
    fnRSSZMCof = std::string(Prior3D_DIR) + "RSSZMCofs.txt";
    fnLAZMPrior = std::string(Prior3D_DIR) + "LAZMPrior.txt";
    fnRadVolPrior = std::string(Prior3D_DIR) + "RadVolPrior.txt";

    typedef float PixelType;
    typedef itk::Image< PixelType, 3 > ImSrcType;
    typedef itk::Image< char, 3 > ImMaskType;
    std::vector<std::string> fn_vec;
    std::stringstream ss;
    ImSrcType::Pointer imSrc, imSrcSmooth;
    ImMaskType::Pointer imLab, imLA;

    double exVolRSS, exVolRSSZM, wCurvatureRSS, wIntensityRSS, wCurvatureRSSZM, wIntensityRSSZM, wRSS;
    int iterRSS, iter;
    std::string strID;

    fn_vec = LASeg::readTextLineToListOfString<char>(fnRSSZMCof.c_str());
    ss.clear();
    ss.str(fn_vec[0]);
    ss >> strID >> wIntensityRSS >> wIntensityRSSZM;

    ss.clear();
    ss.str(fn_vec[1]);
    ss >> strID >> wCurvatureRSS >> wCurvatureRSSZM;

    ss.clear();
    ss.str(fn_vec[2]);
    ss >> strID >> wRSS;

    ss.clear();
    ss.str(fn_vec[3]);
    ss >> strID >> exVolRSS >> exVolRSSZM;

    ss.clear();
    ss.str(fn_vec[4]);
    ss >> strID >> iterRSS >>  iter;

    imSrc = LASeg::readImage<ImSrcType>(fnSrc);
    imLab = LASeg::readImage<ImMaskType>(fnSeed);

    std::cout << "load src and labe\n";

    imSrcSmooth = LASeg::AnisotropicSmoothing<ImSrcType>(imSrc);
    std::cout << "smooth image\n";

    // resample image for 1: obtaining isotropic image 2: reducing computational load
    double spacing[] = {0.5, 0.5, 0.5};
    imSrcSmooth = LASeg::ResampleImage<ImSrcType, ImSrcType>(imSrcSmooth, spacing, 1);
    imLab = LASeg::ResampleImage<ImMaskType, ImMaskType>(imLab, spacing, 0);

//    LASeg::writeImage<ImMaskType>(imLab, "../results/labSmall.nrrd");

    imLA = LASeg::SegWithRSSandZM<PixelType, ImMaskType>(imSrcSmooth, imLab, fnLAZMPrior.c_str(), laSeedRad, wIntensityRSS, wCurvatureRSS, wRSS, exVolRSS, iterRSS, iterRSS);
    imLA = LASeg::SegWithRSSandZM<PixelType, ImMaskType>(imSrcSmooth, imLA, fnLAZMPrior.c_str(), laSeedRad, wIntensityRSSZM, wCurvatureRSSZM, wRSS, exVolRSSZM, iter, 0);
    std::cout << "seg with ZM\n";

//     LASeg::writeImage<ImMaskType>(imLA, "../results/laSmall.nrrd");
     double spacingOri[3];
     spacingOri[0] =  imSrc->GetSpacing()[0];
     spacingOri[1] =  imSrc->GetSpacing()[1];
     spacingOri[2] =  imSrc->GetSpacing()[2];
     imLA = LASeg::ResampleImage<ImMaskType, ImMaskType>(imLA, spacingOri, 0);

    LASeg::writeImage<ImMaskType>(imLA, fnLA);
}

/***************************************************************
        Compute Zernike Moments for binary mask(s)
***************************************************************/
void ComputeZernikeMoments() {

    const char* fnZM = "../LASeg/config/fnZMList.txt";
    std::vector<std::string> fn_vec;
    std::vector<float> zm_mean, zm_vec;
    std::stringstream ss;
    std::string dirSrc, dirSave, preSrc, sufSrc, preSave, sufSave, fnSrc, fnSave, fnSaveMoments, strOrder, strID;
    int order;
    float scaleRecon;

    fn_vec = LASeg::readTextLineToListOfString<char>(fnZM);
    ss.str(fn_vec[0]);
    ss >> strID >> dirSrc;

    ss.clear();
    ss.str(fn_vec[1]);
    ss >> strID >> dirSave;

    ss.clear();
    ss.str(fn_vec[2]);
    ss >> strID >> preSrc;

    ss.clear();
    ss.str(fn_vec[3]);
    ss >> strID >> sufSrc;

    ss.clear();
    ss.str(fn_vec[4]);
    ss >> strID >> preSave;

    ss.clear();
    ss.str(fn_vec[5]);
    ss >> strID >> sufSave;

    ss.clear();
    ss.str(fn_vec[6]);
    ss >> strID >> order;

    ss.clear();
    ss.str(fn_vec[7]);
    ss >> strID >> scaleRecon;

    ss.clear();
    ss.str("");
    ss << order;
    strOrder = ss.str();

    typedef unsigned char TMaskPixelType;
    typedef itk::Image<TMaskPixelType, 3> ImMaskType;
    ImMaskType::Pointer imMask;
    int dim;
    unsigned int count = 0;

    for(unsigned int i = 8; i < fn_vec.size(); i++) {

        count++;

        fnSrc = dirSrc + preSrc + fn_vec[i] + sufSrc + ".nrrd";
//        fnSave = dirSave + preSave + fn_vec[i] + sufSave + ".nrrd";
//        fnSrc = dirSrc + fn_vec[i] + "/" + preSrc + fn_vec[i].substr(1,2) + ".nrrd";

        imMask = LASeg::readImage<ImMaskType>(fnSrc.c_str());

        float* pData = NULL;
        pData = LASeg::ExtractCubicNonZeroRegion<ImMaskType>(imMask, dim);

        ZernikeDescriptor<float, float> zd;
        zd.ComputeZernikeMoments(pData, dim, order);

        fnSave = dirSave + preSave + fn_vec[i] + sufSave + ".txt";
        zd.SaveMoments(fnSave.c_str());

        zd.GetZernikeMoments(zm_vec);

//        double ellapsedTime = (clock() - startingTime) / static_cast<double>(CLOCKS_PER_SEC);
//        std::cout << "runing time: " << ellapsedTime << std::endl;
        if(i == 8) {
            zm_mean.insert(zm_mean.begin(), zm_vec.begin(), zm_vec.end());
        }
        else {
            for(unsigned int j = 0; j < zm_vec.size(); j++) {
                zm_mean[j] += zm_vec[j];
            }
        }

 //        const int dimRecon = int(dim*scaleRecon);
//        typedef std::complex<float>      ComplexT;
//        vector<vector<vector<ComplexT> > > grid(dimRecon, vector<vector<ComplexT> >(dimRecon, vector<ComplexT>(dimRecon)));

//    //    zd.Reconstruct(grid, 0, order, 0, order);
//        zd.Reconstruct(grid);

////        std::string fnSave = "../../Results/Test/zmTest.nrrd";
////        LASeg::WriteArray2Volume<float, ComplexT>(grid, fnSave);

//        // Save results
//        if(dirSave != "") {

//            fnSave = dirSave + preSrc + fn_vec[i] + sufSrc + "_" + strOrder + "_ZM.nrrd";
//            LASeg::WriteArray2Volume<float, ComplexT>(grid, fnSave);

//        }

//        delete [] data;

        // Clean up allocation
        if(pData != NULL) {
            delete []pData;
            pData = NULL;
        }
    }

    // Compute the mean ZM
    for(unsigned int i = 0; i < zm_mean.size(); i++) {
        zm_mean[i] /= count;
    }

    // Save the mean ZM
    fnSave = dirSave + preSave + sufSave + ".txt";
    std::ofstream outfile(fnSave.c_str());
    outfile << order << std::endl;
    for(unsigned int i = 0; i < zm_mean.size(); i += 2) {
        outfile << zm_mean[i] << " " << zm_mean[i + 1] << std::endl;
    }
    outfile.close();
}

///**********************************************************
//            Possion Surface Reconstruction
//***********************************************************/
//template<typename ImMaskType>
//vtkPolyData* PossionSurfaceRecon(typename ImMaskType::Pointer imMask, float depth,  float deciReduction, float deciError, \
//                                 float scale, float solDivide, float isoDivide , float samNode) {

//    typedef itk::Image<float, 3> InteralImageType;

//    typedef itk::BinaryThresholdImageFilter<ImMaskType, InteralImageType >  BinaryTHFilterType;
//    typename BinaryTHFilterType::Pointer binaryTHFilter = BinaryTHFilterType::New();
//    binaryTHFilter->SetInput(imMask);
//    binaryTHFilter->SetInsideValue(1.0f);
//    binaryTHFilter->SetOutsideValue(-1.0f);
//    binaryTHFilter->SetLowerThreshold(1);
//    binaryTHFilter->SetUpperThreshold(1);
//    binaryTHFilter->Update();

////    imMask = binaryTHFilter->GetOutput();

//    /**Convert Itk Image to VtkData*/
//    typedef itk::ImageToVTKImageFilter<InteralImageType> ConnectorType;
//    typename ConnectorType::Pointer connector = ConnectorType::New();
//    connector->SetInput(binaryTHFilter->GetOutput());


//    vtkImageGaussianSmooth *m_VTKGaussianFilter;
//    // Initialize the Gaussian filter
//    m_VTKGaussianFilter = vtkImageGaussianSmooth::New();
//    m_VTKGaussianFilter->ReleaseDataFlagOn();

////    float sigma = 0.8f;
////    const double *spacing = imMask->GetSpacing().GetDataPointer();
////    m_VTKGaussianFilter->SetInput(connector->GetOutput());
////    m_VTKGaussianFilter->SetStandardDeviation(sigma / spacing[0], sigma / spacing[1], sigma / spacing[2]);
////    m_VTKGaussianFilter->SetRadiusFactors(3 * sigma / spacing[0], 3 * sigma / spacing[1], 3 * sigma / spacing[2]);
////    m_VTKGaussianFilter->Update();

////    vtkImageData* imVTKMaskSmooth = m_VTKGaussianFilter->GetOutput();

//    // Create and configure the marching cubes filter
//    vtkMarchingCubes* m_MarchingCubesFilter;
//    m_MarchingCubesFilter = vtkMarchingCubes::New();
//    m_MarchingCubesFilter->SetInput(connector->GetOutput());
//    m_MarchingCubesFilter->ReleaseDataFlagOn();
//    m_MarchingCubesFilter->ComputeScalarsOff();
//    m_MarchingCubesFilter->ComputeGradientsOff();
//    m_MarchingCubesFilter->SetNumberOfContours(1);
//    m_MarchingCubesFilter->ComputeNormalsOn();
//    m_MarchingCubesFilter->SetValue(0,0);
//    m_MarchingCubesFilter->Update();

//    vtkSmartPointer<vtkPoissonReconstruction> poissonFilter =  vtkSmartPointer<vtkPoissonReconstruction>::New();
//    poissonFilter->SetDepth( depth );
//    poissonFilter->SetScale(scale);
//    poissonFilter->SetSolverDivide(solDivide);
//    poissonFilter->SetIsoDivide(isoDivide);
//    poissonFilter->SetSamplesPerNode(samNode);
//    poissonFilter->SetInputConnection(m_MarchingCubesFilter->GetOutputPort());
//    poissonFilter->Update();

//    vtkSmartPointer<vtkDecimatePro> decimatePro = vtkSmartPointer<vtkDecimatePro>::New();
//    decimatePro->SetInput(poissonFilter->GetOutput());
//    decimatePro->SetTargetReduction(deciReduction);
//    decimatePro->SetMaximumError(deciError);
//    decimatePro->SetFeatureAngle(45);
//    decimatePro->PreserveTopologyOn();
//    decimatePro->Update();

//    vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
//    triFilter->SetInput(decimatePro->GetOutput());
//    triFilter->Update();

//    // Adjust the image origin to line up the coordinates difference between ITK and VTK
//    vtkPolyData* surface = triFilter->GetOutput();
//    vtkPoints* points = surface->GetPoints();
//    double* point;
//    for(unsigned int i = 0; i < points->GetNumberOfPoints(); i++) {
//        point = points->GetPoint(i);
//        point[0] = -point[0];
//        point[1] = -point[1];

//        points->SetPoint(i, point);
//    }

//    return triFilter->GetOutput();
//}

/****************************************************************
            Convert VTK to PLY
*****************************************************************/
vtkPolyData* VTK2PolyWithLargestConnectedComponent(vtkPolyData* surface, std::string fnPolySave) {

    // Triangulate VTK data
    vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
    triFilter->SetInput(surface);
    triFilter->Update();

    // Smooth Data
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
    vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(triFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(10);
    smoothFilter->Update();

    // Get Largest Connected Component
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter =
        vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(smoothFilter->GetOutputPort());
    connectivityFilter->SetExtractionModeToLargestRegion();
    connectivityFilter->Update();

    // Save File
    if(fnPolySave != "") {
        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
        writer->SetInput(connectivityFilter->GetOutput());
        writer->SetFileName(fnPolySave.c_str());
        writer->Update();

        return NULL;
    }
    else {
        return connectivityFilter->GetOutput();
    }
}

//template< typename ImMaskType>
//vtkPolyData* ExtractZernikeSurface(typename ImMaskType::Pointer imMask, float depth, std::string fnSave) {

//    // Theshold reconstructed Zernike image
//    typedef itk::OtsuThresholdImageFilter <ImMaskType, ImMaskType> OtsuThresholdImageFilterType;
//    typename OtsuThresholdImageFilterType::Pointer otsuFilter  = OtsuThresholdImageFilterType::New();
//    otsuFilter->SetInput(imMask);
//    otsuFilter->SetInsideValue(0);
//    otsuFilter->SetOutsideValue(1);
//    otsuFilter->Update();


//    std::string fnVol = LASeg::ExtractDirectory(fnSave) + LASeg::ExtractFilenameWithoutExtension(fnSave) + "_Otsu.nrrd";
//    LASeg::writeImage<ImMaskType>(otsuFilter->GetOutput(), fnVol.c_str());

//    vtkPolyData* surface;
//    surface = LASeg::PossionSurfaceRecon<ImMaskType>(otsuFilter->GetOutput(), depth);

////    vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
////    writer->SetInput(surface);
////    writer->SetFileName("../../Results/Test/possion.ply");
////    writer->Update();

//    surface = LASeg::VTK2PolyWithLargestConnectedComponent(surface, fnSave);

//    if(fnSave != "" && surface != NULL) {
//        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
//        writer->SetInput(surface);
//        writer->SetFileName(fnSave.c_str());
//        writer->Update();

//        return NULL;
//    }
//    else {
//        return surface;
//    }
//}

//void ExtractZernikeSurfaces() {

//    const char* fnZM = "../LASeg/config/fnZMSurfaceList.txt";
//    std::vector<std::string> fn_vec;
//    std::string dirSrc, dirSave, sufSrc, sufSave, fnSrc, fnSave, strID;
//    std::stringstream ss;
//    float depth;

//    typedef itk::Image<unsigned char, 3> ImMaskType;
//    ImMaskType::Pointer imMask;

//    fn_vec = LASeg::readTextLineToListOfString<char>(fnZM);
//    ss.str(fn_vec[0]);
//    ss >> strID >> dirSrc;

//    ss.clear();
//    ss.str(fn_vec[1]);
//    ss >> strID >> dirSave;

//    ss.clear();
//    ss.str(fn_vec[2]);
//    ss >> strID >> sufSrc;

//    ss.clear();
//    ss.str(fn_vec[3]);
//    ss >> strID >> sufSave;

//    ss.clear();
//    ss.str(fn_vec[4]);
//    ss >> strID >> depth;

//     for(unsigned int i = 5; i < fn_vec.size(); i++) {

//         fnSrc = dirSrc + fn_vec[i] + sufSrc + ".nrrd";
//         fnSave = dirSave + fn_vec[i] + sufSave + ".vtk";

//         imMask = LASeg::readImage<ImMaskType>(fnSrc.c_str());
//         LASeg::ExtractZernikeSurface<ImMaskType>(imMask, depth, fnSave);

//     }
//}

/*
            Convert VTK to PLY with Largest Connected Component
*/
void VTK2PLYWithLargestLargestConnectedComponent() {

    const char* fnSrcList = "../LASeg/config/fnVTKList.txt";
    std::vector<std::string> fn_vec;
    std::stringstream ss;
    std::string dirSave,dirSrc, fn, fnSave, strID;

//        std::ofstream fout("../Results/Test/VTK2PLYTime.m");
//        fout << "LVVTK2PLYTime = [";

    fn_vec = LASeg::readTextLineToListOfString<char>(fnSrcList);
    ss.clear();
    ss.str(fn_vec[0]);
    ss >> strID >> dirSrc;

    ss.clear();
    ss.str(fn_vec[1]);
    ss >> strID >> dirSave;

    for(unsigned int i = 2; i < fn_vec.size(); i++) {
        fn = dirSrc + fn_vec[i] + ".vtk";

        if(dirSave != "") {
            fnSave = dirSave + fn_vec[i] + ".ply";
            std::cout << "save " + fnSave << endl;
        }
        else
            fnSave = "";

//            double startingTime = clock();

        // Read VTK file
        vtkSmartPointer<vtkPolyDataReader> reader =  vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(fn.c_str());
        reader->Update();

        LASeg::VTK2PolyWithLargestConnectedComponent(reader->GetOutput(), fnSave.c_str());

//            double ellapsedTime = (clock() - startingTime) / static_cast<double>(CLOCKS_PER_SEC);
//            std::cout << "runing time: " << ellapsedTime << std::endl;

//            fout << ellapsedTime << " ";
    }

//        fout << "];\n";
//        fout << "mean(LVVTK2PLYTime(:))";
//        fout.close();
}

///*****************************************
//            Resample Images
//*****************************************/
//void ReSampleImages() {

//    const char* fnResampleList = "../LASeg/config/fnResampleList.txt";
//    std::string fnSave,fnRef, fnSrc;
//    std::string dirRef, dirSrc, dirSave, sufRef, sufSrc, sufSave, strID;
//    std::vector<std::string> fn_vec;
//    std::stringstream ss;
//    float MAXDIS;
//    int flagMaskRef, flagMaskSrc;

//    fn_vec = LASeg::readTextLineToListOfString<char>(fnResampleList);
//    ss.str(fn_vec[0]);
//    ss >> strID >> dirRef;

////        std::cout << "dirRef: " << dirRef << std::endl;

//    ss.clear();
//    ss.str(fn_vec[1]);
//    ss >> strID >> dirSrc;

//    ss.clear();
//    ss.str(fn_vec[2]);
//    ss >> strID >> dirSave;

//    ss.clear();
//    ss.str(fn_vec[3]);
//    ss >> strID >> sufRef;

//    ss.clear();
//    ss.str(fn_vec[4]);
//    ss >> strID >> sufSrc;

//    ss.clear();
//    ss.str(fn_vec[5]);
//    ss >> strID >> sufSave;

//    ss.clear();
//    ss.str(fn_vec[6]);
//    ss >> strID >> flagMaskRef;

//    ss.clear();
//    ss.str(fn_vec[7]);
//    ss >> strID >> flagMaskSrc;

//    ss.clear();
//    ss.str(fn_vec[8]);
//    ss >> strID >> MAXDIS;

//    typedef short RefPixel;
//    if(flagMaskRef)
//        typedef unsigned char RefPixel;
//    typedef itk::Image<RefPixel, 3> ImRefType;

//    typedef short SrcPixel;
//    if(flagMaskSrc)
//        typedef unsigned char SrcPixel;
//    typedef itk::Image<SrcPixel, 3> ImSrcType;
//    ImSrcType::Pointer imResample;

//    for(unsigned int i = 9; i < fn_vec.size(); i++) {
//        fnRef = dirRef + fn_vec[i] + sufRef + ".nrrd";
//        fnSrc = dirSrc + fn_vec[i]  + sufSrc + ".nrrd";
//        fnSave = dirSave + fn_vec[i]  +  sufSave + ".nrrd";

//        ImRefType::Pointer imRef = ImRefType::New();
//        ImSrcType::Pointer imSrc = ImSrcType::New();

//        imRef = LASeg::readImage<ImRefType>(fnRef.c_str());
//        imSrc = LASeg::readImage<ImSrcType>(fnSrc.c_str());
//        imResample = LASeg::ResampleImage<ImSrcType, ImRefType>(imSrc, imRef,MAXDIS);
//        LASeg::writeImage<ImSrcType>(imResample, fnSave.c_str());
//     }
//}

/********************************************************
        Compute Dice Coefficient between 2 volumes
*********************************************************/
float ComputeDiceForVolumes(const char* fnMaskSrc, const char* fnMaskDst, float smoothFactor, std::string measureType,
                                 const char* fnROI, std::string dirSave) {

    std::string fn, fnSave;
    fn = LASeg::ExtractFilenameWithoutExtension(std::string(fnMaskSrc));

    typedef short MaskPixelType;
    typedef float DistPixelType;
    typedef itk::Image<MaskPixelType, 3> ImMaskType;
    typedef itk::Image<DistPixelType, 3> ImDistType;

    ImMaskType::Pointer imMaskSrc, imMaskROI, imMaskDst;
    ImDistType::Pointer imDst;
    imMaskSrc = LASeg::readImage<ImMaskType>(fnMaskSrc);
    imMaskDst = LASeg::readImage<ImMaskType>(fnMaskDst);

//    // Determine ROI for MaskDst
//    if(smoothFactor > vnl_math::eps) {
//        imDst = ComputeDistanceFieldFromMask<ImMaskType, ImDistType>(imMaskSrc);

//        imMaskROI = ImMaskType::New();
//        imMaskROI->SetRegions(imMaskSrc->GetLargestPossibleRegion());
//        imMaskROI->Allocate();
//        imMaskROI->CopyInformation(imMaskSrc);
//        imMaskROI->FillBuffer(0);

//        typedef itk::ImageRegionConstIterator<ImDistType> ImgRegionConstInterator;
//        typedef itk::ImageRegionIterator<ImMaskType> ImgRegionIterator;
//        ImgRegionConstInterator itConst(imDst, imDst->GetLargestPossibleRegion());
//        ImgRegionIterator itROI(imMaskROI,imMaskROI->GetLargestPossibleRegion());

//        double dis;
//        for(itConst.GoToBegin(), itROI.GoToBegin(); !itConst.IsAtEnd(); ++itROI, ++itConst) {
//            dis = itConst.Get();
//            if(dis < smoothFactor) {
//                itROI.Set(1);
//            }
//        }


//        // Crop MaskDst within ROI
//        typedef itk::AndImageFilter<ImMaskType> AndImageFilterType;
//        AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
//        andFilter->SetInput(0, imMaskROI);
//        andFilter->SetInput(1, imMaskDst);
//        andFilter->Update();
//        imMaskDst = andFilter->GetOutput();

//        if(dirSave != "") {
//            fnSave = dirSave + fn + "ManualSmooth.nrrd";
//            LASeg::writeImage<ImMaskType>(imMaskDst, fnSave.c_str());
////            std::cout << "save... " << fnSave << std::endl;

//        }
//    }

    typedef itk::AndImageFilter<ImMaskType> AndImageFilterType;
    if(fnROI != NULL) {
        ImMaskType::Pointer roiImg;
        roiImg = LASeg::readImage<ImMaskType>(fnROI);

        AndImageFilterType::Pointer andFilterSrc = AndImageFilterType::New();
        andFilterSrc->SetInput(0, imMaskSrc);
        andFilterSrc->SetInput(1, roiImg);
        andFilterSrc->Update();
        imMaskSrc = andFilterSrc->GetOutput();

//        std::string fnMaskROI;
//        fnMaskROI = "../Results/tmp/" + LASeg::ExtractFilenameWithoutExtension(fnMaskSrc) + "_MaskSrc_ROI.nrrd";
//        LASeg::writeImage<ImMaskType>(imMaskSrc, fnMaskROI.c_str());

        AndImageFilterType::Pointer andFilterDst = AndImageFilterType::New();
        andFilterDst->SetInput(0, imMaskDst);
        andFilterDst->SetInput(1, roiImg);
        andFilterDst->Update();
        imMaskDst = andFilterDst->GetOutput();

//        fnMaskROI = "../Results/tmp/" + LASeg::ExtractFilenameWithoutExtension(fnMaskSrc) + "_MaskDst_ROI.nrrd";
//        LASeg::writeImage<ImMaskType>(imMaskDst, fnMaskROI.c_str());
    }


    ImMaskType::Pointer andImg;
    AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
    andFilter->SetInput(0, imMaskSrc);
    andFilter->SetInput(1, imMaskDst);
    andFilter->Update();
    andImg = andFilter->GetOutput();

    typedef itk::OrImageFilter<ImMaskType> OrImageFilterType;
    ImMaskType::Pointer orImg;
    OrImageFilterType::Pointer orFilter = OrImageFilterType::New();
    orFilter->SetInput(0, imMaskSrc);
    orFilter->SetInput(1, imMaskDst);
    orFilter->Update();
    orImg = orFilter->GetOutput();

    //    fnSave = dirSave + fn + "_HDist.nrrd";
    //    LASeg::writeImage<ImDistType>(imDst, fnSave.c_str());
    if(dirSave != "") {
        fnSave = dirSave + fn + "_Dice.nrrd";
//        LASeg::writeImage<ImMaskType>(andImg, fnSave.c_str());
    }

    float volAnd, volOr, volSrc, volDst, dice, volOverlap, volDiff, measure;

    volAnd = LASeg::CalVolume<ImMaskType>(andImg);
    volOr = LASeg::CalVolume<ImMaskType>(orImg);
    volSrc = LASeg::CalVolume<ImMaskType>(imMaskSrc);
    volDst = LASeg::CalVolume<ImMaskType>(imMaskDst);

    dice = 2*volAnd/(volSrc + volDst);
    volOverlap = volAnd/volOr;
    volDiff = 2*std::abs(volSrc - volDst)/(volSrc + volDst);

    measure = 0;
    if(measureType == "Dice") {
        measure = dice;
    }
    else if(measureType == "VolOverlap") {
        measure =  volOverlap;
    }
    else if(measureType == "VolDiff") {
        measure = volDiff;
    }

    return measure;
//    return std::min(volSrc, volDst)/std::max(volSrc, volDst);
}

template<typename LabImType>
double CalVolume(typename LabImType::Pointer maskImg) {

    typename LabImType::SpacingType spacing;
    spacing = maskImg->GetSpacing();

    typedef itk::ImageRegionConstIterator<LabImType> ImgRegionConstInterator;
    ImgRegionConstInterator itConst(maskImg, maskImg->GetLargestPossibleRegion());


    double vol;
    vol = 0.0;
    for(itConst.GoToBegin(); !itConst.IsAtEnd(); ++itConst) {
        if(itConst.Get() > 0) {
            vol++;
        }
    }

    vol *= spacing[0]*spacing[1]*spacing[2]/1000.0;

    return vol;
}

/*
        Compute Volume Distribution
*/
void ComputeVolumeDistribution() {

    typedef short MaskPixelType;
    typedef itk::Image<MaskPixelType, 3> ImMaskType;
    ImMaskType::Pointer imMask;

    const char* fnTrainList = "../LASeg/config/fnTrainList.txt";
    std::string dirSrc,dirSave, fnSrc, fnSave, strID;
    std::string preSrc, sufSrc;
    std::vector<std::string> fn_vec;
    std::stringstream ss;
    double vol, mean, sigma;
    int count;

    fn_vec = LASeg::readTextLineToListOfString<char>(fnTrainList);
    ss.str(fn_vec[0]);
    ss >> strID >> dirSrc;

    ss.clear();
    ss.str(fn_vec[1]);
    ss >> strID >> dirSave;

    ss.clear();
    ss.str(fn_vec[2]);
    ss >> strID >> preSrc >> sufSrc;


    fnSave = dirSave + "LAVolPDF.txt";
    std::ofstream fout(fnSave.c_str());
    fout << "vols\n";
    mean = 0.0;
    sigma = 0.0;
    count = 0;
    for(unsigned int i = 3; i < fn_vec.size(); i++) {

        fnSrc = dirSrc + preSrc + fn_vec[i] + sufSrc + ".nrrd";
        imMask = LASeg::readImage<ImMaskType>(fnSrc.c_str());
        vol = LASeg::CalVolume<ImMaskType>(imMask);

        fout << vol << " ";
        mean += vol;
        sigma += vol*vol;
        count++;
    }
    fout << std::endl;

    mean /= count;
    sigma = std::sqrt(sigma/count - mean*mean);

    fout << "mean " << mean << std::endl;
    fout << "sigma " << sigma << std::endl;

    fout.close();
}
}


