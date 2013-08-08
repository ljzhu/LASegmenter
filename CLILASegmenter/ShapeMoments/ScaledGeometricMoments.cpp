/*
                                                                            
                          3D Zernike Moments                                
    Copyright (C) 2003 by Computer Graphics Group, University of Bonn       
           http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/             
                                                                            
Code by Marcin Novotni:     marcin@cs.uni-bonn.de
       
for more information, see the paper:

@inproceedings{novotni-2003-3d,
    author = {M. Novotni and R. Klein},
    title = {3{D} {Z}ernike Descriptors for Content Based Shape Retrieval},
    booktitle = {The 8th ACM Symposium on Solid Modeling and Applications},
    pages = {216--225},
    year = {2003},
    month = {June},
    institution = {Universit\"{a}t Bonn},
    conference = {The 8th ACM Symposium on Solid Modeling and Applications, June 16-20, Seattle, WA}
}        
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/
#ifndef SCALEDGEOMOMENTS_CPP
#define SCALEDGEOMOMENTS_CPP

#include "ScaledGeometricMoments.h"

template<class VoxelT, class MomentT>
ScaledGeometricalMoments<VoxelT,MomentT>::ScaledGeometricalMoments ()
{
}

template<class VoxelT, class MomentT>
ScaledGeometricalMoments<VoxelT,MomentT>::ScaledGeometricalMoments (
        const VoxelT* _voxels,
        int _xDim, int _yDim, int _zDim, 
        double _xCOG, double _yCOG, double _zCOG,
        double _scale,
        int _maxOrder)
{
   Init (_voxels, _xDim, _yDim, _zDim, _xCOG, _yCOG, _zCOG, _scale, _maxOrder);
}


template<class VoxelT, class MomentT>
ScaledGeometricalMoments<VoxelT,MomentT>::ScaledGeometricalMoments (
        const VoxelT* _voxels,
        int _dim,
        double _xCOG, double _yCOG, double _zCOG,
        double _scale,
        int _maxOrder)
{
   Init (_voxels, _dim, _dim, _dim, _xCOG, _yCOG, _zCOG, _scale, _maxOrder);
}


template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>::Init (
        const VoxelT* _voxels,
        int _xDim, int _yDim, int _zDim, 
        double _xCOG, double _yCOG, double _zCOG,
        double _scale,
        int _maxOrder,
        bool bDerivative)
{
    xDim_ = _xDim;
    yDim_ = _yDim;
    zDim_ = _zDim;

    maxOrder_ = _maxOrder;

    size_t totalSize = xDim_ * yDim_ * zDim_;
    voxels_.resize (totalSize);
    for (int i=0; i<totalSize; ++i)
    {
        voxels_[i] = _voxels[i];
    }

    moments_.resize (maxOrder_ + 1);
    dmoments_.resize(maxOrder_ + 1);
    for (int i=0; i<=maxOrder_; ++i)
    {
        moments_[i].resize (maxOrder_ - i + 1);
        dmoments_[i].resize (maxOrder_ - i + 1);
        for (int j=0; j<=maxOrder_ - i; ++j)
        {
            moments_[i][j].resize (maxOrder_ - i - j + 1);
            dmoments_[i][j].resize (maxOrder_ - i - j + 1);

            for(int k = 0; k <= maxOrder_ - i - j; k++) {
                dmoments_[i][j][k].resize(4);
            }
        }
    }

    ComputeSamples(_xCOG, _yCOG, _zCOG, _scale);

    Compute();

    if(bDerivative)
        ComputeNormalizedMomentsDerivative(_scale);
}   

template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>::SetBasicMomentsXYZ(double K0,double K1, double K2, double K3) {

    K0_ = K0;
    K1_ = K1;
    K2_ = K2;
    K3_ = K3;
}


template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>::ComputeSamples (double _xCOG, double _yCOG, double _zCOG, double _scale)
{
    samples_.resize (3);    // 3 dimensions

    int dim[3];
    dim[0] = xDim_;
    dim[1] = yDim_;
    dim[2] = zDim_;

    double min[3];
    min[0] = (-_xCOG) * _scale;
    min[1] = (-_yCOG) * _scale;
    min[2] = (-_zCOG) * _scale;

    for (int i=0; i<3; ++i)
    {
        samples_[i].resize (dim[i]+1);
        for (int j=0; j<=dim[i]; ++j)
        {
            samples_[i][j] = min[i] + j * _scale;
        }
    }      
}


template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>::Compute ()
{
    int arrayDim = zDim_;
    int layerDim = yDim_ * zDim_;

    int diffArrayDim =  zDim_ + 1;
    int diffLayerDim = (yDim_ + 1) * zDim_;
    int diffGridDim  = (xDim_ + 1) * layerDim;

    T1D diffGrid (diffGridDim);
    T1D diffLayer (diffLayerDim);
    T1D diffArray (diffArrayDim);

    T1D layer (layerDim);
    T1D array (arrayDim);
    T   moment;

    typename T1D::iterator iter = voxels_.begin ();
    typename T1D::iterator diffIter = diffGrid.begin ();

    // generate the diff version of the voxel grid in x direction
    for (int x=0; x<layerDim; ++x)
    {
        ComputeDiffFunction (iter, diffIter, xDim_);

        iter += xDim_;
        diffIter += xDim_ + 1;
    }

    for (int i=0; i<=maxOrder_; ++i)
    {
        diffIter = diffGrid.begin ();
        for (int p=0; p<layerDim; ++p)
        {
            // multiply the diff function with the sample values
            T1DIter sampleIter (samples_[0].begin ()); 
            layer[p] = Multiply (diffIter, sampleIter, xDim_ + 1);

            diffIter += xDim_ + 1;
        }              

        iter = layer.begin ();
        diffIter = diffLayer.begin ();
        for (int y=0; y<arrayDim; ++y)
        {
            ComputeDiffFunction (iter, diffIter, yDim_);

            iter += yDim_;
            diffIter += yDim_ + 1;
        }

        for (int j=0; j<maxOrder_+1-i; ++j)
        {
            diffIter = diffLayer.begin ();
            for (int p=0; p<arrayDim; ++p)
            {
                T1DIter sampleIter (samples_[1].begin ()); 
                array[p] = Multiply (diffIter, sampleIter, yDim_ + 1);

                diffIter += yDim_ + 1;
            }

            iter = array.begin ();
            diffIter = diffArray.begin ();
            ComputeDiffFunction (iter, diffIter, zDim_);

            for (int k=0; k<maxOrder_+1-i-j; ++k)
            {
                T1DIter sampleIter (samples_[2].begin ()); 

                moment = Multiply (diffIter, sampleIter, zDim_ + 1);
                moments_[i][j][k] = moment / ((1+i) * (1+j) * (1+k));
            }
        }
    }
}

template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>::ComputeNormalizedMomentsDerivative(double _scale) {


    T eta100, eta010, eta001, eta200, eta020, eta002;
    T R3, K03;

    eta100 = moments_[1][0][0];
    eta010 = moments_[0][1][0];
    eta001 = moments_[0][0][1];
    eta200 = moments_[2][0][0];
    eta020 = moments_[0][2][0];
    eta002 = moments_[0][0][2];
    R3 = std::pow(_scale, 3);
    K03 = std::pow(K0_,3);
    double K0Cx,K0Cy,K0Cz,K00,K1C0, K1C1,K2C0,K2C1,K3C0, K3C1;
    K0Cx = K1_*_scale/(K0_*K0_);
    K0Cy = K2_*_scale/(K0_*K0_);
    K0Cz = K3_*_scale/(K0_*K0_);
    K00 = ( (K1_*eta100 + K2_*eta010 + K3_*eta001)*(_scale*_scale)/K03 \
            -0.5*(eta200 + eta020 + eta002)*(_scale*_scale)/(K0_*K0_));
    K1C0 = _scale/K0_;
    K1C1 = K1_*eta100*R3/K03;
    K2C0 = _scale/K0_;
    K2C1 = K2_*eta010*R3/K03;
    K3C0 = _scale/K0_;
    K3C1 = K3_*eta001*R3/K03;
    for (int i=0; i<= maxOrder_; ++i) {
        for (int j=0; j<maxOrder_+1-i; ++j) {
            for (int k=0; k<maxOrder_+1-i-j; ++k) {

                dmoments_[i][j][k][0] = i*K0Cx*moments_[std::max(i - 1,0)][j][k]  \
                                      + j*K0Cy*moments_[i][std::max(j - 1,0)][k]  \
                                      + k*K0Cz*moments_[i][j][std::max(k - 1,0)]  \
                                      + K00*moments_[i][j][k];

                // dM/dK1
                dmoments_[i][j][k][1] = -(i*K1C0 + K1C1)*moments_[std::max(i - 1,0)][j][k];

                // dM/dK2
                dmoments_[i][j][k][2] = -(j*K2C0 + K2C1)*moments_[i][std::max(j - 1,0)][k];

                // dM/dK3
                dmoments_[i][j][k][3] = -(k*K3C0 + K3C1)*moments_[i][j][std::max(k - 1,0)];

            }
        }
    }
}

template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>::ComputeDiffFunction (T1DIter _iter, T1DIter _diffIter, int _dim)
{
    _diffIter[0] = -_iter[0];
    for (int i=1; i<_dim; ++i)
    {
        _diffIter[i] = _iter[i-1] - _iter[i];
    }
    _diffIter[_dim] = _iter[_dim-1];
}


template<class VoxelT, class MomentT>
MomentT ScaledGeometricalMoments<VoxelT,MomentT>::Multiply (T1DIter _diffIter, T1DIter _sampleIter, int _dim)
{
    T sum (0);
    for (int i=0; i<_dim; ++i)
    {
        _diffIter[i] *= _sampleIter[i];
        sum += _diffIter[i];
    }

    return sum;
}

template<class VoxelT, class MomentT>
MomentT ScaledGeometricalMoments<VoxelT,MomentT>::GetMoment (int _i, int _j, int _k)
{
    return moments_[_i][_j][_k];
}

template<class VoxelT, class MomentT>
void ScaledGeometricalMoments<VoxelT,MomentT>:: GetMomentDerivative(int i, int j, int k, T1D& dMoment) {

    dMoment.resize(4);
    for(unsigned int ii = 0; ii < 4; ii++)
        dMoment[ii] = dmoments_[i][j][k][ii];
}

#endif
