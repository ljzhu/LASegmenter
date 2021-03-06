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

#ifndef SCALEDGEOMETRICALMOMENTS_H
#define SCALEDGEOMETRICALMOMENTS_H

#include <vector>
#include <cmath>
#include <iostream>

using std::vector;

/**
    Class for computing the scaled, pre-integrated geometrical moments.
    These tricks are needed to make the computation numerically stable.
    See the paper for more details.  
    \param VoxelT   type of the voxel values
    \param MomentT  type of the moments -- recommended to be double
 */
template<class VoxelT, class MomentT>
class ScaledGeometricalMoments
{
public:
    // ---- public typedefs ----
    /// the moment type
    typedef MomentT             T;          
    /// vector scalar type
    typedef vector<T>           T1D;        
    /// 2D array scalar type
    typedef vector<T1D>         T2D;        
    /// 3D array scalar type
    typedef vector<T2D>         T3D;
    /// 4D array scalar type
    typedef vector<T3D>         T4D;
    /// vector scalar type
    typedef vector<double>      Double1D;   
    /// vector scalar type
    typedef vector<Double1D>    Double2D;   

    typedef typename T1D::iterator       T1DIter;

    // ----- public methods -----
    
    // ---- construction / init ----
    /// Contructor
    ScaledGeometricalMoments (
        const VoxelT* _voxels,  /**< input voxel grid */
        int _xDim,              /**< x-dimension of the input voxel grid */
        int _yDim,              /**< y-dimension of the input voxel grid */
        int _zDim,              /**< z-dimension of the input voxel grid */
        double _xCOG,           /**< x-coord of the center of gravity */ 
        double _yCOG,           /**< y-coord of the center of gravity */ 
        double _zCOG,           /**< z-coord of the center of gravity */ 
        double _scale,          /**< scaling factor */ 
        int _maxOrder = 1       /**< maximal order to compute moments for */ 
        );

    /// Constructor for equal dimensions for each axis
    ScaledGeometricalMoments (
        const VoxelT* _voxels,  /**< input voxel grid */ 
        int _dim,               /**< the grid is _dim^3 */
        double _xCOG,           /**< x-coord of the center of gravity */ 
        double _yCOG,           /**< y-coord of the center of gravity */ 
        double _zCOG,           /**< z-coord of the center of gravity */ 
        double _scale,          /**< scaling factor */ 
        int _maxOrder = 1       /**< maximal order to compute moments for */
        );

    /// Default constructor
    ScaledGeometricalMoments ();

    /// The init function used by the contructors
    void Init (
        const VoxelT* _voxels,  /**< input voxel grid */
        int _xDim,              /**< x-dimension of the input voxel grid */
        int _yDim,              /**< y-dimension of the input voxel grid */
        int _zDim,              /**< z-dimension of the input voxel grid */
        double _xCOG,           /**< x-coord of the center of gravity */ 
        double _yCOG,           /**< y-coord of the center of gravity */ 
        double _zCOG,           /**< z-coord of the center of gravity */ 
        double _scale,          /**< scaling factor */ 
        int _maxOrder = 1,       /**< maximal order to compute moments for */
        bool bDerivative = false    // compute moments' derivative
        );

    void SetBasicMomentsXYZ(double K0,double K1, double K2, double K3);
    /// Access function
    T GetMoment (
        int _i,                 /**< order along x */         
        int _j,                 /**< order along y */          
        int _k                  /**< order along z */         
        );
    void GetMomentDerivative(int i, int j, int k, T1D& dMoment);

private:
    int xDim_,              // dimensions
        yDim_, 
        zDim_, 
        maxOrder_;          // maximal order of the moments

    T K0_,K1_,K2_,K3_;      // original geometrical moments

    T2D         samples_;   // samples of the scaled and translated grid in x, y, z
    T1D         voxels_;    // array containing the voxel grid
    T3D         moments_;   // array containing the cumulative moments
    T4D         dmoments_;  // derivative of moments with respect to K0,K1,K2,K3

    // ---- private functions ----
    void Compute ();   
    void ComputeSamples (double _xCOG, double _yCOG, double _zCOG, double _scale);
    void ComputeDiffFunction (T1DIter _iter, T1DIter _diffIter, int _dim);

    /**
      Added function to copute the derivatives w.r.t geometrical moments
    */
    void ComputeNormalizedMomentsDerivative(double _scale);
    
    T Multiply (T1DIter _diffIter, T1DIter _sampleIter, int _dim);
                
};

#include "ScaledGeometricMoments.cpp"

#endif

