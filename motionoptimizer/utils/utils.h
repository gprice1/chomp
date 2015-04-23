/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>


#ifdef DEBUG
    #define debug_assert assert
    #define debug std::cout
    #define debug_status( tag, function, point ) std::cout << tag << " - " << function << " - " << point << std::endl
#else 
    #define debug_assert if(0) assert
    #define debug if(0) std::cout
    #define debug_status( tag, function, point ) if(0) std::cout << tag 
#endif

namespace mopt {

/////////////////////typedefs///////////////////////////////

typedef Eigen::MatrixXd MatX;
typedef Eigen::Map<MatX>  MatMap;
typedef const Eigen::Map< const MatX>  ConstMatMap;

//this matrix map is used for subsampling a matrix. 
typedef Eigen::Stride<Eigen::Dynamic,2> SubMatMapStride;
typedef Eigen::Map<MatX, 0, SubMatMapStride > SubMatMap; 
typedef const Eigen::Map<const MatX, 0, SubMatMapStride > ConstSubMatMap;

typedef Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> DynamicStride;
typedef Eigen::Map<MatX, 0, DynamicStride > DynamicMatMap; 
typedef const Eigen::Map<const MatX, 0, DynamicStride > ConstDynamicMatMap;

typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, 
                       Eigen::RowMajor> MatXR;
typedef Eigen::Map<MatXR> MatMapR;
typedef const Eigen::Map<const MatXR> ConstMatMapR;
typedef Eigen::Map<MatXR, 0, Eigen::OuterStride<> > SubMatMapR; 

typedef Eigen::Block<DynamicMatMap, 1, Eigen::Dynamic> Row;
typedef const Eigen::Block<const DynamicMatMap, 1, Eigen::Dynamic> ConstRow;
typedef Eigen::Block<DynamicMatMap, Eigen::Dynamic, 1, Eigen::Dynamic> Col;
typedef const Eigen::Block< const DynamicMatMap, Eigen::Dynamic,
                            1, Eigen::Dynamic> ConstCol;


//Simple enum types.
enum EventType { 
    INIT,
    CHOMP_ITER,
    CHOMP_LOCAL_ITER,
    NLOPT_ITER,
    FINISH,
    TIMEOUT,
};

enum ObjectiveType {
    MINIMIZE_VELOCITY = 0,
    MINIMIZE_ACCELERATION = 1,
};

/////////////////Utility types /////////////////////////
// These classes interface with the motion optimizer, to do something useful


//this is a function to compute the dot product of matrices
template <class Derived1, class Derived2>
static inline double mydot(const Eigen::MatrixBase<Derived1>& a,
                           const Eigen::MatrixBase<Derived2>& b)
{
  return a.cwiseProduct(b).sum();
}


template <class Derived>
inline MatX getPos(const Eigen::MatrixBase<Derived>& x, double h){
    double fac = 1;
    double hn = 1;

    MatX rval = x.row(0);
    for (int i=1; i<x.rows(); ++i) {
      fac *= i;
      hn *= h;
      rval += (hn / fac) * x.row(i);
    }

    return rval;
}


}//namespace

#endif
