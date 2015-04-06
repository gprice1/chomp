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

namespace chomp {

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


/////////////////forward class declarations///////////////
class ChompGradient;
class ChompOptimizerBase;
class ConstraintFactory;
class Constraint;
class HMC;
class Chomp;
class Trajectory;

class OptimizerBase;
class ChompOptimizerBase;
class ChompOptimizer;
class ChompLocalOptimizer;

enum ChompEventType { 
    CHOMP_INIT,
    CHOMP_GLOBAL_ITER,
    CHOMP_LOCAL_ITER,
    NLOPT_ITER,
    CHOMP_FINISH,
    CHOMP_TIMEOUT,
    CHOMP_GOALSET_ITER,
    CHOMP_COVARIANT_ITER

};

enum ChompObjectiveType {
    MINIMIZE_VELOCITY = 0,
    MINIMIZE_ACCELERATION = 1,
};

/////////////////Utility types /////////////////////////
// These classes interface with chomp, to do something useful
class ChompObserver {
  public:
    virtual ~ChompObserver(){}
    virtual int notify(const OptimizerBase& c, 
                       ChompEventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)
    {
        return 0;
    }
};

class DebugChompObserver: public ChompObserver {
  public:
    virtual ~DebugChompObserver(){}
    virtual int notify(const OptimizerBase& c, 
                       ChompEventType e,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)
    {
            
         std::string event_string;

         switch (e) {
            case CHOMP_INIT:
                 event_string = "CHOMP_INIT"; break;
            case CHOMP_GLOBAL_ITER:
                 event_string = "CHOMP_GLOBAL_ITER"; break;
            case CHOMP_LOCAL_ITER:
                 event_string = "CHOMP_LOCAL_ITER"; break;
            case NLOPT_ITER:
                 event_string = "NLOPT_ITER"; break;
            case CHOMP_FINISH:
                 event_string = "CHOMP_FINISH"; break;
            case CHOMP_TIMEOUT:
                 event_string = "CHOMP_TIMEOUT"; break;
            case CHOMP_GOALSET_ITER: 
                 event_string = "CHOMP_GOALSET_ITER"; break;
            case CHOMP_COVARIANT_ITER:
                 event_string = "CHOMP_COVARIANT_ITER"; break;
            default:
                 event_string = "[INVALID]"; break;
         }

         std::cout << "chomp debug: "
              << "event=" << event_string << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10)
              << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10)
              << constraintViolation << "\n";
        if ( e != CHOMP_INIT && e != CHOMP_FINISH && 
            (std::isnan(curObjective) || std::isinf(curObjective) ||
             std::isnan(lastObjective) || std::isinf(lastObjective)) )
        {
            return 1;
        }
        return 0;
    }
};



///////////////////////////////////////////////////////////
/////////////////////Inline Functions//////////////////////

//Copy a MatX into a std::vector of doubles
inline void matToVec( const MatX & mat, std::vector<double> & vec ){
    const double * data = mat.data();
    vec = std::vector<double>(data, data + mat.size());
}

//Copy a vector of doubles into a column-major Matrix.
inline void vecToMat( const std::vector<double> & vec, MatX & mat )
{
    if ( mat.size() == int(vec.size()) ){
        mat = ConstMatMap( vec.data(), mat.rows(), mat.cols() );
    }
    else {
        mat = ConstMatMap( vec.data(), 1, vec.size());
    }
}

//Copy a vector of doubles into a MatX.
inline void vecToMat( const std::vector<double> & vec,
                      int N, int M, MatX & mat)
{
    assert( int( vec.size() ) == N * M );
    mat = ConstMatMap( vec.data(), N, M );
}

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
