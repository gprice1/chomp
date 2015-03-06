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

#ifndef _CLASS_UTILS_H_
#define _CLASS_UTILS_H_

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace chomp {

/////////////////////typedefs///////////////////////////////

typedef Eigen::MatrixXd MatX;
typedef Eigen::Map<MatX>  MatMap;
typedef const Eigen::Map< const MatX>  ConstMatMap;

//this matrix map is used for subsampling a matrix. 
typedef Eigen::Stride<Eigen::Dynamic,2> SubMatMapStride;
typedef Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> DynamicStride;

typedef Eigen::Map<MatX, 0, SubMatMapStride > SubMatMap; 
typedef const Eigen::Map<const MatX, 0, SubMatMapStride > ConstSubMatMap; 

typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, 
                       Eigen::RowMajor> MatXR;
typedef Eigen::Map<MatXR> MatMapR;
typedef const Eigen::Map<const MatXR> ConstMatMapR;
typedef Eigen::Map<MatXR, 0, Eigen::OuterStride<> > SubMatMapR; 

typedef Eigen::Block<MatMap, 1, Eigen::Dynamic> Row;
typedef const Eigen::Block< const MatMap, 1, Eigen::Dynamic> ConstRow;
typedef Eigen::Block<MatMap, Eigen::Dynamic, 1, Eigen::Dynamic> Col;
typedef const Eigen::Block< const MatMap, Eigen::Dynamic,
                            1, Eigen::Dynamic> ConstCol;

typedef Eigen::Block<SubMatMap, 1, Eigen::Dynamic> SubRow;
typedef const Eigen::Block< const SubMatMap, 1, Eigen::Dynamic> ConstSubRow;
typedef Eigen::Block<SubMatMap, Eigen::Dynamic, 1, Eigen::Dynamic> SubCol;
typedef const Eigen::Block< const SubMatMap, Eigen::Dynamic, 1,
                            Eigen::Dynamic> ConstSubCol;


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
    CHOMP_FINISH,
    CHOMP_TIMEOUT,
    CHOMP_GOALSET_ITER,
};


///////////////UTILITY functions/////////////////////////
//prints out a string that corresponds to the variable name
//  of the event type.
const char* eventTypeString(int eventtype) {
    switch (eventtype) {
    case CHOMP_INIT: return "CHOMP_INIT";
    case CHOMP_GLOBAL_ITER: return "CHOMP_GLOBAL_ITER";
    case CHOMP_LOCAL_ITER: return "CHOMP_LOCAL_ITER";
    case CHOMP_FINISH: return "CHOMP_FINISH";
    case CHOMP_TIMEOUT: return "CHOMP_TIMEOUT";
    case CHOMP_GOALSET_ITER: return "CHOMP_GOALSET_ITER";
    default: return "[INVALID]";
    }
}

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
         std::cout << "chomp debug: "
              << "event=" << eventTypeString(e) << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10)
              << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10)
              << constraintViolation << "\n";
        if (std::isnan(curObjective) || std::isinf(curObjective) ||
            std::isnan(lastObjective) || std::isinf(lastObjective)) {
            return 1;
        }
        return 0;
    }
};

enum ChompObjectiveType {
    MINIMIZE_VELOCITY = 0,
    MINIMIZE_ACCELERATION = 1,
};

}//namespace

#endif
