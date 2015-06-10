/*
 Copyright (c) 2008-2014, Matt Zucker and Temple Price
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

#ifndef _SMOOTHNESS_FUNCTION_H_
#define _SMOOTHNESS_FUNCTION_H_

#include "../utils/utils.h"
#include "Trajectory.h"
#include "Metric.h"

namespace mopt {

class SmoothnessFunction {

  private:
    static const char* TAG;
    
    MatX Ax; // A*x of size N-by-M

    MatX g_full; // a working variable for use when calculating
                 // the gradient of a subsampled trajectory
                
    MatX b; // endpoint vectors for this problem of size N-by-M
    
    double c; // c value for objective function

  public:
    
    //prepare the run with the current settings
    void prepareRun(const Trajectory & trajectory,
                    const Metric & metric);
    
    //evaluate the gradient of the objective function
    //  at the current trajectory
    template <class Derived>
    double evaluate( const Trajectory & trajectory,
                     const Metric & metric,
                     const Eigen::MatrixBase<Derived> & g );
    
    // evaluates the objective function for cur. thing.
    double evaluate( const Trajectory & trajectory,
                     const Metric & metric );
    
  private:

    double computeValue( const Trajectory & trajectory ) const;

    //If the trajectory is subsampled, subsample the gradient,
    //because the gradient is computed at a higher resolution
    template <class Derived>
    void subsampleGradient(int N_sub, 
                           const Eigen::MatrixBase<Derived> & g_sub_const);
    
};

//include all of the templated inline functions (or else linking
//  errors occur)
#include "SmoothnessFunction-inl.h"

}//namespace 


#endif
