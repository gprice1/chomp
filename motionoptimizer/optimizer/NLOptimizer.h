/*
* Copyright (c) 2008-2015, Matt Zucker and Temple Price
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

#ifndef _NLOPTIMIZER_H_
#define _NLOPTIMIZER_H_

#include <nlopt.hpp>

#include "OptimizerBase.h"

namespace mopt {

class NLOptimizer : public OptimizerBase{
  
  public:
    
    //the type of algorithm that the optimizer uses.
    nlopt::algorithm algorithm;
    //the type of end that the optimizer comes to.
    nlopt::result result;

    static const char * TAG;

    static const EventType event = NLOPT_ITER;

    NLOptimizer( ProblemDescription & problem,
                 Observer * observer,
                 double obstol = 1e-8,
                 double timeout_seconds = 0,
                 size_t max_iter = size_t(-1)); 

    ~NLOptimizer();

    inline void setAlgorithm( nlopt::algorithm alg ){ algorithm = alg; }

    void solve();

  private:

    void prepareBounds( nlopt::opt & optimizer );
    void prepareCollisionConstraint( nlopt::opt & optimizer );
    void prepareConstraints( nlopt::opt & optimizer );
    
    //a wrapper function for passing the ChompGradient to NLopt.
    static double objectiveFunction(unsigned n,
                                    const double * x,
                                    double* grad,
                                    void *data);
    static void constraintFunction(unsigned constraint_dim,
                                   double * h,
                                   unsigned n_by_m,
                                   const double * x,
                                   double* H,
                                   void *data);
    
    //a wrapper function for passing the collision constraint to NLopt.
    static double collisionConstraintFunction(unsigned n_by_m,
                                              const double * x,
                                              double* h,
                                              void *data);

};


}//namespace 
#endif 
