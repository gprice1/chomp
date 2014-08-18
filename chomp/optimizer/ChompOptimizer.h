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

#ifndef _CHOMP_OPTIMIZER_H_
#define _CHOMP_OPTIMIZER_H_

#include "ChompOptimizerBase.h"

namespace chomp {

class ChompOptimizer : public ChompOptimizerBase {
  public:

    MatX h; // constraint function of size k-by-1
    MatX H; // constraint Jacobian of size k-by-M*N
    
    // working variables
    MatX P, HP, Y, W, delta, delta_trans; 
    
    //A cholesky solver for solving the constraint matrix.
    Eigen::LDLT<MatX> cholSolver;


    ChompOptimizer(
          ConstraintFactory* f,
          const MatX& xi_init, // should be N-by-M
          const MatX& pinit, // q0
          const MatX& pgoal, // q1onst MatX& xi_init, // should be N-by-M
          int nmax,
          double al = 0.1,
          double obstol = 0.01,
          size_t max_global_iter=size_t(-1),
          size_t max_local_iter=size_t(-1),
          double t_total=1.0,
          double timeout_seconds=-1.0,
          bool use_momentum = false);
    
    ~ChompOptimizer(){};
    
  protected: 

    virtual bool iterate(Trajectory & xi);

    virtual void optimize(Trajectory & xi);

    // returns true if performance has converged
    bool goodEnough(double oldObjective, double newObjective);

};

}//Namespace


#endif
