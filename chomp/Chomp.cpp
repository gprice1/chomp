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

#include "Chomp.h"
#include "ConstraintFactory.h"
#include "Constraint.h"
#include "HMC.h"
#include <float.h>
#include <cmath>
#include <iomanip>

#define debug if (0) std::cout
#define debug_assert if (0) assert

namespace chomp {

  Chomp::Chomp(ConstraintFactory* f,
               const MatX& xi_init, // should be N-by-M
               const MatX& pinit, // q0
               const MatX& pgoal, // q1
               int nmax,
               double al,
               double obstol,
               size_t mg,
               size_t ml,
               double tt,
               double timeout_seconds,
               bool use_momentum):
    factory(f),
    observer(0),
    ghelper(0),
    objective_type(MINIMIZE_ACCELERATION),
    maxN(nmax),
    xi(xi_init),
    q0(pinit),
    q1(pgoal),
    alpha(al),
    objRelErrTol(obstol),
    min_global_iter(0),
    max_global_iter(mg),
    min_local_iter(0),
    max_local_iter(ml),
    full_global_at_final(false),
    t_total(tt),
    timeout_seconds( timeout_seconds ),
    didTimeout( false ),
    use_mutex( false ),
    use_goalset( false ),
    use_momentum( use_momentum ),
    hmc( NULL )
  {

    N = xi.rows();
    N_sub = 0;

    M = xi.cols();

    assert(q0.rows() >= 1 && q0.cols() == xi.cols());
    assert(q1.rows() >= 1 && q1.cols() == xi.cols());

    minN = N;
    assert(maxN >= minN);

  }
 

  // calls runChomp and upsamples until N big enough
  // precondition: N <= maxN
  // postcondition: N >= maxN
  void Chomp::solve(bool doGlobalSmoothing, bool doLocalSmoothing) {
  
  
  } 
 
  void Chomp::constrainedUpsampleTo(int Nmax, double htol, double hstep) {

    MatX h, H, delta;
  
    while (N < Nmax) { 

      upsample();
      prepareChomp();

      double hinit = 0, hfinal = 0;
      
      lockTrajectory();
      for (int i=0; i<N; i+=2) {

        Constraint* c = constraints.empty() ? 0 : constraints[i];

        if (!c || !c->numOutputs()) { continue; }
        
        for (int iter=0; ; ++iter) { 
          c->evaluateConstraints(xi.row(i), h, H);
          if (h.rows()) {
            double hn = h.lpNorm<Eigen::Infinity>();
            if (iter == 0) { hinit = std::max(hn, hinit); }
            if (hn < htol) { hfinal = std::max(hn, hfinal); break; }
            delta = H.colPivHouseholderQr().solve(h);
            xi.row(i) -= hstep * delta.transpose();
          }
        }
      
      }
      unlockTrajectory();

      prepareChompIter();
      double f = evaluateObjective();
      if (0) { f = f ? f : f; } // shut up
      //std::cout << "after upsample to " << N << ", objective is " << f << " and ||h|| was reduced from " << hinit << " to " << hfinal << "\n";

    }

  }

void Chomp::useGoalSet( Constraint * goalset ){
      this->goalset = goalset;
      use_goalset = true;
  }


}// namespace

