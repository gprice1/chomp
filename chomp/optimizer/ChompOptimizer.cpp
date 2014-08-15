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

#include "ChompOptimizer.h"
#include "ConstraintFactory.h"
#include "Constraint.h"
#include "HMC.h"
#include <float.h>
#include <cmath>

#define DEBUG_PRINTING 0 
#if DEBUG_PRINTING
    #define debug std::cout
    #define debug_assert assert
#else
    #define debug if (0) std::cout
    #define debug_assert if (0) assert
#endif

namespace chomp {


ChompOptimizer::ChompOptimizer(
               ConstraintFactory* f,
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
    ChompOptimizerBase(f, xi_init, pinit, pgoal, MatX(0,0), MatX(0,0),
                       MINIMIZE_ACCELERATION, tt),
    maxN(nmax),
    xi_sub( NULL, 0, 0, SubMatMapStride(0,2) ),
    alpha(al),
    objRelErrTol(obstol),
    min_global_iter(0),
    max_global_iter(mg),
    min_local_iter(0),
    max_local_iter(ml),
    full_global_at_final(false),
    timeout_seconds( timeout_seconds ),
    didTimeout( false ),
    use_momentum( use_momentum ),
    hmc( NULL )
{
}

virtual void ChompOptimizer::solve( MatX & xi ){
    optimize( xi );
}

virtual void ChompOptimizer::solve( SubMatMap & xi ){
    optimize( xi );
}

// precondition: prepareChomp was called for this resolution level
template <class Derived> 
void ChompOptimizer::prepareIter( Eigen::MatrixBase<Derived> const & xi)
{

    if ( hmc ) {
        //TODO make sure that we are doing Global Chomp,
        //  or add momentum into local chomp
        hmc->iteration( cur_iter, xi, momentum,
                        gradient->getInvAMatrix(),
                        lastObjective );
    }
    
    // If there is a factory, 
    //  get constraints corresponding to the trajectory.
    if ( factory ){
        factory->evaluate(xi, h, H);

        if (h.rows()) {
            hmag = h.lpNorm<Eigen::Infinity>();
        } else {
            hmag = 0;
        }
    }

    gradient->getGradient();
}


// precondition: prepareChomp was called for this resolution level
// updates chomp equation until convergence at the current level
template <class Derived>
void ChompOptimizer::optimize( Eigen::MatrixBase<Derived> const & xi {
    
    if (notify(CHOMP_INIT, 0, lastObjective, -1, hmag)) { 
        global = false;
        local = false;
    }

    if ( timeout_seconds < 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }
    
    if ( hmc ){ 
        use_momentum = true;
        hmc->setupHMC( objective_type, alpha );
        hmc->setupRun();
    }

    if ( use_momentum ){
        momentum.resize( N, M );
        momentum.setZero();
    }
    
    bool not_finished = true;

    //Do the optimization
    prepareIter( xi );
    lastObjective = gradient->evaluateObjective(xi);
    while (not_finished) { not_finished = iterate( xi ); }
    
    notify(CHOMP_FINISH, 0, lastObjective, -1, hmag);
} 
    
template <class Derived>
bool ChompOptimizer::iterate(Eigen::MatrixBase<Derived> const & xi){
    
    debug << "Starting Iteration" << std::endl;

    //perform optimization
    chompGlobal( xi );

    //check and correct the bounds 
    checkBounds( xi );
    
    //increment the iteration.
    cur_iter ++;
    
    //prepare for the next iteration.
    prepareIter( xi );
    
    //check whether optimization is completed.
    return checkFinished( CHOMP_GLOBAL_ITER );
}

// single iteration of chomp
template <class Derived>
void ChompOptimizer::chompGlobal(Eigen::MatrixBase<Derived> const & xi) { 
    
    const MatX& g = gradient->g;
    const MatX& L = gradient->getInvAMatrix();

    //If there are no constraints,
    //  run the update without constraints.
    if (H.rows() == 0) {

        skylineCholSolve(L, g);
      
        //if we are using momentum, add the gradient into the
        //  momentum.
        if ( use_momentum ) {
            momentum += g * alpha;
            const_cast< MatrixBase<Derived>&>(xi) -= momentum;
        }
        else { const_cast< MatrixBase<Derived>&>(xi) -= g * alpha; }

    //chomp update with constraints
    } else {

      P = H.transpose();
      
      // TODO: see if we can make this more efficient?
      for (int i=0; i<P.cols(); i++){ skylineCholSolveMulti(L, P.col(i));}

      debug << "H = \n" << H << "\n";
      debug << "P = \n" << P << "\n";
  
      HP = H*P;

      cholSolver.compute(HP);
      Y = cholSolver.solve(P.transpose());

      debug << "HP*Y = \n" << HP*Y << "\n";
      debug << "P.transpose() = \n" << P.transpose() << "\n";
      debug_assert( P.transpose().isApprox( HP*Y ));

      int newsize = H.cols();
      assert(newsize == N * M);

      assert(g.rows() == N && g.cols() == M);
      
      ConstMatMap g_flat(g.data(), newsize, 1);
      W = (MatX::Identity(newsize,newsize) - H.transpose() * Y)
          * g_flat * alpha;
      skylineCholSolveMulti(L, W);

      Y = cholSolver.solve(h);

      //handle momentum if we need to.
      if (use_momentum){
        MatMap momentum_flat( momentum.data(), newsize, 1);
        momentum_flat += W;
        delta = momentum_flat + P * Y;
      }else {
        delta = W + P * Y;
      }

      assert(delta.rows() == newsize && delta.cols() == 1);

      ConstMatMap delta_rect(delta.data(), N, M);
      
      debug << "delta = \n" << delta << "\n";
      debug << "delta_rect = \n" << delta_rect << "\n";
      
      assert(delta_rect.rows() == N && 
             delta_rect.cols() == M);
      
      const_cast< MatrixBase<Derived>&>(xi) -= delta_rect;
    }
}

}// namespace

