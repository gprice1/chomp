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
#include "../containers/ConstraintFactory.h"
#include "../containers/Constraint.h"
#include "../containers/ChompGradient.h"
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

ChompOptimizer::ChompOptimizer( Trajectory & traj,
                                ConstraintFactory * factory,
                                ChompGradient * gradient,
                                ChompObserver * observer,
                                double obstol,
                                double timeout_seconds,
                                size_t max_iter,
                                const MatX & lower_bounds,
                                const MatX & upper_bounds) : 
    ChompOptimizerBase( traj, factory, gradient, observer,
                        obstol, timeout_seconds, max_iter,
                        lower_bounds, upper_bounds )
{}

// single iteration of chomp
void ChompOptimizer::optimize() { 
    
    const MatX& g = gradient->g;
    const MatX& L = gradient->getInvAMatrix();


    // If there is a factory, 
    //  get constraints corresponding to the trajectory.
    if ( factory ){
        factory->evaluate(h, H);

        if (h.rows()) {
            hmag = h.lpNorm<Eigen::Infinity>();
        } else {
            hmag = 0;
        }
    }

    //If there are no constraints,
    //  run the update without constraints.
    if (H.rows() == 0) {

        skylineCholSolve(L, g);
      
        //if we are using momentum, add the gradient into the
        //  momentum.
        if ( use_momentum ) {
            momentum += g * alpha;
            trajectory.update( momentum );
        }
        else { trajectory.update( g * alpha ) ; }

    //chomp update with constraints
    } else {
      
        const int M = trajectory.M();
        const int N = trajectory.N();

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
               delta_rect.cols() == M );
        
        trajectory.update( delta_rect );
    }
}

}// namespace

