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

namespace chomp {

const char* ChompOptimizer::TAG = "ChompOptimizer";

ChompOptimizer::ChompOptimizer(ProblemDescription & problem,
                               ChompObserver * observer,
                               double obstol,
                               double timeout_seconds,
                               size_t max_iter) : 
    ChompOptimizerBase( problem, observer,
                        obstol, timeout_seconds, max_iter)
{
    event = CHOMP_GLOBAL_ITER;
}

// single iteration of chomp
void ChompOptimizer::optimize() { 
    
    debug_status( TAG, "optimize", "start" );
    
    // If there is a factory, 
    //  get constraints corresponding to the trajectory.
    if ( problem.isConstrained() ){
        constraint_magnitude = problem.evaluateConstraint(h, H);
    }
    
    debug_status( TAG, "optimize", "after getting constraints" );

    //If there are no constraints,
    //  run the update without constraints.
    if (H.rows() == 0) {
        debug_status( TAG, "optimize" , "start unconstrained" );
        if ( !problem.isCovariant() ){
            problem.getMetric().solve( g );
        }
      
        //if we are using momentum, add the gradient into the
        //  momentum.
        if ( use_momentum ) {
            momentum += g * alpha;
            problem.updateTrajectory( momentum );
        }
        else { problem.updateTrajectory( g * alpha ) ; }

        debug_status( TAG, "optimize" , "end unconstrained" );
        
    //chomp update with constraints
    } else {
        debug_status( TAG, "optimize", "simple stuff" );
      
        const int M = problem.M();
        const int N = problem.N();
        
        if (problem.isCovariant() ){
            
            cholSolver.compute( H.transpose()*H );
            
            debug_status( TAG, "optimize", "got P_inv" );

            const int n_by_m = problem.size();
             
            delta = (alpha*(MatX::Identity(n_by_m, n_by_m)
                     - H*cholSolver.solve( H.transpose() ) )
                    * MatMap( g.data(), n_by_m , 1 )
                    + H*cholSolver.solve(h));

        }else {
            debug_status( TAG, "optimize", "constrained case" );

            P = H;
            
            problem.getMetric().solve( MatMap(P.data(), N, M * P.cols()) );
            
            debug_status( TAG, "optimize", "after first skyline" );

            //debug << "H = \n" << H << "\n";
            //debug << "P = \n" << P << "\n";
          
            HP = H.transpose()*P;

            cholSolver.compute(HP);
            Y = cholSolver.solve(P.transpose());

            //debug << "HP*Y = \n" << HP*Y << "\n";
            //debug << "P.transpose() = \n" << P.transpose() << "\n";
            debug_assert( P.transpose().isApprox( HP*Y ));

            int newsize = H.rows();
            
            assert(newsize == N * M);
            assert(g.rows() == N && g.cols() == M);
            
            W = (MatX::Identity(newsize,newsize) - H * Y)
                * MatMap(g.data(), newsize, 1) * alpha;

            problem.getMetric().solve(
                    MatMap(W.data(), N, M * W.cols()) );

            Y = cholSolver.solve(h);

            debug_status( TAG, "optimize", "middle constraint step eval" );
            
            //handle momentum if we need to.
            if (use_momentum){
                MatMap momentum_flat( momentum.data(), newsize, 1);
                momentum_flat += W;
                delta = momentum_flat + P * Y;
            }else {
                delta = W + P * Y;
            }

            assert(delta.rows() == newsize && delta.cols() == 1);

        }
            
        //update the trajectory with the found values.
        problem.updateTrajectory( MatMap(delta.data(), N, M) );

    }

    debug_status( TAG, "optimize", "end" );
}

}// namespace

