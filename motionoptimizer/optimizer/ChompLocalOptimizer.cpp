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

#include "ChompLocalOptimizer.h"
#include <float.h>
#include <cmath>

namespace chomp {

ChompLocalOptimizer::ChompLocalOptimizer(ProblemDescription & problem,
                                        ChompObserver * observer,
                                        double obstol,
                                        double timeout_seconds,
                                        size_t max_iter) : 
    ChompOptimizerBase( problem, observer,
                        obstol, timeout_seconds,
                        max_iter)
{
    event = CHOMP_LOCAL_ITER;
}

// single iteration of local smoothing
//
// precondition: prepareChompIter has been called since the last
// time xi was modified
void ChompLocalOptimizer::optimize( const MatX & g )
{
    
    debug_status( TAG, "optimize", "start" );

    MatX h_t, H_t, P_t, P_t_inv, delta_t;

    constraint_magnitude = 0;
    
    debug_status( TAG, "optimize", "pre-for-loop" );

    for (int t=0; t < problem.N(); ++t){
        
        bool is_constrained = problem.evaluateConstraint( h_t, H_t, t );

        if ( is_constrained ){        
            constraint_magnitude = std::max(constraint_magnitude,
                                            h_t.lpNorm<Eigen::Infinity>());
            
            P_t_inv = ( H_t*H_t.transpose() ).inverse();

            const int M = problem.M();
            problem.updateTrajectory(
                (alpha *(MatX::Identity(M,M) - H_t.transpose()*P_t_inv*H_t )
                  * g.row(t).transpose()
                + H_t.transpose()*P_t_inv*h_t).transpose(),
                t );
        }

        //there are no constraints, so just add the negative gradient
        //  into the trajectory (multiplied by the step size, of course.
        else { 
            debug_status( TAG, "optimize", "gradient" );
            problem.updateTrajectory( alpha * g.row(t), t );
        }
        
    }
    
    debug_status( TAG, "optimize", "end" );
    
}


}// namespace

