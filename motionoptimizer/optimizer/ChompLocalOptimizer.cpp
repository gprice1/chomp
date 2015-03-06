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

ChompLocalOptimizer::ChompLocalOptimizer( Trajectory & traj,
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

// single iteration of local smoothing
//
// precondition: prepareChompIter has been called since the last
// time xi was modified
void ChompLocalOptimizer::optimize( )
{

    debug << "Starting localSmooth" << std::endl;

    MatX h_t, H_t, P_t_inv, delta_t;

    hmag = 0;
    
    const MatX & g = gradient->g;

    for (int t=0; t<trajectory.N(); ++t){

        Constraint* c = factory->constraints.empty() ?
                        NULL : factory->constraints[t];
        
        bool is_constrained = (c && c->numOutputs() > 0);

        //if this timestep could be constrained,
        //  evaluate the constraints
        if (is_constrained) {
            c->evaluateConstraints(trajectory.row(t), h_t, H_t);
            is_constrained = h_t.rows() > 0;
        }
        
        //if there are active constraints this timestep.
        if ( is_constrained ) {

            hmag = std::max(hmag, h_t.lpNorm<Eigen::Infinity>());
            
            P_t_inv = ( H_t*H_t.transpose() ).inverse();

            // transpose g to be a column vector
            delta_t = ( (MatX::Identity(trajectory.N(),trajectory.M() )
                         - H_t.transpose()*P_t_inv*H_t )
                        * g.row(t).transpose() * alpha
                        + H_t.transpose()*P_t_inv*h_t 
                      ).transpose();
        }
        //there are no constraints, so just add the negative gradient
        //  into the trajectory (multiplied by the step size, of course.
        else { delta_t = alpha * g.row(t); }
        
        trajectory.update( delta_t, t );
    }
    
    debug << "Done with localSmooth" << std::endl;
}


}// namespace

