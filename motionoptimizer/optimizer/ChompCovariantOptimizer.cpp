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

#include "ChompCovariantOptimizer.h"
#include <float.h>
#include <cmath>

namespace chomp {

ChompCovariantOptimizer::ChompCovariantOptimizer(
                                        ProblemDescription & problem,
                                        ChompObserver * observer,
                                        double obstol,
                                        double timeout_seconds,
                                        size_t max_iter) : 
    ChompOptimizerBase( problem, observer,
                        obstol, timeout_seconds,
                        max_iter)
{
    event = CHOMP_COVARIANT_ITER;
}

// single iteration of local smoothing
//
// precondition: prepareChompIter has been called since the last
// time xi was modified
void ChompCovariantOptimizer::optimize()
{
    
    debug_status( TAG, "optimize", "start" );

    if( problem.isConstrained() ){

        constraint_magnitude = problem.evaluateConstraint( h, H );
        debug_status( TAG, "optimize", "got_constraint" );
        
        P_inv = ( H.transpose()*H ).inverse();
        
        debug_status( TAG, "optimize", "got P_inv" );

        const int n_by_m = problem.size();
        
        delta = (alpha*(MatX::Identity(n_by_m, n_by_m)
                 - H*P_inv*H.transpose() )
                * MatMap( g.data(), n_by_m , 1 )
                + H*P_inv*h).transpose();
        
        problem.updateTrajectory( MatMap( delta.data(),
                                  problem.N(),
                                  problem.M() ) );
        
        debug_status( TAG, "optimize", "finished constraint eval" );

    }else {

        //if the problem is not constrained, just
        //  update the trajectory with the gradient
        problem.updateTrajectory( alpha * g);
    }
    
    debug_status( TAG, "optimize", "end" );
    
}


}// namespace

