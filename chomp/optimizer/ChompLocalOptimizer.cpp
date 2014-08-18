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
    use_mutex( false ),
    use_goalset( false ),
    use_momentum( use_momentum ),
    hmc( NULL )
{

    N_sub = 0;
    minN = N;
    assert(maxN >= minN);

}

  
void ChompOptimizer::prepareChomp() {
    
    debug << "Preparing Chomp" << std::endl;

    // Subsample
    bool subsample = N > minN && !use_goalset &&
                     !(full_global_at_final && N >= maxN);
    
    if (factory) { factory->getAll( N );}

    //If we are doing goalset chomp, prepare for it.
    if (use_goalset ){ prepareGoalSet(); }

    gradient->prepareRun( N, use_goalset, subsample );

    if( subsample ){
        N_sub = (N+1)/2;
        new (&xi_sub) SubMatMap( xi.data(), N_sub, M,
                                 SubMatMapStride(N,2) );
    }else {
        N_sub = 0;
        if ( use_momentum ){
            momentum.resize( N, M );
            momentum.setZero();
        }
        if( hmc ){ hmc->setupRun(); }
    } 

    debug << "Done Preparing Chomp" << std::endl;
}
  

// precondition: prepareChomp was called for this resolution level
void ChompOptimizer::prepareChompIter() {
    
    if (hmc && !N_sub ) {
        //TODO make sure that we are doing Global Chomp,
        //  or add momentum into local chomp
        hmc->iteration( cur_iter, xi, momentum,
                        gradient->getInvAMatrix(),
                        lastObjective );
    }
    
    debug << "Getting chomp Constraints" << std::endl;
    prepareChompConstraints();
    debug << "Done Getting chomp Constraints" << std::endl;
    
    gradient->getGradient( xi );
    if ( N_sub ){ gradient->getSubsampledGradient( N_sub ); }
}


//TODO make user that this goes in the correct place
void ChompOptimizer::prepareChompConstraints(){

    // compute gradient, constraint & Jacobian and subsample them if
    // needed
    if ( factory ){
        //If we are subsampling, get constraints corresponding
        //  to the subsampled trajectory
        if (N_sub) { factory->evaluate(xi, h_sub, H_sub, 2); }
        //else, get constraints corresponding to the standard trajectory.
        else { factory->evaluate(xi, h, H); }
    }

    if (h.rows()) {
        hmag = h.lpNorm<Eigen::Infinity>();
        debug << "in prepareChompIter, ||h|| = " << hmag << "\n";
    } else if (h_sub.rows()) { 
        hmag = h_sub.lpNorm<Eigen::Infinity>();
        debug << "in prepareChompIter, ||h|| = " << hmag << "\n";
    } else {
        hmag = 0;
    }

}

// precondition: prepareChomp was called for this resolution level
// updates chomp equation until convergence at the current level

void ChompOptimizer::runChomp(bool global, bool local) {
    
    prepareChompIter();
    lastObjective = gradient->evaluateObjective(xi);

    if (notify(CHOMP_INIT, 0, lastObjective, -1, hmag)) { 
        global = false;
        local = false;
    }
    
    //do Global chomp
    cur_iter = 0;
    while (global) { global = iterateChomp( false ); }
    
    //If GoalSet chomp occurred, finish it.
    if (use_goalset){ finishGoalSet(); }
    
    //set the cur_iter to zero, to prepare for 
    //  local chomping.
    cur_iter = 0;
    
    if (full_global_at_final && N >= maxN) { local = false; }

    while (local) { local = iterateChomp( true ); }
    
    //handle subsampled constraint evaluation
    if (factory && N_sub) {
        factory->evaluate(xi, h, H);
        if (h.rows()) {
            hmag = h.lpNorm<Eigen::Infinity>();
        }
    }

    notify(CHOMP_FINISH, 0, lastObjective, -1, hmag);

} 
    
bool ChompOptimizer::iterateChomp( bool local ){
    
    debug << "Starting Iteration" << std::endl;

    ChompEventType event;
    bool not_finished = true;

    //run local or global chomp
    localSmooth();
    checkBounds( xi );
    event = CHOMP_LOCAL_ITER;

    //get the next gradient
    gradient->getGradient( xi );

    
    return not_finished;

}



// calls runChomp and upsamples until N big enough
// precondition: N <= maxN
// postcondition: N >= maxN
void ChompOptimizer::solve(bool doGlobalSmoothing, bool doLocalSmoothing) {
   
    if ( timeout_seconds < 0 ){ canTimeout = false; }
    else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }
    
    if ( hmc ){ 
        use_momentum = true;
        hmc->setupHMC( objective_type, alpha );
    }
    
}

// single iteration of local smoothing
//
// precondition: prepareChompIter has been called since the last
// time xi was modified
void ChompOptimizer::localSmooth() {

    debug << "Starting localSmooth" << std::endl;

    MatX h_t, H_t, P_t, P_t_inv, delta_t;

    hmag = 0;
    
    const MatX & g = gradient->g;

    for (int t=0; t<N; ++t){

        Constraint* c = factory->constraints.empty() ? NULL :
                                                  factory->constraints[t];
        
        bool is_constrained = (c && c->numOutputs() > 0);

        //if this timestep could be constrained,
        //  evaluate the constraints
        if (is_constrained) {
            c->evaluateConstraints(xi.row(t), h_t, H_t);
            is_constrained = h_t.rows() > 0;
        }
        
        //if there are active constraints this timestep.
        if ( is_constrained ) {

            hmag = std::max(hmag, h_t.lpNorm<Eigen::Infinity>());
            
            debug << "ROWS: " << H_t.rows() << " " <<
                      factory->constraints[t]->numOutputs() << "\n";
            debug << "ROWS: " << h_t.rows() << " " <<
                      factory->constraints[t]->numOutputs() << "\n";
    
            P_t = H_t*H_t.transpose();
            P_t_inv = P_t.inverse();

            // transpose g to be a column vector
            delta_t = ( (MatX::Identity(M,M) - H_t.transpose()*P_t_inv*H_t)
                        * g.row(t).transpose() * alpha
                        + H_t.transpose()*P_t_inv*h_t 
                      ).transpose();
        
        }
        //there are no constraints, so just add the negative gradient
        //  into the trajectory (multiplied by the step size, of course.
        else { delta_t = alpha * g.row(t); }
        
        xi.updateTrajectory( delta_t, t, false );
    }
    
    debug << "Done with localSmooth" << std::endl;
}


}// namespace

