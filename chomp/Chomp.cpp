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

#define DEBUG_PRINTING 0 
#if DEBUG_PRINTING
    #define debug std::cout
    #define debug_assert assert
#else
    #define debug if (0) std::cout
    #define debug_assert if (0) assert
#endif

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

 //delete the mutex if one was used.
Chomp::~Chomp(){
    if( use_mutex ){
        pthread_mutex_destroy( &trajectory_mutex );
    }
}
void Chomp::lockTrajectory(){
    if (use_mutex){
        pthread_mutex_lock( &trajectory_mutex );
    }
}
void Chomp::unlockTrajectory(){
    if (use_mutex){
        pthread_mutex_unlock( &trajectory_mutex );
    }
}
void Chomp::initMutex(){
    use_mutex = true;
    pthread_mutex_init( &trajectory_mutex, NULL );
}
  

void Chomp::prepareChomp() {
    
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
void Chomp::prepareChompIter() {
    
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
void Chomp::prepareChompConstraints(){

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

void Chomp::runChomp(bool global, bool local) {
    
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
    
bool Chomp::iterateChomp( bool local ){
    
    debug << "Starting Iteration" << std::endl;

    ChompEventType event;
    bool not_finished = true;

    //run local or global chomp
    if ( local ){
        localSmooth();
        checkBounds( xi );
        event = CHOMP_LOCAL_ITER;

        //get the next gradient
        gradient->getGradient( xi );
    }else{
        chompGlobal();

        //check the bounds 
        if (N_sub == 0 ){ checkBounds( xi ); }
        else { checkBounds( xi_sub ); }

        event = CHOMP_GLOBAL_ITER;

        //prepare for the next iteration.
        prepareChompIter();
    }
    
    cur_iter ++;

    //test for termination conditions
    double curObjective = gradient->evaluateObjective( xi );
    bool greater_than_min = cur_iter >
                            (local ? min_local_iter : min_global_iter);
    bool greater_than_max = cur_iter > 
                            (local ? max_local_iter : max_global_iter);
    
    if (greater_than_max || (
        greater_than_min && goodEnough(lastObjective, curObjective)) ||
        notify(event, cur_iter, curObjective, lastObjective, hmag) )
    {
        not_finished = false;
    } 
    else if ( canTimeout && stop_time < TimeStamp::now() ) {
        not_finished = false;
        
        didTimeout = true;
        notify(CHOMP_TIMEOUT, cur_iter, curObjective, lastObjective, hmag);
    }

    lastObjective = curObjective;

    debug << "Ending Iteration" << std::endl;
    return not_finished;

}



// calls runChomp and upsamples until N big enough
// precondition: N <= maxN
// postcondition: N >= maxN
void Chomp::solve(bool doGlobalSmoothing, bool doLocalSmoothing) {
   
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
    
    //Run Chomp at the current iteration, then upsample, repeat 
    //  until the current trajectory is at the max resolution
    while (1) {
        prepareChomp();
        //run chomp at the current resolution
        runChomp(doGlobalSmoothing, doLocalSmoothing); 

        // if N>=maxN, then we have already upsampled enough
        //    else, we should perform upsampling then 
        //    we will perform chomp on the unsampled trajectory.
        if (N >= maxN) { break; }
        else { upsample(); }
    }
}

// upsamples the trajectory by 2x
void Chomp::upsample() {
  MatX xi_up;
  
  //calls the upsample function from chomputil
  upsampleTrajectory( xi, gradient->q0, gradient->q1, gradient->dt,
            gradient->objective_type, xi_up );

  N = xi_up.rows();

  lockTrajectory();
  xi = xi_up;
  unlockTrajectory();

  h = h_sub = H = H_sub = P = HP = Y = W = delta = MatX();

  N_sub = 0;
}

// single iteration of chomp
void Chomp::chompGlobal() { 
    
    assert(xi.rows() == N && xi.cols() == M);

    // see if we're in our base case (not subsampling)
    bool subsample = N_sub != 0;
    
    const MatX& g = (subsample ? gradient->g_sub : gradient->g );
    const MatX& L = gradient->getInvAMatrix( subsample );

    const MatX& H_which = subsample ? H_sub : H;
    const MatX& h_which = subsample ? h_sub : h;
    const int   N_which = subsample ? N_sub : N;
    
    //If there are no constraints,
    //  run the update without constraints.
    if (H_which.rows() == 0) {

        assert( g.rows() == N_which ); 
      
        skylineCholSolve(L, g);
      
        //if we are using momentum, add the gradient into the
        //  momentum.
        if (!subsample && use_momentum ) {
            momentum += g * alpha;
            updateTrajectory( momentum, false );
        }else {
            updateTrajectory( g * alpha, subsample );
        }

    //chomp update with constraints
    } else {

      P = H_which.transpose();
      
      // TODO: see if we can make this more efficient?
      for (int i=0; i<P.cols(); i++){
        skylineCholSolveMulti(L, P.col(i));
      }

      debug << "H = \n" << H << "\n";
      debug << "P = \n" << P << "\n";
  
      HP = H_which*P;

      cholSolver.compute(HP);
      Y = cholSolver.solve(P.transpose());

      debug << "HP*Y = \n" << HP*Y << "\n";
      debug << "P.transpose() = \n" << P.transpose() << "\n";
      debug_assert( P.transpose().isApprox( HP*Y ));

      int newsize = H_which.cols();
      assert(newsize == N_which * M);

      assert(g.rows() == N_which && g.cols() == M);
      
      Eigen::Map<const MatX> g_flat(g.data(), newsize, 1);
      W = (MatX::Identity(newsize,newsize) - H_which.transpose() * Y)
          * g_flat * alpha;
      skylineCholSolveMulti(L, W);

      Y = cholSolver.solve(h_which);

      //debug_assert( h_which.isApprox(HP*Y) );
      
      //handle momentum if we need to.
      if (!subsample && use_momentum){
        Eigen::Map<MatX> momentum_flat( momentum.data(), newsize, 1);
        momentum_flat += W;
        delta = momentum_flat + P * Y;
      }else {
        delta = W + P * Y;
      }

      assert(delta.rows() == newsize && delta.cols() == 1);

      Eigen::Map<const Eigen::MatrixXd> delta_rect(delta.data(),
                                                   N_which, M);
      
      debug << "delta = \n" << delta << "\n";
      debug << "delta_rect = \n" << delta_rect << "\n";
      
      assert(delta_rect.rows() == N_which && 
             delta_rect.cols() == M);
      
      updateTrajectory( delta_rect, subsample );
    }
}



// single iteration of local smoothing
//
// precondition: prepareChompIter has been called since the last
// time xi was modified
void Chomp::localSmooth() {

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
        
        updateTrajectory( delta_t, t, false );
    }
    
    debug << "Done with localSmooth" << std::endl;
}

// returns true if performance has converged
bool Chomp::goodEnough(double oldObjective, double newObjective )
{
    return (fabs((oldObjective-newObjective)/newObjective)<objRelErrTol);
}

void Chomp::constrainedUpsampleTo(int Nmax,
                                           double htol,
                                           double hstep)
{

    MatX h, H, delta;
  
    while (N < Nmax) { 

      upsample();
      prepareChomp();
      prepareChompIter();

      double hinit = 0, hfinal = 0;
      
      //if there is no factory, or there are no constraints,
      //    do not evaluate the constraints.
      if ( !factory || factory->constraints.empty() ){ continue; }

      for (int i=0; i<N; i+=2) {
    
        Constraint* c = factory->constraints[i];

        if (!c || !c->numOutputs()) { continue; }
        
        for (int iter=0; ; ++iter) { 
          c->evaluateConstraints(xi.row(i), h, H);
          if (h.rows()) {
            double hn = h.lpNorm<Eigen::Infinity>();
            if (iter == 0) { hinit = std::max(hn, hinit); }
            if (hn < htol) { hfinal = std::max(hn, hfinal); break; }
            delta = H.colPivHouseholderQr().solve(h);
            updateTrajectory( hstep * delta.transpose(), i );
          }
        }
      }
    }
}



template <class Derived>
void Chomp::checkBounds( Eigen::MatrixBase<Derived> const & traj ){

    const bool check_upper = (upper_bounds.size() == M);
    const bool check_lower = (lower_bounds.size() == M);

    //if there are bounds to check, check for the violations.

    if ( check_upper || check_lower ){
        
        bool violation;
        bounds_violations.resize( traj.rows(), traj.cols() );

        const int max_checks = 10;
        int count = 0;

        do{
            count ++;
            double max_violation = 0.0;
            std::pair< int, int > max_index;

            for ( int j = 0; j < traj.cols(); j ++ ) {
                const double upper = (check_upper ? upper_bounds(j) : 0 );
                const double lower = (check_upper ? lower_bounds(j) : 0 );
                
                for( int i = 0; i < traj.rows(); i ++ ){

                    //is there a violation of a lower bound?
                    if ( check_lower && traj(i,j) < lower ){
                        const double magnitude = traj(i,j) - lower;
                        bounds_violations(i,j) = magnitude; 
                        
                        //save the max magnitude, and its index.
                        if ( -magnitude > max_violation ){
                            max_violation = -magnitude;
                            max_index.first = i;
                            max_index.second = j;
                        }
                    //is there a violation of an upper bound?
                    }else if ( check_upper && traj(i,j) > upper ){
                        const double magnitude = traj(i,j) - upper; 
                        bounds_violations(i,j) = magnitude; 

                        //save the max magnitude, and its index.
                        if ( magnitude > max_violation ){
                            max_violation = magnitude;
                            max_index.first = i;
                            max_index.second = j;
                        }
                    }else {
                        bounds_violations(i,j) = 0;
                    }
                }
            }
            

            //There are violations in the largest violation does not
            //  have a magnitude of zero.
            violation = ( max_violation > 0.0 );
            
            if( violation ){
                //smooth out the bounds violation matrix. With
                //  the appropriate smoothing matrix.
                if ( bounds_violations.rows() == traj.rows() ){
                    skylineCholSolve( gradient->L, bounds_violations );
                }else {
                    assert( gradient->L_sub.rows() == traj.rows() );
                    skylineCholSolve( gradient->L_sub, bounds_violations );
                }

                //scale the bounds_violation matrix so that it sets the
                //  largest violation to zero.
                double current_mag = bounds_violations( max_index.first,
                                                        max_index.second );
                const double scale = (current_mag > 0 ?
                                      max_violation/current_mag :
                                     -max_violation/current_mag  );
                const_cast<Eigen::MatrixBase<Derived>&>(traj) -= 
                                             bounds_violations * scale;
            }
        //continue to check and fix limit violations as long as they exist.
        }while( violation && count < max_checks );
    }
}

////////////////GOAL SET FUNCTIONS//////////////////////////////////

void Chomp::useGoalSet( Constraint * goalset ){
      this->goalset = goalset;
      use_goalset = true;
}

void Chomp::prepareGoalSet(){
    
    //do not subsample if doing goalset run.
    N_sub = 0;

    //resize xi, and add q1 into it.
    xi.conservativeResize( xi.rows() + 1, xi.cols() );
    xi.row( xi.rows() - 1 ) = gradient->q1;
    
    //set N to the current size of xi.
    N = xi.rows();

    //add the goal constraint to the constraints vector.
    factory->constraints.push_back( goalset );
}

void Chomp::finishGoalSet(){
    
    debug << "Finishing goal set" << std::endl;

    use_goalset = false;
    
    //copy the last state in the trajectory back into q1
    gradient->q1 = xi.row( xi.rows() - 1 );

    //resize xi, keeping old values, and get rid of the 
    //  final state. And set N to the correct trajectory size
    xi.conservativeResize( xi.rows() -1, xi.cols() );
    N = xi.rows();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory->constraints.pop_back();
    
    //call prepare chomp to reset important stuff.
    prepareChomp();
}

}// namespace

