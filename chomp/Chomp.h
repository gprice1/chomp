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

#ifndef _CHOMP_H_
#define _CHOMP_H_

#include "chomputil.h"

#include <vector>
#include <pthread.h>
#include "../mzcommon/TimeUtil.h"
#include "ChompGradient.h"

namespace chomp {

class Chomp {
  public:

    ConstraintFactory* factory;        

    ChompObserver* observer;

    ChompGradient* gradient;

    ChompObjectiveType objective_type;
    
    int M; // degrees of freedom

    // the actual desired number of timesteps
    int maxN; 

    // the base (minimum) number of timesteps
    int minN;

    int N; // number of timesteps
    int N_sub; // number of timesteps for subsampled trajectory
    


    // current trajectory of size N-by-M
    MatX xi;
    SubMatMap xi_sub; // current trajectory of size N_sub-by-M

    MatX h; // constraint function of size k-by-1
    MatX h_sub; // constraint function of size k_sub-by-1

    MatX H; // constraint Jacobian of size k-by-M*N
    MatX H_sub; // constraint Jacobian of size k_sub-by-1

    double hmag; // inf. norm magnitude of constraint violation

    // working variables
    MatX P, HP, Y, W, delta, delta_trans; 
    
    double alpha;       // the gradient step size
    double objRelErrTol; //Objective function value relative to the
                         //previous value

    // last_objective : save the last objective, for use with rejection.
    double lastObjective;
    
    size_t cur_iter; // the number of the current iteration of chomp
    //chomp will not stop global chomp until it has exceeded 
    //  min_global_iter, and must terminate the run at the current 
    //  resolution if the number of iterations reaches max_global_iter
    size_t min_global_iter, max_global_iter;
    //chomp will not stop local chomp until it has exceeded 
    //  min_local_iter, and must terminate the run at the current 
    //  resolution if the number of iterations reaches max_local_iter
    size_t min_local_iter, max_local_iter;

    bool full_global_at_final; //perform an iteration of global chomp
                               //on the whole trajectory at the end?


    
    //timeout_seconds : the amount of time from the start of chomp
    //                  to a forced timeout.
    //canTimeout : is timing out a possible termination condition?
    // didTimeout : has chomp timed out ?
    // stop_time : the time at which chomp will timeout.
    double timeout_seconds;
    bool canTimeout, didTimeout;
    TimeStamp stop_time;

    //A mutex for locking the trajectory when updates are being made.
    //  Only needed if a concurrent thread wants data out of 
    //  chomp.
    pthread_mutex_t trajectory_mutex;
    bool use_mutex;
    
    //A cholesky solver for solving the constraint matrix.
    Eigen::LDLT<MatX> cholSolver;

    std::vector<Constraint*> constraints; // vector of size N
    
    //used for goal set chomp.
    Constraint * goalset;
    bool use_goalset;
     
    bool use_momentum;
    MatX momentum;
    
    //an HMC object for performing the Hamiltonian Monte Carlo method
    HMC * hmc;

    Chomp(ConstraintFactory* f,
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
    
    //delete the mutex if one was used.
    ~Chomp();

    //methods for using a mutex 
    void lockTrajectory();
    void unlockTrajectory();
    void initMutex();
    
    //clear the constraint vector of constraints.
    void clearConstraints();

    //prepares chomp to be run at a resolution level
    void prepareChomp();    

    // precondition: prepareChomp was called for this resolution level
    void prepareChompIter();
    
    //Prepare the constraint matrix 
    void prepareChompConstraints();
    
    //Do one iteration of chomp, locally or globally
    bool iterateChomp( bool global );

    // precondition: prepareChomp was called for this resolution level
    // updates chomp equation until convergence at the current level
    void runChomp(bool global, bool local);

    // calls runChomp and upsamples until N big enough
    // precondition: N <= maxN
    // postcondition: N >= maxN
    void solve(bool doGlobalSmoothing, bool doLocalSmoothing);

    // get the tick, respecting endpoint repetition
    MatX getTickBorderRepeat(int tick) const;

    // upsamples the trajectory by 2x
    void upsample();

    // single iteration of chomp
    void chompGlobal();
    
    // single iteration of local smoothing
    //
    // precondition: prepareChompIter has been called since the last
    // time xi was modified
    void localSmooth();

    // upsamples trajectory, projecting onto constraint for each new
    // trajectory element.
    void constrainedUpsampleTo(int Nmax, double htol, double hstep=0.5);
    
    // call the observer if there is one
    int notify(ChompEventType event,
               size_t iter,
               double curObjective, 
               double lastObjective,
               double constraintViolation) const;

    //Give a goal set in the form of a constraint for chomp to use on the
    //  first resolution level.
    void useGoalSet( Constraint * goalset );
    
    //Call before running goal set chomp
    void prepareGoalSet();

    //call after running goal set chomp.
    void finishGoalSet(); 

    // returns true if performance has converged
    bool goodEnough(double oldObjective, double newObjective);

    //updates the trajectory via a matrix delta. Delta
    // must be the same size and shape as the trajectory,
    //  or the subsampled trajectory
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta,
                           bool subsample );

    //updates the trajectory at the given row, by the vector delta,
    //  delta should have the same number of columns as the trajectory.
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta,
                           int row, bool subsample );

};

template <class Derived>
inline void Chomp::updateTrajectory( 
                              const Eigen::MatrixBase<Derived> & delta,
                              bool subsample )
{
    lockTrajectory();
    if ( subsample ){ xi_sub -= delta; }
    else{ xi -= delta; }
    unlockTrajectory();
}

template <class Derived>
inline void Chomp::updateTrajectory( 
                              const Eigen::MatrixBase<Derived> & delta,
                              int index, bool subsample )
{
    lockTrajectory();
    if ( subsample ){ xi_sub.row( index ) -= delta; }
    else{ xi.row( index ) -= delta; }
    unlockTrajectory();
}



}//Namespace


#endif
