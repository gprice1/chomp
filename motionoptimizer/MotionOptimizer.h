/*
* Copyright (c) 2008-2015, Matt Zucker and Temple Price
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
#ifndef _MOTION_OPTIMIZER_H_
#define _MOTION_OPTIMIZER_H_

#include "utils/utils.h"
#include "utils/Observer.h"

#include "containers/ProblemDescription.h"

#include "containers/ConstraintFactory.h"
#include "containers/Constraint.h"

#include "optimizer/ChompLocalOptimizer.h"
#include "optimizer/ChompOptimizer.h"
#include "optimizer/TestOptimizer.h"

#ifdef NLOPT_FOUND
    #include "optimizer/NLOptimizer.h"
#endif


namespace mopt {

/**
 * \enum OptimizationAlgorithm
 * An enum that tells the MotionOptimizer class the type of 
 * algorithm that will be used for optimization.
 * All algorithms from the NLOPT optimization package are appended
 * with '_NLOPT'
 */    
enum OptimizationAlgorithm {
    LOCAL_CHOMP, ///< The local chomp algorithm, optimizes with gradient
                 /// descent without using the metric
    CHOMP,       ///< the standard CHOMP algorithm
    TEST,        ///< An algorithm to test the interface between NLOPT and
                 /// the problem description
    MMA_NLOPT,   ///< The nlopt MMA algorithm
    CCSAQ_NLOPT, ///< The nlopt MMA algorithm
    SLSQP_NLOPT, ///< The nlopt MMA algorithm
    LBFGS_NLOPT, ///< The nlopt MMA algorithm
    TNEWTON_PRECOND_RESTART_NLOPT, ///< The nlopt Truncated Newton
                                   /// algorithm with LBFGS preconditioner
                                   /// and restarting
    TNEWTON_RESTART_NLOPT, ///< The nlopt Truncated Newton algorithm with 
                           /// restarting
    TNEWTON_PRECOND_NLOPT, ///< The nlopt Truncated Newton algorithm with 
                           /// and LBFGS preconditioner
    TNEWTON_NLOPT, ///< The nlopt Truncated Newton algorithm
    VAR1_NLOPT, ///< The nlopt VAR1 algorithm
    VAR2_NLOPT, ///< The nlopt VAR2 algorithm
    NONE        /// < a placeholder for no algorithm
};


#ifdef NLOPT_FOUND
/**
 * Converts between an nlopt::algorithm type and a
 *  mopt::OptimizationAlgorithm type. 
 */
nlopt::algorithm getNLoptAlgorithm( OptimizationAlgorithm alg );
#endif

/**
 * Takes an OptimizationAlgorithm enum and turns it into
 * a string. Useful for printing information
 * \param alg the algorithm to be stringified
 */
std::string algorithmToString( OptimizationAlgorithm alg );

/**
 * takes a string and turns it into an OptimizationAlgorithm
 * \param str the string that will be turned into an OptimizationAlgorithm
 * \sa OptimizationAlgorithm
 */
OptimizationAlgorithm algorithmFromString( const std::string & str );

/**
 * \class MotionOptimizer
 * This is the main class that a user will interface with to
 * do motion optimization. It is used by creating an optimizer,
 * passing in arguments, and then calling solve()
 */
class MotionOptimizer {
    

  private:
    
    /**
     * The description of the problem to be optimized.
     * This variable is passed on to the optimizer and is
     * the interface between the MotionOptimizer and the 
     * optimizer.
     */
    ProblemDescription problem;
    
    /** 
     * A pointer to an observer, if there is no observer for
     * the given problem, then this variable is NULL.
     */
    Observer * observer;

    /** The max and min values of the length of the trajectory.
     * If N_max is greater than N_min, then the MotionOptimizer
     * will use upsampling to eventually get a trajectory >= N_max
     */
    int N_max, N_min;

    /** 
     * If this is set to true, then the algorithm will perform an
     * optimization at full resolution as the last stage of optimzation.
     * This variable is set to false as a default.
     */ 
    bool full_global_at_final;

    /** 
     * If this is set to false, then there will be no subsampling
     * after stages of upsampling.
     * By default, this variable is set to true.
     */
    bool do_subsample;

    double obstol;
    double timeout_seconds, alpha;
    size_t max_iterations;

    OptimizationAlgorithm algorithm1, algorithm2;

    const static char* TAG;

  public:
    //constructor.
    MotionOptimizer( Observer * observer = NULL,
                     double obstol = 1e-8,
                     double timeout_seconds = 0,
                     size_t max_iter = size_t(-1),
                     const MatX & lower_bounds=MatX(0,0),
                     const MatX & upper_bounds=MatX(0,0),
                     OptimizationAlgorithm algorithm1 = LBFGS_NLOPT,
                     OptimizationAlgorithm algorithm2 = LBFGS_NLOPT,
                     int N_max = 0);

    void solve();
    
  private:
    //sets up the factory, gradient, and optimizer for the current
    //  resolution.
    void optimize( OptimizerBase * optimizer, bool subsample = false);
    OptimizerBase * getOptimizer( OptimizationAlgorithm algorithm );
    
  public:

    //Functions for setting the upper and lower bounds of MotionOptimizer.
    void setLowerBounds( const MatX & lower );
    void setLowerBounds( const std::vector<double> & lower);
    void setLowerBounds( const double * lower, int M);

    void setUpperBounds(const MatX & upper );
    void setUpperBounds(const std::vector<double> & upper );
    void setUpperBounds(const double * upper, int M );

    void setBounds( const MatX & lower, const MatX & upper );
    void setBounds( const std::vector<double> & lower,
                    const std::vector<double> & upper );
    void setBounds( const double * lower,
                    const double * upper,
                    int M );

    //add a constraint starting at start_time, and
    //  ending at end_time
    void addConstraint( Constraint * c,
                        double start_time,
                        double end_time);

    //simple getters and setters, all of which are inline
    //  NOTE: these are all implemented in MotionOptimizer-inl.h
    void setNMax( int n_max );
    int  getNMax() const;
    
    void setGoalset( Constraint * goal);
    const Constraint * getGoalset() const;

    void   setTimeoutSeconds( double s );
    double getTimeoutSeconds() const;

    void setMaxIterations( size_t max );
    size_t getMaxIterations() const;

    void setFunctionTolerance( double tol );
    double getFunctionTolerance() const;
    
    void setAlgorithm(OptimizationAlgorithm a1, OptimizationAlgorithm a2=NONE);
    void setAlgorithm1(OptimizationAlgorithm a1);
    void setAlgorithm2(OptimizationAlgorithm a2);
    
    void setAlgorithm( const std::string & alg1,
                       const std::string & alg2="NONE");
    void setAlgorithm1( const std::string & alg_string );
    void setAlgorithm2( const std::string & alg_string );
    
    OptimizationAlgorithm getAlgorithm() const;
    OptimizationAlgorithm getAlgorithm1() const;
    OptimizationAlgorithm getAlgorithm2() const;

    void dontSubsample();
    void doSubsample();
    void setSubsample( bool subsample );
    
    void setAlpha( double a );
    double getAlpha() const;

    void doFullGlobalAtFinal();
    void dontFullGlobalAtFinal();
    bool getFullGlobalAtFinal() const;

    Trajectory & getTrajectory();
    const Trajectory & getTrajectory() const;
    void setTrajectory( const Trajectory & trajectory );

    void setCollisionFunction( CollisionFunction * coll_func);
    const CollisionFunction * getCollisionFunction() const;

    void setObserver( Observer * obs );
    Observer * getObserver();
    const Observer * getObserver() const;

    void doCovariantOptimization();
    void dontCovariantOptimization();
    void setCovariantOptimization( bool covariant );
    bool isCovariantOptimization() const;

    void doCollisionConstraint();
    void dontCollisionConstraint();
    void setCollisionConstraint( bool do_collision_constraint );
    bool isCollisionConstraint() const ;
    
    //TODO implement.
    void doMomentum(){}
    void dontMomentum(){}
    void setMomentum( bool momentum ){}

    void setHMC( double lambda ){}
    double getHMC() const { return 0.0; }

};

}// namespace

#endif 
