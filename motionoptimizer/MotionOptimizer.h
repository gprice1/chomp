
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
enum OptimizationAlgorithm {
    LOCAL_CHOMP,
    CHOMP,
    TEST,
    MMA_NLOPT,
    CCSAQ_NLOPT,
    SLSQP_NLOPT,
    LBFGS_NLOPT,
    TNEWTON_PRECOND_RESTART_NLOPT,
    TNEWTON_RESTART_NLOPT,
    TNEWTON_PRECOND_NLOPT,
    TNEWTON_NLOPT,
    VAR1_NLOPT,
    VAR2_NLOPT,
    NONE
};


#ifdef NLOPT_FOUND
nlopt::algorithm getNLoptAlgorithm( OptimizationAlgorithm alg );
#endif

std::string algorithmToString( OptimizationAlgorithm alg );

OptimizationAlgorithm algorithmFromString( const std::string & str );

class MotionOptimizer {
    

  private:

    ProblemDescription problem;
    
    Observer * observer;

    int N_max, N_min;

    bool full_global_at_final, do_subsample;

    double obstol, timeout_seconds, alpha;
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
    
    void setAlgorithm(OptimizationAlgorithm a);
    OptimizationAlgorithm getAlgorithm() const;

    void setAlgorithm1(OptimizationAlgorithm a1);
    void setAlgorithm2(OptimizationAlgorithm a2);
    
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

};

#include "MotionOptimizer-inl.h"

}// namespace

#endif 
