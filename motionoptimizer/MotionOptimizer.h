
#ifndef _MOTION_OPTIMIZER_H_
#define _MOTION_OPTIMIZER_H_

#include "utils/utils.h"

#include "containers/ProblemDescription.h"
#include "containers/Trajectory.h"
#include "containers/ChompGradient.h"

#include "containers/ConstraintFactory.h"
#include "containers/Constraint.h"

#include "optimizer/ChompLocalOptimizer.h"
#include "optimizer/ChompOptimizer.h"


#ifdef NLOPT_FOUND
    #include <nlopt.hpp>
    #include "optimizer/NLOptimizer.h"
#endif


namespace chomp {
enum OptimizationAlgorithm {
    LOCAL_CHOMP,
    GLOBAL_CHOMP,
    MMA_NLOPT,
    CCSAQ_NLOPT,
    SLSQP_NLOPT,
    LBFGS_NLOPT,
    TNEWTON_PRECOND_RESTART_NLOPT,
    TNEWTON_RESTART_NLOPT,
    TNEWTON_NLOPT,
    VAR1_NLOPT,
    VAR2_NLOPT,
    NONE
};

nlopt::algorithm getNLoptAlgorithm( OptimizationAlgorithm alg );

class MotionOptimizer {
    

  private:

    ProblemDescription problem;
    
    ChompObserver * observer;

    int N_max, N_min;

    bool full_global_at_final, do_subsample;

    double obstol, timeout_seconds, alpha;
    size_t max_iterations;

    OptimizationAlgorithm algorithm1, algorithm2;

    const static char* TAG;
    
  public:
    //constructor.
    MotionOptimizer( ChompObserver * observer = NULL,
                     double obstol = 1e-8,
                     double timeout_seconds = 0,
                     size_t max_iter = size_t(-1),
                     const MatX & lower_bounds=MatX(0,0),
                     const MatX & upper_bounds=MatX(0,0),
                     OptimizationAlgorithm algorithm1 = VAR2_NLOPT,
                     OptimizationAlgorithm algorithm2 = VAR2_NLOPT,
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
    inline void setNMax( int n_max ){ N_max = n_max; }
    inline int  getNMax( ){ return N_max; }

    inline void setGoalset( Constraint * goal){ problem.setGoalset(goal);}
    inline const Constraint * getGoalset() const 
                              {return problem.getGoalset();}

    inline void   setTimeoutSeconds( double s ){ timeout_seconds = s; }
    inline double getTimeoutSeconds(){ return timeout_seconds; }

    inline void setMaxIterations( size_t max ){ max_iterations = max; }
    inline size_t getMaxIterations(){ return max_iterations; }

    inline void setFunctionTolerance( double tol ){ obstol = tol; }
    inline double getFunctionTolerance(){ return obstol; }
    
    inline void setAlgorithm(OptimizationAlgorithm a){ algorithm1 = a; }
    inline OptimizationAlgorithm getAlgorithm(){ return algorithm1; }

    inline void setAlgorithm1(OptimizationAlgorithm a1){ algorithm1 = a1; }
    inline void setAlgorithm2(OptimizationAlgorithm a2){ algorithm2 = a2; }
    
    inline OptimizationAlgorithm getAlgorithm1(){ return algorithm1; }
    inline OptimizationAlgorithm getAlgorithm2(){ return algorithm2; }

    inline void dontSubsample(){ do_subsample = false; }
    inline void doSubsample(){ do_subsample = true; }
    inline void setSubsample( bool subsample ){ do_subsample = subsample;}
    
    inline void setAlpha( double a ){ alpha = a; }
    inline double getAlpha(){ return alpha; }

    inline void doFullGlobalAtFinal(){ full_global_at_final = true; }
    inline void dontFullGlobalAtFinal(){ full_global_at_final = false; }
    inline bool getFullGlobalAtFinal(){ return full_global_at_final; }

    inline Trajectory & getTrajectory(){ return problem.trajectory; }
    inline const Trajectory & getTrajectory() const 
                              { return problem.getTrajectory();}
    inline void setTrajectory( const Trajectory & trajectory )
                { problem.trajectory = trajectory; }

    inline void setGradientHelper(ChompGradientHelper * helper)
                        { problem.gradient.setGradientHelper(helper); }
    
    inline void setObserver( ChompObserver * obs ){ observer = obs; }
    inline ChompObserver * getObserver(){ return observer; }
    inline const ChompObserver * getObserver() const { return observer; }

    inline void doCovariantOptimization()
                { problem.doCovariantOptimization(); }
    inline void dontCovariantOptimization()
                { problem.dontCovariantOptimization(); }
    inline bool isCovariantOptimization()
                { return problem.isCovariantOptimization(); }


};

}// namespace

#endif 
