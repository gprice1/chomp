
#ifndef _MOTION_OPTIMIZER_H_
#define _MOTION_OPTIMIZER_H_

#include "utils/utils.h"

#include "containers/Trajectory.h"
#include "containers/ChompGradient.h"

#include "containers/ConstraintFactory.h"
#include "containers/Constraint.h"

#include "optimizer/ChompLocalOptimizer.h"
#include "optimizer/ChompOptimizer.h"

#if NLOPT_FOUND
    #include "optimizer/NLOptimizer.h"
#endif


namespace chomp {
enum OptimizationAlgorithm {
    CHOMP,
    THE_OTHER
};

class MotionOptimizer {

  public:

    Trajectory trajectory; //the trajectory
    ChompGradient gradient;
    ConstraintFactory factory;
    
    ChompObserver * observer;

    int N_max, N_min;

    bool use_goalset, full_global_at_final;

    Constraint * goalset;

    double obstol, timeout_seconds;
    size_t max_iterations;

    MatX lower_bounds, upper_bounds;
    
    OptimizationAlgorithm algorithm;

    const static char* TAG;
    
    //constructor.
    MotionOptimizer( ChompObserver * observer = NULL,
                     double obstol = 1e-8,
                     double timeout_seconds = 0,
                     size_t max_iter = size_t(-1),
                     const MatX & lower_bounds=MatX(0,0),
                     const MatX & upper_bounds=MatX(0,0),
                     OptimizationAlgorithm algorithm = CHOMP,
                     int N_max = 0);

    void solve();
    
    //sets up the factory, gradient, and optimizer for the current
    //  resolution.
    void optimize();

    //setbounds;

    //setTime.

    //setGoalSet
    void setGoalset( Constraint * goalset );

  private: 
    void prepareGoalSet();
    void finishGoalSet();

    
  public:

    //Functions for setting the upper and lower bounds of MotionOptimizer.
    void setLowerBounds( const MatX & lower );
    void setLowerBounds( const std::vector<double> & lower);
    void setLowerBounds( const double * lower);

    void setUpperBounds(const MatX & upper );
    void setUpperBounds(const std::vector<double> & upper );
    void setUpperBounds(const double * upper );

    void setBounds( const MatX & lower, const MatX & upper );
    void setBounds( const std::vector<double> & lower,
                    const std::vector<double> & upper );
    void setBounds( const double * lower,
                            const double * upper );

    //simple getters and setters
    inline void setNMax( int n_max ){ N_max = n_max; }
    inline int  getNMax( ){ return N_max; }

    inline void   setTimeoutSeconds( double s ){ timeout_seconds = s; }
    inline double getTimeoutSeconds(){ return timeout_seconds; }

    inline void setMaxIterations( size_t max ){ max_iterations = max; }
    inline size_t getMaxIterations(){ return max_iterations; }

    inline void setFunctionTolerance( double tol ){ obstol = tol; }
    inline double getFunctionTolerance(){ return obstol; }

    inline void setAlgorithm( OptimizationAlgorithm a ){ algorithm = a; }
    inline OptimizationAlgorithm getAlgorithm(){ return algorithm; }

};

}// namespace

#endif 
