
#ifndef _MOTION_OPTIMIZER_H_
#define _MOTION_OPTIMIZER_H_

#include "utils/class_utils.h"
#include "utils/function_utils.h"

#include "containers/Trajectory.h"
#include "containers/ChompGradient.h"

#include "constraint/ConstraintFactory.h"
#include "constraint/Constraint.h"

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
 
    ConstraintFactory * factory;
    ChompObserver * observer;

    int N_max, N_min, N_sub;

    bool use_goalset, full_global_at_final;

    Constraint * goalset;

    MatX upper_bounds, lower_bounds;
    
    OptimizationAlgorithm algorithm;

    double obstol, timeout_seconds;
    int max_iterations;

    //constructor.
    MotionOptimizer( ConstraintFactory * factory = NULL,
                     ChompObserver * observer = NULL,
                     double obstol = 1e-8,
                     double timeout_seconds = 0,
                     size_t max_iter = size_t(-1),
                     const MatX & lower_bounds=MatX(0,0),
                     const MatX & upper_bounds=MatX(0,0),
                     OptimizationAlgorithm algorithm = CHOMP);

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

    inline void setAlgorithm(OptimizationAlgorithm alg ){ algorithm = alg;}

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



};

}// namespace

#endif 
