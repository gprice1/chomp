
#ifndef _MOTION_OPTIMIZER_H_
#define _MOTION_OPTIMIZER_H_

#include "class_utils.h"
#include "Trajectory.h"

namespace chomp {

class MotionOptimizer {

  public:
 
    ConstraintFactory * factory;
    ChompGradient * gradient;

    Trajectory trajectory; //the trajectory
    
    int N_max, N_min, N_sub;

    bool use_goalset, full_global_at_final;

    Constraint * goalset;

    MatX upper_bounds, lower_bounds;
    
    //constructor.
    MotionOptimizer( 



         );

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
    virtual void setLowerBounds( const MatX & lower );
    virtual void setLowerBounds( const std::vector<double> & lower);
    virtual void setLowerBounds( const double * lower);

    virtual void setUpperBounds(const MatX & upper );
    virtual void setUpperBounds(const std::vector<double> & upper );
    virtual void setUpperBounds(const double * upper );

    virtual void setBounds( const MatX & lower, const MatX & upper );
    virtual void setBounds( const std::vector<double> & lower,
                            const std::vector<double> & upper );
    virtual void setBounds( const double * lower,
                            const double * upper );



};

}// namespace

#endif 
