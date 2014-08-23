#ifndef _OPTIMIZER_BASE_H_
#define _OPTIMIZER_BASE_H_

#include "chomputil.h"

namespace chomp{

class OptimizerBase{
    
  public:
    Trajectory & xi;

    ConstraintFactory * factory;
    ChompGradient * gradient;
    ChompObserver * observer;

    double obstol, timeout_seconds;
    double last_objective, current_objective;
    size_t max_iter;

    //the current trajectory of size N*M.
    MatX lower_bounds, upper_bounds;
    
    OptimizerBase( Trajectory & xi,
                   ConstraintFactory * factory,
                   ChompGradient * gradient,
                   ChompObserver * observer,
                   double obstol = 1e-8,
                   double timeout_seconds = 0,
                   size_t max_iter = size_t(-1),
                   const MatX & lower_bounds=MatX(0,0),
                   const MatX & upper_bounds=MatX(0,0));

    virtual ~OptimizerBase();

    virtual void solve()=0;

    //notify the observer
    int notify(ChompEventType event,
               size_t iter,
               double curObjective,
               double lastObjective,
               double constraintViolation) const;
    


};

}//namespace

#endif
