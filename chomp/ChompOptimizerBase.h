#ifndef _CHOMP_OPTIMIZER_BASE_H_
#define _CHOMP_OPTIMIZER_BASE_H_

#include "chomputil.h"

namespace chomp{

class ChompOptimizerBase{
    
  public:
    ConstraintFactory * factory;
    ChompGradient * gradient;
    ChompObserver * observer;
    ChompObjectiveType objective_type;

    //N: the current number of waypoints in the trajectory,
    //      not including the endpoints.
    //M: the degrees of freedom.
    int N, M;
    
    //the current trajectory of size N*M.
    MatX xi;
    
    
    ChompOptimizerBase( ConstraintFactory * f,
                        const MatX & xi,
                        const MatX & pinit,
                        const MatX & pgoal,
                      ChompObjectiveType object_type=MINIMIZE_ACCELERATION,
                        double total_time=1.0);

    ~ChompOptimizerBase();

    virtual void solve( bool global=true, bool local=true)=0;

    //notify the observer
    int notify(ChompEventType event,
               size_t iter,
               double curObjective,
               double lastObjective,
               double constraintViolation) const;

};

}//namespace

#endif
