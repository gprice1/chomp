
#include "OptimizerBase.h"

namespace chomp {

OptimizerBase::OptimizerBase( Trajectory & traj, 
                              ConstraintFactory * f,
                              ChompGradient * g,
                              ChompObserver * o,
                              double obstol,
                              double timeout,
                              size_t max_iter,
                              const MatX & lower,
                              const MatX & upper):
    xi( traj ),
    factory( f ),
    gradient( g ), 
    observer( o ),
    obstol( obstol ),
    timeout_seconds( timeout ),
    max_iter( max_iter ),
    lower_bounds( lower ), upper_bounds( upper )
{
}

int OptimizerBase::notify(ChompEventType event,
                          size_t iter,
                          double curObjective,
                          double lastObjective,
                          double constraintViolation) const
{
    if (observer) {
        return observer->notify(*this, event, iter, 
                                curObjective, lastObjective,
                                constraintViolation);
    } else {
        return 0;
    }
}



}//namespace
