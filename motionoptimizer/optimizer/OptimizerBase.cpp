
#include "OptimizerBase.h"

namespace mopt {

OptimizerBase::OptimizerBase(ProblemDescription & problem,
                             Observer * observer,
                             double obstol,
                             double timeout,
                             size_t max_iter):
    problem( problem ),
    observer( observer ),
    obstol( obstol ),
    timeout_seconds( timeout ),
    max_iter( max_iter ),
    last_objective( HUGE_VAL ),
    current_objective( HUGE_VAL ),
    constraint_magnitude( HUGE_VAL ),
    current_iteration( 0 )
{
}

int OptimizerBase::notify(EventType event) const
{
    if (observer) {
        return observer->notify(*this, event, current_iteration, 
                                current_objective, last_objective,
                                constraint_magnitude);
    } else {
        return 0;
    }
}



}//namespace
