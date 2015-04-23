#ifndef _OPTIMIZER_BASE_H_
#define _OPTIMIZER_BASE_H_

#include "../utils/utils.h"
#include "../utils/Observer.h"
#include "../containers/ProblemDescription.h"

namespace mopt {

class OptimizerBase{
   
    //this allows the observer to peek at the 
    //  problem member variables
    friend class Observer;
    
  public:
    
    ProblemDescription & problem;
    Observer * observer;

    double obstol, timeout_seconds;
    size_t max_iter;
    
    double last_objective, current_objective, constraint_magnitude;
    size_t current_iteration;

    
    OptimizerBase( ProblemDescription & problem,
                   Observer * observer = NULL,
                   double obstol = 1e-8,
                   double timeout_seconds = 0,
                   size_t max_iter = size_t(-1));

    virtual ~OptimizerBase(){}

    virtual void solve()=0;

    //notify the observer
    int notify(EventType event) const;
    
};

}//namespace

#endif
