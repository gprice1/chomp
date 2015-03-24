#ifndef _OPTIMIZER_BASE_H_
#define _OPTIMIZER_BASE_H_

#include "../utils/utils.h"
#include "../containers/ProblemDescription.h"

namespace chomp{

class OptimizerBase{
   
    //this allows the observer to peak at the 
    //  problem member function
    friend class ChompObserver;
    
  public:
    
    ProblemDescription & problem;
    ChompObserver * observer;

    double obstol, timeout_seconds;
    size_t max_iter;
    
    double last_objective, current_objective, constraint_magnitude;
    size_t current_iteration;

    
    OptimizerBase( ProblemDescription & problem,
                   ChompObserver * observer = NULL,
                   double obstol = 1e-8,
                   double timeout_seconds = 0,
                   size_t max_iter = size_t(-1));

    virtual ~OptimizerBase(){}

    virtual void solve()=0;

    //notify the observer
    int notify(ChompEventType event) const;
    
};

}//namespace

#endif
