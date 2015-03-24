#ifndef _CHOMP_OPTIMIZER_BASE_H_
#define _CHOMP_OPTIMIZER_BASE_H_

#include "mzcommon/TimeUtil.h"
#include "OptimizerBase.h"

namespace chomp{

class ChompOptimizerBase : public OptimizerBase{
    
  public:
    
    //either global or local iteration, depending
    //on the type of the optimization, this should be assigned in 
    // inheriting class constructor    
    ChompEventType event; 

    static const char* TAG;

    double alpha;       // the gradient step size

    MatX bounds_violations;
    
    //timeout_seconds : the amount of time from the start of chomp
    //                  to a forced timeout.
    //canTimeout : is timing out a possible termination condition?
    // didTimeout : has chomp timed out ?
    // stop_time : the time at which chomp will timeout.
    bool canTimeout, didTimeout;
    TimeStamp stop_time;

    //TODO set to zero.
    size_t min_iter;
   
    bool use_momentum;
    MatX g, momentum;
    
    //an HMC object for performing the Hamiltonian Monte Carlo method
    HMC * hmc;
    
    ChompOptimizerBase(ProblemDescription & problem,
                       ChompObserver * observer = NULL,
                       double obstol = 1e-8,
                       double timeout_seconds = 0,
                       size_t max_iter = size_t(-1)); 

    virtual ~ChompOptimizerBase(){};
    
    void solve();

    inline void setAlpha( double a ){ alpha = a; }
    inline double getAlpha( ){ return alpha ; }

  protected:

    virtual void optimize()=0;

  private:
    
    //Checks the bounds of chomp, and smoothly pushes the trajectory
    //  back into the bounds.
    void checkBounds();
    
    //called for every interation. It returns true, if 
    //  the optimization is not finished.
    bool iterate();

    //check if chomp is finished.
    bool checkFinished(ChompEventType event);
    
    // returns true if performance has converged
    bool goodEnough(double oldObjective, double newObjective );

};


}//namespace

#endif
