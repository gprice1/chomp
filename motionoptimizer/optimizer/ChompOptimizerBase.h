#ifndef _CHOMP_OPTIMIZER_BASE_H_
#define _CHOMP_OPTIMIZER_BASE_H_

#include "../utils/function_utils.h"
#include "../utils/class_utils.h"
#include "mzcommon/TimeUtil.h"
#include "OptimizerBase.h"

namespace chomp{

class ChompOptimizerBase : public OptimizerBase{
    
  public:

    double alpha;       // the gradient step size

    MatX bounds_violations;
    
    //TODO set hmag to zero in the initializer.
    // The magnitude of the constraint violations
    double hmag;
    
    //the previous and current objective function values.
    double last_objective, current_objective;
    
    //timeout_seconds : the amount of time from the start of chomp
    //                  to a forced timeout.
    //canTimeout : is timing out a possible termination condition?
    // didTimeout : has chomp timed out ?
    // stop_time : the time at which chomp will timeout.
    bool canTimeout, didTimeout;
    TimeStamp stop_time;

    //TODO set to zero.
    size_t curr_iter, min_iter;
   
    bool use_momentum;
    MatX momentum;
    
    //an HMC object for performing the Hamiltonian Monte Carlo method
    HMC * hmc;
    
    ChompEventType event; //either global or local iteration, depending
                          //on the type of the optimization.

    ChompOptimizerBase(  Trajectory & traj,
                         ConstraintFactory * factory,
                         ChompGradient * gradient,
                         ChompObserver * observer,
                         double obstol = 1e-8,
                         double timeout_seconds = 0,
                         size_t max_iter = size_t(-1),
                         const MatX & lower_bounds=MatX(0,0),
                         const MatX & upper_bounds=MatX(0,0)); 

    virtual ~ChompOptimizerBase(){};
    
    void solve();

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
