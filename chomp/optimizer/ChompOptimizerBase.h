#ifndef _CHOMP_OPTIMIZER_BASE_H_
#define _CHOMP_OPTIMIZER_BASE_H_

#include "chomputil.h"
#include "TimeUtil.h"

namespace chomp{

class ChompOptimizerBase : public OptimizerBase{
    
  public:

    double alpha;       // the gradient step size
    double objRelErrTol; //Objective function value relative to the
                         //previous value

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
    double timeout_seconds;
    bool canTimeout, didTimeout;
    TimeStamp stop_time;

    //TODO set to zero.
    size_t curr_iter, min_iter;
   
    bool use_momentum;
    MatX momentum;
    
    //an HMC object for performing the Hamiltonian Monte Carlo method
    HMC * hmc;

    ChompOptimizerBase( ConstraintFactory * f,
                        const MatX & xi,
                        const MatX & pinit,
                        const MatX & pgoal,
                        const MatX & lower_bounds=MatX(0,0),
                        const MatX & upper_bounds=MatX(0,0),
                      ChompObjectiveType object_type=MINIMIZE_ACCELERATION,
                        double total_time=1.0);

    virtual ~ChompOptimizerBase();

    virtual void solve( Trajectory & xi );

  protected:
    bool iterate( Trajectory & xi );
    virtual void prepareIter( const Trajectory & xi)=0;

    //notify the observer
    int notify(ChompEventType event,
               size_t iter,
               double curObjective,
               double lastObjective,
               double constraintViolation) const;
    
    //check if chomp is finished.
    bool checkFinished(ChompEventType event);
    
    // returns true if performance has converged
    bool goodEnough(double oldObjective, double newObjective );


};


}//namespace

#endif
