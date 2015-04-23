#ifndef _TEST_OPTIMIZER_H_
#define _TEST_OPTIMIZER_H_

#include "mzcommon/TimeUtil.h"
#include "OptimizerBase.h"

namespace mopt {

class TestOptimizer : public OptimizerBase{
    
  public:
    
    static const EventType event = CHOMP_ITER; 

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
   
    MatX P, HP, Y, W, delta;
    
    double * x_data, * g_data, *H_data, *h_data;
    MatMap x, g, H, h;
    
    //A cholesky solver for solving the constraint matrix.
    Eigen::LDLT<MatX> cholSolver;
    
    MatX violations;
    
    TestOptimizer(ProblemDescription & problem,
                   Observer * observer = NULL,
                   double obstol = 1e-8,
                   double timeout_seconds = 0,
                   size_t max_iter = size_t(-1)); 

    virtual ~TestOptimizer(){
        if ( x_data ){ delete x_data; }
        if ( g_data ){ delete g_data; }
        if ( h_data ){ delete h_data; }
        if ( H_data ){ delete H_data; }
    };
    
    void solve();

    inline void setAlpha( double a ){ alpha = a; }
    inline double getAlpha( ){ return alpha ; }


    
  private:
    //updates the trajectory
    void optimize();

    //Checks the bounds of chomp, and smoothly pushes the trajectory
    //  back into the bounds.
    void checkBounds();
    
    //called for every interation. It returns true, if 
    //  the optimization is not finished.
    bool iterate();

    //check if chomp is finished.
    bool checkFinished(EventType event);
    
    // returns true if performance has converged
    bool goodEnough(double oldObjective, double newObjective );

};


}//namespace

#endif
