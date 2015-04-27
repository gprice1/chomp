#ifndef _STOMP_OPTIMIZER_H_
#define _STOMP_OPTIMIZER_H_

#include "mzcommon/TimeUtil.h"
#include "OptimizerBase.h"

namespace mopt {

class StompOptimizer : public OptimizerBase{
    
  //private member variables
  private:
    //on the type of the optimization, this should be assigned in 
    // inheriting class constructor    
    static const EventType event = STOMP_ITER; 

    static const char* TAG;
    
    const int K;
    double metric_scale;
    
    MatX noisy_samples;
    MatX timestep_cost;
    MatX maximums, minimums;
    MatX q0, q1, q2, cspace_vel, cspace_accel;

  public:

    StompOptimizer(ProblemDescription & problem,
                   Observer * observer = NULL,
                   double obstol = 1e-8,
                   double timeout_seconds = 0,
                   size_t max_iter = size_t(-1)); 

    virtual ~StompOptimizer(){};

    void solve();

  private:
    void updateTrajectory();

    void generateNoisyTrajectories();

    void calculatePerTimestepCost();

    double calculateCollisionCost( const MatX & state );

    double calculateConstraintCost( const MatX & state, int timestep);

};

} //namespace 

#endif
