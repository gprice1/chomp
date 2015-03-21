

#ifndef _NLOPTIMIZER_H_
#define _NLOPTIMIZER_H_

#include <nlopt.hpp>

#include "OptimizerBase.h"

namespace chomp {

class NLOptimizer : public OptimizerBase{
  
  public:
    
    //the type of algorithm that the optimizer uses.
    nlopt::algorithm algorithm;
    //the type of end that the optimizer comes to.
    nlopt::result result;

    int iteration;

    static const ChompEventType event = NLOPT_ITER;

    NLOptimizer( Trajectory & traj,
                 ConstraintFactory * factory,
                 ChompGradient * gradient,
                 ChompObserver * observer,
                 double obstol = 1e-8,
                 double timeout_seconds = 0,
                 size_t max_iter = size_t(-1),
                 const MatX & lower_bounds=MatX(0,0),
                 const MatX & upper_bounds=MatX(0,0)); 

    ~NLOptimizer();

    inline void setAlgorithm( nlopt::algorithm alg ){ algorithm = alg; }

    void solve();

  private:

    void giveBoundsToNLopt( nlopt::opt & optimizer );

    void prepareNLoptConstraints( nlopt::opt & optimizer );
    
    void copyNRows( const MatX & original_bounds, 
                    std::vector<double> & result);

    //a wrapper function for passing the ChompGradient to NLopt.
    static double objectiveFunction(unsigned n,
                                    const double * x,
                                    double* grad,
                                    void *data);

};


}//namespace chomp
#endif 
