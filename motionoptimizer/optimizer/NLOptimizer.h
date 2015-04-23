

#ifndef _NLOPTIMIZER_H_
#define _NLOPTIMIZER_H_

#include <nlopt.hpp>

#include "OptimizerBase.h"

namespace mopt {

class NLOptimizer : public OptimizerBase{
  
  public:
    
    //the type of algorithm that the optimizer uses.
    nlopt::algorithm algorithm;
    //the type of end that the optimizer comes to.
    nlopt::result result;

    static const char * TAG;

    static const EventType event = NLOPT_ITER;

    NLOptimizer( ProblemDescription & problem,
                 Observer * observer,
                 double obstol = 1e-8,
                 double timeout_seconds = 0,
                 size_t max_iter = size_t(-1)); 

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
    static void constraintFunction(unsigned constraint_dim,
                                   double * h,
                                   unsigned n_by_m,
                                   const double * x,
                                   double* H,
                                   void *data);

};


}//namespace 
#endif 
