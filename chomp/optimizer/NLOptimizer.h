

#ifndef _CHOMP_NLOPT_H_
#define _CHOMP_NLOPT_H_

#include <nlopt.hpp>

#include "ChompOptimizerBase.h"

namespace chomp {

class ChompNLopt : public ChompOptimizerBase{
  
  public:
    
    //the type of algorithm that the optimizer uses.
    nlopt::algorithm algorithm;
    //the type of end that the optimizer comes to.
    nlopt::result result;
    
    ChompNLopt(ConstraintFactory * factory,
               const MatX& xi_init,
               const MatX& pinit,
               const MatX& pgoal,
               double obstol = 0.01,
               int max_iter=0,
               int N_max = 0,
               ChompObjectiveType obj_t = MINIMIZE_ACCELERATION,
               const MatX & lower_bounds = MatX(0,0),
               const MatX & upper_bounds = MatX(0,0));

    ~ChompNLopt();

    void solve(Trajectory & xi);

  private:

    void giveBoundsToNLopt( nlopt::opt & optimizer );

    void prepareNLoptConstraints( nlopt::opt & optimizer );
    
    void copyNRows( const MatX & original_bounds, 
                    std::vector<double> & result);
  
};


}//namespace chomp
#endif 
