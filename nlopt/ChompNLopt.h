

#ifndef _CHOMP_NLOPT_H_
#define _CHOMP_NLOPT_H_

#include <nlopt.hpp>

#include "ChompOptimizerBase.h"

namespace chomp {

class ChompNLopt : public ChompOptimizerBase{
  
  public:
    int N_max, max_iter;
    double obstol, objective_value;
    
    nlopt::opt * optimizer;
    nlopt::algorithm algorithm;
    nlopt::result result;
    
    std::vector< double > lower, upper;

    

    ChompNLopt(const MatX& xi_init,
               const MatX& pinit,
               const MatX& pgoal,
               double obstol = 0.01,
               int max_iter=0,
               int N_max = 0,
               ChompObjectiveType obj_t = MINIMIZE_ACCELERATION,
               const std::vector<double> & upper=std::vector<double>(),
               const std::vector<double> & lower=std::vector<double>());
    
    ~ChompNLopt();

    void solve(bool global=true, bool local=true);

    void setLowerBounds( const std::vector<double> & lower );
    void setUpperBounds( const std::vector<double> & upper );

  private:
    double optimize();

    void giveBoundsToNLopt();


    void copyNRows( const std::vector<double> & original, 
                    std::vector<double> & result);

};


}//namespace chomp
#endif 
