

#ifndef _CHOMP_NLOPT_H_
#define _CHOMP_NLOPT_H_

#include "ChompGradient.h"
#include <nlopt.hpp>

namespace chomp {

class ChompNLopt{
  
  public:
    int N, M, N_max;
    int max_iter;
    double obstol, objective_value;

    ChompGradient * gradient;
    ChompObserver* observer;
    
    nlopt::opt * optimizer;

    MatX xi;
    nlopt::algorithm algorithm;

    std::vector< double > lower, upper;

    nlopt::result result;
    

    ChompNLopt(const MatX& xi_init,
               const MatX& pinit,
               const MatX& pgoal,
               double obstol = 0.01,
               int max_iter=0,
               int N_max = 0,
               const std::vector<double> & upper=std::vector<double>(),
               const std::vector<double> & lower=std::vector<double>());
    
    ~ChompNLopt();

    void solve();

  private:
    double optimize();

    void setLowerBounds( const std::vector<double> & lower );
    void setUpperBounds( const std::vector<double> & upper );

    void copyNRows( const std::vector<double> & original, 
                    std::vector<double> & result);

};


}//namespace chomp
#endif 
