#ifndef _CHOMP_OPTIMIZER_BASE_H_
#define _CHOMP_OPTIMIZER_BASE_H_

#include "chomputil.h"

namespace chomp{

class ChompOptimizerBase{
    
  public:
    ConstraintFactory * factory;
    ChompGradient * gradient;
    ChompObserver * observer;
    ChompObjectiveType objective_type;
    
    //N: the current number of waypoints in the trajectory,
    //      not including the endpoints.
    //M: the degrees of freedom.
    int N, M;

    //the current trajectory of size N*M.
    MatX xi;
    MatX lower_bounds, upper_bounds;
    
    ChompOptimizerBase( ConstraintFactory * f,
                        const MatX & xi,
                        const MatX & pinit,
                        const MatX & pgoal,
                        const MatX & lower_bounds=MatX(0,0),
                        const MatX & upper_bounds=MatX(0,0),
                      ChompObjectiveType object_type=MINIMIZE_ACCELERATION,
                        double total_time=1.0);

    virtual ~ChompOptimizerBase();

    virtual void solve( bool global=true, bool local=true)=0;

    //notify the observer
    int notify(ChompEventType event,
               size_t iter,
               double curObjective,
               double lastObjective,
               double constraintViolation) const;
    
    //Functions for setting the upper and lower bounds of chomp.
    virtual void setLowerBounds( const MatX & lower );
    virtual void setLowerBounds( const std::vector<double> & lower);

    virtual void setUpperBounds(const MatX & upper );
    virtual void setUpperBounds(const std::vector<double> & upper );

    virtual void setBounds( const MatX & lower, const MatX & upper );
    virtual void setBounds( const std::vector<double> & lower,
                            const std::vector<double> & upper );
    

};

}//namespace

#endif
