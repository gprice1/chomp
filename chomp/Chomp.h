
#ifndef _CHOMP_H_
#define _CHOMP_H_

#include "chomputil.h"

namespace chomp {

class Chomp {

  public:
    ConstraintFactory * factory;
    ChompGradient * gradient;

    Trajectory xi; //the trajectory
    
    //constructor.
    Chomp( 



         );

    solve();
    
    //sets up the factory, gradient, and optimizer for the current
    //  resolution.
    optimize();

    //setbounds;

    //setTime.

    //setGoalSet
    void setGoalset( Constraint * goalset );
    void prepareGoalSet();
    void finishGoalSet();


    //Functions for setting the upper and lower bounds of chomp.
    virtual void setLowerBounds( const MatX & lower );
    virtual void setLowerBounds( const std::vector<double> & lower);
    virtual void setLowerBounds( const double * lower);

    virtual void setUpperBounds(const MatX & upper );
    virtual void setUpperBounds(const std::vector<double> & upper );
    virtual void setUpperBounds(const double * upper );

    virtual void setBounds( const MatX & lower, const MatX & upper );
    virtual void setBounds( const std::vector<double> & lower,
                            const std::vector<double> & upper );
    virtual void setBounds( const double * lower,
                            const double * upper );

};

}// namespace

#endif 
