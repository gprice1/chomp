
#ifndef _CHOMP_H_
#define _CHOMP_H_

#include "chomputil.h"

namespace chomp {

class Chomp {

    ConstraintFactory * factory;
    ChompGradient * gradient;
    ChompOptimizerBase * optimizer;

    MatX upper_bounds, lower_bounds;
    MatX xi; //the trajectory

    solve( bool global = true, bool local = false );
    
    //sets up the factory, gradient, and optimizer for the current
    //  resolution.
    optimize( bool global, bool local );

    //setbounds;

    //setTime.

    //setGoalSet
    void setGoalset( Constraint * goalset );
    void prepareGoalSet();
    void finishGoalSet();
    
    void constrainedUpsampleTo(int Nmax, double htol, double hstep);

    //
};

}// namespace

#endif 
