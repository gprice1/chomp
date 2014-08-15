
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
    optimize();

    //setbounds;

    //setTime.

    //setGoalSet

    //
};

}// namespace

#endif 
