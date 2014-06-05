#include "Pipe.h"

using namespace chomp;

Constraint * Plan::getConstraint(size_t t, size_t total){

    for ( size_t i = 0; i < phases.size() ; i ++ ){
        if ( phases[i].get_start_index() >= t &&
             phases[i].get_end_index() < t ){
            return phases[i].get_traj_constraint() ;
        }
        else if (phases[i].get_end_index() == t ){
            return phases[i].get_goal_constraint();
        }
    }
    
    std::cerr << "\nReached end of getConstraint " 
              << "without finding valid constraint";

    return NULL;
}

void Pipe::solve( bool doGlobalSmooth, bool doLocalSmooth ){
    
    MatX xi;
    pathConstructor_ptr->createPath( plan_ptr, xi );

    Chomp chomper( plan_ptr, xi,
                   plan_ptr->get_start_state(), 
                   plan_ptr->get_end_state(),
                   Nmax, alpha, errorTol );
}




