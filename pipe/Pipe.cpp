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


virtual void PathConstructor::createPath( Plan * plan, MatX & trajectory ){
    trajectory << plan->start_state;
    for ( size_t i = 0; i < plan->phases.size(); i ++ ){
        
        //get the current phase
        Phase & current_phase = plan->phases[i];
        
        //set the beginning of the phase index:
        current_phase.start_index = trajectory.cols();

        //initialize the path as an empty matrix
        MatX path;
        
        //the end of the trajectory matrix is the last planned state,
        //  so it is the start state of the next phase.
        const MatX & start = trajectory.col( trajectory.cols() - 1 );
        assert( start.rows() == plan->start_state.rows() );

        planPhase( current_phase, start, path );

        //make sure that the path is not empty
        // TODO : this should be removed in the release mode
        assert( path.cols() != 0 );
        //make sure that the dimensions are correct, before they are joined
        assert( trajectory.rows() == path.rows() );
        
        trajectory << trajectory, path;
        current_phase.end_index = trajectory.cols() - 1;
    }
}


void Pipe::solve( bool doGlobalSmooth, bool doLocalSmooth ){
    
    MatX xi;
    pathConstructor_ptr->createPath( plan_ptr, xi );
    size_t Nmax = upsample * xi.rows();

    Chomp chomper( plan_ptr, xi,
                   plan_ptr->start_state, 
                   plan_ptr->end_state,
                   Nmax, alpha, errorTol );

    chomper.setDoDeleteConstraints( false );
    chomper.solve( doGlobalSmooth, doLocalSmooth );

}




