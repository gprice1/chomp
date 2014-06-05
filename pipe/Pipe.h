#include "chomp/ConstraintFactory.h"
#include "chomp/Chomp.h"
#include "chomp/Constraint.h"
#include <stdio.h>

using namespace chomp;

class Phase{
  public:
    size_t start_index, end_index;
  private: 
    size_t _id;
    Constraint* trajectory_constraint_ptr;
    Constraint* goal_constraint_ptr;


  public:
    Phase( size_t id, Constraint * trajectory_constraint,
                       Constraint * goal_constraint ):
           _id(id), 
           trajectory_constraint_ptr( trajectory_constraint),
           goal_constraint_ptr( goal_constraint ) 
    {
    }
    

    size_t get_start_index(){ return start_index; }
    size_t get_end_index(){return end_index; }
    size_t get_id(){ return _id; }
    Constraint* get_traj_constraint(){return trajectory_constraint_ptr;}
    Constraint* get_goal_constraint(){return goal_constraint_ptr;}
    
};

class Plan : public ConstraintFactory{

  public:

    Plan(){}
    ~Plan(){}
    
    //holds all of the phases
    std::vector<Phase> phases;
    MatX start_state;
    MatX end_state;

    
    //adding phases
    void addPhase( Constraint * trajectory_constraint,
                   Constraint * goal_constraint ){
        Phase newPhase( phases.size(), trajectory_constraint,
                        goal_constraint);
        phases.push_back( newPhase );
    }

    Constraint* getConstraint( size_t t, size_t total);

        
};


class PathConstructor{ 
    public:
    //Trajectory matrix should be an n*m matrix that can be submitted to
    //chomp as the inital trajectory.
    //It is not necessary to make a create path function, because 
    // the default just calls plan phase for each of the phases, and
    // then chains the trajectories for each path into a complete trajectory
    virtual void createPath( Plan * plan, MatX & trajectory);
    
    //Plans a trajectory for a single phase. 
    virtual void planPhase( Phase & phase,
                            const MatX & start_state,
                            MatX & trajectory ) = 0;

};


class Pipe{

  public: 

    /////////////////member variables
    PathConstructor * pathConstructor_ptr;
    Plan * plan_ptr;
    size_t upsample;
    double alpha, errorTol;


    ////////////////Public Member Functions//////////////////
    Pipe(PathConstructor* pathConstructor,
             Plan * plan,
             size_t upsample,
             double alpha = 0.01,
             double errorTol= 0.01) : 
        pathConstructor_ptr( pathConstructor),
        plan_ptr(plan),
        upsample(upsample),
        alpha(alpha),
        errorTol(errorTol)
    {
    }

    ~Pipe(){}

    
    void solve( bool doGlobalSmooth, bool doLocalSmooth );
};
