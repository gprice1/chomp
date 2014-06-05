#include "chomp/ConstraintFactory.h"
#include "chomp/Chomp.h"
#include "chomp/Constraint.h"
#include <stdio.h>

using namespace chomp;

class Phase{
  private: 
    size_t _id;
    size_t start_index, end_index;
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
  private: 
     MatX _start_state;
     MatX _end_state;

  public:

    Plan(){}
    ~Plan(){}
    
    //holds all of the phases
    std::vector<Phase> phases;
    
    
    //adding phases
    void addphase( Constraint * trajectory_constraint,
                   Constraint * goal_constraint ){
        Phase newPhase( phases.size(), trajectory_constraint,
                        goal_constraint);
        phases.push_back( newPhase );
    }

    //getters and setters for the start and end states
    void set_start_state( MatX& start_state )
                        { _start_state = start_state;}
    void set_end_state( MatX& end_state )
                        { _end_state = end_state;}
    MatX get_start_state(){ return _start_state; }
    MatX get_end_state(){ return _end_state;}


    Constraint* getConstraint( size_t t, size_t total);

        
};


class PathConstructor{ 
    public:
    //Trajectory matrix should be an n*m matrix that can be submitted to
    //chomp as the inital trajectory.

    virtual void createPath( const Plan * plan, MatX & trajectory) = 0;
};


class Pipe{

  public: 

    /////////////////member variables
    PathConstructor * pathConstructor_ptr;
    Plan * plan_ptr;
    int Nmax;
    double alpha, errorTol;


    ////////////////Public Member Functions//////////////////
    Pipe(PathConstructor* pathConstructor,
             Plan * plan,
             size_t Nmax,
             double alpha = 0.01,
             double errorTol= 0.01) : 
        pathConstructor_ptr( pathConstructor),
        plan_ptr(plan),
        Nmax(Nmax),
        alpha(alpha),
        errorTol(errorTol)
    {
    }

    ~Pipe(){}

    
    void solve( bool doGlobalSmooth, bool doLocalSmooth );
};
