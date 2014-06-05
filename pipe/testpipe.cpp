#include "Pipe.h"

class RRTConstructor : public PathConstructor{ 

    virtual void planPhase( Phase & phase,
                            const MatX & start_state,
                            MatX & trajectory ){
    }

    virtual 
};

class LineConstructor : public PathConstructor{
    
  public: 
    // the length of the path
    int N;

    //a basic constructor
    LineConstructor( int N ) : N(N){}


    virtual void createPath( Plan * plan, MatX & trajectory){
        trajectory.resize( N, 2 );
        for (int i = 0; i < N; i++ ){
            trajectory.row(i) = (i+1)
                               * (plan->start_state - plan->end_state)
                               / (N+1)
                               + plan->start_state;
        }
        
        int step = N/plan->phases.size() 
        for (int i = 0; i < plan->phases.size(); i ++ ){
            plan->phases.start_index = 
    }

    virtual void planPhase( Phase & phase, const MatX & start_state,
                            MatX & trajectory )
    {
        std::cerr << "\nLineConstructor::planPhase should never be called"
                  << ", but it was called";
        exit(1);
    }
                        


};


//this is the circle constraint class from the thing.
class CircleConstraint : public Constraint {
public:

    virtual size_t numOutputs(){
        return 1;
    }

    virtual void evaluateConstraints(const MatX& qt,
                                           MatX& h, 
                                           MatX& H){

        assert((qt.cols()==1 && qt.rows()==2) ||
               (qt.cols()==2 && qt.rows()==1  ));

        h.conservativeResize(1,1);
        H.conservativeResize(1,2);
        
        h(0) = mydot(qt,qt) - 1; 
        
        for(int i=0; i<2; i++){
            H(i) = 2*qt(i);
        }
    }
};

void generatePlan( Plan & plan ){ 

    //setup the vectors for the first goal position
    std::vector<double> consValue1, consValue2, consValue3;
    std::vector<size_t> consIndex;

    consIndex.push_back(0);
    consIndex.push_back(1);
    
    //setup the positions of the goal vectors.
    consValue1.push_back(-1 );
    consValue1.push_back( 0 );

    consValue2.push_back( 1 );
    consValue2.push_back( 0 );
    
    consValue3.push_back(  2 );
    consValue3.push_back( -2 );

    //from the goal poses, create the constraints
    ConstantConstraint* goal1 = new ConstantConstraint( consIndex, consValue1 );
    ConstantConstraint* goal2 = new ConstantConstraint( consIndex, consValue2 );
    ConstantConstraint* goal3 = new ConstantConstraint( consIndex, consValue3 );

    //create the constraints for the trajectories
    CircleConstraint* circle = new CircleConstraint();
    NullConstraint* null = new NullConstraint(); 

    //Add the phases to the plan
    plan.addPhase( null,   goal1 );
    plan.addPhase( circle, goal2 );
    plan.addPhase( null,   goal3 );


    //set the start and end states 
    plan.start_state.resize(1,2);
    plan.end_state.resize(1,2);

    plan.start_state << -2,  2;
    plan.end_state   <<  2, -2;


}

int main(){
    
    size_t upsample = 2;

    bool doGlobalSmooth = true;
    bool doLocalSmooth = true;
    double alpha = 0.05;
    double errorTol = 1e-7;
    
    PathConstructor * path = new LineConstructor();
    //PathConstructor * path = new RRTConstructor();

    //construct the plan
    Plan plan;
    generatePlan( plan );

    //create the pipe
    Pipe pipe( path, &plan, upsample, alpha, errorTol );
    pipe.solve( doGlobalSmooth, doLocalSmooth );


}
