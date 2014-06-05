#include "Pipe.h"
class RRTConstructor : public PathConstructor{ 


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
    
    h(0) = mydot(qt,qt) - 4; 
    
    for(int i=0; i<2; i++){
      H(i) = 2*qt(i);
    }
  }
};



int main(){

    PathConstructor * path = NULL;

    //construct the plan
    Plan plan;

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
    
    consValue3.push_back( 2 );
    consValue3.push_back( 2 );

    
    ConstantConstraint goal1( consIndex, consValue1 ),
                       goal2( consIndex, consValue1 ),
                       goal3( consIndex, consValue1 );

    CircleConstraint circle;

    NullConstraint null; 


    plan.addPhase( null&, goal1& );
    plan.addPhase( circle&, goal2& );
    plan.addPhase( null&, goal3& );

    Pipe pipe( path, *plan, 50, 0.1, 0.001 );

    delete plan;

}
