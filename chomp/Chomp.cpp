
#include "Chomp.h"


namespace chomp {

void Chomp::solve(){

    //optimize at the current resolution
    optimize();

    //if the current resolution is not the final resolution,
    //  upsample and then optimize. Repeat until the final resolution
    //  is reached or exceeded.
    while ( xi.rows() < N_max ){

        //upsample the trajectory and prepare for the next
        //  stage of optimization.
        xi.upsample( gradient->q0, gradient->q1, gradient->dt,
                     gradient->objective_type );

        optimize();
    }
}


void ChompOptimizer::optimize(){
    
    if ( use_goalset ){ prepareGoalSet(); }
    gradient->prepareRun( N, use_goalset, subsample );

    // Subsample
    bool subsample = N > minN && !use_goalset &&
                     !(full_global_at_final && N >= maxN);
    if( subsample ){ xi.subsample(); }

}

void ChompOptimizer::useGoalSet( Constraint * goalset ){
      this->goalset = goalset;
      use_goalset = true;
}

void ChompOptimizer::prepareGoalSet(){
    
    //do not subsample if doing goalset run.
    N_sub = 0;

    //resize xi, and add q1 into it.
    xi.conservativeResize( xi.rows() + 1, xi.cols() );
    xi.row( xi.rows() - 1 ) = gradient->q1;
    
    //set N to the current size of xi.
    N = xi.rows();

    //add the goal constraint to the constraints vector.
    factory->constraints.push_back( goalset );
}

void ChompOptimizer::finishGoalSet(){
    
    use_goalset = false;
    
    //copy the last state in the trajectory back into q1
    gradient->q1 = xi.row( xi.rows() - 1 );

    //resize xi, keeping old values, and get rid of the 
    //  final state. And set N to the correct trajectory size
    xi.conservativeResize( xi.rows() -1, xi.cols() );

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory->constraints.pop_back();

}


void ChompOptimizerBase::setLowerBounds( const MatX & lower )
{
    assert( lower.size() == M );
    lower_bounds = lower;
}
void ChompOptimizerBase::setLowerBounds( const std::vector<double> & lower)
{
    assert( lower.size() == size_t(M) );
    lower_bounds = ConstMatMap(lower.data(), M, 1 );
}
void ChompOptimizerBase::setLowerBounds( const double * lower)
{
    lower_bounds = ConstMatMap(lower, M, 1 );
}


void ChompOptimizerBase::setUpperBounds(const MatX & upper )
{
    assert( upper.size() == M );
    upper_bounds = upper;
}
void ChompOptimizerBase::setUpperBounds( const std::vector<double> & upper)
{
    assert( upper.size() == size_t(M) );
    upper_bounds = ConstMatMap(upper.data(), M, 1 );
}
void ChompOptimizerBase::setUpperBounds( const double * upper)
{
    upper_bounds = ConstMatMap(upper, M, 1 );
}


void ChompOptimizerBase::setBounds( const MatX & lower, const MatX & upper )
{
    setLowerBounds( lower );
    setUpperBounds( upper );
}
void ChompOptimizerBase::setBounds( const std::vector<double> & lower, 
                                    const std::vector<double> & upper){
    setLowerBounds( lower );
    setUpperBounds( upper );
}
void ChompOptimizerBase::setBounds( const double* lower, 
                                    const double* upper){
    setLowerBounds( lower );
    setUpperBounds( upper );
}






}//namespace
