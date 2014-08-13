
#include "ChompOptimizerBase.h"
#include "ChompGradient.h"
#include "Constraint.h"

namespace chomp {

ChompOptimizerBase::ChompOptimizerBase( ConstraintFactory * f,
                                        const MatX & xinit,
                                        const MatX & pinit,
                                        const MatX & pgoal,
                                        const MatX & lower_bounds,
                                        const MatX & upper_bounds,
                                        ChompObjectiveType object_type,
                                        double total_time ) :
    factory( f ), 
    observer( NULL ),
    objective_type( object_type ),
    N(xinit.rows()), M(xinit.cols()),
    xi( xinit ),
    lower_bounds( lower_bounds ), upper_bounds( upper_bounds )
{

    assert( pinit.size() == M );
    assert( pgoal.size() == M );
    
    //either the bounds vector is equivalent to the number of
    //  DOFs, or it is empty.
    assert( lower_bounds.size() == 0 || lower_bounds.size() == size_t(M) ); 
    assert( upper_bounds.size() == 0 || upper_bounds.size() == size_t(M) ); 

    gradient = new ChompGradient( pinit, pgoal, object_type, total_time );
}

ChompOptimizerBase::~ChompOptimizerBase(){
    if (gradient){ delete gradient; }
}

int ChompOptimizerBase::notify(ChompEventType event,
                               size_t iter,
                               double curObjective,
                               double lastObjective,
                               double constraintViolation) const
{
    if (observer) {
        return observer->notify(*this, event, iter, 
                                curObjective, lastObjective,
                                constraintViolation);
    } else {
        return 0;
    }
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



}//namespace
