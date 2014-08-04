
#include "ChompOptimizerBase.h"
#include "ChompGradient.h"

namespace chomp {

ChompOptimizerBase::ChompOptimizerBase( ConstraintFactory * f,
                                        const MatX & xinit,
                                        const MatX & pinit,
                                        const MatX & pgoal,
                                        ChompObjectiveType object_type,
                                        double total_time ) :
    factory( f ), 
    observer( NULL ),
    objective_type( object_type ),
    N(xinit.rows()), M(xinit.cols()),
    xi( xinit )
{
    std::cout << "xi_init: N: " << xinit.rows()
              << "\tM: " << xinit.cols() << "\n" << xinit <<std::endl;

    assert( pinit.size() == M );
    assert( pgoal.size() == M );
    
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

}//namespace
