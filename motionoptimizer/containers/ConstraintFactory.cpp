/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#include "ConstraintFactory.h"
#include "Constraint.h"

namespace chomp {

ConstraintFactory::ConstraintFactory( Trajectory & trajectory ):
    trajectory( trajectory ),
    interval_is_sorted( false )
{
}

ConstraintFactory::~ConstraintFactory(){
}

//TODO this may not be needed anymore
void ConstraintFactory::clearConstraints() {

    for ( std::vector<Constraint*>::iterator it = constraints.begin(); 
        it != constraints.end();
        it ++ )
    {
        if (*it){ delete *it; }
    }
    constraints.clear();
}

void ConstraintFactory::sortIntervals(){

    std::sort( constraint_intervals.begin(),
               constraint_intervals.end(),
               ConstraintInterval::compare );

    interval_is_sorted = true;

}

void ConstraintFactory::addConstraint( Constraint * constraint, 
                                       double start,
                                       double stop )
{
    constraint_intervals.resize( constraint_intervals.size() + 1 );
    constraint_intervals.back().start = start;
    constraint_intervals.back().stop = stop;
    constraint_intervals.back().constraint = constraint;
    interval_is_sorted = false;
}

//TODO make this able to take more constraint
Constraint* ConstraintFactory::getConstraint(size_t t, size_t total){
    
    if ( constraint_intervals.size() == 0 ){ return NULL; }

    //if the interval vector is not sorted, sort it.
    if ( !interval_is_sorted ){ sortIntervals(); }
    
    double time = double(t) / double(total);

    //iterate over all of the constriant intervals;
    for (std::vector<ConstraintInterval>::const_iterator it
            = constraint_intervals.begin();
         it != constraint_intervals.end();
         ++ it )
    {
        if ( it->start <= time ){
            if ( time <= it->stop ){ return it->constraint; }
            else { return NULL; }
        }
    }

    return NULL;
}

void ConstraintFactory::getAll(size_t total) {

    if ( constraint_intervals.size() == 0 ){ return; }
    
    //fill up the constraints vector.
    constraints.assign( total, NULL );
    
    //if the interval vector is not sorted, sort it.
    if ( !interval_is_sorted ){ sortIntervals(); }
    
    //iterate over all of the constriant intervals;
    for (std::vector<ConstraintInterval>::const_iterator it
            = constraint_intervals.begin();
         it != constraint_intervals.end();
         ++ it )
    {
        const int start_point = int( it->start * total );
        const int end_point = int( it->stop * total );
        
        for ( int i = start_point; i <= end_point; i++ ){
            constraints[i] = it->constraint;
        }
    }
}


size_t ConstraintFactory::numOutput()
{

    if ( constraint_intervals.size() == 0 ){ return 0; }
    
    //compute the total dimensionality of the constraints
    int constraint_dims = 0;
    for ( std::vector<Constraint*>::iterator it = constraints.begin();
        it != constraints.end();
        it ++ )
    {
        //if there is a constraint, get its dimensions.
        if ( *it ) { constraint_dims += (*it)->numOutputs(); }
    }

    return constraint_dims;
}

void ConstraintFactory::evaluate( unsigned constraint_dim,
                                double* result,
                                unsigned n_by_m,
                                const double * x,
                                double* grad )
{
    assert( int( n_by_m ) == trajectory.size() );
    assert( constraint_dim == numOutput() );

    //put the data into matrix maps.
    trajectory.setData( x );
    MatMap h_total(result, constraint_dim, 1);

    if ( grad == NULL ){
        evaluate( h_total );
    }else{
        MatMap H_total( grad, n_by_m, constraint_dim);

        //evaluate the data.
        evaluate( h_total, H_total);
    }
}

} //namespace
