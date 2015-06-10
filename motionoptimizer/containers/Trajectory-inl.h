/*
* Copyright (c) 2008-2015, Matt Zucker and Temple Price
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
    
template <class Derived>
void Trajectory::update(const Eigen::MatrixBase<Derived> & delta)
{
    xi -= delta;
}

template <class Derived>
void Trajectory::update(const Eigen::MatrixBase<Derived> & delta, int index)
{
    xi.row( index ) -= delta;
}

////////////////////TEMPLATE FUNCTIONS////
//initialize the data from a preexistant matrix
template <class Derived>
void Trajectory::initializeData( const Eigen::MatrixBase<Derived> & traj )
{

    q0 = traj.row(0);
    q1 = traj.row( traj.rows() - 1 );

    const int n = traj.rows() - 2;
    const int m = traj.cols();

    //if there is previously existant data,
    //  delete it.
    if (cached_data){ 
        delete cached_data;
        cached_data = NULL;
    } else if (data){
        delete data;
    }

    //allocate the data
    data = new double[ n * m ];

    //remap the data
    remapXi( n, n, m );

    //copy over the data from the trajectory
    full_xi = traj.block( 1, 0, n, m );

    dt = total_time / xi.rows()+1;
    
}

//definitions of the inline functions.
template <class Derived>
void Trajectory::initializeData( const Eigen::MatrixBase<Derived> & pinit,
                                 const Eigen::MatrixBase<Derived> & pgoal,
                                 int n )
{
    debug_status( TAG, "initializeData", "start" );

    q0 = pinit;
    q1 = pgoal;

    debug_status( TAG, "initializeData", "early middle" );
    //if there is previously existant data,
    //  delete it.
    if (cached_data){ 
        delete cached_data;
        cached_data = NULL;
    } else if (data){
        delete data;
    }
    
    debug_status( TAG, "initializeData", "middle" );
    
    const int m = q0.size();
    //allocate the data
    data = new double[ n * m ];

    //remap the data
    remapXi( n, n, m );

    debug_status( TAG, "initializeData", "before create initial" );
    
    //copy over the data from the trajectory
    createInitialTrajectory();

    dt = total_time / xi.rows()+1;
    
    debug_status( TAG, "initializeData", "end" );
}





