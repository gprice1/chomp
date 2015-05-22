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





