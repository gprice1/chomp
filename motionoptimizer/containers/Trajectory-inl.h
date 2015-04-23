

//lots of simple inline getters and setters.
inline double * Trajectory::getData(){ return data; }
inline const double * Trajectory::getData() const { return data; } 

inline const DynamicMatMap & Trajectory::getTraj()const { return xi; }

inline double Trajectory::getDt() const { return dt; }

inline bool Trajectory::isSubsampled() const { return is_subsampled; }

inline double & Trajectory::operator ()( int i, int j ){
    return xi( i, j );
}
inline double const & Trajectory::operator()(int i, int j)const {
    return xi( i, j );
}

inline int Trajectory::N() const { return xi.rows(); }
inline int Trajectory::M() const { return xi.cols(); }

inline int Trajectory::fullN() const { return full_xi.rows(); }
inline int Trajectory::fullM() const { return full_xi.cols(); }

inline int Trajectory::rows() const { return xi.rows();}
inline int Trajectory::cols() const { return xi.cols();}

inline int Trajectory::size() const { return xi.size();}

inline const MatX & Trajectory::getQ0() const { return q0; }
inline const MatX & Trajectory::getStartPoint() const { return q0; }

inline const MatX & Trajectory::getQ1() const { return q1; };
inline const MatX & Trajectory::getEndPoint() const { return q1; }

//get individual matrix blocks or rows
inline Row Trajectory::row( int i ) { return xi.row(i); }
inline ConstRow Trajectory::row( int i ) const{ return xi.row(i); }

inline Col Trajectory::col( int i ){  return xi.col(i); }
inline ConstCol Trajectory::col( int i ) const { return xi.col(i); }

inline const DynamicMatMap & Trajectory::getXi() const { return xi; } 
inline const MatMap & Trajectory::getFullXi() const { return full_xi; } 

inline void Trajectory::setObjectiveType( ObjectiveType otype)
{
    objective_type = otype;
}
inline ObjectiveType Trajectory::getObjectiveType() const
{
    return objective_type;
}


inline void Trajectory::remapXi( int n, int full_n, int m ){
    
    //find out which inner stride to use.
    //  if full_n is not equal to n, then we are subsampling,
    //  and the correct stride is 2, otherwise, the correct stride is 1.
    const int inner_stride = full_n == n ? 1 : 2;
    
    new (&xi) DynamicMatMap( data, n, m, 
                             DynamicStride(full_n, inner_stride) );
    new (&full_xi) MatMap( data, full_n, m);
}

////////////////////MORE COMPLICATED INLINE AND TEMPLATE FUNCTIONS////
//initialize the data from a preexistant matrix
template <class Derived>
inline void Trajectory::initializeData( 
                               const Eigen::MatrixBase<Derived> & traj )
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
inline void Trajectory::initializeData( 
                               const Eigen::MatrixBase<Derived> & pinit,
                               const Eigen::MatrixBase<Derived> & pgoal,
                               int n )
{
    q0 = pinit;
    q1 = pgoal;

    const int m = q0.size();

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
    createInitialTrajectory();

    dt = total_time / xi.rows()+1;
    
}


inline MatX Trajectory::getTick(int tick) const
{
    //if the tick is negative, get a state that falls off the
    //  the edge of the trajectory
    if (tick < 0) { return getPos(q0, (tick+1)*dt); }
    
    //if the tick is larger than the number of states,
    //  get a state that falls off the positive edge of the
    //  trajectory
    if (tick >= N() ) { return getPos(q1, (tick-N())*dt);}

    //if the tick corresponds to a state in the trajectory,
    //  return the corresponding state.
    return xi.row( tick );
} 


template <class Derived>
inline void Trajectory::update(const Eigen::MatrixBase<Derived> & delta)
{
    xi -= delta;
}

template <class Derived>
inline void Trajectory::update(const Eigen::MatrixBase<Derived> & delta,
                               int index)
{
    xi.row( index ) -= delta;
}


