

#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include "chomputil.h"

namespace chomp{

class Trajectory {
  public:
    MatX xi;
    SubMatMap xi_sub;

    MatX q0, q1;
    ChompObjectiveType objective_type;

  private:
    double dt, total_time;
    bool do_subsample;

  public:
    Trajectory();
    
    Trajectory( const MatX & q0, const MatX & q1, int N,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory( const std::vector<double> & q0,
                const std::vector<double> & q1,
                int N,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory(const MatX & xi,
               const MatX & q0,
               const MatX & q1,
               ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);

    //Initialization from a whole prexisting trajectory.
    Trajectory( const MatX & traj,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory( const std::vector<double> & traj,
                int M,
                ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
                double t_total=1.0);
    Trajectory(const std::vector<std::vector<double> > & traj,
               ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);
    Trajectory(const double * traj,
               int M, int N,
               ChompObjectiveType o_type=MINIMIZE_ACCELERATION, 
               double t_total=1.0);

    ~Trajectory();

    MatX const & getTraj()const { return xi; }
    MatX       & getTraj(){ return xi; }
    void getTraj( std::vector<double> & vec ) const ;
    void getTraj( double * array ) const ;

    void getState( std::vector<double> & vec, int i ) const ;
    void getState( double * array, int i ) const ;
    
    inline double getDt(){ return dt; }

    void subsample();
    void endSubsample();

    bool isSubsampled(){ return do_subsample; }
    
    //simple element accessors. 
    double & operator ()( int i, int j );
    double const & operator()(int i, int j ) const;
    int N() const; 
    int M() const;
    int rows() const;
    int cols() const;
    int size() const;

    virtual inline MatX const & getXi() const { return xi; } 
    virtual inline MatX & getXi() { return xi; } 
    
    virtual inline SubMatMap const & getXiSub() const { return xi_sub; }
    virtual inline SubMatMap & getXiSub() { return xi_sub; }
    
    //upsample the trajectory by two times.
    void upsample();
    //upsample the trajectory until it is greater than Nmax.
    void upsampleTo( int Nmax );

    void constrainedUpsample( ConstraintFactory * factory,
                              double htol,
                              double hstep);
    void constrainedUpsampleTo( ConstraintFactory * factory,
                                double htol,
                                double hstep,
                                int Nmax );

    //initialize from a matrix.
    virtual void initialize(const MatX & traj);
    virtual void initialize(const std::vector<double> & traj, int M);
    virtual void initialize(const std::vector<std::vector<double> > & traj);
    virtual void initialize(const double * traj, int M, int N);
    
    //initialize from a two endpoints.
    virtual void initialize(const MatX & q0, const MatX & q1, int N);
    virtual void initialize( const MatX & q0, const MatX & q1 );

    virtual void initialize( const std::vector<double> & q0,
                           const std::vector<double> & q1 );
    virtual void initialize( const std::vector<double> & q0,
                           const std::vector<double> & q1,
                           int N);
    virtual void initialize( const double * q0, const double * q1);
    virtual void initialize( const double * q0, const double * q1, int N);

    //updates the trajectory via a matrix delta. Delta
    // must be the same size and shape as the trajectory,
    //  or the subsampled trajectory
    template <class Derived>
    void update( const Eigen::MatrixBase<Derived> & delta );

    //updates the trajectory at the given row, by the vector delta,
    //  delta should have the same number of columns as the trajectory.
    template <class Derived>
    void update( const Eigen::MatrixBase<Derived> & delta,
                           int row);

    inline MatX getTick(int tick);

  private:
    void createInitialTrajectory(); 
       
};

//definitions of the inline functions.

inline double & Trajectory::operator ()( int i, int j ){
    if (do_subsample){ return xi_sub( i, j );}
    return xi( i, j );
}
inline double const & Trajectory::operator()(int i, int j)const {
    if (do_subsample){ return xi_sub( i, j );}
    return xi( i, j );
}

inline int Trajectory::N() const { return xi.rows(); }
inline int Trajectory::M() const { return xi.cols(); }

inline int Trajectory::rows() const { 
    if ( do_subsample ){ return xi_sub.rows(); }
    return xi.rows();
}
inline int Trajectory::cols() const {        
    if ( do_subsample ){ return xi_sub.cols(); }
    return xi.cols();
}

inline int Trajectory::size() const {         
    if ( do_subsample ){ return xi_sub.size(); }
    return xi.size();
}

inline MatX Trajectory::getTick(int tick)
{
    //if the tick is negative, get a state that falls off the
    //  the edge of the trajectory
    if (tick < 0) { return getPos(q0, (tick+1)*dt); }
    
    //if the tick is larger than the number of states,
    //  get a state that falls off the positive edge of the
    //  trajectory
    const int rows = ( do_subsample ? xi_sub.rows() : xi.rows() );
    if (tick >= rows ) { return getPos(q1, (tick-rows))*dt;}

    //if the tick corresponds to a state in the trajectory,
    //  return the corresponding state.
    if( do_subsample ){ return xi_sub.row( tick ); }
    return xi.row( tick );
} 


template <class Derived>
inline void Trajectory::update( 
                              const Eigen::MatrixBase<Derived> & delta)
{
    if ( do_subsample ){ xi_sub -= delta; }
    else{ xi -= delta; }
}

template <class Derived>
inline void Trajectory::update( 
                              const Eigen::MatrixBase<Derived> & delta,
                              int index)
{
    if ( do_subsample ){ xi_sub.row( index ) -= delta; }
    else{ xi.row( index ) -= delta; }
}

}//namespace

#endif 
