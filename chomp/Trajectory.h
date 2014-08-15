
#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include <pthread>
#include "chomputil.h"

namespace chomp{

class Trajectory {
  public:
    MatX xi;
    SubMatMap xi_sub;
  
  private:

    //A mutex for locking the trajectory when updates are being made.
    //  Only needed if a concurrent thread wants data out of 
    //  chomp.
    pthread_mutex_t trajectory_mutex;
    bool use_mutex;

    bool do_subsample;

  public:
    Trajectory( const MatX & q0, const MatX & q1, N );
    Trajectory( const std::vector<double> & q0,
                const std::vector<double> & q1,
                int N);
    Trajectory(const MatX & xi);
    Trajectory();

    ~Trajectory();

    MatX const & getTraj(){ return xi; }
    MatX       & getTraj(){ return xi; }
    void getTraj( std::vector<double> & vec );
    
    void subsample();
    void endSubsample();
    
    inline Eigen::ConstRowXpr row( int i ) const { return xi.row( i ); }
    inline Eigen::ConstColXpr col( int i ) const { return xi.col( i ); }

    inline Eigen::RowXpr row( int i ) const { return xi.row( i ); }
    inline Eigen::ColXpr col( int i ) const { return xi.col( i ); }

    inline double & operator ()( int i, int j ){ return xi( i, j ); }
    inline double & operator ()( int i ){ return xi( i ); }
    
    inline int rows(){ return xi.rows(); }
    inline int cols(){ return xi.cols(); }

    void initMutex(){
        use_mutex = true;
        pthread_mutex_init( &trajectory_mutex, NULL );
    }

    inline void lock(){
        if (use_mutex){ pthread_mutex_lock( &trajectory_mutex );}
    }
    inline void unlock(){
        if (use_mutex){ pthread_mutex_unlock( &trajectory_mutex );}
    }
    
    void upsampleTrajectory(const MatX & q0,
                            const MatX & q1,
                            double dt,
                            ChompObjectiveType objective_type);

    virtual void initTraj(const MatX & q0, const MatX & q1, int N);
    virtual void initTraj( const MatX & q0, const MatX & q1 );

    virtual void initTraj( const std::vector<double> & q0,
                           const std::vector<double> & q1 );
    virtual void initTraj( const std::vector<double> & q0,
                           const std::vector<double> & q1,
                           int N);

    virtual void initTraj( const double * q0, const double * q1);
    virtual void initTraj( const double * q0, const double * q1, int N);

    //updates the trajectory via a matrix delta. Delta
    // must be the same size and shape as the trajectory,
    //  or the subsampled trajectory
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta );

    //updates the trajectory at the given row, by the vector delta,
    //  delta should have the same number of columns as the trajectory.
    template <class Derived>
    void updateTrajectory( const Eigen::MatrixBase<Derived> & delta,
                           int row);

  private:
    template <class Derived> 
    void createInitialTrajectory( const Eigen::MatrixBase<Derived> & q0,
                                  const Eigen::MatrixBase<Derived> & q1); 
    
};

template <class Derived>
inline void ChompOptimizer::updateTrajectory( 
                              const Eigen::MatrixBase<Derived> & delta)
{
    lockTrajectory();
    if ( do_subsample ){ xi_sub -= delta; }
    else{ xi -= delta; }
    unlockTrajectory();
}

template <class Derived>
inline void ChompOptimizer::updateTrajectory( 
                              const Eigen::MatrixBase<Derived> & delta,
                              int index)
{
    lock();
    if ( do_subsample ){ xi_sub.row( index ) -= delta; }
    else{ xi.row( index ) -= delta; }
    unlock();
}

}//namespace

#endif 
