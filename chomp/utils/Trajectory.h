

#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#define USE_PTHREAD 0

#if USE_PTHREAD
    #include <pthread.h>
#endif

#include "chomputil.h"

namespace chomp{

class Trajectory {
  public:
    MatX xi;
    SubMatMap xi_sub;

  private:
    bool do_subsample;

  public:
    Trajectory( const MatX & q0, const MatX & q1, int N );
    Trajectory( const std::vector<double> & q0,
                const std::vector<double> & q1,
                int N);
    Trajectory(const MatX & xi);
    Trajectory();

    ~Trajectory();

    MatX const & getTraj()const { return xi; }
    MatX       & getTraj(){ return xi; }
    void getTraj( std::vector<double> & vec ) const ;
    void getTraj( double * array ) const ;

    void getState( std::vector<double> & vec, int i ) const ;
    void getState( double * array, int i ) const ;

    void subsample();
    void endSubsample();

    bool isSubsampled(){ return do_subsample; }
    
    inline double & operator ()( int i, int j ){
        if (do_subsample){ return xi_sub( i, j );}
        return xi( i, j );
    }
    inline double & operator ()( int i ){ 
        if (do_subsample){ return xi_sub( i );}
        return xi( i );
    }
    
    inline double const & operator()(int i, int j)const {
        if (do_subsample){ return xi_sub( i, j );}
        return xi( i, j );
    }
    inline double const & operator()( int i ) const {         
        if (do_subsample){ return xi_sub( i );}
        return xi( i );
    }

    inline int N() const { 
        if ( do_subsample ){ return xi_sub.rows(); }
        return xi.rows();
    }
    inline int M() const {        
        if ( do_subsample ){ return xi_sub.cols(); }
        return xi.cols();
    }

    inline int rows() const { 
        if ( do_subsample ){ return xi_sub.rows(); }
        return xi.rows();
    }
    inline int cols() const {        
        if ( do_subsample ){ return xi_sub.cols(); }
        return xi.cols();
    }

    inline int size() const {         
        if ( do_subsample ){ return xi_sub.size(); }
        return xi.size();
    }
    
    virtual inline MatX const & getXi() const { return xi; } 
    virtual inline MatX & getXi() { return xi; } 
    
    virtual inline SubMatMap const & getXiSub() const { return xi_sub; }
    virtual inline SubMatMap & getXiSub() { return xi_sub; }
    
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
    if ( do_subsample ){ xi_sub -= delta; }
    else{ xi -= delta; }
}

template <class Derived>
inline void ChompOptimizer::updateTrajectory( 
                              const Eigen::MatrixBase<Derived> & delta,
                              int index)
{
    if ( do_subsample ){ xi_sub.row( index ) -= delta; }
    else{ xi.row( index ) -= delta; }
}

}//namespace

#endif 
