
#include "Trajectory.h"

namespace chomp {

Trajectory::Trajectory( const MatX & q0, const MatX & q1, N ):
    xi_sub( NULL, 0, 0 ),
    use_mutex(false),
    subsample(false)
{
    initTraj( q0, q1, N );
}

Trajectory::Trajectory( const std::vector<double> & q0,
                        const std::vector<double> & q1,
                        int N):
    xi_sub( NULL, 0, 0 ),
    use_mutex(false),
    subsample(false)
{
    initTraj( q0, q1, N );
}

//row major matrix as a vector
Trajectory::Trajectory( const std::vector<double> & vec, int N, int M):
    xi_sub( NULL, 0, 0 ),
    use_mutex(false),
    subsample(false)
{
    xi = MatMapR( vec.data, N, M );
}


Trajectory::Trajectory(const MatX & xi) : 
    xi( xi ),
    xi_sub( NULL, 0, 0 ),
    use_mutex(false),
    subsample(false)
{
}

Trajectory::Trajectory() : 
    xi_sub( NULL, 0, 0 ),
    use_mutex(false),
    subsample(false)
{
}

Trajectory::~Trajectory(){
    if( use_mutex ){ pthread_mutex_destroy( &trajectory_mutex ); }
}

void Trajectory::subsample(){
    do_subsample = true;
    new (&xi_sub) SubMatMap( xi, xi.rows()/2 + 1, xi.cols() );
}

void Trajectory::endSubsample(){
    do_subsample = false;
}

//copy the trajecotyr into the vector in row_major format.
void Trajectory::getTraj( std::vector<double> & vec )
{
    vec.resize( xi.rows() * xi.cols() );
    MatMapR( vec.data(), xi.rows(), xi.cols() ) = xi;
}

//copy the trajecotyr into the vector in row_major format.
void Trajectory::getState( std::vector<double> & vec, int i )
{
    vec.resize( xi.cols() );
    MatMap( vec.data(), 1, xi.cols() ) = xi.row( i );
}


void Trajectory::upsample(const MatX & q0,
                          const MatX & q1,
                          double dt,
                          ChompObjectiveType objective_type)
{

  const int N = xi.rows();
  const int M = xi.cols();
  const int N_up = 2*N+1; // e.g. size 3 goes to size 7
  
  MatX xi_up(N_up, M);

  // q0    d0    d1    d2    q1   with n = 3
  // q0 u0 u1 u2 u3 u4 u5 u6 q1   with n = 7
  //
  // u0 = 0.5*(q0 + d0)
  // u1 = d0
  // u2 = 0.5*(d0 + d1)
  // u3 = d1 
  // u4 = 0.5*(d1 + d2)
  // u5 = d2
  // u6 = 0.5*(d2 + q1)

  for (int t=0; t<N_up; ++t) { // t is timestep in new (upsampled) regime

    if (t % 2 == 0) {

      assert(t == N_up-1 || (t/2) < xi.rows());
      assert(t < xi_up.rows());

      if (objective_type == MINIMIZE_VELOCITY) {

        MatX qneg1 = getTickBorderRepeat(t/2-1, xi, q0, q1, dt);
        MatX qpos1 = getTickBorderRepeat(t/2,   xi, q0, q1, dt);
        xi_up.row(t) = 0.5 * (qneg1 + qpos1);

      } else { 

        MatX qneg3 = getTickBorderRepeat(t/2-2, xi, q0, q1, dt);
        MatX qneg1 = getTickBorderRepeat(t/2-1, xi, q0, q1, dt);
        MatX qpos1 = getTickBorderRepeat(t/2,   xi, q0, q1, dt);
        MatX qpos3 = getTickBorderRepeat(t/2+1, xi, q0, q1, dt);

        const double c3 = -1.0/160;
        const double c1 = 81.0/160;

        xi_up.row(t) = c3*qneg3 + c1*qneg1 + c1*qpos1 + c3*qpos3;

      }


    } else {
      xi_up.row(t) = xi.row(t/2);
    }

  }

  xi = xi_up;

}


void Trajectory::initTraj(const MatX & q0, const MatX & q1, int N)
{
    xi.resize( N, q0.size() );
    createInitialTraj( q0, q1 );
}
void Trajectory::initTraj(const MatX & q0, const MatX & q1 )
{
    createInitialTraj( q0, q1 );
}


void Trajectory::initTraj( const std::vector<double> & q0,
                           const std::vector<double> & q1 )
{
    ConstMatMap pinit( q0.data(), q0.size(), 1 );
    ConstMatMap pgoal( q1.data(), q1.size(), 1 );
    createInitialTraj( pinit, pgoal );
}

void Trajectory::initTraj( const std::vector<double> & q0,
                           const std::vector<double> & q1,
                           int N)
{
    xi.resize( N, q0.size() );
    initTraj( q0, q1 );
}
void Trajectory::initTraj( const double * q0, const double * q1, int M)
{
    ConstMatMap pinit( q0, M, 1 );
    ConstMatMap pgoal( q1, M, 1 );

    initTraj( pinit, pgoal );
}

void Trajectory::initTraj( const double * q0, const double * q1,
                           int M, int N)
{
    xi.resize( N, M );
    initTraj( q0, q1, M );
}


template <class Derived> 
void Trajectory::createInitialTrajectory(
                              const Eigen::MatrixBase<Derived> & q0,
                              const Eigen::MatrixBase<Derived> & q1); 
{
    assert( q0.size() == q1.size() );
    assert( xi.cols() == q0.size() );
    
    const int N = xi.rows();

    //if the goal is Minimize Velocity, simply linearly interpolate
    //  from the start to goal state.
    for ( int i=0; i < N; i ++ ){

        double t = double(i+1) / double(N+1);
        xi.row( i ) = q0 + (q1-q0)*t;
    }
}


}//namespace
