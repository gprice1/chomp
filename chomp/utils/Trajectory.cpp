
#include "Trajectory.h"

namespace chomp {

Trajectory::Trajectory( const MatX & q0, const MatX & q1, int N,
                        ChompObjectiveType o_type, 
                        t_total ):
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( q0, q1, N );
}

Trajectory::Trajectory( const std::vector<double> & pinit,
                        const std::vector<double> & pgoal,
                        int N, ChompObjectiveType o_type, 
                        t_total ):
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( pinit, pgoal, N );
}

//row major matrix as a vector
Trajectory::Trajectory( const std::vector<double> & vec, int M,
                        ChompObjectiveType o_type, 
                        t_total ):
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( vec, M );
}


Trajectory::Trajectory(const MatX & xinit,
                       const MatX & pinit,
                       const MatX & qgoal,
                       ChompObjectiveType o_type, 
                       t_total ) : 
    xi( xinit ),
    xi_sub( NULL, 0, 0 ),
    q0( pinit ), q1( pgoal ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    dt = total_time / xi.rows()+1;
    
}
//Initialization from a whole prexisting trajectory.
Trajectory::Trajectory( const MatX & traj,
                        ChompObjectiveType o_type, 
                        double t_total) :
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( traj );
}

Trajectory::Trajectory( const std::vector<double> & traj,
                        int M,
                        ChompObjectiveType o_type, 
                        double t_total) : 
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( traj , M);
}

Trajectory::Trajectory(const std::vector<std::vector<double> > & traj,
                       ChompObjectiveType o_type, 
                       double t_total) :
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( traj );
}

Trajectory::Trajectory(const double * traj,
                       int M, int N,
                       ChompObjectiveType o_type, 
                       double t_total) :
    xi_sub( NULL, 0, 0 ),
    objective_type( o_type ),
    total_time( t_total ),
    subsample(false)
{
    initialize( traj, M, N );
}

Trajectory::Trajectory() : 
    xi_sub( NULL, 0, 0 ),
    objective_type( MINIMIZE_ACCELERATION ),
    total_time( 1.0 ),
    subsample(false)
{
}

void Trajectory::subsample(){
    do_subsample = true;
    new (&xi_sub) SubMatMap( xi, xi.rows()/2 + 1, xi.cols() );
}

void Trajectory::endSubsample(){
    do_subsample = false;
}

//copy the trajecotyr into the vector in row_major format.
void Trajectory::getTraj( std::vector<double> & vec ) const 
{
    vec.resize( xi.size() );
    MatMapR( vec.data(), xi.rows(), xi.cols() ) = xi;
}

//copy the trajecotyr into the vector in row_major format.
void Trajectory::getTraj( double * array ) const 
{
    MatMapR( array, xi.rows(), xi.cols() ) = xi;
}

//copy the trajecotyr into the vector in row_major format.
void Trajectory::getState( std::vector<double> & vec, int i ) const 
{
    vec.resize( xi.cols() );
    MatMap( vec.data(), 1, xi.cols() ) = xi.row( i );
}

//copy the trajecotyr into the vector in row_major format.
void Trajectory::getState( double * array, int i ) const 
{
    MatMap( array, 1, xi.cols() ) = xi.row( i );
}

void Trajectory::upsample()
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
    dt = total_time / xi.rows()+1;
}

void Trajectory::upsampleTo( int Nmax ){
    while ( xi.rows < Nmax ){ upsample; }
}

void Trajectory::constrainedUpsample(ConstraintFactory * factory,
                                     double htol,
                                     double hstep)
{
    upsample();
    
    //if there is no factory, or there are no constraints,
    //    do not evaluate the constraints.
    if ( factory ){ factory->getAll( xi.rows() ); }
    else{ continue; }
    if ( factory->constraints.empty() ){ continue; }


    MatX h, H, delta;
    double hinit = 0, hfinal = 0;

    for (int i=0; i < xi.rows(); i+=2) {
  
        Constraint* c = factory->constraints[i];

        if (!c || !c->numOutputs()) { continue; }
      
        for (int iter=0; ; ++iter) { 
            c->evaluateConstraints(xi.row(i), h, H);
            if (h.rows()) {
                double hn = h.lpNorm<Eigen::Infinity>();
                if (iter == 0) { hinit = std::max(hn, hinit); }
                if (hn < htol) { hfinal = std::max(hn, hfinal); break; }
                delta = H.colPivHouseholderQr().solve(h);
                updateTrajectory( hstep * delta.transpose(), i );
            }
        }
    }
}

void Trajectory::constraintUpsampleTo( ConstraintFactory * factory,
                                       double htol,
                                       double hstep,
                                       int Nmax)
{
    while (xi.rows() < Nmax ){
        constrainedUpsample( factory, htol, hstep );
    }
}

void Trajectory::initialize(const MatX & traj)
{
    q0 = traj.row( 0 );
    q1 = traj.row( traj.rows() - 1 );
    xi = traj.block( 1, 0, traj.rows()-2, traj.cols() );

    dt = total_time / xi.rows()+1;
}

void Trajectory::initialize(const std::vector<double> & traj, int M)
{
    int N = traj.size() / M;
    q0 = MatMapR( traj.data(), 1, M );
    q1 = MatMapR( traj.data() + M*(N-1), 1, M );
    xi = MatMapR( traj.data()+M, N-2, M );

    dt = total_time / xi.rows()+1;
}

void Trajectory::initialize(const std::vector< std::vector<double> > & traj)
{
    const int N = traj.size();
    const int M = traj[0].size();

    q0 = MatMapR( traj[0].data(), 1, M );
    q1 = MatMapR( traj.back().data(), 1, M );
    
    xi.resize( N-2, M );
    for( int i = 0; i < N-2; i ++ ){
        assert( traj[i+1].size() == M );
        xi.row( i ) = MatMap( traj[i+1].data(), 1 , M );
    }

    dt = total_time / xi.rows()+1;
    
}

void Trajectory::initialize(const double * traj, int M, int N)
{
    int N = vec.size() / M;
    q0 = MatMapR( traj, 1, M );
    q1 = MatMapR( traj + M*(N-1), 1, M );
    xi = MatMapR( traj+M, N-2, M );

    dt = total_time / xi.rows()+1;
    
}


void Trajectory::initialize(const MatX & pinit, const MatX & pgoal, int N)
{
    xi.resize( N, pinit.size() );
    initialize( pinit, pgoal );
}

void Trajectory::initialize(const MatX & pinit, const MatX & pgoal )
{
    q0 = pinit;
    q1 = pgoal;

    dt = total_time / xi.rows()+1;

    createInitialTraj();
}


void Trajectory::initialize( const std::vector<double> & pinit,
                             const std::vector<double> & pgoal )
{
    q0 = ConstMatMap( pinit.data(), 1, pinit.size() );
    q1 = ConstMatMap( pgoal.data(), 1, pgoal.size() );

    dt = total_time / xi.rows()+1;
    assert( xi.rows() > 0 );
    createInitialTraj();
}

void Trajectory::initialize( const std::vector<double> & pinit,
                             const std::vector<double> & pgoal,
                             int N)
{
    xi.resize( N, pinit.size() );
    initialize( pinit, pgoal );
}

void Trajectory::initialize( const double * pinit, 
                             const double * pgoal,
                             int M)
{
    q0 = ConstMatMap( pinit, 1, M );
    q1 = ConstMatMap( pgoal, 1, M );
    
    assert( xi.rows() > 0 );

    dt = total_time / xi.rows()+1;
    createInitialTraj();
}

void Trajectory::initialize( const double * q0,
                             const double * q1,
                             int M, int N)
{
    xi.resize( N, M );
    initialize( q0, q1, M );
}


void Trajectory::createInitialTrajectory(); 
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
