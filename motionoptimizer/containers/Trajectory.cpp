#include "Trajectory.h"
#include "ConstraintFactory.h"
#include "Constraint.h"

#include <sstream>

namespace chomp {

Trajectory::Trajectory( const MatX & q0, const MatX & q1, int N,
                        ChompObjectiveType o_type, 
                        double t_total ):
    xi( NULL, 0, 0 ),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( q0, q1, N );
}

Trajectory::Trajectory( const std::vector<double> & pinit,
                        const std::vector<double> & pgoal,
                        int N, ChompObjectiveType o_type, 
                        double t_total ):
    xi( NULL, 0, 0 ),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( pinit, pgoal, N );
}

//row major matrix as a vector
Trajectory::Trajectory( const std::vector<double> & vec,
                        int M,
                        ChompObjectiveType o_type, 
                        double t_total ):
    xi( NULL, 0, 0 ),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( vec, M );
}


Trajectory::Trajectory(const MatX & xinit,
                       const MatX & pinit,
                       const MatX & pgoal,
                       ChompObjectiveType o_type, 
                       double t_total ) : 
    xi( NULL, 0, 0 ),    
    q0( pinit ), q1( pgoal ),
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    
    const int N = xinit.rows();
    const int M = xinit.cols();

    data = new double[ N*M  ];
    new (&xi) MatMap( data, N, M);

    xi = xinit;
    
    dt = total_time / xi.rows()+1;
    
}
//Initialization from a whole prexisting trajectory.
Trajectory::Trajectory( const MatX & traj,
                        ChompObjectiveType o_type, 
                        double t_total) :
    xi( NULL, 0, 0 ),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( traj );
}


Trajectory::Trajectory(const std::vector<std::vector<double> > & traj,
                       ChompObjectiveType o_type, 
                       double t_total) :
    xi( NULL, 0, 0 ),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( traj );
}

Trajectory::Trajectory(const double * traj,
                       int M, int N,
                       ChompObjectiveType o_type, 
                       double t_total) :
    xi( NULL, 0, 0 ),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( traj, M, N );
}

Trajectory::Trajectory() : 
    xi( NULL, 0, 0 ),    
    objective_type( MINIMIZE_ACCELERATION ),
    total_time( 1.0 ),
    is_subsampled(false)
{
}


// increase the size of the trajectory by one waypoint, and
//  add in q1 as the last point of the trajectory.
void Trajectory::startGoalSet(){
    
    const int N = xi.rows();
    const int M = xi.cols();

    //copy the old data vector to a new data vector.
    double * new_data = new double [ size() + M ];
    
    //copy the old data over to the new data vector
    MatMap( new_data, N+1, M ).block(0,0, N, M) = xi;

    //delete the old data, and replace it with new_data.
    delete data;
    data = new_data;
    
    //set the xi matrix map.
    new (&xi) MatMap( data, N+1, M );
    xi.row( N ) = q1;
}


void Trajectory::endGoalSet(){
    const int N = xi.rows();
    const int M = xi.cols();

    q1 = xi.row( N-1 ); 
        
    //set the xi matrix map.
    new (&xi) MatMap( data, N-1, M );

    //there is no need to reallocate the array,
    //  because there is already enough room for it.
}


//This is dangerous, and it may be a bad idea.
void Trajectory::setData( const double * new_data ){
    
    data = const_cast<double*>(new_data);

    const int N = xi.rows();
    const int M = xi.cols();
    
    new (&xi) MatMap( data, N, M );
}

void Trajectory::setData( double * new_data ){
    
    data = new_data;

    const int N = xi.rows();
    const int M = xi.cols();
    
    new (&xi) MatMap( data, N, M );
}

void Trajectory::copyToData( const std::vector<double> & vec ){
    std::copy( vec.begin(), vec.end(), data );
}
void Trajectory::copyToData( const double * new_data ){
    std::copy( new_data, new_data + size(), data );
}


void Trajectory::copyDataTo( double * copy) const {
    std::copy( data, data + size(), copy );
}
//This is dangerous, and it may be a bad idea.
void Trajectory::copyDataTo( std::vector<double> & vec ) const {
    vec.resize( size() );
    std::copy( data, data + size(), vec.begin() );
}


void Trajectory::subsample(){

    assert( !is_subsampled );
    is_subsampled = true;

    const int n = N();
    const int n_sub = (n - 1)/2;
    const int m = M();

    sampled_data = new double [n_sub * m];
    MatMap sampled_xi( sampled_data, n_sub, m );

    //copy over the sampled data
    for ( int i = 0; i < n_sub; i ++ ){
        sampled_xi.row( i ) = xi.row( i*2 + 1 );
    }

    std::swap( sampled_data, data );
    new (&xi) MatMap( data, n_sub, m );
}

void Trajectory::endSubsample(){
    assert( is_subsampled );
    is_subsampled = false;

    const int n_sub = N();
    const int n = n_sub*2 + 1;
    const int m = M();

    MatMap sampled_xi( sampled_data, n, m );

    //copy over the newly calculated data. That is
    //  housed in xi, into the sized n sampled matrix.
    for ( int i = 0; i < n_sub; i ++ ){
        sampled_xi.row( i*2 + 1 ) = xi.row( i );
    }

    delete data;
    data = sampled_data;
    sampled_data = NULL;
    
    new (&xi) MatMap( data, n, m );

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
    
    double * upsampled_data = new double [N_up*M];
    MatMap xi_up(upsampled_data, N_up, M);

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

                MatX qneg1 = getTick(t/2-1);
                MatX qpos1 = getTick(t/2);
                xi_up.row(t) = 0.5 * (qneg1 + qpos1);

            } else { 

                MatX qneg3 = getTick(t/2-2);
                MatX qneg1 = getTick(t/2-1);
                MatX qpos1 = getTick(t/2  );
                MatX qpos3 = getTick(t/2+1);

                const double c3 = -1.0/160;
                const double c1 = 81.0/160;

                xi_up.row(t) = c3*qneg3 + c1*qneg1 + c1*qpos1 + c3*qpos3;

            }

        } else {
            xi_up.row(t) = xi.row(t/2);
        }

    }
    

    //clean up the old data array.
    delete data;
    
    //create the new matrix map, and set the data variable to the
    //  new data array. 
    new (&xi) MatMap( upsampled_data, N_up, M);
    data = upsampled_data;


    dt = total_time / xi.rows()+1;
}

void Trajectory::upsampleTo( int Nmax ){
    while ( xi.rows() < Nmax ){ upsample(); }
}

void Trajectory::constrainedUpsample(ConstraintFactory * factory,
                                     double htol,
                                     double hstep)
{
    upsample();
    
    //if there is no factory, or there are no constraints,
    //    do not evaluate the constraints.
    if ( factory && !factory->constraints.empty()  ){
        factory->getAll( xi.rows() ); 
    } else{ 
        return; 
    }


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
                update( hstep * delta.transpose(), i );
            }
        }
    }
}

void Trajectory::constrainedUpsampleTo( ConstraintFactory * factory,
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

    data = new double[ (traj.rows()-2) * traj.cols() ];
    new (&xi) MatMap( data, traj.rows()-2, traj.cols() );

    xi = traj.block( 1, 0, traj.rows()-2, traj.cols() );

    dt = total_time / xi.rows()+1;
}

void Trajectory::initialize(const std::vector<double> & traj, int M)
{
    int N = traj.size() / M;
    q0 = ConstMatMap( traj.data(), 1, M );
    q1 = ConstMatMap( traj.data() + M*(N-1), 1, M );

    data = new double[ traj.size()-2 ];
    new (&xi) MatMap( data, N-2, M );

    xi = ConstMatMapR( traj.data()+M, N-2, M );

    dt = total_time / xi.rows()+1;
}

void Trajectory::initialize(const std::vector< std::vector<double> > & traj)
{
    const int N = traj.size();
    const int M = traj[0].size();

    q0 = ConstMatMap( traj[0].data(), 1, M );
    q1 = ConstMatMap( traj.back().data(), 1, M );
    
    data = new double[ (N-2)*M  ];
    new (&xi) MatMap( data, N-2, M);

    for( int i = 0; i < N-2; i ++ ){
        assert( traj[i+1].size() == size_t(M) );
        xi.row( i ) = ConstMatMap( traj[i+1].data(), 1 , M );
    }

    dt = total_time / xi.rows()+1;
    
}

void Trajectory::initialize(const double * traj, int M, int N)
{
    q0 = ConstMatMap( traj, 1, M );
    q1 = ConstMatMap( traj + M*(N-1), 1, M );

    data = new double[ (N-2)*M  ];
    new (&xi) MatMap( data, N-2, M);

    xi = ConstMatMapR( traj+M, N-2, M );

    dt = total_time / xi.rows()+1;
    
}


void Trajectory::initialize(const MatX & pinit, const MatX & pgoal, int N)
{
    
    const int M = pgoal.size();
    
    data = new double[ N*M  ];
    new (&xi) MatMap( data, N, M);

    initialize( pinit, pgoal );
}

void Trajectory::initialize(const MatX & pinit, const MatX & pgoal )
{
    q0 = pinit;
    q1 = pgoal;

    dt = total_time / xi.rows()+1;

    createInitialTrajectory();
}

void Trajectory::initialize( const std::vector<double> & pinit,
                             const std::vector<double> & pgoal,
                             int N)
{
    
    const int M = pinit.size();
    data = new double[ N*M  ];
    new (&xi) MatMap( data, N, M);

    dt = total_time / xi.rows()+1;
    
    initialize( pinit, pgoal );
}

void Trajectory::initialize( const std::vector<double> & pinit,
                             const std::vector<double> & pgoal )
{

    assert( xi.rows() != 0 );
    
    q0 = ConstMatMap( pinit.data(), 1, pinit.size() );
    q1 = ConstMatMap( pgoal.data(), 1, pgoal.size() );

    
    createInitialTrajectory();
}



void Trajectory::initialize( const double * q0,
                             const double * q1,
                             int M, int N)
{
    
    data = new double[ N*M ];
    new (&xi) MatMap( data, N, M);
    dt = total_time / xi.rows()+1;

    initialize( q0, q1, M );
}

void Trajectory::initialize( const double * pinit, 
                             const double * pgoal,
                             int M)
{

    assert( xi.rows() > 0 );

    q0 = ConstMatMap( pinit, 1, M );
    q1 = ConstMatMap( pgoal, 1, M );
    
    createInitialTrajectory();
}


void Trajectory::createInitialTrajectory()
{
    assert( q0.size() == q1.size() );
    assert( xi.cols() == q0.size() );
    
    const int N = xi.rows();

    //simply linearly interpolate from the start to goal state.
    for ( int i=0; i < N; i ++ ){

        double t = double(i+1) / double(N+1);
        xi.row( i ) = q0 + (q1-q0)*t;
    }
}

std::string Trajectory::toString() const {
    
    std::ostringstream s;
    Eigen::IOFormat print_format(Eigen::StreamPrecision,
                                 Eigen::DontAlignCols,
                                 ", ", ", ", "", "", " = ", "");
    
    s << "[ q0" << q0.format( print_format ) << "\n";

    for ( int i = 0; i < N() ; i ++ ){
        for ( int j = 0; j < M() ; j ++ ){
            s << xi( i, j ) << " ";
        }
        s << "\n";
    }
    
    s << "q1" << q1.format( print_format ) <<  "]";

    return s.str();

}
void Trajectory::print() const{
    std::cout << toString();

}


}//namespace
