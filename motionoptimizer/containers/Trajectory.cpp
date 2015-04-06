#include "Trajectory.h"
#include "ConstraintFactory.h"
#include "Constraint.h"

#include <sstream>

namespace chomp {

const char* Trajectory::TAG = "Trajectory";


Trajectory::Trajectory( const MatX & q0, const MatX & q1, int N,
                        ChompObjectiveType o_type, 
                        double t_total ):
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),
    full_xi(NULL, 0, 0),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    debug_status( TAG, "construction", "start");
    initialize( q0, q1, N );
    debug_status( TAG, "construction", "end");
}

Trajectory::Trajectory( const std::vector<double> & pinit,
                        const std::vector<double> & pgoal,
                        int N, ChompObjectiveType o_type, 
                        double t_total ):
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
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
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
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
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
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
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( traj );
}


Trajectory::Trajectory(const std::vector<std::vector<double> > & traj,
                       ChompObjectiveType o_type, 
                       double t_total) :
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
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
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
    objective_type( o_type ),
    total_time( t_total ),
    is_subsampled(false)
{
    initialize( traj, M, N );
}

Trajectory::Trajectory() : 
    data(NULL),
    cached_data(NULL),
    xi(NULL, 0, 0, DynamicStride(1,1) ),    
    full_xi(NULL, 0, 0),    
    objective_type( MINIMIZE_ACCELERATION ),
    total_time( 1.0 ),
    is_subsampled(false)
{
}


Trajectory::~Trajectory() { 
    if (cached_data ){ delete cached_data; }
    else {delete data; };
}


//assignment operator
Trajectory & Trajectory::operator= (const Trajectory & other)
{
    if (this != &other) // protect against invalid self-assignment
    {
        
        //delete the old data.
        if (this->data){ 
            delete data;
            this->data = NULL;
        }

        //if there is data in the other.data array, copy it over
        //  and create the xi matrix
        if (other.data){ 
            this->data = new double[other.size()];
            std::copy( other.data, 
                       other.data + other.size(),
                       this->data );

            int inner = other.getXi().innerStride();
            int outer = other.getXi().outerStride();

            new ( &(this->full_xi) ) MatMap( this->data,
                                             other.fullN(),
                                             other.fullM() );
            new ( &(this->xi) ) DynamicMatMap( 
                                this->data,
                                other.N(), other.M(),
                                DynamicStride( inner, outer ) );
        }


        this->q0 = other.q0;
        this->q1 = other.q1;

        this->dt = other.dt;
        this->total_time = other.total_time;
        this->is_subsampled = other.is_subsampled;
        
        this->objective_type = other.objective_type;

        this->cached_data = NULL;

    }
    
    //return the trajectory we just created.
    return *this;
}


// increase the size of the trajectory by one waypoint, and
//  add in q1 as the last point of the trajectory.
void Trajectory::startGoalset(){
    
    const int N = xi.rows();
    const int M = xi.cols();

    //copy the old data vector to a new data vector.
    double * new_data = new double [ size() + M ];

    //copy the old data over to the new data vector
    MatMap( new_data, N+1, M ).block(0,0, N, M) = full_xi;

    //delete the old data, and replace it with new_data.
    delete data;
    data = new_data;

    //set the xi matrix map.
    new (&xi) DynamicMatMap( data, N+1, M , DynamicStride(N,1));
    new (&full_xi) MatMap( data, N+1, M );
    
    xi.row( N ) = q1;
}


void Trajectory::endGoalset(){
    
    const int N = xi.rows();
    const int M = xi.cols();

    q1 = xi.row( N-1 ); 
        
    //copy the old data vector to a new data vector.
    double * new_data = new double [ size() - M ];
    new (&full_xi) MatMap( new_data, N-1, M );

    full_xi = xi.block(0,0, N-1, M);
    
    delete data;
    data = new_data;
    
    //set the xi matrix map.
    new (&xi) DynamicMatMap( data, N-1, M, DynamicStride(N-1,1));
}

void Trajectory::restoreData(){
    data = cached_data;
    cached_data = NULL;
    remapXi( N(), fullN(), M() );
}


//This is dangerous, and it may be a bad idea.
void Trajectory::setData( const double * new_data ){
    
    //if we are subsampling, just copy the data over
    if (is_subsampled){
        xi = ConstMatMap( new_data, N(), M() );

    //if we are not subsampling, this are more complicated,
    //  because we want to avoid copying all of the data,
    //  in this case, we cache the data, if the data has not 
    //  already been cached, and we set the data pointer
    //  to the new data, making sure to change the mat mappings 
    //  as well
    } else {
        
        if (cached_data == NULL ){ cached_data = data;}
        data = const_cast<double*>(new_data);

        remapXi( N(), N(), M() );
    }
}

void Trajectory::setData( double * new_data ){
    
    //if we are subsampling, just copy the data over
    if (is_subsampled){
        xi = MatMap( new_data, N(), M() );
        
    //if we are not subsampling, this are more complicated,
    //  because we want to avoid copying all of the data,
    //  in this case, we cache the data, if the data has not 
    //  already been cached, and we set the data pointer
    //  to the new data, making sure to change the mat mappings 
    //  as well
    } else {
        
        if (cached_data == NULL ){ cached_data = data;}
        
        remapXi( N(), N(), M() );
    }
}

void Trajectory::copyToData( const std::vector<double> & vec )
{
    if ( cached_data ){ restoreData(); }
    xi = ConstMatMap( vec.data(), N(), M() );
}
void Trajectory::copyToData( const double * new_data )
{
    if ( cached_data ){ restoreData(); }
    xi = ConstMatMap( new_data, N(), M() );
}

void Trajectory::copyDataTo( double * copy) const 
{
    MatMap( copy, N(), M() ) = xi;
}

void Trajectory::copyDataTo( std::vector<double> & vec ) const
{
    vec.resize( size() );
    MatMap( vec.data(), N(), M() ) = xi;
}


void Trajectory::subsample(){
    
    debug_status( TAG, "subsample", "start" );
    
    assert( !cached_data );
    assert( !is_subsampled );
    
    is_subsampled = true;

    const int n = full_xi.rows();
    const int n_sub = n/2 + 1;
    const int m = M();

    //re-map the xi map
    new (&xi) DynamicMatMap( data, n_sub, m, DynamicStride( n, 2 ));

    debug_status( TAG, "subsample", "end" );
}

void Trajectory::endSubsample(){
    assert( is_subsampled );
    is_subsampled = false;

    const int n = full_xi.rows();
    const int m = M();

    new (&xi) DynamicMatMap( data, n, m, DynamicStride(n,1) );

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
    
    //TODO this should be an error. After a set of calls to
    //  setData, the data should have been restored.
    //  
    assert( !cached_data );
    
    debug_status( TAG, "upsample", "start" );
    const int N = xi.rows();
    const int M = xi.cols();
    const int N_up = 2*N+1; // e.g. size 3 goes to size 7
    
    double * upsampled_data = new double [N_up*M];
    new (&full_xi) MatMap(upsampled_data, N_up, M);

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
            assert(t < full_xi.rows());

            if (objective_type == MINIMIZE_VELOCITY) {

                MatX qneg1 = getTick(t/2-1);
                MatX qpos1 = getTick(t/2);
                full_xi.row(t) = 0.5 * (qneg1 + qpos1);

            } else { 

                MatX qneg3 = getTick(t/2-2);
                MatX qneg1 = getTick(t/2-1);
                MatX qpos1 = getTick(t/2  );
                MatX qpos3 = getTick(t/2+1);

                const double c3 = -1.0/160;
                const double c1 = 81.0/160;

                full_xi.row(t) = c3*qneg3 + c1*qneg1 + c1*qpos1 + c3*qpos3;

            }

        } else {
            full_xi.row(t) = xi.row(t/2);
        }
    }
    

    //clean up the old data array.
    delete data;
    data = upsampled_data;
    
    //create the new matrix map, and set the data variable to the
    //  new data array. 
    new (&xi) DynamicMatMap( data, N_up, M, DynamicStride( N_up, 1 ));
    
    dt = total_time / xi.rows()+1;

    debug_status( TAG, "upsample", "end" );
}

void Trajectory::upsampleTo( int Nmax ){
    while ( xi.rows() < Nmax ){ upsample(); }
}

void Trajectory::constrainedUpsample(ConstraintFactory & factory,
                                     double htol,
                                     double hstep)
{
    upsample();
    
    //if there is no factory, or there are no constraints,
    //    do not evaluate the constraints.
    if ( !factory.empty()  ){ factory.getAll( xi.rows() ); }
    else{ return; }


    MatX h, H, delta;
    double hinit = 0, hfinal = 0;

    for (int i=0; i < xi.rows(); i+=2) {
  
        Constraint* c = factory.getConstraint(i);

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

void Trajectory::constrainedUpsampleTo(ConstraintFactory & factory,
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
    initializeData( traj );
}

void Trajectory::initialize(const std::vector<double> & traj, int N)
{
    initializeData( ConstMatMapR( traj.data(), N, traj.size()/N ) );
}

void Trajectory::initialize(const std::vector< std::vector<double> > & traj)
{
    const int N = traj.size()-2;
    const int M = traj[0].size();

    q0 = ConstMatMap( traj[0].data(), 1, M );
    q1 = ConstMatMap( traj.back().data(), 1, M );
    
    data = new double[ N*M ];
    remapXi( N, N, M );

    for( int i = 0; i < N; i ++ ){
        assert( traj[i+1].size() == size_t(M) );
        xi.row( i ) = ConstMatMap( traj[i+1].data(), 1 , M );
    }

    dt = total_time / xi.rows()+1;
    
}

void Trajectory::initialize(const double * traj, int N, int M)
{
    initializeData( ConstMatMapR( traj, N, M ) );
}

void Trajectory::initialize(const MatX & pinit, const MatX & pgoal, int N)
{
    initializeData( pinit, pgoal, N );
}


void Trajectory::initialize( const std::vector<double> & pinit,
                             const std::vector<double> & pgoal,
                             int N)
{
    assert( pinit.size() == pgoal.size() );

    const int M = pinit.size();
    
    initializeData( ConstMatMap( pinit.data(), 1, M ),
                    ConstMatMap( pgoal.data(), 1, M ),
                    N );

}


void Trajectory::initialize( const double * pinit,
                             const double * pgoal,
                             int N, int M)
{
    initializeData( ConstMatMap( pinit, 1, M ),
                    ConstMatMap( pgoal, 1, M ), N );
}


void Trajectory::createInitialTrajectory()
{

    debug_status( TAG, "createInitialTrajectory", "start");
    assert( q0.size() == q1.size() );
    assert( xi.cols() == q0.size() );
    
    const int N = xi.rows();

    //simply linearly interpolate from the start to goal state.
    for ( int i=0; i < N; i ++ ){

        double t = double(i+1) / double(N+1);
        full_xi.row( i ) = q0 + (q1-q0)*t;
    }
    

    debug_status( TAG, "createInitialTrajectory", "start");
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

void Trajectory::getNonCovariantTrajectory( const Metric & metric,
                                            Trajectory & other ) const
{
    debug_status( TAG, "getNonCovariantTrajectory", "start" );
    
    copyDataToOther( other );
    metric.multiplyLowerInverseTranspose( other.full_xi );

    debug_status( TAG, "getNonCovariantTrajectory", "end" );
}

void Trajectory::getCovariantTrajectory( const Metric & metric,
                                         Trajectory & other ) const
{
    debug_status( TAG, "getCovariantTrajectory", "start" );

    copyDataToOther( other );
    metric.multiplyLowerTranspose( other.full_xi );

    debug_status( TAG, "getCovariantTrajectory", "end" );
}

void Trajectory::copyDataToOther( Trajectory & other ) const{
    
    debug_status( TAG, "copyDataToOther", "start" );
    //set up the shape of the matrix
    //
    if ( other.fullN() != this->fullN() ){

        if (other.data){  delete other.data; }

        const int fullN = this->fullN();
        const int N = this->N();
        const int M = this->M();

        other.data = new double [ fullN * M ];
        other.remapXi( N, fullN, M );

        other.dt = this->dt;
        other.total_time = this->total_time;

    }else if ( other.N() != this->N() ){
        other.remapXi( this->N(), this->fullN(), this->M());
    }

    other.is_subsampled = this->is_subsampled;
    other.objective_type = this->objective_type;

    debug_status( TAG, "copyDataToOther", "beforeAssign" );
    
    //copy over the data.
    other.full_xi = this->full_xi;

    debug_status( TAG, "copyDataToOther", "end" );
}


}//namespace
