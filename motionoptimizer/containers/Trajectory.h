

#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include "../utils/utils.h"

namespace chomp{

class Trajectory {
  public:

    double * data, * sampled_data, * cached_data;
    MatMap xi, sampled_xi;

    MatX q0, q1;

    //TODO : remove this from the trajectory. This should be in
    //  the gradient.
    ChompObjectiveType objective_type;
    
    static const char* TAG;
  private:
    double dt, total_time;
    bool is_subsampled;

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
    Trajectory & operator= (const Trajectory & other);
    
    void startGoalSet();
    void endGoalSet();

    //restores the data to the original after 
    //  it had been replaced in a call to setData.
    void restoreData();
    
    void setData( const double * new_data );
    void setData( double * new_data );

    double * getData(){ return data; }
    const double * getData() const { return data; } 
    
    void copyToData( const std::vector<double> & vec );
    void copyToData( const double * new_data );
    
    void copyDataTo( double * copy ) const ;
    void copyDataTo( std::vector<double> & vec ) const;

    MatMap const & getTraj()const { return xi; }
    MatMap       & getTraj()      { return xi; }
    void getTraj( std::vector<double> & vec ) const ;
    void getTraj( double * array ) const ;

    void getState( std::vector<double> & vec, int i ) const ;
    void getState( double * array, int i ) const ;
    
    inline double getDt() const { return dt; }

    void subsample();
    void endSubsample();

    bool isSubsampled() const { return is_subsampled; }
    
    inline double & operator ()( int i, int j ){
        return xi( i, j );
    }
    inline double const & operator()(int i, int j)const {
        return xi( i, j );
    }

    inline int N() const { return xi.rows(); }
    inline int M() const { return xi.cols(); }

    inline int sampledN() const { return sampled_xi.rows(); }

    inline int rows() const { return xi.rows();}
    inline int cols() const { return xi.cols();}

    inline int size() const { return xi.size();}

    inline const MatX & getQ0() const { return q0; }
    inline const MatX & getStart() const { return q0; }

    inline const MatX & getQ1() const { return q1; };
    inline const MatX & getEnd() const { return q1; }

    //get individual matrix blocks or rows
    inline Row row( int i ) { return xi.row(i); }
    inline ConstRow row( int i ) const{ return xi.row(i); }

    inline Col col( int i ){  return xi.col(i); }
    inline ConstCol col( int i ) const { return xi.col(i); }

    inline MatMap const & getXi() const { return xi; } 
    inline MatMap & getXi() { return xi; } 

    inline MatMap const & getSampledXi() const { return sampled_xi; } 
    inline MatMap & getSampledXi() { return sampled_xi; } 

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
    void initialize(const MatX & traj);
    void initialize(const std::vector<double> & traj, int M);
    void initialize(const std::vector<std::vector<double> > & traj);
    void initialize(const double * traj, int M, int N);
    
    //initialize from two endpoints.
    void initialize(const MatX & q0, const MatX & q1, int N);
    void initialize( const MatX & q0, const MatX & q1 );

    void initialize( const std::vector<double> & q0,
                   const std::vector<double> & q1 );
    void initialize( const std::vector<double> & q0,
                   const std::vector<double> & q1,
                   int N);
    void initialize( const double * q0, const double * q1, int N);
    void initialize( const double * q0, const double * q1,
                             int M, int N);

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

    inline MatX getTick(int tick) const;

    std::string toString() const;
    void print() const;

  private:
    void createInitialTrajectory(); 
       
};

//definitions of the inline functions.



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
    if (is_subsampled)
    { 
        SubMatMap( sampled_data, N(), M(),
                  SubMatMapStride( sampledN(), 2 ) ) = xi;
    }
}

template <class Derived>
inline void Trajectory::update(const Eigen::MatrixBase<Derived> & delta,
                               int index)
{
    xi.row( index ) -= delta;
    if (is_subsampled){
        SubMatMap( sampled_data, N(), M(),
                   SubMatMapStride( sampledN(),2 ) ).row(index)
            = xi.row(index);
    }
}

}//namespace

#endif 
