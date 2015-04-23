//THis is the implementation file of the templated functions
//  in the constraint factory


inline void ConstraintFactory::addGoalset( Constraint * goalset )
{
    constraints.push_back( goalset );
}
inline void ConstraintFactory::removeGoalset()
{
    constraints.resize( constraints.size() - 1 );
}

inline bool ConstraintFactory::empty() const 
{ 
    return constraint_intervals.empty();
}

inline Constraint* ConstraintFactory::getConstraint( size_t t ) const
{
    return constraints[t];
}

inline const std::vector<Constraint*> & ConstraintFactory::getConstraints()
    const
{
    return constraints;
}

inline int ConstraintFactory::numOutput(){ return constraint_dims; }

//THe definitions of the templated functions
template <class Derived1, class Derived2>
double ConstraintFactory::evaluate(
        const Trajectory & trajectory,
        const Eigen::MatrixBase<Derived1> & h_tot_const,
        const Eigen::MatrixBase<Derived2> & H_tot_const)
{
    
    debug_status( TAG, "evaluate", "start" );

    if (empty()){ return 0; }
    
    Eigen::MatrixBase<Derived1>& h_tot = 
        const_cast<Eigen::MatrixBase<Derived1>&>(h_tot_const);

    Eigen::MatrixBase<Derived2>& H_tot = 
        const_cast<Eigen::MatrixBase<Derived2>&>(H_tot_const);
    
    int M = trajectory.cols();
    int N = trajectory.rows();
    H_tot.setZero();    
    
    debug_assert(size_t( N ) == constraints.size());

    MatX h, H; // individual constraints at time t
    
    debug_status( TAG, "evaluate", "end" );
    

    for (int i=0, row=0; size_t(i) < constraints.size(); ++i) {

        Constraint* c = constraints[i];

        //get individual h and H from qt
        if (!c || c->numOutputs()==0){ continue; }

        c->evaluateConstraints( trajectory.row(i), h, H);

        //stick all the h and H into h_tot and H_tot
        debug_assert(H.cols() == M);
        debug_assert(h.cols() == 1);
        debug_assert(size_t( H.rows()) == constraints[i]->numOutputs());
        debug_assert(H.rows() == h.rows());
        
        for (int r = 0; r < h.rows(); r++, row++){
            debug_assert( h_tot.rows() > row );
            debug_assert( h_tot.cols() == h.size() );
            h_tot(row) = h(r);

            for (int j = 0; j<M; j++){
                debug_assert( H_tot.cols() > row   );
                debug_assert( H_tot.rows() > j*N + i );
              
                H_tot(j*N + i, row) = H(r,j);
            }
        }
    }  

    debug_status( TAG, "evaluate", "end" );
    
    //TODO return the magnitude of the constraint violations
    return h.lpNorm<Eigen::Infinity>();
}


template <class Derived>
double ConstraintFactory::evaluate( 
        const Trajectory & trajectory,
        const Eigen::MatrixBase<Derived> & h_tot_const)
{
    if (empty()){ return 0; }
    
    Eigen::MatrixBase<Derived>& h_tot = 
        const_cast<Eigen::MatrixBase<Derived>&>(h_tot_const);

    assert(size_t(trajectory.rows()) == constraints.size());

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints

    for (size_t i=0; i<constraints.size(); ++i) {

        Constraint* c = constraints[i];

            //get individual h and H from qt
        if (!c || c->numOutputs()==0){ continue; }

        c->evaluateConstraints( trajectory.row(i), h, H);
      
        for (size_t r = 0; r<(size_t)h.rows(); r++, row++){
            h_tot(row) = h(r);
        }
    }  
    
    //TODO return the magnitude of the constraint violations
    return h.lpNorm<Eigen::Infinity>();

}


