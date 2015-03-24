

template <class Derived>
inline void ChompGradient::subsampleGradient(int N_sub, 
        const Eigen::MatrixBase<Derived> & g_sub_const)
{   
    debug_status( TAG, "getSubsampledGradient", "start" );
    
    //cast away the const to edit the g_sub matrix
    Eigen::MatrixBase<Derived>& g_sub = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_sub_const);

    
    for ( int i = 0; i < g_sub.rows(); i ++ ){
        g_sub.row( i ) = g_full.row( i * 2 );
    }

    debug_status( TAG, "getSubsampledGradient", "end" );
    
}

template <class Derived>
inline void ChompGradient::evaluate(
                           const Trajectory & trajectory,
                           const Eigen::MatrixBase<Derived> & g)
{
    debug_status( TAG, "evaluate", "start" );

    if ( trajectory.isSubsampled() ){
        g_full.resize( trajectory.fullN(), trajectory.M() );

        //compute the smoothness gradient on a full-size trajectory.
        evaluateSmoothness( trajectory, g_full );

        //There is no reason to compute collision gradients
        //  on a subsampled trajectory, because even if there are
        //  obstacles, a subsampled matrix doesn't have the
        //  freedom of motion to avoid them.

        //subsample the full-sized gradient, and store the result in g
        subsampleGradient(trajectory.N(), g );
        
    } else {

        //TODO : the evaluateCollision thing needs to be working
        evaluateSmoothness( trajectory, g );
        //evaluateCollision( trajectory, g );
    }
        
    debug_status( TAG, "evaluate", "end" );
}

template <class Derived>
inline void ChompGradient::evaluateSmoothness(
                    const Trajectory & trajectory,
                    const Eigen::MatrixBase<Derived> & g_const)
{
    debug_status( TAG, "computeSmoothnessGradient", "start" );
    
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& grad = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);

    //Performs the operation: A * x.
    //  (fill the matrix Ax, with the results.
    //
    const MatMap & xi = trajectory.getFullXi();
    
    if( use_goalset ){
        diagMul(coeffs, coeffs_goalset, xi, Ax);
    } else { 
        diagMul(coeffs, xi, Ax);
    }
    
    //add in the b matrix to get the contribution from the
    //  endpoints, and set this equal to the gradient.
    grad = Ax + b; 

    debug_status( TAG, "computeSmoothnessGradient", "end" );
}

template <class Derived>
inline void ChompGradient::evaluateCollision( const Trajectory & trajectory,
                        const Eigen::MatrixBase<Derived> & g )
{

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    g.setZero();
    
    if (ghelper) {
        fextra = ghelper->addToGradient(trajectory, g);
    } else {
        fextra = 0;
    }
}

