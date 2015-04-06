

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
                           const Eigen::MatrixBase<Derived> & g,
                           const Trajectory * covariant_trajectory )
{
    debug_status( TAG, "evaluate", "start" );
    
    const bool is_covariant = ( covariant_trajectory != NULL );
    
    if ( trajectory.isSubsampled() ){
        g_full.resize( trajectory.fullN(), trajectory.M() );

        //compute the smoothness gradient on a full-size trajectory.
        if ( is_covariant ){
            evaluateSmoothness( *covariant_trajectory, g_full, true);
        }else {
            evaluateSmoothness( trajectory, g_full, false );
        }

        //There is no reason to compute collision gradients
        //  on a subsampled trajectory, because even if there are
        //  obstacles, a subsampled matrix doesn't have the
        //  freedom of motion to avoid them.

        //subsample the full-sized gradient, and store the result in g
        subsampleGradient(trajectory.N(), g );
        
    } else {
        if ( is_covariant ){
            evaluateSmoothness( *covariant_trajectory, g, true);
        }else {
            evaluateSmoothness( trajectory, g, false );
        }
        
        evaluateCollision( trajectory, g, is_covariant );
    }
        
    debug_status( TAG, "evaluate", "end" );
}

template <class Derived>
inline void ChompGradient::evaluateSmoothness(
                    const Trajectory & trajectory,
                    const Eigen::MatrixBase<Derived> & g_const,
                    bool is_covariant)
{
    debug_status( TAG, "computeSmoothnessGradient", "start" );
    
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& g = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);
    
    if (is_covariant){
        g = trajectory.getFullXi() + b;
    }else {
        
        //Performs the operation: A * x.
        //  (fill the matrix Ax, with the results.
        metric.multiply( trajectory.getFullXi(), Ax);
        
        //add in the b matrix to get the contribution from the
        //  endpoints, and set this equal to the gradient.
        g = Ax + b; 
    }

    debug_status( TAG, "computeSmoothnessGradient", "end" );
}

template <class Derived>
inline void ChompGradient::evaluateCollision(
                        const Trajectory & trajectory,
                        const Eigen::MatrixBase<Derived> & g_const,
                        bool is_covariant )
{

    //If there is a gradient helper, add in the contribution from
    //  that source, and set the fextra variable to the cost
    //  associated with the additional gradient.
    if (ghelper) {
        //TODO : this method seems hacky, and is unfortunately slow.
        MatX g_full = MatX::Zero( g_const.rows(), g_const.cols() );
        fextra = ghelper->addToGradient(trajectory, g_full);

            //cast away the const-ness of g_const
        Eigen::MatrixBase<Derived>& g = 
            const_cast<Eigen::MatrixBase<Derived>&>(g_const);
        
        //if the problem is covariant, get the covariant gradient step
        if ( is_covariant ) { metric.multiplyLowerInverse( g_full ); }
        g += g_full;
        
    } else {
        fextra = 0;
    }
}

