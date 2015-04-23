

template <class Derived>
inline void SmoothnessFunction::subsampleGradient(int N_sub, 
        const Eigen::MatrixBase<Derived> & g_sub_const)
{   
    debug_status( TAG, "getSubsampledSmoothnessFunction", "start" );
    
    //cast away the const to edit the g_sub matrix
    Eigen::MatrixBase<Derived>& g_sub = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_sub_const);

    //Subsample the g_full matrix. This uses the matrix map method
    //  to subsample the matrix, basically, we make a Matrix Map that
    //  make g_full look like the subsampled matrix that we want,
    //  and then we set g_sub to be equal to that image.
    g_sub = DynamicMatMap( g_full.data(),
                           g_sub.rows(), g_sub.cols(), 
                           DynamicStride(g_full.rows(), 2) );

    debug_status( TAG, "getSubsampledSmoothnessFunction", "end" );
    
}


template <class Derived>
inline double SmoothnessFunction::evaluate(
                    const Trajectory & trajectory,
                    const Metric & metric,
                    const Eigen::MatrixBase<Derived> & g_const)
{
    debug_status( TAG, "evaluateGradient", "start" );
    
    //cast away the const-ness of g_const
    Eigen::MatrixBase<Derived>& g = 
        const_cast<Eigen::MatrixBase<Derived>&>(g_const);
    
            
    //Performs the operation: A * x.
    //  (fill the matrix Ax, with the results.
    metric.multiply( trajectory.getFullXi(), Ax);
    
    //add in the b matrix to get the contribution from the
    //  endpoints, and set this equal to the gradient.
    if (trajectory.isSubsampled() ){
        g_full = Ax + b; 
        subsampleGradient(trajectory.N(), g );
    }else {
        g = Ax + b; 
    }
    
    debug_status( TAG, "evaluateGradient", "end" );
    
    //return the computed smoothness value of the trajectory
    return computeValue( trajectory );
    
}


