
#include "Chomp.h"


namespace chomp {

void Chomp::solve( bool global = true, bool local = false ){

    //optimize at the current resolution
    optimize(global, local);

    //if the current resolution is not the final resolution,
    //  upsample and then optimize. Repeat until the final resolution
    //  is reached or exceeded.
    while ( xi.rows() < N_max ){

        //upsample the trajectory and prepare for the next
        //  stage of optimization.
        MatX xi_up;
        upsampleTrajectory( xi, gradient->q0,
                                gradient->q1,
                                gradient->dt,
                                gradient->objective_type, xi_up );
        xi = xi_up;

        optimize(global, local);
    }
}


void ChompOptimizer::optimize( bool global, bool local ){
    
    if ( use_goalset ){ prepareGoalSet(); }
    gradient->prepareRun( N, use_goalset, subsample );

    // Subsample
    bool subsample = N > minN && !use_goalset &&
                     !(full_global_at_final && N >= maxN);
    if( subsample ){
        N_sub = (N+1)/2;
        new (&xi_sub) SubMatMap( xi.data(), N_sub, M,
                                 SubMatMapStride(N,2) );
    }


}
void ChompOptimizer::useGoalSet( Constraint * goalset ){
      this->goalset = goalset;
      use_goalset = true;



}

void ChompOptimizer::prepareGoalSet(){
    
    //do not subsample if doing goalset run.
    N_sub = 0;

    //resize xi, and add q1 into it.
    xi.conservativeResize( xi.rows() + 1, xi.cols() );
    xi.row( xi.rows() - 1 ) = gradient->q1;
    
    //set N to the current size of xi.
    N = xi.rows();

    //add the goal constraint to the constraints vector.
    factory->constraints.push_back( goalset );
}

void ChompOptimizer::finishGoalSet(){
    
    use_goalset = false;
    
    //copy the last state in the trajectory back into q1
    gradient->q1 = xi.row( xi.rows() - 1 );

    //resize xi, keeping old values, and get rid of the 
    //  final state. And set N to the correct trajectory size
    xi.conservativeResize( xi.rows() -1, xi.cols() );

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    factory->constraints.pop_back();

}

void Chomp::constrainedUpsampleTo(int Nmax, double htol, double hstep)
{

    MatX h, H, delta;
  
    while (N < Nmax) { 

      upsample();
      prepareChomp();
      prepareChompIter();

      double hinit = 0, hfinal = 0;
      
      //if there is no factory, or there are no constraints,
      //    do not evaluate the constraints.
      if ( !factory || factory->constraints.empty() ){ continue; }

      for (int i=0; i<N; i+=2) {
    
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
}





}//namespace
