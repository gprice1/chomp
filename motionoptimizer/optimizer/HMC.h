#ifndef _HMC_H_
#define _HMC_H_

#include "../utils/utils.h"
#include "../containers/Metric.h"

namespace chomp{

class HMC{

  private:
    
    // hmc_lambda : the parameter that determines the frequency, and
    //              magnitude of random resampling
    double lambda, alpha, previous_energy;

    // doNotReject : if true, the energy of the system will not
    //               be evaluated, and low energy systems will not be
    //               rejected.
    bool doNotReject;

    // hmc_resample_iter : the next iteration that resampling will be
    //                     done on
    size_t resample_iter;
    
    //old_data : the previous xi, saved from the last resample iteration,
    //          if rejection is on, this will be restored if the, 
    //          energy of the trajectory does not increase.
    //old_momentum : like old_xi, this is a saved previous state,
    //               for use if the current trajectory is rejected.
    MatX old_momentum;
    double * old_data;
    
    static const std::string TAG;

  //PUBLIC MEMBER FUNCIONS
  public:


    HMC( double lambda=0.02, bool doNotReject=true );

    ~HMC();    
    

    //SHould be called each iteration of the optimizer,
    //  if it is on a resample iteration,
    //  it will resample the momentum.
    void iterate( size_t current_iteration,
                  double lastObjective
                  const Metric & metric,
                  Trajectory & traj,
                  MatX & momentum);
    
    //setup the random seed for HMC
    void setSeed(unsigned long seed=0);

  private: 
    
    
    //checks the current HMC iteration, and rejects it if the 
    //  energy of the system is too low.
    bool checkForRejection( Trajectory & traj,
                            MatX & momentum,
                            double lastObjective);
    
    //samples a random momentum from the probability dist given by:
    //  exp( -0.5*xAx )
    void getRandomMomentum( const Metric & metric,
                            size_t current_iteration,
                            MatX & momentum );
        
};

}//namespace

#endif
