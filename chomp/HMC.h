#ifndef _HMC_H_
#define _HMC_H_

#include "Chomp.h"

namespace chomp{

enum ChompObjectiveType;

class HMC{

  public:
    
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
    
    // old_xi : the previous xi, saved from the last resample iteration,
    //          if rejection is on, this will be restored if the, 
    //          energy of the trajectory does not increase.
    // old_momentum_delta : like old_xi, this is a saved previous state,
    //                      for use if the current trajectory is rejected.
    MatX old_xi, old_momentum;
   
    
    //the smoothing operator to be applied to the random momentum vector,
    //  and the size of the smoothing operator.
    double smoothing_operator[3];
    int n_operator;

  //PUBLIC MEMBER FUNCIONS


    HMC( double lambda=0.02, bool doNotReject=true ) :
        lambda( lambda ), doNotReject( doNotReject ){}

    //setup the random seed for HMC
    void setSeed(unsigned long seed=0);
    
    //resamples the momentum.
    void iteration( size_t cur_iteration, MatX & xi, MatX & momentum,
                    const MatX L, double lastObjective);
    
    //checks the current HMC iteration, and rejects it if the 
    //  energy of the system is too low.
    bool checkForRejection( MatX & xi, MatX & momentum,
                            double lastObjective);
    
    //samples a random momentum from the probability dist given by:
    //  exp( -0.5*xAx )
    void getRandomMomentum( MatX & momentum, size_t cur_iteration );
    
    //called in prepareChomp. Sets up a run of HMC.
    void setupRun();

    void setupHMC( ChompObjectiveType objective_type, double chomp_alpha );

};

}//namespace

#endif
