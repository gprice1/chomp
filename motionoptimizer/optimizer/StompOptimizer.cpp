
#include "StompOptimizer.h"

namespace mopt{

StompOptimizer::StompOptimizer(ProblemDescription & problem,
                               Observer * observer = NULL,
                               double obstol = 1e-8,
                               double timeout_seconds,
                               size_t max_iter) :
    OptimizerBase( problem, observer, obstol, timeout_seconds, max_iter ),
    K(10),i
    update( problem.N(), problem.M() ),
    maximums( problem.N(), 1),
    minimums( problem.N(), 1),
    q0(1, problem.M()), 
    q1(1, problem.M()), 
    q2(1, problem.M()),
    cspace_vel(1, problem.M()),
    cspace_accel(1, problem.M())
{
}

void StompOptimizer::solve()
{
    
    generateNoisyTrajectories();

    calculatePerTimestepCost();

    updateTrajectory();
}

void StompOptimizer::updateTrajectory()
{

    //TODO this complex series of matrix operations could use some 
    //  speedening up.
    maximums = timestep_cost.rowwise().maxCoeff();
    minimums = timestep_cost.rowwise().minCoeff();

    timestep_cost = Eigen::pow( M_E, 
                     (timestep_cost.colwise() - minimums).array()/
                     (maximums - minimums).array() );

    timestep_cost = timestep_cost.array().colwise() / 
                    timestep_cost.rowwise().sum();
    
    noisy_samples.leftCols( problem.M() ).colwise().array()
        *= timestep_cost.col(0);

    for( int i = 0; i < K; ++i ){
        noisy_samples.leftCols( problem.M() ) += 
            noisy_samples.block( 0, K*i, problem.N(), problem.M()
            ).colwise().array() * timestep_cost.col(i);
    }

    metric.solve( noisy_samples.leftCols( problem.M() ) );
    problem.updateTrajectory( - noisy_samples.leftCols( problem.M() ) );
    
}


void StompOptimizer::generateNoisyTrajectories()
{
    noisy_samples.resize( problem.N(), problem.M() * K );
    problem.getMetric().sampleNormalDistribution( standard_deviation, 
                                                  noisy_samples );
    noisy_samples *= metric_scale;
    
}

void StompOptimizer::calculatePerTimestepCost(int i)
{
    
    q1 = trajectory.getQ0;
    q2 = noisy_samples.block( 0, K*i, 1, problem.M() ) +
       problem.getTrajectory.row(0);

    for (int j = 0; j < problem.N(); ++j ){
        q0 = q1;
        q1 = q2;

        if (i == problem.N()-1 ){
            q2 = problem.getTrajectory().getQ1();
        }else {
            q2 = noisy_samples.block( j+1, K*i, 1, problem.M() ) +
                 problem.getTrajectory.row(j+1);
        }

        cspace_vel = q2 - q1;
        cspace_accel = q1*2 - q2 - q0;
        
        double collision_cost = calculateCollisionCosts( q1 );
        double constraint_cost = calculateConstraintCost( q1, j );
        timestep_cost(j, i) = collision_cost + 
                              torque_cost +
                              constraint_cost;
    }
}

double calculateTrajectoryCost(){
    
}

void StompOptimizer::calculateCollisionCost( const MatX & state )
{
    

}

void StompOptimizer::calculateConstraintCost( const MatX & state, 
                                              int timestep)
{


}



}//namespace
