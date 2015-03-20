//THis is the implementation file of the templated functions
//  in the constraint factory


//THe definitions of the templated functions
template <class Derived1, class Derived2>
void ConstraintFactory::evaluate(
        const Eigen::MatrixBase<Derived1> & h_tot_const,
        const Eigen::MatrixBase<Derived2> & H_tot_const)
{
    
    Eigen::MatrixBase<Derived1>& h_tot = 
        const_cast<Eigen::MatrixBase<Derived1>&>(h_tot_const);

    Eigen::MatrixBase<Derived2>& H_tot = 
        const_cast<Eigen::MatrixBase<Derived2>&>(H_tot_const);

    //TODO - Make the trajectory aware of subsampled matrices.
    size_t DoF = trajectory.cols();
    size_t N = trajectory.rows();
    H_tot.setZero();    

    assert(N == constraints.size());

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints

    for (size_t i=0; i<constraints.size(); ++i) {

        Constraint* c = constraints[i];

        //get individual h and H from qt
        if (!c || c->numOutputs()==0){ continue; }

        c->evaluateConstraints( trajectory.row(i), h, H);

        //stick all the h and H into h_tot and H_tot
        assert(H.cols() == (int)DoF);
        assert(h.cols() == 1);
        assert(H.rows() == (int)constraints[i]->numOutputs());
        assert(H.rows() == h.rows());
        

        for (size_t r = 0; r<(size_t)h.rows(); r++, row++){
            h_tot(row) = h(r);

            for (size_t j = 0; j<DoF; j++){
                //TODO change to debug assert
                assert( H_tot.cols() > (int)row   );
                assert( H_tot.rows() > (int)(j*N + i) );
              
                H_tot(j*N + i, row) = H(r,j);
            }
        }
    }  
}


template <class Derived>
void ConstraintFactory::evaluate( 
        const Eigen::MatrixBase<Derived> & h_tot_const)
{
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
}


