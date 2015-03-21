void ConstraintFactory::evaluate( MatX& h_tot, 
                                MatX& H_tot, 
                                int step)
{

    size_t DoF = xi.cols();

    assert(size_t(xi.rows()) == constraints.size());

    size_t numCons = 0;

    size_t timesteps = 0;

    for (size_t t=0; t<constraints.size(); t+=step) {
        Constraint* c = constraints[t];
        if (c) {
            numCons += c->numOutputs();
        }
        ++timesteps;
    }

    if (size_t(h_tot.rows()) != numCons || size_t(h_tot.cols()) != 1) {
        h_tot.resize(numCons,1);
    }

    if (   size_t(H_tot.rows()) != numCons 
        || size_t(H_tot.cols()) != DoF*timesteps)
    {
        H_tot.resize(numCons,DoF*timesteps);
    }

    H_tot.setZero(); 

    MatX h, H; // individual constraints at time t
    size_t row = 0; // starting row for the current bundle o constraints

    for (size_t t=0, i=0; t<constraints.size(); t+=step, ++i) {

        Constraint* c = constraints[t];

        //get individual h and H from qt
        if (!c || c->numOutputs()==0){ continue; }
        c->evaluateConstraints(xi.row(t), h, H);

        if (h.rows() == 0) { 
            assert(H.rows() == 0);
            continue;
        }

        //stick all the h and H into h_tot and H_tot
        assert((size_t)H.cols()==DoF);
        assert(h.cols() == 1);
        assert((size_t)H.rows()<=constraints[t]->numOutputs());
        assert(H.rows() == h.rows());

        for (size_t r = 0; r<(size_t)h.rows(); r++){
            h_tot(row) = h(r);
            for (size_t j = 0; j<DoF; j++){
                assert( row < size_t(H_tot.rows()) );
                assert( j*timesteps+i < size_t(H_tot.cols()) );
                assert( j < size_t(H.cols()) );
                H_tot(row, j*timesteps+i) = H(r,j);
            }
          row++;
        }

    }  

    assert( row <= (size_t)numCons );
    if ( row < numCons ) {
        h_tot.conservativeResize(row, 1);
        H_tot.conservativeResize(row, DoF*timesteps);
    }

    assert(h_tot.rows() == H_tot.rows());
    assert(h_tot.cols() == 1);
    assert(size_t(H_tot.cols()) == DoF*timesteps);

}

