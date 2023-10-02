 
  void createGlobalMatrices(const double timeStep, const double _alpha, const double _beta)
  {
   
    /*************************create the M, D, K matrices from alpha, beta, poisson ratio, and Young's modulus. Afterward create the matrix "A" with the given timeStep that is the left hand side of the entire system.
     *********/
    
    double lambda = (poissonRatio * youngModulus) / ((1 + poissonRatio) * (1 - 2 * poissonRatio));
    double mu = youngModulus / (2 * (1 + poissonRatio));

    SparseMatrix<double> C = SparseMatrix<double>(6,6);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            C.coeffRef(i, j) += lambda;
        }
        C.coeffRef(i, i) += 2 * mu;
        C.coeffRef(i + 3, i + 3) += mu;
    }

    SparseMatrix<double> D_local = SparseMatrix<double>(6, 9);
    D_local.insert(0,0) = 1;
    D_local.insert(1,4) = 1;
    D_local.insert(2,8) = 1;

    D_local.insert(3, 1) = 0.5;
    D_local.insert(3, 3) = 0.5;
    D_local.insert(4, 5) = 0.5;
    D_local.insert(4, 7) = 0.5;
    D_local.insert(5, 2) = 0.5;
    D_local.insert(5, 6) = 0.5;

    SparseMatrix<double> J_e = SparseMatrix<double>(9, 12);
    SparseMatrix<double> G_e = SparseMatrix<double>(3, 4);
    SparseMatrix<double> B_e;
    MatrixXd P_e = MatrixXd(4, 4);
    MatrixXd P_e_inverse = MatrixXd(4, 4);
    SparseMatrix<double> KTet = SparseMatrix<double>(12, 12);
    SparseMatrix<double> KPrime = SparseMatrix<double>(12 * tetVolumes.size(), 12 * tetVolumes.size());

    vector<Triplet<double>> tripletListK;
    tripletListK.reserve(tetVolumes.size() * 12 * 12);
    for (int i = 0; i < tetVolumes.size(); i++)
    {
        G_e = SparseMatrix<double>(3, 4);
        P_e = MatrixXd(4,4);
        J_e = SparseMatrix<double>(9, 12);
        for (int j = 0; j < 4; j++)
        {
            P_e(j, 0) = 1;
            P_e(j, 1) = origPositions(T.coeff(i, j) * 3);
            P_e(j, 2) = origPositions(T.coeff(i, j) * 3 + 1);
            P_e(j, 3) = origPositions(T.coeff(i, j) * 3 + 2);
        }
        P_e_inverse = P_e.inverse();
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                if (P_e_inverse(j + 1, k) != 0)
                    G_e.insert(j, k) = P_e_inverse(j + 1, k);
            }
        }

        for (int k = 0; k < 3; k++)
        {
            for (int l = 0; l < 4; l++)
            {
                J_e.coeffRef(k, l ) = G_e.coeffRef(k,l);
                J_e.coeffRef(k + 3, l + 4) = G_e.coeffRef(k, l);
                J_e.coeffRef(k + 6, l + 8) = G_e.coeffRef(k, l);
            }
        }
        B_e = D_local * J_e;
        KTet = tetVolumes(i) * B_e.transpose() * C * B_e;
        for (int j = 0; j < 12; j++)
        {
            for (int k = 0; k < 12; k++)
            {
                tripletListK.push_back(Triplet<double>(12 * i + j, 12 * i + k, KTet.coeffRef(j, k)));
            }
        }
    }
    // KPrime is supposed to be 12|T| X 12|T|
    KPrime.setFromTriplets(tripletListK.begin(), tripletListK.end());

    //Q coordinate dominant
    SparseMatrix<double> Q = SparseMatrix<double>(12 * tetVolumes.size(), voronoiVolumes.size() * 3);
    int temp;
    for (int i = 0; i < tetVolumes.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            temp = T.coeff(i, j); // Get which vertex is the j-th vertex of the i-th tet
            Q.insert(12 * i + j, temp * 3) = 1;
            Q.insert(12 * i + j + 4, temp * 3 + 1) = 1;
            Q.insert(12 * i + j + 8, temp * 3 + 2) = 1;
        }
    }

    K = Q.transpose() * KPrime * Q; // positions in K will be relevant to a certain vertex

    vector<Triplet<double>> tripletList;
    M = SparseMatrix<double>(invMasses.size() * 3, invMasses.size() * 3);
    tripletList.reserve(invMasses.size() * 3);
    for (int i = 0; i < invMasses.size(); i++)
    {
        tripletList.push_back(Triplet<double>(3 * i, 3 * i, 1 / invMasses[i])); // or maybe 1/4 * density * voronoiVolumes[i]
        tripletList.push_back(Triplet<double>(3 * i + 1, 3 * i + 1, 1 / invMasses[i])); // or maybe 1/4 * density * voronoiVolumes[i]
        tripletList.push_back(Triplet<double>(3 * i + 2, 3 * i + 2, 1 / invMasses[i])); // or maybe 1/4 * density * voronoiVolumes[i]
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());

    D = _alpha * M + _beta * K;

    A=M+D*timeStep+K*(timeStep*timeStep);
    
    //Should currently fail since A is empty
    if (ASolver==NULL)
      ASolver=new SimplicialLLT<SparseMatrix<double>>();
    ASolver->analyzePattern(A);
    ASolver->factorize(A);
    
  }
  
  
  //performing the integration step of the soft body.
  void integrateVelocity(double timeStep){
    
    if (isFixed)
      return;
    
    VectorXd fext = VectorXd(invMasses.size() * 3);
    for (int i = 0; i < invMasses.size(); i++)
    {
        fext(i * 3) = 0;
        fext(i * 3 + 1) = -9.8 / invMasses[i];
        fext(i * 3 + 2) = 0;
    }

    VectorXd rhs = M * currVelocities - timeStep * (K * (currPositions - origPositions) - fext);
    currVelocities=ASolver->solve(rhs);
  }
  
