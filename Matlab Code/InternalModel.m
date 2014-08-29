function dVars = InternalModel( t, vars )

	global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcrit
    global V2 V3 J ICa IK Ileak h_inf m_inf sigma
    
    X = vars(1);
    Y = vars(2);

    dX =  J(X, Y) - K.*X - phi.*U;
    dY = -J(X, Y);
    
    dVars = [dX dY]';
