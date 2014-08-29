function dVars = SingleCellModel( t, vars )

	global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcrit
    global V2 V3 J ICa IK Ileak h_inf m_inf sigma
    
    V = vars(1);
    X = vars(2);
    Y = vars(3);

    dV = -1/C * (ICa(V) + IK(X,V) + Ileak(V));
    dX =  J(X, Y) - K.*X - phi.*ICa(V);
    dY = -J(X, Y);


    dVars = [dV dX dY]';
