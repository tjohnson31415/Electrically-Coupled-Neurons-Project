function dVars = NetworkModel_strong( t, vars )

	global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcrit
    global V2 V3 J ICa IK Ileak h_inf m_inf sigma
    
    len = length(vars);
    
    V = vars(1);
    Xs = vars(2:2:len);
    Ys = vars(3:2:len);
 
    dVars = zeros(size(vars));
    
    IKavg = mean( IK(Xs, V) );
    dV = -1/C * (ICa(V) + IKavg + Ileak(V));
    
    dXs =  J(Xs, Ys) - K.*Xs - phi.*ICa(V);
    dYs = -J(Xs, Ys);

    dVars(1)       = dV;
    dVars(2:2:len) = dXs;
    dVars(3:2:len) = dYs;
