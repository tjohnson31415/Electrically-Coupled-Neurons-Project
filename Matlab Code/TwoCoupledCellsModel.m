function dVars = TwoCoupledCellsModel( t, vars )

	global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcrit
    global V2 V3 J ICa IK Ileak h_inf m_inf sigma
    global Icoupled gcoupled
    
    gcoupled = 1e4;
    Icoupled  = @(V1, V2) gcoupled * (V1 - V2);
    
    Vcell1 = vars(1);
    Xcell1 = vars(2);
    Ycell1 = vars(3);
    
    Vcell2 = vars(4);
    Xcell2 = vars(5);
    Ycell2 = vars(6);

    dVcell1 = -1/C * (ICa(Vcell1) + IK(Xcell1,Vcell1) + Ileak(Vcell1) + Icoupled(Vcell1, Vcell2));
    dXcell1 =  J(Xcell1, Ycell1) - K.*Xcell1 - phi.*ICa(Vcell1);
    dYcell1 = -J(Xcell1, Ycell1);
    
    dVcell2 = -1/C * (ICa(Vcell2) + IK(Xcell2,Vcell2) + Ileak(Vcell2) + Icoupled(Vcell2, Vcell1));
    dXcell2 =  J(Xcell2, Ycell2) - K.*Xcell2 - phi.*ICa(Vcell2);
    dYcell2 = -J(Xcell2, Ycell2);


    dVars = [dVcell1 dXcell1 dYcell1 dVcell2 dXcell2 dYcell2]';
