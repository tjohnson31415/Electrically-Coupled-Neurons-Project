function dVars = SingleCellModel_Shunt( t, vars )

	global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcirt
    global V2 V3 J ICa IK Ileak h_inf m_inf sigma
    
    if isempty(whos('global','gshunt'))
        % Reversal Potential at the rest potential
        Vshunt = -59;
        % High shunt conductance
        gshunt = 2e4;

        Ishunt = @(V) gshunt * (V - Vshunt);
    else
        global Vshunt gshunt Ishunt
    end
    
    V = vars(1);
    X = vars(2);
    Y = vars(3);

    dV = -1/C * (ICa(V) + IK(X,V) + Ileak(V) + Ishunt(V) );
    dX =  J(X, Y) - K.*X - phi.*ICa(V);
    dY = -J(X, Y);


    dVars = [dV dX dY]';
