function dVars = NetworkModel( t, vars )

	global K K2 K4 K5 Vm2 Vm3 phi U Vm Tm Vh Th VCa gCa VK gK Vleak gleak C beta Xcrit
    global V2 V3 J ICa IK Ileak h_inf m_inf sigma
    
    len = length(vars);
    Ncells = len/3;
    
    Vs = vars(1:3:len);
    Xs = vars(2:3:len);
    Ys = vars(3:3:len);
    
    gcoupled = @(i,j) (4 + 2.*i + 2.*j).*1e3;
    if Ncells == 2
        gcoupled = @(i,j) 1e4+(i+j)*eps; % got an error without the last bit
    end
    gij = bsxfun(gcoupled, 1:Ncells, (1:Ncells)');
    Vij = bsxfun(@minus, Vs, Vs');
    
    Icoupling = sum( gij .* Vij, 2);
 
    dVars = zeros(size(vars));
    
    dVs = -1/C * (ICa(Vs) + IK(Xs, Vs) + Ileak(Vs) + Icoupling);
    
    dXs =  J(Xs, Ys) - K.*Xs - phi.*ICa(Vs);
    dYs = -J(Xs, Ys);

    dVars(1:3:len) = dVs;
    dVars(2:3:len) = dXs;
    dVars(3:3:len) = dYs;
