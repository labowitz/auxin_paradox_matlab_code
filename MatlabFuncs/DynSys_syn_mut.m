function dydt = DynSys_syn_mut(t, y, theta, phi_s, B, AP)
    dydt = zeros(3,1);

    auxin = y(1);
    N = y(2); 
    Nmut = y(3);
    
    NormSys = DynSys_syn(t, [auxin, N + Nmut], theta, phi_s, B, AP);
    MutSys = DynSys_syn(t, [0, N + Nmut], theta, phi_s, B, AP);
    
    dydt(1) = NormSys(1);
    dydt(2) = NormSys(2) .* N ./ (N + Nmut);
    dydt(3) = MutSys(2) .* Nmut ./ (N + Nmut);