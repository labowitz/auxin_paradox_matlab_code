function dydt = DynSys_mut(t, y, theta, B, AP)
    dydt = zeros(3,1);

    auxin = y(1);
    N = y(2); 
    Nmut = y(3);
    
    NormSys1 = DynSys(t, [auxin, N + Nmut], theta, B, AP);
    NormSys2 = DynSys(t, [auxin, N ], theta, B, AP);
    MutSys = DynSys(t, [0, Nmut], theta, B, AP);
    
    dydt(1) = NormSys1(1);
    dydt(2) = NormSys2(2) ;
    dydt(3) = MutSys(2) ;