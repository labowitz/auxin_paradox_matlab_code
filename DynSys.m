function dydt = DynSys(t, y, theta, B, AP, in)
    dydt = zeros(2,1);

    N_cap = 1;
    
    Nu1 = 100 ./ (N_cap * 24); %Ncap cells secrete 100uM auxin a day 
    Delta1 = log(2) / 24; %24 hour half life

    auxin = y(1);
    N = y(2);
    
    eta_A = 1/0.76;

    rate_func2 = @(x,xdata)([
        ((x(1) + x(4)) .* transReduc(x, eta_A, xdata)) - x(4);
        x(1)- x(2) .* apopActi(x, eta_A, xdata)
        ]);
    
    rate_comb = @(x, xdata)(sum(rate_func2(x, xdata), 1) - x(1));
    
    theta(5) = theta(5) * (B / 50);
    theta(3) = theta(3) * (AP / 20).^(1/theta(8));
    
    dydt(1) = Nu1.*N - Delta1.*auxin;
    dydt(2) = rate_comb(theta, auxin) .* N .* (1 - N / N_cap);