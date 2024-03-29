% Differential equastions of the dynmaic system
% Blast and AP1903 are input as ug/ml and uM
function dydt = DynSys_syn(t, y, theta, phi_s, B, AP)
    run('FuncGenerator.m');
    
    dydt = zeros(2,1);
    N_cap = 1;
    nu = 2.69; %Correction for generalized logistic dynmaics
    
    Nu1 = 4.18; %Fitted data
    Delta1 = log(2) / 24; %24 hour half life

    auxin = y(1);
    N = y(2);
    
    % Factor in the AP1903 and blasticidin concentration into the
    % parameters
    theta(5) = theta(5) * (B / 50);
    theta(3) = theta(3) * (AP / 20).^(1/theta(8));
    
    dydt(1) = Nu1.*N - Delta1.*auxin;
    if rate_syn(theta, phi_s, auxin) > 0
        dydt(2) = rate_syn(theta, phi_s, auxin) .* N .* (1 - (N / N_cap)^nu);
    else
        dydt(2) = rate_syn(theta, phi_s, auxin) .* N;
    end
    