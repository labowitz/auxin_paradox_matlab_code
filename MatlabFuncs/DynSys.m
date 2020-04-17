% Differential equastions of the dynmaic system
% Blast and AP1903 are input as ug/ml and uM
function dydt = DynSys(t, y, theta, B, AP)
    run('FuncGenerator.m');
    
    dydt = zeros(2,1);

    N_cap = 1;
    nu = 3.14; %Correction for generalized logistic dynmaics
    
    Nu1 = 100 ./ (N_cap * 24); %Ncap cells secrete 100uM auxin a day 
    Delta1 = log(2) / 24; %24 hour half life

    auxin = y(1);
    N = y(2);
    
    % Factor in the AP1903 and blasticidin concentration into the
    % parameters
    theta(5) = theta(5) * (B / 50);
    theta(3) = theta(3) * (AP / 20).^(1/theta(8));
    
    dydt(1) = Nu1.*N - Delta1.*auxin;
    dydt(2) = rate_comb(theta, auxin) .* N .* (1 - (N / N_cap)^nu);