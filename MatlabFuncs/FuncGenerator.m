eta_A = 1/0.76;

rate_func2 = @(x,xdata)([
    ((x(1) + x(4)) .* transReduc(x, eta_A, xdata)) - x(4);
    x(1) - x(2) .* apopActi(x, eta_A, xdata)
    ]);

rate_comb = @(x, xdata)(sum(rate_func2(x, xdata), 1) - x(1));

syn = @(x, phi_s, xdata)(phi_s .* (1 - transReduc(x, eta_A, xdata)) .* apopActi(x, eta_A, xdata));
rate_syn = @(x, phi_s, xdata)(sum(rate_func2(x, xdata), 1) - x(1) - syn(x, phi_s, xdata));