function rate = growthRate_mv(param, A, B, AP)
    eta_A = 1/0.76;
    x = param;
    x(4) = x(4) / 50 * AP;
    x(6) = x(6) / 20 * B;
    rate = ((x(1) + x(5)) .* transReduc(x,A)) - x(5) - x(3) .* apopActi(x, xdata);
end