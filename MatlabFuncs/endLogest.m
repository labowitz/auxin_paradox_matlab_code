function Y = endLogest(r, Y0, K, t, nu)
    Q = (K / Y0) ^ nu - 1;
    Y = K ./ (1 + Q .* exp(-r * nu .* t)) .^ (1/nu);
end
