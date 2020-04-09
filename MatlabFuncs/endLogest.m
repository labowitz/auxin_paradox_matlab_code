function Y = endLogest(r, Y0, K, t)
    Y = K * Y0 .* exp(r .* t) ./ (K + Y0 .*  exp(r .* t));
end
