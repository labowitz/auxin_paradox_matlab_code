function r = revLogest(Y, Y0, K, t, nu)
    r = log(((K ./ Y0).^ nu - 1) ./ ((K ./ Y).^ nu - 1)) / t / nu;
end