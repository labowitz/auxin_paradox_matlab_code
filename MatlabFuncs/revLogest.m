function r = revLogest(Y, Y0, K, t)
    r = log(Y .* (K - Y0)  ./ Y0 ./ (K - Y)) / t;
end