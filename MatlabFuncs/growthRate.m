function rate = growthRate(param, xdata)
    x = param;
    rate = ((x(1) + x(5)) .* transReduc(x,xdata)) - x(5) - x(3) .* apopActi(x, xdata);
end