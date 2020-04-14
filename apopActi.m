function acti = apopActi(param, xdata)
    x = param;
    acti = (x(4)./x(1))./((x(4)./x(1)) + (x(2).*xdata + 1));
end