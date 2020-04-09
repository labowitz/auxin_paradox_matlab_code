function acti = apopActi(param, xdata)
    x = param;
    acti = (x(4)./x(1)).^2./((x(4)./x(1)).^2 + (x(2).*xdata + 1).^2);
end