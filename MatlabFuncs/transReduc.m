function reduc = transReduc(param, xdata)
    x = param;
    reduc = ((x(2).*xdata +1).^2 + (x(7)./x(1)).^2)./((x(6) + 1).*(x(2).*xdata +1).^2 + (x(7)./x(1)).^2);
end