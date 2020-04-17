function reduc = transReduc(param, eta_A, xdata)
    x = param;
    reduc = ((eta_A.*xdata + 1).^x(7) + (x(6)./x(1)).^x(7))./((x(5) + 1).*(eta_A.*xdata +1).^x(7) + (x(6)./x(1)).^x(7));
end