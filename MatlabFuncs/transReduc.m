function reduc = transReduc(param, xdata)
    x = param;
%     reduc = ((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./((x(6) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4);
    reduc = ((x(2).*xdata +1) + (x(7)./x(1)))./((x(6) + 1).*(x(2).*xdata +1) + (x(7)./x(1)));
end