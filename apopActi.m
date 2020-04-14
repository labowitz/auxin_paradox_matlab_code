function acti = apopActi(param, eta_A, xdata)
    x = param;
    acti = (x(3)./x(1)).^x(8)./((x(3)./x(1)).^x(8) + (eta_A.*xdata + 1).^x(8));
%     acti = (x(4)./x(1))./((x(4)./x(1)) + (eta_A.*xdata + 1));
end