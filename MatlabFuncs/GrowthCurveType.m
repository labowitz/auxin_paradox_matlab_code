function grType = GrowthCurveType(rateFunc, A_cap)
    testAs = linspace(0, A_cap, 1e4);
    GRs = rateFunc(testAs);
    
    signSwitches = [];
    prevGR = GRs(1);
    for i = 2:numel(GRs)
        if GRs(i) * prevGR < 0
            signSwitches = [signSwitches, testAs(i)];
        elseif GRs(i) * prevGR == 0
            error("Exact zero occored at %0.2e", testAs(i));
        end
            
        prevGR = GRs(i);
    end
    
    if GRs(1) < 0
        switch length(signSwitches)
            case 0
                grType = -1;
            case 1
                grType = 2;
            case 2
                grType = 3;
            otherwise
%                 disp(signSwitches);
                grType = -2;
%                 warning("Started negative and switched %d times", length(signSwitches));
        end
            
    elseif GRs(1) > 0
        switch length(signSwitches)
            case 0
                grType = 1;
            case 1
                grType = 0;
            otherwise
%                 disp(signSwitches);
                grType = -2;
%                 warning("Started positive and switched %d times", length(signSwitches));
        end
        
    else
        error("Exact zero occored at %0.2e", testAs(1));
    end
    
    