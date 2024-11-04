if sum(H)==2
    if H(1)==1 && H(2) == 1 %% vp and vs
        d = [d(1); d(2)];
        s = [s(1); s(2)];
    elseif H(1)==1 && H(3) == 1 %% vp and rho
        d = [d(1); d(3)];
        s = [s(1); s(3)];
    elseif H(2)==1 && H(3) == 1 %% vs and rho   
        d = [d(2); d(3)];
        s = [s(2); s(3)];
    else
        error('What?')
    end
elseif sum(H)==1
    if H(1)==1 % vp
        d = d(1); s = s(1);
    elseif H(2)==1 % vs
        d = d(2); s = s(2);
    elseif H(3) == 1 % rho
        d = d(3); s = s(3);
    else
        error('What?')
    end
elseif sum(H)==3
    %% use all three data points
else
    error('What?')
end