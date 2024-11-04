function [logpi,dM] = myLogPi(x,lb,ub,d,s,H)
n = length(x);
dM = myBerry(x,H);
if sum(x>lb)==n && sum(x<=ub)==n
    if sum(H)==2 && H(1)==1 && H(2)==1 %% vp, vs
        if dM(1)>dM(2)
            logpi = -0.5*norm(s.\(d-dM)).^2;
        else
            logpi = -inf;
        end
    elseif sum(H)==3 %% vp, vs, rho
        if dM(1)>dM(2) 
            logpi = -0.5*norm(s.\(d-dM)).^2;
            % fprintf('Satisfying all constraints \n')
        else
            logpi = -inf;
        end
    elseif sum(H)==2 && H(1)==1 || H(2)==1 %% vp OR vs and rho
        if dM(2) < 3100 && dM(2) > 2500
            logpi = -0.5*norm(s.\(d-dM)).^2;
        else
            logpi = -inf;
        end
    elseif sum(H)==1 && H(3)==1 %% rho
        if dM < 3100 && dM > 2500
            logpi = -0.5*norm(s.\(d-dM)).^2;
        else
            logpi = -inf;
        end
    else
        logpi = -0.5*norm(s.\(d-dM)).^2;
    end
else
    logpi = -inf;
end
end