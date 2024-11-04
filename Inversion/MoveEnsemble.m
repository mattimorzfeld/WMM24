function [cX,D,clogpi,accMoves] = MoveEnsemble(cX,a,logpi,nData)
Ne = size(cX,2);
n = size(cX,1);
z = sampleG(Ne,a);
accMoves = 0;
clogpi = zeros(Ne,1);
D = zeros(nData,Ne);
for jj=1:Ne
    x = cX(:,jj);
    pind = findPartner(jj,Ne);
    xj = cX(:,pind);
    % stretch move
    y = xj + z(jj)*(x-xj);
    % accept/reject
    [logpiy,dy] = logpi(y);
    [logpix,dx] = logpi(x);
    logalpha = (n-1)*log(z(jj))+logpiy-logpix;
    if rand<min(1,exp(logalpha)) % accept
        cX(:,jj) = y;
        D(:,jj) = dy;
        accMoves = accMoves + 1;
    else
        D(:,jj) = dx;
    end
    clogpi(jj) = -logpi(cX(:,jj));
end