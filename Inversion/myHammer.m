function [X, D, LogPi, AccRatio] = myHammer(Nsteps,Xo,a,logpi,H)
nData = sum(H);
n = size(Xo,1);
Ne = size(Xo,2);
X = zeros(n,Nsteps,Ne);
D = zeros(nData,Nsteps,Ne);
LogPi = zeros(Ne,Nsteps);
%% initialization
for kk=1:Ne
    X(:,1,kk) = Xo(:,kk);
end

accMovesTotal = 0;
%% move NSteps
for Step=1:Nsteps-1
    cX = squeeze(X(:,Step,:));
    [cX,DX,clogpi,accMoves] = MoveEnsemble(cX,a,logpi,nData);
    X(:,Step+1,:) = cX;
    D(:,Step+1,:) = DX;
    if sum(isnan(cX))>0
        break
    end
    LogPi(:,Step+1)= clogpi; 
    accMovesTotal = accMovesTotal + accMoves;
    if ~mod(Step,1e2)
        fprintf('Acc. ratio at step %g/%g: %g\r',Step,Nsteps,accMovesTotal/Step/Ne)
    end
end
AccRatio = accMovesTotal/Nsteps/Ne;