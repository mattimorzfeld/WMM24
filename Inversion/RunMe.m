clearvars
close all
clc
Colors = brewermap(8,'Dark2');

%% Data
%% ------------------------------------------------------------------------
d =    [4.1; 2.5;   2537];  % data: vp, vs, rho_m
s =    [0.2; 0.3;    167];  % standard deviation of the noise

H = [1;1;1]; % select which data to use
nData = sum(H);
AssignData
%% ------------------------------------------------------------------------


%% Model
%% ------------------------------------------------------------------------
% Parameter bounds
%   asp  phi   Water content k mu rho_min
lb = [0.0  0.0    0 75.6 25.6 2680]';     % lower bound
ub = [1    0.50   1 80   40   2900]';     % upper bound
n = length(ub);                           % number of unknown parameters
%% ------------------------------------------------------------------------
 

%% Set up the MCMC and inversion
%% ------------------------------------------------------------------------
% Make a function for the log posterior
logpi=@(x)myLogPi(x,lb,ub,d,s,H);

%% Run the emcee hammer
Ne = 3*n;
Xo = zeros(n,Ne);
numModels = 0;
% initialize with parameters that satisfy all constraints
go = 1;
counter = 0;
while go == 1
    xo = lb+(ub-lb).*rand(n,1);
    if isfinite(myLogPi(xo,lb,ub,d,s,H))
        counter = counter+1;
        Xo(:,counter) = xo;
        if counter>=Ne
            go = 0;
        end
    end
end
% Cold start
Nsteps = 5e2;
[X, ~, LogPi,~]=myHammer(Nsteps,Xo,2.6,logpi,H);
Xrs = reshape(X,[n,size(X,2)*size(X,3)]);
LogPirs = reshape(LogPi,[1,size(LogPi,1)*size(LogPi,2)]);
RMSE = sqrt(2*LogPirs/length(d));
Xrs = Xrs(:,RMSE<3);
Xo = Xrs(:,randi(length(Xrs),Ne,1));

% Warm start
Nsteps = 5e4;
[X, D, LogPi, AccRatio]=myHammer(Nsteps,Xo,2.6,logpi,nData);

% MCMC aftermath
BurnIn = 1e3;
X = X(:,BurnIn:end,:);
D = D(:,BurnIn:end,:);
LogPi = LogPi(:,BurnIn:end);
Xrs = reshape(X,[n,size(X,2)*size(X,3)]);
Drs = reshape(D,[nData,size(D,2)*size(D,3)]);
LogPirs = reshape(LogPi,[1,size(LogPi,1)*size(LogPi,2)]);
RMSE = sqrt(2*LogPirs/length(d));
%% ------------------------------------------------------------------------

%% Plots
%% ------------------------------------------------------------------------
PlotScript
%% ------------------------------------------------------------------------

%% Display results
%% --------------------------------------------
m = mean(Xrs,2);
sp = std(Xrs,[],2);

disp(' '), disp(' ')
fprintf('Asp. ratio: %g +/- %g \n',m(1),sp(1))
fprintf('Porosity: %g +/- %g \n',m(2),sp(2))
fprintf('Water: %g +/- %g \n',m(3),sp(3))
fprintf('K: %g +/- %g \n',m(4),sp(4))
fprintf('mu: %g +/- %g \n',m(5),sp(5))
fprintf('rho_min: %g +/- %g \n',m(6),sp(6))
%% --------------------------------------------

% save('Results_weighted.mat')
