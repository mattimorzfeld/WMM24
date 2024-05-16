function out = myBerry(theta,H)
asp = [1 theta(1)];
x_phi = theta(2);
rock_vol = 1-theta(2);
x = [rock_vol, x_phi]; 
rock_density = theta(6).* rock_vol; % density of solid phase: mineral density of basalt is used -- 2900
gas_density = 0.020.* x_phi;        % density of fluid phase: gas density is used -- 0.020
rhob1 = rock_density + gas_density; % bulk density
P_water = theta(3);                 % percentage of water in pore space
k  = [theta(4)*1e9 0];
mu = [theta(5)*1e9 0];

[~,~,vp,vs,rhob,~] = berryscm(k,mu,asp,x,rhob1,P_water);



if sum(H) == 3
    out = [[vp; vs]./1e3; rhob];
elseif sum(H)==2
    if H(1)==1 && H(2) == 1 %% vp and vs
         out = [vp; vs]./1e3;
    elseif H(1)==1 && H(3) == 1 %% vp and rho
        out = [vp/1e3; rhob];
    elseif H(2)==1 && H(3) == 1 %% vs and rho   
        out = [vs/1e3; rhob];
    else
        error('What?')
    end
elseif sum(H)==1
    if H(1)==1
        out = vp/1e3;
    elseif H(2)==1
       out = vs/1e3;
    elseif H(3) == 1
        out = rhob;
    else
        error('What?')
    end
end