
% First, make sure that the ESC toolbox functions are in the path
addpath readonly
load ./readonly/CellModel.mat % load ESC model

% GRADED FUNCTION (do not modify this line)

% function [pChg,pDis] = HPPCpower(z0,T,dT,eta,ns,np,model,limits)
%
% z0  -  the SOC to use when computing discharge and charge resistances RDis and RChg, 
%       and when computing available power based on SOC limits; the same value is 
%       used for every cell in the battery pack
% T   - the temperature to use when computing discharge and charge resistances RDis and RChg, 
%       and when computing available power; the same value is used for every cell in the battery pack
% dT  - the pulse duration to use when computing discharge and charge resistances RDis and RChg;
%       note that the pulse magnitude should use a 10C rate (just like example code from lesson 5.3.3)
% eta - the coulombic efficiency to use when computing available power based on SOC limits
% ns  - the number of cells in series in the battery pack
% np  - the number of cells in parallel in the battery pack
% model - an ESC model type
% limits - design limits on SOC, voltage, current, and power (see code for how these are stored)
%
% pChg - your computed value for charge power (W)
% pDis - your computed value for discharge power (W)
function [pChg,pDis] = HPPCpower(z0,T,dT,eta,ns,np,model,limits)

  % First, unpack the design limits from the "limits" structure.
  % These limits have the same meaning as in the example code for Lesson 5.3.4
  zmin = limits.zMin; zmax = limits.zMax; % Retrieve SOC limits [unitless]
  vmin = limits.vMin; vmax = limits.vMax; % Retrieve voltage limits [V]
  imin = limits.iMin; imax = limits.iMax; % Retrieve current limits [A]
  pmin = limits.pMin; pmax = limits.pMax; % Retrieve design power limits [W]
  
  % BEGIN MODIFYING CODE AFTER THIS

  % NOTE: Resistance calculated using a 10C dis/charge pulse for dT samples
  % You will need to modify this code to work correctly with the specific input
  % parameter list to this function: z0, T, dT ...
  % Note that rChg and rDis are calculated correctly for the default input parameters
  % but are not calculated correctly for non-default inputs. You will need to change
  % this code to calculate rChg and rDis correctly for arbitrary function inputs
  
  Q = getParamESC('QParam',T,model); 
  iChgPulse = 10*Q*[zeros(dT/2,1); -ones(dT,1); zeros(dT/2,1)];  % [A] charge pulse
  iDisPulse = 10*Q*[zeros(dT/2,1);  ones(dT,1); zeros(dT/2,1)];  % [A] discharge pulse
  [vk,~,~,~,~] = simCell(iChgPulse,T,model,1,z0,0,0);
  rChg  = abs((max(vk)-vk(1))/min(iChgPulse));
  [vk,~,~,~,~] = simCell(iDisPulse,T,model,1,z0,0,0);
  rDis  = abs((min(vk)-vk(1))/max(iDisPulse));

  % Now, compute pDis and pChg using rChg and rDis from above, and the equations
  % from the notes. Be sure to incorporate z0, T, dT, eta, ns, np, and the limits
  % correctly (The example code from Lesson 5.3.4 does not implement all of this
  % functionality! You will need to study Lessons 5.3.2 and 5.3.3 to see which
  % equations need to be implemented.)
 % pDis = 1;  % You will need to change this to compute it correctly
 % pChg = -1; % You will need to change this to compute it correctly

soc = z0 ;

% HPPC Power Estimation: Truth
OCV      = OCVfromSOCtemp(soc,T,model);
iDisMaxV = (OCV-vmin)/rDis;
iDisMaxZ = (soc - zmin)*3600*Q/dT;
iDisMax  = max(0,min([iDisMaxV;iDisMaxZ;imax*ones(size(soc))]));
pDisMax  = min(vmin*iDisMax,pmax*ones(size(soc)));
iChgMinV = (OCV-vmax)/rChg;
iChgMinZ = (soc - zmax)*3600*Q/eta/dT;
iChgMin  = max([iChgMinV;iChgMinZ;imin*ones(size(soc))]);
pChgMin  = min(0,max(vmax*iChgMin,pmin*ones(size(soc))));
pDis = pDisMax;
pChg = pChgMin;


end

% This code tests your HPPCpower function using default input values
% You should verify that your code operates for reasonable non-default values as well
% You will be graded on how closely your function results agree with mine
default.z0 = 0.5;
default.T = 25;
default.dT = 10;
default.eta = 1;
default.ns = 1;
default.np = 1;
limits.zMin = 0.1;
limits.zMax = 0.9;
limits.vMin = 2.8;
limits.vMax = 4.3;
limits.iMin = -200;
limits.iMax = 350;
limits.pMin = -1000;
limits.pMax =  1000;
default.limits = limits;
[pChg,pDis] = HPPCpower(default.z0,default.T,default.dT,default.eta,default.ns,default.np,model,default.limits)
% Note that the correct answer for the default set of limits is:
% pChg = -385.18
% pDis = 885.55
% It may also be helpful to know that the correct values for rChg and rDis for the default set of parameters is
% rChg = 3.6787 mOhm
% rDis = 3.7009 mOhm
