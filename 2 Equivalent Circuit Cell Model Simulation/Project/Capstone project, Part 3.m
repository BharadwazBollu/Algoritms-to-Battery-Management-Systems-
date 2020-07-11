
% Initialize workspace, load the E2 circuit model as well as the E2 dynamic data
addpath readonly
load readonly/E2model.mat; % load parameter values already created for the E2 cell
load readonly/E2_DYN_P25.mat; % load raw test data for the E2 cell at 25 degC

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (2000:deltaT:2500); % select short segment to speed up simulation
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;

% GRADED FUNCTION (do not modify this line)

% function [vpack,vcell,icell,zcell,qcell,rcell] = simSCMTemp(Ns,Np,current,temp,deltaT,model)
%
% Simulate series-connected-module packs (cells are connected in series
% to make modules; these modules are connected in parallel to make packs).
% Note: This function must work for models having a single R-C pair. It is
% not required to work for models having multiple R-C pairs.
%
% Ns - number of cells connected in series to make a module
% Np - number of modules connected in parallel in each pack
% current - battery pack current, where (+) is discharge. Size is N x 1.
% temp  - temperature (degC). Size is N x 1 (all cells at same temperature).
% deltaT = sampling interval in data (s)
% model - standard model structure
% delta - variability vector: variable initial SOC if delta(1)==1; 
%                         variable total capacity if delta(2)==1;
%                         variable series resistance if delta(3)==1;
%
% vpack - battery pack voltage. Size is N x 1.
% vcell - individual cell voltages. Size is N x Ns x Np
% icell - individual cell currents. Size is N x Ns x Np
% zcell - individual cell states of charge. Size is N x Ns x Np
% qcell - individual cell capacities. Size is N x Ns x Np
% rcell - individual cell series resistances. Size is N x Ns x Np

function [vpack,vcell,icell,zcell,qcell,rcell] = simSCMTemp(Ns,Np,current,temp,deltaT,model,delta) %#ok<FNDEF>

  % Force current to be column vector in case user entered data incorrectly
  current = current(:); N = length(current); temp = temp(:);
  % Initialize function outputs
  vpack = zeros(N,1); vcell = zeros(N,Ns,Np); icell = zeros(N,Ns,Np);
  zcell = zeros(N,Ns,Np); qcell = zeros(N,Ns,Np); rcell = zeros(N,Ns,Np);
  
  % Do some error checking on the function inputs
  Nr = length(getParamESC('RCParam',25,model)); % number of R-C pairs.
  if Nr ~= 1,
    error('This code does not work for models having multiple R-C pairs.');
  end
  if length(temp) ~= N,
    error('Input "temp" vector not the correct dimension.');
  end
  
  % Initialize states for ESC cell model
  if delta(1),
    z = reshape(linspace(0.3,0.7,Ns*Np),Ns,Np); % Different initial SOCs
  else
    z = 0.5*ones(Ns,Np);
  end
  irc = zeros(Ns,Np);
  h   = zeros(Ns,Np);

  % BEGIN MODIFYING CODE AFTER THIS
  for k = 1:N,
  
  T = temp(k); % sample code uses only single temperature -- you will need to change this!
  
  % The code reproduced below as a starting point for you is modified from
  % "simSCM.m".
  %
  % Get model parameters from model structure -- notice that these retreive
  % parameter values for only a single temperature. You will need to change
  % this to load parameter values for all temperatures in temp.
  
  % Default initialization for cells within the pack
  q  = getParamESC('QParam',T,model)*ones(Ns,Np); 
  rc = exp(-deltaT./abs(getParamESC('RCParam',T,model)))'*ones(Ns,Np);
  r  = (getParamESC('RParam',T,model))';
  m  = getParamESC('MParam',T,model)*ones(Ns,Np);
  g  = getParamESC('GParam',T,model)*ones(Ns,Np);
  r0 = getParamESC('R0Param',T,model)*ones(Ns,Np); 
  rt = 0.000125; % 125 microOhm resistance for each tab

  % How to modify capacity at any temperature if that input option is set... 
  % Don't change this functionality since the grader assumes it.
  if delta(2), 
    q = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*q; 
  end
  % How to modify resistances at any temperature if that input option is set...
  % Don't change this functionality since the grader assumes it.
  if delta(3),
    r0 = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*r0; 
  end
  r0 = r0 + 2*rt; % add tab resistance to cell resistance

  % Okay... now to simulate pack performance using ESC cell model.
  
    v = OCVfromSOCtemp(z,T,model); % get OCV for each cell: Ns * Np matrix
    v = v + m.*h - r.*irc; % add in capacitor voltages and hysteresis

    V = (sum(sum(v,1)./sum(r0,1),2)-current(k))./sum(1./sum(r0,1),2); % Bus V
    ik = (sum(v,1)-repmat(V,1,Np))./sum(r0,1); % 1*Np cell currents
    ik = repmat(ik,Ns,1); % Ns*Np cell currents

    z = z - (1/3600)*ik./q;  % Update each cell SOC
    irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
    fac = exp(-abs(g.*ik)./(3600*q));
    h = fac.*h + (fac-1).*sign(ik); % Update hysteresis voltages

    vpack(k)     = V; % Store pack voltage
    vcell(k,:,:) = v - ik.*r0; % Store cell voltages
    zcell(k,:,:) = z; % Store cell SOCs
    icell(k,:,:) = ik; % Store cell currents
    qcell(k,:,:) = q; % Store cell capacities
    rcell(k,:,:) = r0; % Store cell resistances
  end % for k

  % FINISH MODIFYING CODE BEFORE THIS

end
% END GRADED FUNCTION

% Execute simCell to determine voltage and other internal states/variables

temp = 25*ones(size(current)); % for now, use constant 25 degC temperature.
% temp = linspace(25,45,length(current)); % uncomment to simulate temperature ramp 

Ns = 2; Np = 2; 
[vpack,vcell,icell,zcell,qcell,rcell] = simSCMTemp(Ns,Np,current,temp,deltaT,model,[1,0,0]);

% The grader will input custom pack configurations, pack current, temperature, and combination 
% of delta inputs and compare output to expected output. You will not be told what these custom 
% inputs will be, but you can still test basic functionality and visualize how the cell behavior 
% changes for different temperatures

% For example, for visualization purposes, plot 

% Plot the individual SOC vs. time for all cells in all 
% series PCMs. There is one subplot for each PCM.
t = (0:(length(zcell(:,:,1))-1))/60; 
xplots = round(1.0*ceil(sqrt(Ns))); yplots = ceil(Ns/xplots); means = [];
for k = 1:Np,
  zr=squeeze(100*zcell(:,:,k));
  subplot(yplots,xplots,k); plot(t,zr); axis([0 ceil(max(t)) 0 100]);
  title(sprintf('Cells in SCM %d',k)); 
  ylabel('SOC (%)'); xlabel('Time (min)'); 
end
