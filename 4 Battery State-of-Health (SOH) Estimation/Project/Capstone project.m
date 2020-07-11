
% Plot the clean and noisy SOC used by this project
load readonly/shortDataset.mat
plot(0:199999,100*noisyData,'b'); hold on; plot(0:199999,100*cleanData,'c','linewidth',2);
xlabel('Time (seconds of operation while BMS operational)');
ylabel('State of charge (%)');
title('Noisy SOC estimates used for total-capacity updates (zoom for detail using "xlim")')
legend('Noise-added SOC data','Clean (noise-free) SOC data','location','southwest')
% for example, uncomment the next line to show more detail
% xlim([0 20000])

% Plot the true capacity versus time
load readonly/Qtrue.mat
plot(Qtrue);
xlabel('Time (seconds of operation while BMS operational)');
ylabel('True cell capacity (Ah)');
title('Capacity fading over time for the simulation');

% Determine vectors of total-capacity updates (x,y) with corresponding variances.
% zhatk - roughly 5 million values of noisy SOC versus time
% etiak - roughly 5 million values of noisy coulombic efficiency times measured current
% dz - threshold value in SOC that must be exceeded before we create a new output data pair
% SigmaZ - variance of each SOC measurement error
% Sigmai - variance of each current-sensor measurement error
%
% xi - vector of change in SOC, one value for every total-capacity update
% yi - vector of accumulated ampere-seconds over that interval for every update
% SigmaXi - vector of uncertainty in change in SOC
% SigmaYi - vector of uncertainty in accumulated current
% k - vector of absolute time of beginning of interval used in measurement update
function [xi,yi,SigmaXi,SigmaYi,k] = processData(zhatk,etaik,dz,SigmaZ,Sigmai)
  xi = []; yi = []; SigmaXi = []; SigmaYi = []; k = [];
  ind1 = 1;
  while 1,
    % find next time absolute change in SOC is greater than dz (note: ind2 >= 2)
    ind2 = find(abs(zhatk(ind1:end) - zhatk(ind1)) > dz,1,'first');
    % Next two lines force xi to be somewhat random in the neighborhood of absolute change of dz 
    % Otherwise, xi are biased! (i.e. abs(xi)>dz always, so noise is
    % one-sided and math for xLS breaks down)
    ind2 = ind2 + randi(10); if ind1+ind2-1>length(zhatk), ind2 = length(zhatk) - ind1 + 1; end
    if isempty(ind2), break; end % no more changes greater than dz, so return to calling function
    xi = [xi; zhatk(ind1) - zhatk(ind1+ind2-1)];  % add this "x" value to xi
    yi = [yi; sum(etaik(ind1:ind1+ind2-2))/3600]; % add this "y" value to yi
    SigmaXi = [SigmaXi; 2*SigmaZ];                % add 2*variance of one SOC estimate
    SigmaYi = [SigmaYi; (ind2-1)*Sigmai];         % add accumulated current-sensor variance
    k = [k; ind1];                                % add absolute time of start of interval
    ind1 = ind1+ind2;                             % move starting time for next interval
  end
end

addpath readonly % make sure xLSalgos.m is in the path
load readonly/Qdata.mat % load the pre-processed {xi, yi, SigmaXi, SigmaYi, k} data

% Evaluate xLS algorithms for a set of tuning factors, compute RMS errors, plot results
% dz - index of threshold to use. For example, for threshold of 0.15 use dz = 15
% gamma - forgetting factor, <= 1.0
% method - xLS method number = 1 for WLS, 2 for WTLS, 3 for TLS, 4 for AWTLS
% Qdata - compiled dataset of (xi,yi) etc. loaded from Qdata.mat
% Qtrue - true capacity, loaded from Qtrue.mat
%
% rmsErr - root-mean-squared total-capacity estimation error for this tuning set
function rmsErr = computeResults(dz,gamma,method,Qdata,Qtrue)
  % compute capacity estimates for this dataset, using exact initialization
  [Qhat,SigmaQ]=xLSalgos(Qdata(dz).xi,Qdata(dz).yi,Qdata(dz).SigmaXi,Qdata(dz).SigmaYi,gamma,8,1e-2);
  dataLen = length(Qtrue); % used later on -- number of samples in simulation, nearly 5 million

  % First, note that Qhat updates only every time Qdata(dz).k changes -- it stays constant between xLS updates
  % So, we need to replicate Qhat estimates from one value of "k" until the next.
  Qest = repelems(Qhat(:,method),[1:length(Qdata(dz).k);diff([Qdata(dz).k; length(Qtrue)+1])']);
  Qerr = Qtrue - Qest';
  rmsErr = sqrt(mean(Qerr.^2));

  % Plot results with 3-sigma bounds
  hold on; % use "stairs" to extend estimates until next update automatically
  stairs([Qdata(dz).k; dataLen],[Qhat(:,method); Qhat(end,method)],'b','linewidth',3); % WLS
  % Plot true capacity
  plot(1:dataLen,Qtrue,'k-','linewidth',1);
  % Plot bounds
  stairs([Qdata(dz).k; dataLen],[Qhat(:,method)+3*sqrt(SigmaQ(:,method)); ...
                                 Qhat(end,method)+3*sqrt(SigmaQ(end,method))],'b--','linewidth',0.5);
  stairs([Qdata(dz).k; dataLen],[Qhat(:,method)-3*sqrt(SigmaQ(:,method)); ...
                                 Qhat(end,method)-3*sqrt(SigmaQ(end,method))],'b--','linewidth',0.5);
  
  switch method,
    case 1, title('Capacity estimates, bounds: WLS'); 
    case 2, title('Capacity estimates, bounds: WTLS'); 
    case 3, title('Capacity estimates, bounds: TLS'); 
    case 4, title('Capacity estimates, bounds: AWTLS'); 
  end
  xlabel('Data sample number'); ylabel('Capacity estimate (Ah)');
  legend('Capacity estimate','True capacity','Confidence bounds on estimate')
end  

dz = 15;     % look at xLS performance for threshold = dz/100 = 0.15
gamma = 1.0; % look at xLS performance for forgetting factor gamma = 1.0
method = 2;  % method = 1 for WLS, 2 for WTLS, 3 for TLS, 4 for AWTLS
rmsErr = computeResults(dz,gamma,method,Qdata,Qtrue)
ylim([7 9])

% GRADED FUNCTION (do not modify this line)

% function [dz, gamma] = tunexLS(method)
%
% method - the method to tune: 1 for WLS, 2 for WTLS, 3 for TLS, 4 for AWTLS
%
% dz - your tuning value for "best" difference in SOC threshold for each xLS algorithm
% gamma - your tuning value for "best" forgetting factor for each xLS algorithm

function [dz, gamma] = tunexLS(method)

  % BEGIN MODIFYING CODE AFTER THIS
  switch(method)
    case 1, % for the WLS method ... you are not required to tune these values
      dz = 15;      % A sample value... you are not required to change it
      gamma = 1.0;  % A sample value... you are not required to change it
    case 2, % for the WTLS method ... you are required to tune these values
      dz = 15;      % This is a sample value. You will need to change it.
      gamma = 1.0;  % This is a sample value. You will need to change it.
    case 3, % for the TLS method ... you are required to tune these values
      dz = 38;      % This is a sample value. You will need to change it.
      gamma = 1.0;  % This is a sample value. You will need to change it.
    case 4, % for the AWTLS method ... you are required to tune these values
      dz = 15;      % This is a sample value. You will need to change it.
      gamma = 1.0;  % This is a sample value. You will need to change it.
  end
end  

method = 3; % for example... but change this to test different methods
[dz,gamma] = tunexLS(method);
rmsErr = computeResults(dz,gamma,method,Qdata,Qtrue)

% compute estimate of grade for this tuning...
gradingTable = [...
    0.00    0.0600    0.0790    0.0600
    0.00    0.0570    0.0760    0.0570
    0.00    0.0550    0.0720    0.0550
    0.00    0.0520    0.0700    0.0520
    0.00    0.0510    0.0670    0.0510
    0.00    0.0490    0.0610    0.0490
    0.00    0.0470    0.0570    0.0470
    0.00    0.0410    0.0510    0.0410
    0.00    0.0360    0.0460    0.0360
    0.00    0.0330    0.0330    0.0330];
grade = find(rmsErr < gradingTable(:,method),1,'last');
if isempty(grade), grade = 0; end
fprintf('Assuming that your 3-sigma bounds are okay, you would receive a grade of %d for this method.\n',grade)
