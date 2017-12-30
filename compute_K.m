% Compute the K quantity defined in Niels' whitepaper, based on simulation
% results from LRd myocyte model. K should remain constant for certain
% kinds of stimulus protocols. 

clear variables;

datasourcefolder = 'lrddata_fromfixedpointsearch_highres/';
% This folder is the source for simulations
% where ICs were chosen as fixed points

bcls = 1000; 

simoutfname = ['lrddata_1cell_b' num2str(bcls(1))];

eval(['load ' datasourcefolder simoutfname ' Y data subdiv_per_cyc;']);

% define membrane capacitance and ionic charges 
%C = 1; %uF/cm^2
Cm = data.F*data.AF; % membrane capacitance, uF 
qca = 2;
qna = 1;
qk = 1; 

% Define "a" vector
avec = [Cm zeros(1,6) -data.F*qca*data.vmyo -data.F*qna*data.vmyo -data.F*qk*data.vmyo ... 
    -data.F*qca*data.vjsr -data.F*qca*data.vnsr zeros(1,5)]; 

% Compute K 
K = avec*Y; 
time = (1:size(Y,2))*bcls(1)/subdiv_per_cyc; 

figure
plot(time, K)
title('K vs. time, BCL = 1000 ms')
ylabel('K, nF')
xlabel('time, ms')

% Next steps: do a step-by-step accounting of currents, preferrably
% for a non-steady-state trajectory. Check if adjusting the Na-K pump
% strength fixes the Na & K drift you saw during pacedowns at 900ms. 