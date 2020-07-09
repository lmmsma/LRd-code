% Compute the K quantity defined in Niels' whitepaper, based on simulation
% results from LRd myocyte model. K should remain constant for certain
% kinds of stimulus protocols. 
% Modified to check all fixed points and compare with corresponding value
% from Livshitz-provided IC. 

clear variables;

% datasourcefolder = 'lrddata_fromfixedpointsearch_highres/';
% This folder is the source for simulations
% where ICs were chosen as fixed points
fixedpointfolder = 'fixedpoints/';%This folder is the source for fixed points
datasourcefolder = 'lrddata_fromfixedpointsearch121718/';%This folder is the source simulation runs where IC = fixed point

constantsLRd_strand % load settings and parameters for model

eval(['load ' fixedpointfolder 'compiled_fp selected_bcls_for_fps'])

bcls = selected_bcls_for_fps; 

% define membrane capacitance and ionic charges 
%C = 1; %uF/cm^2
Cm = data.F*data.AF; % membrane capacitance, uF 
qca = 2;
qna = 1;
qk = 1; 

% Define "a" vector
avec = [Cm zeros(1,6) -data.F*qca*data.vmyo -data.F*qna*data.vmyo -data.F*qk*data.vmyo ... 
    -data.F*qca*data.vjsr -data.F*qca*data.vnsr zeros(1,5)]; 

K = cell(length(bcls),1); 
time = cell(length(bcls),1);

h=figure; 
hold on;

for ii=1:length(bcls)
    simoutfname = ['lrddata_1cell_b' num2str(bcls(ii))]
    eval(['load ' datasourcefolder simoutfname ' Y subdiv_per_cyc']);

    % Compute K
    K{ii} = avec*Y;
    time{ii} = (1:size(Y,2))*bcls(ii)/subdiv_per_cyc;
    figure(h)
    plot(time{ii}, K{ii})
end

% Add corresponding K value for Livshitz's IC
load ..\Livshitz_LRd_Strand_2009\initial_alternans
Xsize=1;Xmesh=50;%% 1 cm=100 cells
x = linspace(0,Xsize,Xmesh);
yinit = ppval(cs_0,x(1));
K_livshitz = avec(1)*yinit(1,1) + avec(2:end)*yinit(2:17,1); 

figure(h) 
p1=plot(0,K_livshitz,'ro');
legend(p1,'Livshitz IC') 
title('K vs. time, all BCLs')
ylabel('K, nF')
xlabel('time, ms')

% Next steps: do a step-by-step accounting of currents, preferrably
% for a non-steady-state trajectory. Check if adjusting the Na-K pump
% strength fixes the Na & K drift you saw during pacedowns at 900ms. 