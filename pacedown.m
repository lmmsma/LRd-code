3

%LMM: perform pacedown protocol on LRD model
clear variables;

folder = ['lrddata/']; % save data here

constantsLRd_strand % load settings and parameters for model

% The variable data.stimflag is defined in constantsLRd_strand. 
% If data.stimflag is 0, apply biphasic stimuli to the cell through
% V dynamics. Otherwise, apply monophasic stimuli through K+-dynamics.
if data.stimflag 
    systemselect = 'kmonovsolem12' 
    %The systemselect string is used in certain filenames
else
    systemselect = 'solem12'
end

% BCL steps loosely based on Koller/Riccio/Gilmour dynamic protocol:
%bcls = [1000 90];
bcls = [1000:-50:400 390:-10:70];% 69:-1:50]; %full pacedown; cycle lengths in ms
%bcls = [50:-1:40]; %cycle lengths in ms

% Load initial condition, if there happens to be a preferred one for the 
% starting bcl value 
if bcls(1) == 1000
    % Filenames and label settings based on fixed point type:
    if strcmp(systemselect, 'solem12')
        unpertfilename = 'b1000fsolem12_fwde_shift0_newpulse';
    elseif strcmp(systemselect,'kmonovsolem12')
        unpertfilename = 'b1000fkmonovsolem12_fwde_shift0_newpulse_relpert1e-7';
    end
    % Load the unperturbed fixed point
    eval(['load ' unpertfilename ' ' systemselect])
    yinit = eval(systemselect);
elseif bcls(1) == 400
    % Filenames and label settings based on fixed point type:
    if strcmp(systemselect, 'solem12')
        unpertfilename = 'b400fsolem12_fwde_shift0_newpulse';
    elseif strcmp(systemselect,'kmonovsolem12')
        unpertfilename = 'initial_state_400ms';
    end
    % Load the unperturbed fixed point
    eval(['load ' unpertfilename ' ' 'yin'])
    yinit = eval('yin');

elseif bcls(1) == 900
    if data.stimflag
        eval(['load ' folder 'pacedownKstim1000_900_pace30000ms_samp0p5ms/lrddata_1cell_b900 Y'])
    else
        eval(['load ' folder 'pacedownVstim1000_900_pace30000ms_samp0p5ms/lrddata_1cell_b900 Y'])
    end
    yinit = Y(:,end); % load final condition of previous recording
elseif bcls(1) == 50
    if data.stimflag
        eval(['load ' folder 'pacedownKstim1000_50_pace30000ms_samp0p5ms/lrddata_1cell_b50 Y'])
    else
        eval(['load ' folder 'pacedownVstim1000_50_pace30000ms_samp0p5ms/lrddata_1cell_b50 Y'])
    end
    yinit = Y(:,end); % load final condition of previous recording
end


% if bcls(1) ~= 1000
%     bflag = input({'Are you sure you want to use the BCL = 1000 ms', ...
%         'fixed point as the initial condition?', 'Enter 1 to proceed, 0 to quit: '});
%     if ~bflag
%         return;
%     end
% end


% For now, try fixed-time pacedown (whole number of cycles closest to
% 30sec)

pacetimeperbcl = 30*1000; % ms; apply fixed-BCL stimuli for this number of ms

ncycs = round(pacetimeperbcl./bcls); % approx number of cycles per pacetime interval
% (pacing will be applied for the whole number of bcls that is closest to 
% pacetimeperbcl)

pacedownstarttimer = tic; 
for i = 1:length(bcls)
    bcl = bcls(i);
    % print current BCL to screen 
    disp(['BCL = ' num2str(bcl) ' ms'])
    ncyc = ncycs(i);
    
    %    subdiv_per_cyc = bcl/data.dt;
    %    subdiv_per_cyc = bcl;
    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    %subdiv_per_cyc = 5*bcl;
    
    % Run simulation
    yend = lrd_p2p(yinit, bcl, ncyc, subdiv_per_cyc);
    
    % Next IC = final condition of previous run
    yinit = yend;
end

eval(['save ' folder 'pacedownsettings data bcls ncycs pacetimeperbcl'])

% Time elapsed during pacedown 
pacedownstoptimer = toc(pacedownstarttimer)

% To produce a time-series plot of the entire pacedown (V vs. t), 
% uncomment the following line: 
plotpacedowntimeseries 