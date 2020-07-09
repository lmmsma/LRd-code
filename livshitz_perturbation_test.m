%LMM: Modify pacedown code to try to reproduce Livshitz's dynamic
%perturbation test from his 2009 paper. I could have included this test in
%perturbation_test.m, but it's easier to modify this script.

clear variables;

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');

%folder = ['lrddata/']; % save data here
fixedpointfolder = 'fixedpoints/';%This folder is the source for fixed points
pertfolder = 'perturbation_tests/';%Save perturbation outputs here
modelinputfolder = 'lrddata/'; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script. It is also the default location where
% simulation output is stored.

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

%bcls = [1000:-50:400 390:-10:70];% 69:-1:50]; %full pacedown; cycle lengths in ms
bcls = [1000 1000]; %cycle lengths in ms. Duplicate BCL in order to apply a perturbation partway through

% % Load initial condition, if there happens to be a preferred one for the
% % starting bcl value
% if bcls(1) == 1000
%     % Filenames and label settings based on fixed point type:
%     if strcmp(systemselect, 'solem12')
%         unpertfilename = 'b1000fsolem12_fwde_shift0_newpulse';
%     elseif strcmp(systemselect,'kmonovsolem12')
%         unpertfilename = 'b1000fkmonovsolem12_fwde_shift0_newpulse_relpert1e-7';
%     end
%     % Load the unperturbed fixed point
%     eval(['load ' unpertfilename ' ' systemselect])
%     yinit = eval(systemselect);
% elseif bcls(1) == 400
%     % Filenames and label settings based on fixed point type:
%     if strcmp(systemselect, 'solem12')
%         unpertfilename = 'b400fsolem12_fwde_shift0_newpulse';
%     elseif strcmp(systemselect,'kmonovsolem12')
%         unpertfilename = 'initial_state_400ms';
%     end
%     % Load the unperturbed fixed point
%     eval(['load ' unpertfilename ' ' 'yin'])
%     yinit = eval('yin');
%
% elseif bcls(1) == 900
%     if data.stimflag
%         eval(['load ' folder 'pacedownKstim1000_900_pace30000ms_samp0p5ms/lrddata_1cell_b900 Y'])
%     else
%         eval(['load ' folder 'pacedownVstim1000_900_pace30000ms_samp0p5ms/lrddata_1cell_b900 Y'])
%     end
%     yinit = Y(:,end); % load final condition of previous recording
% elseif bcls(1) == 50
%     if data.stimflag
%         eval(['load ' folder 'pacedownKstim1000_50_pace30000ms_samp0p5ms/lrddata_1cell_b50 Y'])
%     else
%         eval(['load ' folder 'pacedownVstim1000_50_pace30000ms_samp0p5ms/lrddata_1cell_b50 Y'])
%     end
%     yinit = Y(:,end); % load final condition of previous recording
% end
%

% Unlike pacedown code, load previously-computed fixed point for first BCL in
% the list
% load fixed points
eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps'])
bclfpindices = 1:length(selected_bcls_for_fps);
bclind = bclfpindices(selected_bcls_for_fps == bcls(1));
yinit = fp_found(:, bclind);

%pacetimeperbcl = 30*1000; % ms; apply fixed-BCL stimuli for this number of ms

%ncycs = round(pacetimeperbcl./bcls); % approx number of cycles per pacetime interval
% (pacing will be applied for the whole number of bcls that is closest to
% pacetimeperbcl)
%ncycs = [9 round(2*60*1000/bcls(2))]; % Apply perturbation at 10th cycle,
% then run for 2min
%ncycs = [9 round(5*60*1000/bcls(2))]; % Apply perturbation at 10th cycle
ncycs = [9 round(10*60*1000/bcls(2))]; % Apply perturbation at 10th cycle

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
    
    % default filename for saved data
    simoutfname = ['lrddata_1cell_b' num2str(bcl) ];
    % Copy model output to perturbation folder, and append an index to
    % filename so that runs at same bcl won't be overwritten
    eval(['movefile ' modelinputfolder simoutfname '.mat ' pertfolder [simoutfname '_' num2str(i)] '.mat'])
    
    % Next IC = final condition of previous run
    yinit = yend;
    if i == 1 % Apply a perturbation to Na and K only at the 10th cycle.
        % Na is the 9th state variable (decrease by 1.5mM)
        % K is the 10th state variable (increase by 1.5mM)
                yinit(9) = yinit(9) - 1.5;
                yinit(10) = yinit(10) + 1.5;
%         yinit(9) = yinit(9) - 2; % In the text, Livshitz says 2mM, but in the figures, the perturbation sizes look more like 1.5mM
%         yinit(10) = yinit(10) + 2;
    end
end



eval(['save ' pertfolder 'livshitzpertsettings data bcls ncycs'])

% Time elapsed during pacedown
pacedownstoptimer = toc(pacedownstarttimer)

% To produce a time-series plot of the entire pacedown (V vs. t),
% uncomment the following line:
%plotpacedowntimeseries
if data.stimflag % 0 if through V, otherwise through K+
    stimtitlestr = 'K+';
else
    stimtitlestr = 'V';
    
end

for i = 9:10%1:17
    h(i) = figure;
    hold on;
    grid on; 
end

numtimes = zeros(1,length(bcls)+1);
for ibcl=1:length(bcls)
    bcl = bcls(ibcl);
    ncyc = ncycs(ibcl);
    
    fname = [pertfolder '\lrddata_1cell_b' num2str(bcl) '_' num2str(ibcl)]; %  simulation data was saved in this file
    load(fname);
    
    %    indexend = (ncyc*bcl/outstep)+1;
    numtimes(ibcl+1) = size(Y,2);
    for i = 9:10%1:17
        figure(h(i));
        plot(((1:numtimes(ibcl+1)) + numtimes(ibcl))*(bcl/subdiv_per_cyc),Y(i,:),'.-');
        %plot(bcl*ncyc*(ibcl-1)+(1:indexend)*outstep,Y(i,:),'.-');
        xlabel('ms')
        %    ylabel(statenames(i,:))
        ylabel(statenames_latex(i,:),'Interpreter','latex')
        
        title(['BCL = ' num2str(bcl) ' ms, stimulated through ' stimtitlestr])
    end
end