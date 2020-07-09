% Load obsv and ctrb files for individual shift values, then average over
% shift values.
% L. Munoz, 2019

clear variables;

% I can't find a way to change LaTeX-interpreted labels back to the default
% font (Helvetica), so I can try changing everything else to the LaTeX
% font, cmr12. Update: this doesn't work either since the Matlab pdf
% printer doesn't interpret the font correctly. Set defaults globally
% instead.
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

fs = 20; % fontsize
lw = 2; % linewidth for thresholds
lwV = 1; % linewidth for V plots
set(0,'defaultaxesfontsize',fs);

% Choose default or adjusted parameter set:
param = 'def';
%param = 'adj';

allshiftstrings = {'_shift0p4Vnormdepol', '_shift0p6Vnormdepol', '_shift0p8Vnormdepol', '_shift1Vnormdepol', '_shift0p8Vnormrepol', '_shift0p6Vnormrepol', '_shift0p4Vnormrepol', '_shift0p2Vnormrepol', '_shift0p001Vnormrepol', ''};
numshifts = length(allshiftstrings);

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex_nooverwrite = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$'); 
% the statenames appear to get overwritten by later load commands, so use
% the version above instead

ocfolder = ['ocvalues' allshiftstrings{1} '/']; %folder where obsv and ctrb values will be saved
eval(['load ' ocfolder param '/ocfile nbcls']); % load number of BCLs

eigthreshvals = [0.9 0.5]; % list all lower bounds used for modal averages
numeigthresh = length(eigthreshvals);
hoavg1 = figure;
hold on;
hoavg2 = figure;
hold on;
hoavg3 = figure;
hold on;
hoavg4 = figure;
hold on;
hcavg1 = figure;
hold on;

for eigthreshctr = 1:numeigthresh
    eigthr = eigthreshvals(eigthreshctr); 
    eigthreshstr = num2str(eigthr,'%1.2f');
    temp_obsv_avgd_over_bcls_modes_shifts_sc  = zeros(length(statenames),1);
    temp_ctrb_avgd_over_bcls_modes_shifts_sc  = zeros(length(statenames),1);
    
    temp_obsv_avgd_over_modes_shifts_sc  = zeros(length(statenames),nbcls);
    temp_ctrb_avgd_over_modes_shifts_sc  = zeros(length(statenames),nbcls);
    
    for shiftctr = 1:numshifts
        shiftstring = allshiftstrings{shiftctr}
        ocfolder = ['ocvalues' shiftstring '/']; %folder where obsv and ctrb values will be saved
        if eigthr == 0.9 % filename should have suffix for eigthr ~= 0.9
            %        disp('Error: check eigthresh value, which should be 0.9.')
            %        pause
            eval(['load ' ocfolder param '/ocfile *']);
        else
            eval(['load ' ocfolder param '/ocfile_b1to46_eigthresh' eigthreshstr(1) 'p' eigthreshstr(3:end) ' *']); %this assumes the default ocfile used eigthresh = 0.9
        end
        temp_obsv_avgd_over_bcls_modes_shifts_sc  = temp_obsv_avgd_over_bcls_modes_shifts_sc + obsv_avgd_over_bcls_and_modes_sc;
        temp_ctrb_avgd_over_bcls_modes_shifts_sc  = temp_ctrb_avgd_over_bcls_modes_shifts_sc + ctrb_avgd_over_bcls_and_modes_sc;
        temp_obsv_avgd_over_modes_shifts_sc  = temp_obsv_avgd_over_modes_shifts_sc + obsv_avgd_over_modes_sc;
        temp_ctrb_avgd_over_modes_shifts_sc  = temp_ctrb_avgd_over_modes_shifts_sc + ctrb_avgd_over_modes_sc;
    end
    
    obsv_avgd_over_bcls_modes_shifts_sc = temp_obsv_avgd_over_bcls_modes_shifts_sc/numshifts;
    ctrb_avgd_over_bcls_modes_shifts_sc = temp_ctrb_avgd_over_bcls_modes_shifts_sc/numshifts;
    
    obsv_avgd_over_modes_shifts_sc = temp_obsv_avgd_over_modes_shifts_sc/numshifts;
    ctrb_avgd_over_modes_shifts_sc = temp_ctrb_avgd_over_modes_shifts_sc/numshifts;
    
    % Select variables that are more likely to be measurable
    selected_meas_indices = [1 8:10];
    
    % Observability
    % Plot avg obsv for individual measurement variables
    
    figure(hoavg1)
    set(gca,'ColorOrderIndex',1) % Reset color order, otherwise subsequent rounds of plots will pick up where previous colors left off
    if eigthr == 0.9
        plot(selected_bcls_for_fps,obsv_avgd_over_modes_shifts_sc)
    else
        plot(selected_bcls_for_fps,obsv_avgd_over_modes_shifts_sc,'--')
    end
    ylabel('$|\cos \phi_{\,i}|_{k\breve{V}-av}$')
    xlabel('$T$, ms')
    
    % Only show more measurable variables
    figure(hoavg2)
    set(gca,'ColorOrderIndex',1) % Reset color order, otherwise subsequent rounds of plots will pick up where previous colors left off
    if eigthr == 0.9
        p=plot(selected_bcls_for_fps,obsv_avgd_over_modes_shifts_sc(selected_meas_indices,:), 'LineWidth',lwV);
    else
        plot(selected_bcls_for_fps,obsv_avgd_over_modes_shifts_sc(selected_meas_indices,:), '--', 'LineWidth',lwV)
    end
    
    % Average over all measurements
    figure(hoavg3)
    hold on;
    if eigthr == 0.9
        plot(selected_bcls_for_fps,mean(obsv_avgd_over_modes_shifts_sc))
    else
        plot(selected_bcls_for_fps,mean(obsv_avgd_over_modes_shifts_sc),'--')
    end
    ylabel('$|\cos \phi_{\,i}|_{k\breve{V}i-av}$')
    xlabel('$T$, ms')
    title('Obsv. averaged over all measurements') 
    
    % Average over all measurements that seem more accessible (V and intracellular concentrations):
    figure(hoavg4)
    hold on;
    if eigthr == 0.9
        plot(selected_bcls_for_fps,mean(obsv_avgd_over_modes_shifts_sc(selected_meas_indices,:)))
    else
        plot(selected_bcls_for_fps,mean(obsv_avgd_over_modes_shifts_sc(selected_meas_indices,:)),'--')
    end
    ylabel('$|\cos \phi_{\,i}|_{k\breve{V}i-av}$')
    xlabel('$T$, ms')
    title('Obsv. averaged over V and intrac. ionic conc. measurements') 
    
    % Controllability
    % Average over all inputs:
    figure(hcavg1)
    hold on;
    if eigthr == 0.9
        plot(selected_bcls_for_fps,mean(ctrb_avgd_over_modes_shifts_sc))
    else
        plot(selected_bcls_for_fps,mean(ctrb_avgd_over_modes_shifts_sc),'--')
    end
    xlabel('$T$, ms')
    
    % Sort in descending order of obsv or ctrb
    disp('Scaled observability, averaged over bcls, modes above threshold, and shifts:')
    [sortedobsvbms_sc,obsvsortindexbms_sc] = sort(abs(obsv_avgd_over_bcls_modes_shifts_sc),'descend');
    statenames_latex_nooverwrite(obsvsortindexbms_sc,:)
    obsv_avgd_over_bcls_modes_shifts_sc(obsvsortindexbms_sc)
    
    % Print LaTeX tables
    obsv_avgd_over_bcls_modes_shifts_sc_table = cell(numstate,2);
    obsv_avgd_over_bcls_modes_shifts_sc_table(:,1) = cellstr(statenames_latex_nooverwrite(obsvsortindexbms_sc,:));
    obsv_avgd_over_bcls_modes_shifts_sc_table(:,2) = cellstr(num2str(log10(obsv_avgd_over_bcls_modes_shifts_sc(obsvsortindexbms_sc)),'%f'));
    latextableinput.data = obsv_avgd_over_bcls_modes_shifts_sc_table;
    latexTable(latextableinput)
    
    disp('Scaled controllability, averaged over bcls, modes above threshold, and shifts:')
    [sortedctrbbms_sc,ctrbsortindexbms_sc] = sort(abs(ctrb_avgd_over_bcls_modes_shifts_sc),'descend');
    statenames_latex_nooverwrite(ctrbsortindexbms_sc,:)
    ctrb_avgd_over_bcls_modes_shifts_sc(ctrbsortindexbms_sc)
    
    ctrb_avgd_over_bcls_modes_shifts_sc_table = cell(numstate,2);
    ctrb_avgd_over_bcls_modes_shifts_sc_table(:,1) = cellstr(statenames_latex_nooverwrite(ctrbsortindexbms_sc,:));
    ctrb_avgd_over_bcls_modes_shifts_sc_table(:,2) = cellstr(num2str(log10(ctrb_avgd_over_bcls_modes_shifts_sc(ctrbsortindexbms_sc)),'%f'));
    latextableinput.data = ctrb_avgd_over_bcls_modes_shifts_sc_table;
    latexTable(latextableinput)
    
end

figure(hoavg2)
ylabel('$|\cos \phi_{\,i}|_{k\breve{V}-av}$')
xlabel('$T$, ms')
leghandle = legend(p,statenames_latex_nooverwrite(selected_meas_indices,:));
axis([70 1000 -0.03 1.03])
if param == 'def'
    text(750,0.4, 'default ','fontsize',fs)
    figname = ['avgobsv_vs_bcl_def'];
elseif param == 'adj'
    text(750,0.4, 'adjusted ','fontsize',fs)
    figname = ['avgobsv_vs_bcl_adj'];
end
set(leghandle,'Position',[0.6364    0.5729    0.2434   0.2874])
%         set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.55    0.79])
%         set(gca,'Position', [0.165    0.1618    0.7551    0.7632])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,figname,'-dpdf')
saveas(fig,figname)
