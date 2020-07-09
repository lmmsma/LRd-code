% LMM: Use my version of the LRd code to plot APD vs. BCL and APD vs. DI 
% using simple thresholding method. This file computes APDs and DIs for the
% last several cycles at each BCL setting and adds them to the appropriate
% plots. 

clear variables;

% I can't find a way to change LaTeX-interpreted labels back to the default
% font (Helvetica), so I can try changing everything else to the LaTeX
% font by setting the LaTeX interpreter as default.
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth', 2)

plottsflagV = 0; % if nonzero, produce a V vs time timeseries plot for each BCL
plottsflagCa = 0; %if nonzero, produce a Ca vs Time timeseries plot for each BCL

equalaxesflag = 0; % Set to nonzero for equal scaling on plot axes
fs = 20; % fontsize
if equalaxesflag
    fs = 24; % fontsize
end
markersize = 10; 
%folder = ['lrddata/pacedownKstim1000_50_pace30000ms_samp0p5ms/'];
%folder = ['Stored_Runs/Default-400_70/'];
%folder = ['lrddata/'];
folder = ['lrddata_pacedown_def/'];

% Load settings recorded from pacedown simulation run
load([folder 'pacedownsettings']) % contains data bcls ncycs pacetimeperbcl
datapacedown = data;

if data.stimflag % 0 if through V, otherwise through K+
    stimtitlestr = 'K+';
else
    stimtitlestr = 'V';
end

% If bcls is defined below, it overwrites the list that was recorded 
% during the actual pacedown. Leave bcls undefined here to compute and 
% plot APDs and DIs for all bcls. If you only want to plot values for 
% certain BCLs, you can list them below: 
%bcls = [250 200 190 180];
%bcls = [1000 190 54];
%bcls = [1000:-50:200 190:-10:70 69:-1:50]; %cycle lengths in ms
%bcls = [50:-1:43]; %cycle lengths in ms
nbcls = length(bcls);
%ncyc = 4; %number of cycles to run (starting from fixed point for that BCL)
%drepol_threshold = -40; %depol and repol threshold, mV
drepol_threshold = -75; %depol and repol threshold, mV
interpflag = 1; % if nonzero,use linear interpolation to improve estimates
%of depol and repol times

allapd = cell(nbcls,1);
alldi = cell(nbcls,1);

for i = 1:nbcls
    bcl = bcls(i);
    ncyc = ncycs(i);
    
    fname = [folder 'lrddata_1cell_b' num2str(bcl)]; 
    % load simulation data for this bcl 
    load(fname);
    
    if ~isequal(datapacedown.stimflag,data.stimflag) % Check for mismatch between stim methods
        return;
    end
    
    % Ending index
    indexend = (ncyc*bcl/outstep)+1;
    
    indexstart_all = [];
    % Extract starting indices of last 4 BCLs
    if ncyc >= 1
        inminus1bcl = (ncyc-1)*bcl/outstep + 1;
        indexstart_all = [inminus1bcl indexstart_all];
        if ncyc >= 2
            inminus2bcl = (ncyc-2)*bcl/outstep + 1;
            indexstart_all = [inminus2bcl indexstart_all];
            if ncyc >= 3
                inminus3bcl = (ncyc-3)*bcl/outstep + 1;
                indexstart_all = [inminus3bcl indexstart_all];
                if ncyc >= 4
                    inminus4bcl = (ncyc-4)*bcl/outstep + 1;
                    indexstart_all = [inminus4bcl indexstart_all];
                end
            end
        end
    end
    
    % lowest starting index for this BCL 
    indexstart = min(indexstart_all);
    
    V = Y(1,indexstart:indexend); % membrane potential over last ncyc cycles
    
    Ca = Y(8,indexstart:indexend); % Calcium concentration over last ncyc cylces
    
    % time range for this part of the pacedown
    time = pacetimeperbcl*(i-1)+(indexstart:indexend)*outstep;
    
    % Search for rising edges
    i_depol = find(V(1:end-1)<=drepol_threshold & V(2:end)>drepol_threshold);
    
    % Search for falling edges
    i_repol = find(V(1:end-1)>drepol_threshold & V(2:end)<=drepol_threshold);
    
    % Since APD and DI are being computed from (typically downsampled) recorded
    % values, include correction factors from linear interpolation. These
    % factors will disappear if depol/repol times coincide with the chosen
    % drepol_threshold.
    % If true repol time is trt, Vrt = V(trt) = drepol_threshold. Let tr =
    % ts(i_repol) and trp1 = tr(i_repol+1). Then the slope of the line
    % connecting V(trp1) and V(tr) is m = (V(trp1)-V(tr))/(trp1-tr). Solve for
    % trt from V(trt)-V(tr)= m(trt-tr). trt = tr +
    % (trp1-tr)*(V(trt)-V(tr))/(V(trp1)-V(tr)).
    % Similarly, for the depolarization, m = (V(tdp1)-V(td))/(tdp1-td),
    % V(tdt)-V(td)= m(tdt-td), and tdt = td +
    % (tdp1-td)*(V(tdt)-V(td))/(V(tdp1)-V(td))
    
    if interpflag
        dfract = (drepol_threshold-V(i_depol))./(V(i_depol+1)-V(i_depol));
        rfract = (drepol_threshold-V(i_repol))./(V(i_repol+1)-V(i_repol));
        % Compute depol/repol times with interpolation
        dtemp = time(i_depol)+dfract.*(time(i_depol+1)-time(i_depol)); % interpolated depol time
        rtemp = time(i_repol)+rfract.*(time(i_repol+1)-time(i_repol)); % interpolated repol time
        %rtemp = asamp*(i_repol+rfract);
    else
        % Compute depol/repol times with no interpolation
        dtemp = time(i_depol);
        rtemp = time(i_repol);
    end
    
    % Compute APDs for ith BCL 
    allapd{i} = rtemp-dtemp;
    % Compute DIs for ith BCL 
    alldi{i} = dtemp(2:end)-rtemp(1:end-1);
    
    % Warning: the script assumes that V starts at rest, so include a check...
    if dtemp(1) > rtemp(1)
        return; % exit if first depol comes after first repol
    end
    
    if plottsflagV % Create time-series plots for V vs T
        figure
        hold on;
        plot(time, V)
        % Mark depol and repol points on plot
        if interpflag
            plot(dtemp,drepol_threshold,'g*')
            plot(rtemp,drepol_threshold,'rs')
        else
            plot(time(i_depol),V(i_depol),'g*')
            plot(time(i_repol),V(i_repol),'rs')
        end
       xlabel('time, ms')
       ylabel('V, mV')
%       title(['BCL = ' num2str(bcl) ' ms, stimulated through ' stimtitlestr])
       title(['BCL = ' num2str(bcl) ' ms'])
    end
    
    if plottsflagCa % Create time-series plots for Ca vs T
        figure
        hold on;
        plot(time, Ca)
       xlabel('time, ms')
       ylabel('Ca, mmol/L')
%       title(['BCL = ' num2str(bcl) ' ms, stimulated through ' stimtitlestr])
       title(['BCL = ' num2str(bcl) ' ms'])
    end
end

% Plot APDs vs. BCLs 
h=figure;
hold on;
for i = 1:nbcls
    plot(bcls(i),allapd{i},'b.', 'MarkerSize', markersize);
end
%xlabel(['BCL, ms'])
xlabel(['$T$, ms'])
ylabel(['APD, ms'])
%title('APD vs. BCL: Default Parameters')
if equalaxesflag
axis equal
set(gcf,'units','normalized','outerposition',[0.01 0.35 0.8 0.48])
%https://www.mathworks.com/matlabcentral/answers/352024-programmatically-performing-expand-axes-to-fill-figure
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1)+0.05, InSet(2), 0.92-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
end
text(700,100, 'default ','fontsize',fs)
xticks([0 200 400 600 800 1000])
axis([0 1050 30 180])
set(gca,'FontSize',fs)
%grid
saveas(gcf,'apdvsbcldef')
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print('apdvsbcldef','-depsc','-r0')
%print('apdvsbcldef','-dpdf','-r0')
print(fig,'apdvsbcldef','-dpdf')

% Plot APDs vs. DIs 
h2=figure;
hold on;
for i = 1:nbcls
    plot(alldi{i},allapd{i}(2:end),'b.', 'MarkerSize', markersize);
end
xlabel(['DI, ms'])
ylabel(['APD, ms'])
%title('APD vs. DI: Default Parameters')
if equalaxesflag
axis equal
set(gcf,'units','normalized','outerposition',[0.01 0.35 0.8 0.48])
%https://www.mathworks.com/matlabcentral/answers/352024-programmatically-performing-expand-axes-to-fill-figure
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1)+0.05, InSet(2), 0.92-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
end
text(700,100, 'default ','fontsize',fs)
xticks([0 200 400 600 800 1000])
axis([-50 1000 30 180])
set(gca,'FontSize',fs)
%grid
saveas(gcf,'apdvsdidef')
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print('apdvsdidef','-depsc','-r0')
%print('apdvsdidef','-dpdf','-r0')
print(fig,'apdvsdidef','-dpdf')



