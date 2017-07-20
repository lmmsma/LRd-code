

clear variables;

plottsflagV = 0; % if nonzero, produce a V vs time timeseries plot for each BCL
plottsflagCa = 0; %if nonzero, produce a Ca vs Time timeseries plot for each BCL

%folder = ['Stored_Runs/Default-400_70/'];
folder = ['Stored_Runs/Default-1000_70/'];

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
drepol_threshold = -40; %depol and repol threshold, mV
interpflag = 1; % if nonzero,use linear interpolation to improve estimates
%of depol and repol times

allapd = cell(nbcls,1);
alldi = cell(nbcls,1);

i = find(bcls == input('Enter BCL\n'))
bcl = bcls(i);
ncyc = ncycs(i);

fname = [folder 'lrddata_1cell_b' num2str(bcl)]; 
% load simulation data for this bcl 
load(fname);

if ~isequal(datapacedown.stimflag,data.stimflag) % Check for mismatch between stim methods
    return;
end

V = Y(1,:); % membrane potential over last ncyc cycles

Ca = Y(8,:); % Calcium concentration over last ncyc cylces

% time range for this part of the pacedown
time = pacetimeperbcl*(i-1)+((1:ncyc*bcl/.5)+1)*.5;

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

l= length(allapd{i});
apd_diff=allapd{i}(2:l)-allapd{i}(1:l-1)
% figure
% hold on;
% for j=1:l-1
%     plot(j, apd_diff(j), 'b*')
% end
% 
% xlabel('Cycle Number')
% ylabel('APD Diff, ms')
% title(['BCL = ' num2str(bcl) ' ms APD Difference'])
% grid
% hold off;

figure
hold on;
for j=1:l-1
    scatter(j, log10(abs(apd_diff(j))), 'b*')
end

xlabel('Cycle Number')
ylabel('log(APD Diff), ms')
title(['BCL = ' num2str(bcl) ' ms APD Difference Semilog'])
grid