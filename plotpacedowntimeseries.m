% LNM: Generate timeseries plot of entire pacedown and compute beat-to-beat 
% errors in state vectors for final period-1 and period-2 intervals
% (small values indicate that system is close to steady state) 

clear variables;

%ms = 10; fs = 14; % marker size and font size
statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');

%folder = ['lrddata/'];
folder = ['C:\Users\laura\Google Drive\REU 2017\Data\400-70 ms taufc1 = .1576\']

load([folder 'pacedownsettings']) % contains data bcls ncycs pacetimeperbcl
datapacedown = data;

if data.stimflag % 0 if through V, otherwise through K+
    stimtitlestr = 'K+';
else
    stimtitlestr = 'V';
end


for i = 1:17
    h(i) = figure;
    hold on;
end

for ibcl=1:length(bcls)
    bcl = bcls(ibcl);
    ncyc = ncycs(ibcl);
    
    fname = [folder 'lrddata_1cell_b' num2str(bcl)]; % save simulation data in this file
    load(fname);
    
    if ~isequal(datapacedown.stimflag,data.stimflag) % Check for mismatch between stim methods
        return;
    end
    
    indexend = (ncyc*bcl/outstep)+1;
    
    % Compute cycle-to-cycle errors
    if ncyc >= 1
        inminus1bcl = (ncyc-1)*bcl/outstep + 1;
        error_p2p = Y(:,indexend) - Y(:,inminus1bcl);
        normep2p = norm(error_p2p);
        disp(['For BCL = ' num2str(bcl) ', 1-cycle error is ' num2str(normep2p)])
        if ncyc > 1
            inminus2bcl = (ncyc-2)*bcl/outstep + 1;
            error_p2p2bcl = Y(:,indexend) - Y(:,inminus2bcl);
            normep2p2bcl = norm(error_p2p2bcl);
            disp(['For BCL = ' num2str(bcl) ', 2-cycle error is ' num2str(normep2p2bcl)])
        end
    end
    
    for i = 1:17
        figure(h(i));
        plot(pacetimeperbcl*(ibcl-1)+(1:indexend)*outstep,Y(i,:),'.-');
    end
end

for i = 1:17
    figure(h(i));
    xlabel('ms')
    ylabel(statenames(i,:))
    if length(bcls) == 1
        title(['BCL = ' num2str(bcl) ' ms, stimulated through ' stimtitlestr])
    else
        title([num2str(pacetimeperbcl) ' ms pacedown, stimulated through ' stimtitlestr])
    end
%   figsuffix = ['b' num2str(bcl)];
%    savefig(h(i),[folder deblank(statenames(i,:)) figsuffix '.fig'])
%    saveas(h(i),[folder deblank(statenames(i,:)) figsuffix '.pdf'])
    pause
end