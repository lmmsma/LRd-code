% Using unshifted fixed points as ICs, plot all APs on the same plot
% March 2019: Modified to pre-compute measurement timings based on
% normalized voltage thresholds and save timings to folders. Previously
% this was done in fp_compile2_shiftedtraj.m but it's easier to perform a
% visual check on the values here.

clear variables;
fs = 20; % fontsize
lw = 2; % linewidth for thresholds
lwV = 1; % linewidth for V plots
set(0,'defaultaxesfontsize',fs);
% I can't find a way to change LaTeX-interpreted labels back to the default
% font (Helvetica), so I can try changing everything else to the LaTeX
% font, cmr12. Update: this doesn't work either since the Matlab pdf
% printer doesn't interpret the font correctly.
%set(0,'defaultAxesFontName', 'cmr12')
%set(0,'defaultTextFontName', 'cmr12')
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth', lw)

%%folder = ['lrddata_fromfixedpointsearch121718/']; %This folder is the source for fixed points
%folder = ['lrddata_fromfixedpointsearch_highres011719/'];% source for
%1-cyc trajectories
folder = ['lrddata_fromfixedpointsearch_highres_4cyc/']; % source for
%4-cyc trajectories

% list all shift folder name extensions
%shiftstrings = {'','_shift1Vnormdepol'};
%shiftstrings = {'_shift1Vnormdepol'};
%shiftstrings = {'_shift0p2Vnormdepol'};
%Vthreshs = 0.2;
%shiftstrings = {'_shift0p4Vnormdepol','_shift0p6Vnormdepol','_shift0p8Vnormdepol','_shift0p2Vnormrepol','_shift0p4Vnormrepol','_shift0p6Vnormrepol','_shift0p8Vnormrepol','_shift0p001Vnormrepol'};
%Vthreshs = [0.4 0.6 0.8 0.2 0.4 0.6 0.8 0.001]; % Threshold values in normalized membrane potential units.
shiftstrings = {'_shift0p4Vnormdepol','_shift0p6Vnormdepol','_shift0p8Vnormdepol','_shift1Vnormdepol','_shift0p2Vnormrepol','_shift0p4Vnormrepol','_shift0p6Vnormrepol','_shift0p8Vnormrepol','_shift0p001Vnormrepol'};
%shiftstrings = {''};
Vthreshs = [0.4 0.6 0.8 1.0 0.2 0.4 0.6 0.8 0.001]; % Threshold values in normalized membrane potential units.
% This vector must align with shiftstrings

bcls = [1000:-50:400 390:-10:70];%Vector of bcls in ms. Make sure this matches values in fp_compile2*.m files.

nbcls = length(bcls);

for shiftctr = 1:length(shiftstrings)
    shiftstring = shiftstrings{shiftctr}
    
    fixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '/']; %alternative location if using fp_compile2.m to produce high-resolution trajectories
    % Create folder if it doesn't yet exist
    if ~exist(fixedpointfolder,'dir')
        mkdir(fixedpointfolder)
    end
    
    % Save crossing times (or leave as zeros) based on file name
    repolflag = strfind(shiftstring,'repol');
    depolflag = strfind(shiftstring,'depol');
    
    % Put this here since I initially didn't save all values in
    % fp_compile2_shiftedtraj.m.
    i_repol = zeros(1,length(bcls));
    %offset = 7; % this a voltage threshold in mV, not an offset in ms
    % Now add i_depol:
    i_depol = zeros(1,length(bcls));
    
    %    Vthresh = 1; % threshold value for normalized voltages
    %    Vthresh = 0.2; % threshold value for normalized voltages
    Vthresh = Vthreshs(shiftctr); % threshold value for normalized voltages
    
    h=figure
    hold on;
    
    hnorm = figure
    hold on;
    
    if shiftctr == 1 % Produce a figure showing thresholds for the paper
        hnorm_allthresh = figure;
        hold on;
    end
    
    for i = 1:nbcls
        bcl = bcls(i);
        
        fname = [folder 'lrddata_1cell_b' num2str(bcl)];
        % load simulation data for this bcl
        load(fname);
        V = Y(1,:); % membrane potential
        time = (1:length(V))*outstep;
        
        % borrow falling-edge detection from plot_restitution.m:
        % Search for falling edge crossing of threshold
        %    i_repol(i) = find(V(1:end-1)>offset & V(2:end)<=offset);
        
        ampl = abs(max(V)-min(V));
        Vmin = min(V);
        
        Vnorm = (V-Vmin)/ampl; % normalized membrane potential
        
        
        % borrow falling-edge detection from plot_restitution.m:
        % Search for falling and rising edge crossings of threshold
        idtemp1 = find(Vnorm(1:end-1)<Vthresh & Vnorm(2:end)>=Vthresh);
        irtemp1 = find(Vnorm(1:end-1)>Vthresh & Vnorm(2:end)<=Vthresh);
        
        if ~isempty(idtemp1) & depolflag
            idtemp = idtemp1(1); 
            i_depol(i) = idtemp;
        end
        if ~isempty(irtemp1) & repolflag
            irtemp = irtemp1(1); 
            i_repol(i) = irtemp;
        end
        
        figure(h)
        plot(time, V)
        if i_depol(i)
            plot(time(i_depol(i)),V(i_depol(i)),'rx')
        elseif i_repol(i)
            plot(time(i_repol(i)),V(i_repol(i)),'rx')
        end
        
        figure(hnorm)
        plot(time/bcl, (V-Vmin)/ampl)
        if i_depol(i)
            plot(time(i_depol(i))/bcl,Vnorm(i_depol(i)),'rx')
        elseif i_repol(i)
            plot(time(i_repol(i))/bcl,Vnorm(i_repol(i)),'rx')
        end
        
        if shiftctr == 1
            figure(hnorm_allthresh)
            if bcl==1000
                p(1)=plot(time(1:1:end)/bcl, (V(1:1:end)-Vmin)/ampl,'k-','linewidth',5*lwV)
            elseif bcl == 500
                p(2)=plot(time(1:1:end)/bcl, (V(1:1:end)-Vmin)/ampl,'g-','linewidth',5*lwV)
            elseif bcl == 70
                p(3)=plot(time(1:1:end)/bcl, (V(1:1:end)-Vmin)/ampl,'m-','linewidth',5*lwV)
            else
                plot(time/bcl, (V-Vmin)/ampl,'linewidth',lwV)
            end
            
        end
    end
    figure(h)
    xlabel('time, ms')
    ylabel('V, mV')
    
    figure(hnorm)
    xlabel('time/BCL')
    ylabel('(V-minV)/(maxV-minV), dimensionless')
    
    if shiftctr == 1
        figure(hnorm_allthresh)
        plot([-0.025 0.04],[0 0],'b:','linewidth',lw)
        % add labels to orient graph
        for kk = 1:length(shiftstrings)
            if contains(shiftstrings{kk},'depol')
                plot([-0.025 0.04],[Vthreshs(kk) Vthreshs(kk)],'b:','linewidth',lw)
            elseif Vthreshs(kk)==0.8
                plot([0.055 0.35],[Vthreshs(kk) Vthreshs(kk)],'r-.','linewidth',lw)
            elseif Vthreshs(kk)==0.6
                plot([0.1 0.6],[Vthreshs(kk) Vthreshs(kk)],'r-.','linewidth',lw)
            elseif Vthreshs(kk)==0.4
                plot([0.1 0.7],[Vthreshs(kk) Vthreshs(kk)],'r-.','linewidth',lw)
            elseif Vthreshs(kk)==0.2
                plot([0.1 0.8],[Vthreshs(kk) Vthreshs(kk)],'r-.','linewidth',lw)
            elseif Vthreshs(kk)==0.001
                plot([0.1 1.02],[Vthreshs(kk) Vthreshs(kk)],'r-.','linewidth',lw)
            end
        end
        axis([-0.1 1.1 -0.1 1.1])
%        set(gca,'position',[.16 .18 .761 .7769])
        set(gcf,'position',[277 583 1400 395]) % pos on 
        set(gca,'position',[0.08 0.2 0.87 0.75])
        xlabel('t/T, dimensionless')
        yl=ylabel('$\breve{V}(t)$, dimensionless','Interpreter','latex');
        legend(p,'$T$ = 1000 ms', '$T$ = 500 ms', '$T$ = 70 ms')
        figname = 'Vnorm_4cyc';
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,figname,'-dpdf')
        saveas(fig,figname)
    end
    
    pause
    disp('Press any key to continue')
    
    %    save([fixedpointfolder 'threshold_times.mat'],'i_depol','i_repol');
end

% Be careful with this one. Just put it here since I forgot to compute it
% as an array earlier:
%save('lrddata_fromfixedpointsearch_highres_shift7mVrepol/compiled_fp.mat','i_repol','-append')