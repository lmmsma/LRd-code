% Plot results of test_nonlinear_luen.m more clearly.

clear variables

fs = 24; % fontsize
lw = 1.5; % linewidth for thresholds
ms = 10; % markersize
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

largeampindices = [1 8:12 17];
smallampindices = [2:7 13:16];

%statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
%statename_symbols = char('v-',  'x-', 'o-', '*-', 'v-', 's-', 'p-',    '^-',               '>-',          '<-',       'p-',                'h-','^-','d-','+-','<-','o-');
%statename_symbols = char('-',  '-', '-', '-', '-', '-', '-',    '-',               '-',          '-',       '-',                '-','-','--','--','--','-');
statename_symbols = char('v-',  'x-', 'o-', '*-', 'v-', 's-', 'p-',    '^-',               '>-',          '<-',       'p-',                'h-','^-','d--','+--','<--','o-');

pertfolder = 'D:/nonlinear_luen_tests/';%Save perturbation outputs here
%loadfilename = 'nonl_luen_results_b200_m0_shift0p2Vnormrepol_ncyc30_pfl4_pdirrand_psz-5.mat';
%loadfilename = 'nonl_luen_results_b200_m1_shift0p2Vnormrepol_ncyc30_pfl4_pdirrand_psz-5.mat';
loadfilename = 'nonl_luen_results_b200_m10_shift0p2Vnormrepol_ncyc30_pfl4_pdirrand_psz-5.mat';
eval(['load ' pertfolder loadfilename]);

%     figure
%     hold on;
%     plot((1:numtimes)*(bcl/subdiv_per_cyc), Smat*trajectorydiff)
%     xlabel('time, ms')
%     ylabel('estimation errors, nondimensionalized')
%     title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex');
%
% Some error quantities may not have been saved in certain earlier runs. If so, recompute
% them here:
if ~exist('trajectorydiffrpd','var')
    Yold = Y;
    eval(['load ' fixedpointfolder simoutfname ' Y;']);
    Yfp_traj = Y;
    %relative "percent" (really just fractional) difference error
    trajectorydiffrpd = 2*(Yold-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1)))./(abs(Yold)+abs(Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1))));
end
if ~exist('trajectorydiff','var')
    Yold = Y;
    eval(['load ' fixedpointfolder simoutfname ' Y;']);
    Yfp_traj = Y;
    trajectorydiff = Yold-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1));
end

figure
hold on;
plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiffrpd)
xlabel('time, ms')
%ylabel('RPD estimation errors')
ylabel('$E_{RD}$')
title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex');

% Plot abs val RPD errors 
figure
subplot(2,1,1)
hold on;
for ii=1:length(largeampindices)
%    p1(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiffrpd(largeampindices(ii),:),statename_symbols(largeampindices(ii),:),'markersize',ms);
    p1(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), abs(trajectorydiffrpd(largeampindices(ii),:)),statename_symbols(largeampindices(ii),:),'markersize',ms);
    p1(ii).MarkerIndices = 1:round(bcl/data.dt):numtimes; 
end
set(gca, 'YScale', 'log')
ylabel('$|E_{RD}|$')
%legend(p1,statenames_latex(largeampindices,:))
if param == 'def'
%    text(200,80, 'default ','fontsize',fs)
    figname = ['nonl_est_err_def_b'  num2str(bcl) '_m' num2str(measurementindex) shiftstring];
elseif param == 'adj'
%    text(200,80, 'adjusted ','fontsize',fs)
    figname = ['nonl_est_err_adj_b'  num2str(bcl) '_m' measlabel shiftstring];
end
%title(['$T = $' num2str(bcl) ' ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 
axis([0 6000 10^-16.5 10^0.5])
yticks([10^-15 10^-10 10^-5 10^0])

subplot(2,1,2)
hold on;
for ii=1:length(smallampindices)
%    p2(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiffrpd(smallampindices(ii),:),statename_symbols(smallampindices(ii),:),'markersize',ms);
    p2(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), abs(trajectorydiffrpd(smallampindices(ii),:)),statename_symbols(smallampindices(ii),:),'markersize',ms);
    p2(ii).MarkerIndices = 1:round(bcl/data.dt):numtimes; 
end
set(gca, 'YScale', 'log')
ylabel('$|E_{RD}|$')
%legend(statenames_latex(smallampindices,:))
xlabel('time, ms')
axis([0 6000 10^-16.5 10^0.5])
yticks([10^-15 10^-10 10^-5 10^0])

% %        set(gca,'position',[.16 .18 .761 .7769])
set(gcf,'position',[106.3333  120.3333  894.0000  737.3333])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,figname,'-dpdf')
saveas(fig,figname)

legendflag = 1; 
% Plot abs val RPD errors, zoomed in 
figure
subplot(2,1,1)
hold on;
for ii=1:length(largeampindices)
    p1(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), abs(trajectorydiffrpd(largeampindices(ii),:)),statename_symbols(largeampindices(ii),:),'markersize',ms);
    p1(ii).MarkerIndices = 1:round(0.5*bcl/data.dt):numtimes; 
end
set(gca, 'YScale', 'log')
ylabel('$|E_{RD}|$')
%legend(p1,statenames_latex(largeampindices,:),'location','EastOutside')
if legendflag
    l1=legend(p1,statenames_latex(largeampindices,:)); 
end
if param == 'def'
%    text(200,80, 'default ','fontsize',fs)
    figname = ['nonl_est_err_def_b'  num2str(bcl) '_m' num2str(measurementindex) shiftstring '_zoom'];
elseif param == 'adj'
%    text(200,80, 'adjusted ','fontsize',fs)
    figname = ['nonl_est_err_adj_b'  num2str(bcl) '_m' num2str(measurementindex) shiftstring '_zoom'];
end
%title(['$T = $' num2str(bcl) ' ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 
%axis([5000 6300 10^-16 10^0])
axis([5200 6200 10^-16.5 10^0.5])
yticks([10^-15 10^-10 10^-5 10^0])

subplot(2,1,2)
hold on;
for ii=1:length(smallampindices)
%    p2(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiffrpd(smallampindices(ii),:),statename_symbols(smallampindices(ii),:),'markersize',ms);
    p2(ii) = plot((1:numtimes)*(bcl/subdiv_per_cyc), abs(trajectorydiffrpd(smallampindices(ii),:)),statename_symbols(smallampindices(ii),:),'markersize',ms);
    p2(ii).MarkerIndices = 1:round(0.5*bcl/data.dt):numtimes; 
end
set(gca, 'YScale', 'log')
ylabel('$|E_{RD}|$')
%legend(statenames_latex(smallampindices,:),'location','EastOutside')
if legendflag
    l2=legend(statenames_latex(smallampindices,:));
end
xlabel('time, ms')
%axis([5000 6300 10^-16 10^1])
axis([5200 6200 10^-16.5 10^0.5])
yticks([10^-15 10^-10 10^-5 10^0])
set(gcf,'position',[106.3333  120.3333  694.0000  737.3333])
if legendflag
    set(l1,'position',[0.775    0.595    0.2189    0.3376])
    set(l2,'position',[0.83   0.0460    0.1274    0.4802]);
end
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,figname,'-dpdf')
saveas(fig,figname)

% Compute error norms and average absolute error over last 5 cycles
errorendncyc = ncyc; % ending cycle index over which to compute error norm
errorstartindex = round(bcl*(errorendncyc-5)/data.dt)+1; %1; % Enter a starting index for error norm computation
errorendindex = round(bcl*errorendncyc/data.dt)+1;
error_nondim_norm = norm(Smat*trajectorydiff(:,errorstartindex:errorendindex))
error_nondim_mean = mean(mean(abs(Smat*trajectorydiff(:,errorstartindex:errorendindex))))
errorrpd_norm = norm(trajectorydiffrpd(:,errorstartindex:errorendindex))
errorrpd_mean = mean(mean(abs(trajectorydiffrpd(:,errorstartindex:errorendindex))))
