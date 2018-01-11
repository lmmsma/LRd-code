% Plot high-resolution period-1 trajectories 

clear variables;

linewidth = 1; 
statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');

datasourcefolder = 'lrddata_fromfixedpointsearch_highres/';
% This folder is the source for simulations
% where ICs were chosen as fixed points

bcls = 1000;
plotvarindices = [8]; % Choose state variable indicies of variables to plot 

simoutfname = ['lrddata_1cell_b' num2str(bcls(1))];
eval(['load ' datasourcefolder simoutfname ' Y data subdiv_per_cyc;']);
time = (1:size(Y,2))*bcls(1)/subdiv_per_cyc;

figure
plot(time/1000,Y(plotvarindices,:),'Linewidth',linewidth)
title([statenames_latex(plotvarindices,:) ' vs. time'],'Interpreter','latex')
xlabel('time (sec)')
ylabel('mmol/L') 
axis([0 1 -0.07 0.07])