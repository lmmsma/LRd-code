% Jan 2020: Plot fixed points vs. period for selected AP phases. Based off
% plotV_allbcls. 

clear variables;
fs = 24; % fontsize
lw = 2; % linewidth for thresholds
lwV = 1; % linewidth for V plots
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

%folder = ['lrddata_fromfixedpointsearch121718/']; %This folder is the source for fixed points
folder = ['lrddata_fromfixedpointsearch_highres011719/'];%alternative location if using fp_compile2.m to produce high-resolution trajectories
%folder = ['lrddata/'];
%folder = ['lrddata_fromfixedpointsearch_highres_shift7mVrepol/'];

% list all shift folder name extensions
%
%shiftstrings = {'_shift0p4Vnormdepol','_shift0p6Vnormdepol','_shift0p8Vnormdepol','_shift1Vnormdepol','_shift0p2Vnormrepol','_shift0p4Vnormrepol','_shift0p6Vnormrepol','_shift0p8Vnormrepol','_shift0p001Vnormrepol'};
%Vthreshs = [0.4 0.6 0.8 1.0 0.2 0.4 0.6 0.8 0.001]; % Threshold values in normalized membrane potential units.
%shiftaxislabels = {'$0.4 \uparrow$', '$0.6 \uparrow$', '$0.8 \uparrow$', '$1.0 \;\;\;$', '$0.8 \downarrow$', '$0.6 \downarrow$', '$0.4 \downarrow$', '$0.2 \downarrow$', '$0.001 \downarrow$', '$0  \;\;\;$'};
%
%shiftstrings = {'','_shift0p8Vnormrepol'};
shiftstrings = {''};
%shiftstrings = {'_shift0p8Vnormrepol'};
%shiftstrings = {'_shift0p2Vnormrepol'};
%Vthreshs = [0]; % Threshold values in normalized membrane potential units.
%shiftaxislabels = {'$0.4 \uparrow$', '$0.6 \uparrow$', '$0.8 \uparrow$', '$1.0 \;\;\;$', '$0.8 \downarrow$', '$0.6 \downarrow$', '$0.4 \downarrow$', '$0.2 \downarrow$', '$0.001 \downarrow$', '$0  \;\;\;$'};

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
%statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
%statenames_latex = char('$V$, mV','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{IT}$, mmol/L','$[Na^+]_I$, mmol/L','$[K^+]_I$, mmol/L','$[Ca^{2+}]_{JT}$, mmol/L','$[Ca^{2+}]_N, mmol/L$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$, mmol/L/ms');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
statename_symbols = char('v',  'x', 'o', '*', 'v', 's', 'p',    '^',               '>',          '<',       'p',                'h','^','d','+','<','o');

largeampindices = [1 8:12 17];
smallampindices = [2:7 13:16];

param = 'def';


% bcls = [1000:-50:400 390:-10:70];%Vector of bcls in ms. Make sure this matches values in fp_compile2*.m files.
% 
% nbcls = length(bcls);

for shiftctr = 1:length(shiftstrings)
    shiftstring = shiftstrings{shiftctr}
    if isempty(shiftstring) 
        fixedpointfolder = ['lrddata_fromfixedpointsearch_highres011719/'];
    else
        fixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '/']; %alternative location if using fp_compile2.m to produce high-resolution trajectories
    end
    %Vthresh = Vthreshs(shiftctr); % threshold value for normalized voltages
    
    load([fixedpointfolder 'compiled_fp.mat'], 'selected_bcls_for_fps','fp_found')

    h=figure
    subplot(2,1,1)
    hold on;
    for ii=1:length(largeampindices)
        p1(ii) = plot(selected_bcls_for_fps,fp_found(largeampindices(ii),:),statename_symbols(largeampindices(ii),:),'markersize',ms);
    end
    %ylim([-10 150])
    ylabel('Fixed-point value')
    set(gca,'position',[ 0.1300    0.5973    0.67    0.3277])
%    legend(p1,statenames_latex(largeampindices,:))
    lh1=legend(p1,statenames_latex(largeampindices,:)); 
    set(lh1,'position',[0.8240    0.5921    0.1699    0.3376]) 
    if param == 'def'
    %    text(200,80, 'default ','fontsize',fs)
        figname = ['fixedptsdef' shiftstring];
    elseif param == 'adj'
        text(200,80, 'adjusted ','fontsize',fs)
        figname = ['fixedptsadj' shiftstring];
    end

    
    subplot(2,1,2)
    hold on;    
      for ii=1:length(smallampindices)
        p2(ii) = plot(selected_bcls_for_fps,fp_found(smallampindices(ii),:),statename_symbols(smallampindices(ii),:),'markersize',ms);
      end
   % ylim([-.25 1.25])
    ylabel('Fixed-point value')
    set(gca,'position',[ 0.1300    0.1235   0.67    0.3277])

%    legend(statenames_latex(smallampindices,:))
%    legend(p2, statenames_latex(smallampindices,:),'location','EastOutside')
    lh2 = legend(p2, statenames_latex(smallampindices,:)); 
    set(lh2,'position',[0.86    0.0179    0.0989    0.4802]) 
    xlabel('$T$, ms')
    set(gcf,'position',[106.3333  120.3333  894.0000  737.3333])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)

end

