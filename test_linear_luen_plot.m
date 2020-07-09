% Plot results of test_linear_luen.m.

clear variables

combineplots = 1; % set to 0 for two separate V and K plots

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

lu_folder = 'linear_luen_tests/';

pertsize = 10^-5;
bcllist = [100 200 500 1000];
numsimsteps = 30;

measurementindex = 1;

if combineplots
    % Plot V-meas cases on one plot
    shiftstring = '_shift1Vnormdepol'
    bcl = 100;
    
    hV=figure;
    hold on;
    
    pctr = 1;
    
    loadfilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
    eval(['load ' lu_folder loadfilename ' *']);
    ph(pctr) = plot(1:numsimsteps, mean(abs(err_scaled)));
    legendstr = ['$T = \;\;$ ' num2str(bcl) ' ms, $\breve{V} = 1.0 \;\;\;$, meas. ' deblank(statenames_latex(measurementindex,:))];
    
    % errorendindex = numsimsteps; % ending cycle index over which to compute error norm
    % errorstartindex = errorendindex-4; %1; % Enter a starting index for error norm computation
    % error_lin_nondim_mean = mean(mean(abs(err_scaled(:,errorstartindex:errorendindex))))
    
    % Produce V-meas plots for other shift value
    
    shiftstring = '_shift0p2Vnormrepol'
    
    for bcl = bcllist
        pctr = pctr + 1;
        loadfilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
        eval(['load ' lu_folder loadfilename ' *']);
        ph(pctr) = plot(1:numsimsteps, mean(abs(err_scaled)));
        if bcl < 1000
            legendstr = char(legendstr, ['$T = \;\;$  ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$, meas. ' deblank(statenames_latex(measurementindex,:))]);
        else
            legendstr = char(legendstr, ['$T =$ ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$, meas. ' deblank(statenames_latex(measurementindex,:))]);
        end
    end
    figure(hV)
    %set(gca, 'YScale', 'log')
    %ylabel('$\arrowvert \bar{e}_j \arrowvert_{i-av}$')
    %xlabel('$j$, cycle index')
    %ylim([10^-8 10^1.5])
    %legend(strjust(legendstr))  %strjust (default right justify) doesn't seem
    %to work the way I want, whitespaces aren't preserved as expected
    %legend(legendstr)
    %title(['Measured variable: ' deblank(statenames_latex(measurementindex, :)) ])
    %title(['$T = $' num2str(bcl) 'ms, Meas = ' statenames_latex(measurementindex, :) ', ' shiftstring(7:end) ', ' param],'Interpreter','Latex');
    % set(gcf,'position',[106.3333  120.3333  946.0000  737.3333])
    % % if legendflag
    % %     set(l1,'position',[0.775    0.595    0.2189    0.3376])
    % %     set(l2,'position',[0.83   0.0460    0.1274    0.4802]);
    % % end
    % if param == 'def'
    % %    text(200,80, 'default ','fontsize',fs)
    %     figname = ['lin_est_err_def_m' num2str(measurementindex)];
    % elseif param == 'adj'
    % %    text(200,80, 'adjusted ','fontsize',fs)
    %     figname = ['lin_est_err_adj_m' num2str(measurementindex)];
    % end
    % fig = gcf;
    % fig.PaperPositionMode = 'auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    % % Adjust legend manually before saving
    % print(fig,figname,'-dpdf')
    % saveas(fig,figname)
    co=get(gca,'ColorOrder');
    
    
    
    % Add meas=K plots
    
    measurementindex = 10;
    pctr = 1;
    
    for bcl = bcllist
        pctr = pctr + 1;
        loadfilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
        eval(['load ' lu_folder loadfilename ' *']);
        ph(pctr) = plot(1:numsimsteps, mean(abs(err_scaled)), '--', 'Color',co(pctr,:));
        if bcl < 1000
            legendstr = char(legendstr, ['$T = \;\;$  ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$, meas. ' deblank(statenames_latex(measurementindex,:))]);
        else
            legendstr = char(legendstr, ['$T =$ ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$, meas. ' deblank(statenames_latex(measurementindex,:))]);
        end
    end
    figure(hV)
    set(gca, 'YScale', 'log')
    ylabel('$\arrowvert \bar{e}_j \arrowvert_{i-av}$')
    xlabel('$j$, cycle index')
    ylim([10^-8 10^1.5])
    %legend(legendstr,'location','NorthEastOutside')
    %legend(legendstr,'location','EastOutside')
    lh=legend(legendstr);
    set(lh,'box','off')
    set(gca,'position',[0.100    0.1100    0.5    0.8150])
    set(gcf,'position',1000*[0.0950    0.0010    1.3460    0.8873])
    if param == 'def'
        figname = ['lin_est_err_def_m1_and_10'];
    elseif param == 'adj'
        figname = ['lin_est_err_adj_m1_and_10'];
    end
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    % Adjust legend manually before saving
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    
else % two separate plots
    % Plot V-meas cases on one plot
    shiftstring = '_shift1Vnormdepol'
    bcl = 100;
    
    hV=figure;
    hold on;
    
    pctr = 1;
    
    loadfilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
    eval(['load ' lu_folder loadfilename ' *']);
    ph(pctr) = plot(1:numsimsteps, mean(abs(err_scaled)));
    legendstr = ['$T = \;\;$ ' num2str(bcl) ' ms, $\breve{V} = 1.0 \;\;\;$'];
    
    % errorendindex = numsimsteps; % ending cycle index over which to compute error norm
    % errorstartindex = errorendindex-4; %1; % Enter a starting index for error norm computation
    % error_lin_nondim_mean = mean(mean(abs(err_scaled(:,errorstartindex:errorendindex))))
    
    % Produce V-meas plots for other shift value
    
    shiftstring = '_shift0p2Vnormrepol'
    
    for bcl = bcllist
        pctr = pctr + 1;
        loadfilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
        eval(['load ' lu_folder loadfilename ' *']);
        ph(pctr) = plot(1:numsimsteps, mean(abs(err_scaled)));
        if bcl < 1000
            legendstr = char(legendstr, ['$T = \;\;$  ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$']);
        else
            legendstr = char(legendstr, ['$T =$ ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$']);
        end
    end
    figure(hV)
    set(gca, 'YScale', 'log')
    ylabel('$\arrowvert \bar{e}_j \arrowvert_{i-av}$')
    xlabel('$j$, cycle index')
    ylim([10^-8 10^1.5])
    %legend(strjust(legendstr))  %strjust (default right justify) doesn't seem
    %to work the way I want, whitespaces aren't preserved as expected
    legend(legendstr)
    %title(['Measured variable: ' deblank(statenames_latex(measurementindex, :)) ])
    %title(['$T = $' num2str(bcl) 'ms, Meas = ' statenames_latex(measurementindex, :) ', ' shiftstring(7:end) ', ' param],'Interpreter','Latex');
    set(gcf,'position',[106.3333  120.3333  946.0000  737.3333])
    % if legendflag
    %     set(l1,'position',[0.775    0.595    0.2189    0.3376])
    %     set(l2,'position',[0.83   0.0460    0.1274    0.4802]);
    % end
    if param == 'def'
        %    text(200,80, 'default ','fontsize',fs)
        figname = ['lin_est_err_def_m' num2str(measurementindex)];
    elseif param == 'adj'
        %    text(200,80, 'adjusted ','fontsize',fs)
        figname = ['lin_est_err_adj_m' num2str(measurementindex)];
    end
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    % Adjust legend manually before saving
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    co=get(gca,'ColorOrder');
    
    
    
    % Produce separate plot for meas=K
    hK=figure;
    hold on;
    set(gca,'ColorOrder',co(2:end,:));
    
    measurementindex = 10;
    pctr = 1;
    
    for bcl = bcllist
        pctr = pctr + 1;
        loadfilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
        eval(['load ' lu_folder loadfilename ' *']);
        ph(pctr) = plot(1:numsimsteps, mean(abs(err_scaled)));
        if bcl == 100
            legendstr2 = char(['$T = \;\;$  ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$']);
        elseif bcl < 1000
            legendstr2 = char(legendstr2, ['$T = \;\;$  ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$']);
        else
            legendstr2 = char(legendstr2, ['$T =$ ' num2str(bcl) ' ms, $\breve{V} = 0.2 \downarrow$']);
        end
    end
    figure(hK)
    set(gca, 'YScale', 'log')
    ylabel('$\arrowvert \bar{e}_j \arrowvert_{i-av}$')
    xlabel('$j$, cycle index')
    ylim([10^-8 10^1.5])
    legend(legendstr2)
    set(gcf,'position',[106.3333  120.3333  946.0000  737.3333])
    if param == 'def'
        figname = ['lin_est_err_def_m' num2str(measurementindex)];
    elseif param == 'adj'
        figname = ['lin_est_err_adj_m' num2str(measurementindex)];
    end
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    % Adjust legend manually before saving
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
end