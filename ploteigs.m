%This script is to be run after computeeig.m
%It takes in the calculated eigenvalues and plots either their a)Magnitudes
%vs. BCL or b)their Real vs. Imaginary parts
% Script developed by Anthony and Ryan. Modified by Laura to include
% epsilon values used in Jacobian computation, along with time shifts

clear variables;

% I can't find a way to change LaTeX-interpreted labels back to the default
% font (Helvetica), so I can try changing everything else to the LaTeX
% font, cmr12. Update: this doesn't work either since the Matlab pdf
% printer doesn't interpret the font correctly. Set defaults globally
% instead.
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

logepsln = -5; % This is the log10 of the epsilon value used in Jacobian computation.

%shiftstrings = {'','_shift0p75ms','_shift6p5ms','_shift7mVrepol','_shift-50mVrepol'};
%shiftstrings = {'','_shift0p2Vnormdepol', '_shift0p4Vnormdepol', '_shift0p6Vnormdepol', '_shift0p8Vnormdepol', '_shift1Vnormdepol', '_shift0p8Vnormrepol', '_shift0p6Vnormrepol', '_shift0p4Vnormrepol', '_shift0p2Vnormrepol', '_shift0p001Vnormrepol'};
shiftstrings = {'_shift0p4Vnormdepol', '_shift0p6Vnormdepol', '_shift0p8Vnormdepol', '_shift1Vnormdepol', '_shift0p8Vnormrepol', '_shift0p6Vnormrepol', '_shift0p4Vnormrepol', '_shift0p2Vnormrepol', '_shift0p001Vnormrepol', ''};
%shiftstrings = {''}
%shiftstrings = {'_shift0p001Vnormrepol'}
%shiftstrings = {'_shift0p2Vnormrepol'}
% The following ticks and labels must match the ordering in shiftstrings
%shiftaxisticks = [0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2-0.001]-1;
shiftaxislabels = {'$0.4 \uparrow$', '$0.6 \uparrow$', '$0.8 \uparrow$', '$1.0 \;\;\;$', '$0.8 \downarrow$', '$0.6 \downarrow$', '$0.4 \downarrow$', '$0.2 \downarrow$', '$0.001 \downarrow$', '$0  \;\;\;$'};
%shiftaxisticks = [0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2-0.001 2]-1;
shiftaxisticks = [0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2-0.001 2.1]-1;
%shiftaxislabels = {'0  ', '0.4 \uparrow', '0.6 \uparrow', '0.8 \uparrow', '1.0 \uparrow', '0.8 \downarrow', '0.6 \downarrow', '0.4 \downarrow', '0.2 \downarrow', '0.001 \downarrow'};
xlold = {'70','240','400','600','800','1000'}; % default BCL labels
xtold = [70 240 400 600 800 1000]; % default BCL ticks
% Surf plots require a lot of fine adjustments to labels to get them to
% line up with anything in a reasonable way. If the plot ordering (or
% values) change, offsets below will need to be adjusted. In general, they should
% be close to half-widths of the relevant x or y-axis interval.
surfyaxisoffset = 0.1;
surfyaxisoffset_nearrest = 0.05;
surfxaxisoffset_shortbcl = 5;
surfxaxisoffset_longbcl = 25;

individualfigflag = 0; % nonzero to plot one set of figures per shift. Otherwise, put all values on one plot.
eveccompplotflag = 0; % nonzero to plot eigenvector components vs. BCL and shift index
surfplotflag = 1; % nonzero to plot obsv mags vs. BCL and offset
%approxshifttimes = [0 0.75 6.5 70 100]; % in ms
% paramflag = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
% if paramflag
%     param= 'adj';
% else
%     param = 'def';
% end
param = 'def';
%param = 'adj';

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
statename_symbols = char('v',  'rh', '+', '.', 'rx', '>', 'x',    '<',               's',          'p',       'h',                'yx','kx','b+','g+','r+','^');
numstate = length(statenames);
eigsymbols = char('gs','b+','rx','v','^'); % for plotting eigenvalues in complex plane
%selected_bcls_eig_complex = [1000 600 400 200 70]; % bcls to use in complex plane plot
selected_bcls_eig_complex = [1000 200 70]; % bcls to use in complex plane plot

fontsize = 28; % fontsize
ma = 80; % marker area for scatter plot
mksz = 10; % marker size
lw = 2; % linewidth

% load approximate "state normalization" scaling matrix
load b1000fsolem12variable_amplitudes varamp
% Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
Smat = diag(1./varamp); % scaling matrix
Smatinv = inv(Smat); % only need to compute once

numstate = size(Smat,1);
% input matrix for all possible individual inputs
B = eye(numstate);

% output matrix for all possible individual measurements
C = eye(numstate);

umax = diag(ones(1,size(B,2))); % This is just a placeholder.
% For inputs, also need to know approximate maximum values of each element
% of input vector. Note that u is a deviational quantity that may or
% may not attain the size of the stimulus.
% Each umax diagonal element should probably be the size of
% deviational max for integrated system.
Bs = Smat*B*umax; % scaled B matrix
Cs = C*Smatinv; % scaled C matrix

kc = input('Enter a measurement index from 1 to 17: ');
kb = input('Enter a control input index from 1 to 17: ');

if eveccompplotflag
    hv1alt = figure;
    hv1alt_sc = figure;
    hv1m1 = figure;
    hv1m1_sc = figure;
    hv1m2 = figure;
    hv1m2_sc = figure;
    %hv2m2 = figure;
end

eigfolder = ['eigenvalues' '/' param];
eval(['load ' eigfolder '/eigfile' num2str(logepsln) ' selected_bcls_for_fps'])

nbcls = length(selected_bcls_for_fps);
nshiftstrings = length(shiftstrings);
allmaxalteiginds = zeros(nshiftstrings,nbcls);
allmaxalteiginds_sc = zeros(nshiftstrings,nbcls);
allmaxuniteiginds = zeros(nshiftstrings,nbcls);
allmaxnearuniteiginds = zeros(nshiftstrings,nbcls);
allcdotvmags_valt_sc = zeros(nshiftstrings,nbcls);
allobsvmags_valt_sc = zeros(nshiftstrings,nbcls);
allobsvmags_vm1_sc = zeros(nshiftstrings,nbcls);
allobsvmags_vm2_sc = zeros(nshiftstrings,nbcls);
meanmaxabsaltmode_sc = zeros(1,nshiftstrings);
% Not sure how to initialize a structure array. The following casuses
% errors:
%legendobjs = zeros(nshiftstrings,nbcls);
%legendobjs_sc = zeros(nshiftstrings,nbcls);

bclinds = 1:nbcls;

for shiftctr = 1:nshiftstrings%nshiftstrings%
    %shiftstring = '' % if using default shift of data.dt
    %shiftstring = '_shift0p75ms'
    %shiftstring = '_shift6p5ms'
    %shiftstring = '_shift7mVrepol'
    %shiftstring = '_shift-50mVrepol'
    shiftstring = shiftstrings{shiftctr}
    
    %eigfolder = ['eigenvalues/' param];
    eigfolder = ['eigenvalues' shiftstring '/' param];
    eval(['load ' eigfolder '/eigfile' num2str(logepsln) ' *'])
    % Compute direction of eigenvector of largest eigenvalue with negative real part
    % (alternans eigenvalue).
    % Choose largest alternans eigenvalue, though it's
    % hard to tell which alt eig is which below ~100ms.
    maxabsalteigcomponent = zeros(1,nbcls);
    maxalteigcomponentind = zeros(1,nbcls);
    maxalteig = NaN*ones(1,nbcls);
    maxabsalteigcomponent_sc = zeros(1,nbcls);
    maxalteigcomponentind_sc = zeros(1,nbcls);
    maxalteig_sc = NaN*ones(1,nbcls);
    for ii=1:nbcls
        % find alternans eigenvalues
        [alteigrow, alteigcol] = find(real(alleigs{ii})<0);
        if ~isempty(alteigrow)
            % find index of largest alternans eigenvalue
            [maxalteig(ii), maxalteigind] = max(abs(alleigs{ii}(alteigrow)));
            % store alternans eigenvalue (without taking abs val)
            maxalteig(ii) = alleigs{ii}(alteigrow(maxalteigind));
            % normalized alternans eigenvector (return abs val of each
            % component)
            absvofalteig = abs(allv{ii}(:, alteigrow(maxalteigind)));
            % find largest component of eigenvector
            [maxabsalteigcomponent(ii),maxalteigcomponentind(ii)] = max(absvofalteig);
            % normalized alternans eigenvector for scaled system
            valt_sc = Smat*allv{ii}(:, alteigrow(maxalteigind))/norm(Smat*allv{ii}(:, alteigrow(maxalteigind)));
            % return abs val of each component
            absvofalteig_sc = abs(valt_sc);
            % find largest component of eigenvector
            [maxabsalteigcomponent_sc(ii),maxalteigcomponentind_sc(ii)] = max(absvofalteig_sc);
            % compute observability measure (don't divide by C mag, but do
            % divide by vmag)
            allcdotvmags_valt_sc(shiftctr,ii) = abs(Cs(kc,:)*valt_sc);
            % Nondim obsv measure:
            allobsvmags_valt_sc(shiftctr,ii) = abs(Cs(kc,:)*valt_sc)/norm(Cs(kc,:));
            
        end
    end
    
    max(maxabsalteigcomponent_sc)
    min(maxabsalteigcomponent_sc)
    statenames(maxalteigcomponentind_sc,:)
    
    
    % Plot max component of eigenvalue nearest unit circle
    %    eigtol1 = 1e-7; % how close eigenvalue should be to 1.0
    eigtol1 = 1e-6; % how close eigenvalue should be to 1.0
    % Other mode near unit circle ranges from 0.99035 to 0.99938, approx
    maxabsuniteigcomponent = zeros(1,nbcls);
    maxuniteigcomponentind = zeros(1,nbcls);
    maxabsuniteigcomponent_sc = zeros(1,nbcls);
    maxuniteigcomponentind_sc = zeros(1,nbcls);
    maxabsnearuniteigcomponent = zeros(1,nbcls);
    maxnearuniteigcomponentind = zeros(1,nbcls);
    maxabsnearuniteigcomponent_sc = zeros(1,nbcls);
    maxnearuniteigcomponentind_sc = zeros(1,nbcls);
    sortedabsnearuniteigcomponents = zeros(numstate,nbcls);
    sortedabsnearuniteigcomponentind = zeros(numstate,nbcls);
    sortedabsnearuniteigcomponents_sc = zeros(numstate,nbcls);
    sortedabsnearuniteigcomponentind_sc = zeros(numstate,nbcls);
    unityeigs = zeros(1,nbcls);
    nearunityeigs = zeros(1,nbcls);
    for ii=1:nbcls
        % find index of unity eigenvalue
        [uniteigrow, uniteigcol] = find(alleigsabs{ii} > (1-eigtol1) & alleigsabs{ii} < (1+eigtol1));
        % store unity eigenvalue
        unityeigs(ii) = alleigs{ii}(uniteigrow,uniteigcol);
        % find index of near-unity eigenvalue
        %        [nearuniteigrow, nearuniteigcol] = find(alleigs{ii} > 0.9903 & alleigs{ii} < 0.9994);
        [nearuniteigrow, nearuniteigcol] = find(alleigs{ii} > 0.9903 & alleigs{ii} < 0.9995);
        % store near-unity eigenvalue
        nearunityeigs(ii) = alleigs{ii}(nearuniteigrow,nearuniteigcol);
        % find largest element of eigenvector belonging to unity eigenvalue
        absvofuniteig = abs(allv{ii}(:, uniteigrow));
        % normalized m1 eigenvector for scaled system
        vm1_sc = Smat*allv{ii}(:, uniteigrow)/norm(Smat*allv{ii}(:, uniteigrow));
        absvofuniteig_sc = abs(vm1_sc);
        [maxabsuniteigcomponent(ii),maxuniteigcomponentind(ii)] = max(absvofuniteig);
        [maxabsuniteigcomponent_sc(ii),maxuniteigcomponentind_sc(ii)] = max(absvofuniteig_sc);
        % Nondim obsv measure:
        allobsvmags_vm1_sc(shiftctr,ii) = abs(Cs(kc,:)*vm1_sc)/norm(Cs(kc,:));
        
        % Perform calculations for near-unity mode:
        absvofnearuniteig = abs(allv{ii}(:, nearuniteigrow));
        % normalized m2 eigenvector for scaled system
        vm2_sc = Smat*allv{ii}(:, nearuniteigrow)/norm(Smat*allv{ii}(:, nearuniteigrow));
        absvofnearuniteig_sc = abs(vm2_sc);
        [maxabsnearuniteigcomponent(ii),maxnearuniteigcomponentind(ii)] = max(absvofnearuniteig);
        [maxabsnearuniteigcomponent_sc(ii),maxnearuniteigcomponentind_sc(ii)] = max(absvofnearuniteig_sc);
        [sortedabsnearuniteigcomponents(:,ii), sortedabsnearuniteigcomponentind(:,ii)] = sort(abs(allv{ii}(:, nearuniteigrow)),'descend','ComparisonMethod','abs');
        [sortedabsnearuniteigcomponents_sc(:,ii), sortedabsnearuniteigcomponentind_sc(:,ii)] = sort(abs(vm2_sc),'descend','ComparisonMethod','abs');
        % Nondim obsv measure:
        allobsvmags_vm2_sc(shiftctr,ii) = abs(Cs(kc,:)*vm2_sc)/norm(Cs(kc,:));
    end
    
    max(maxabsuniteigcomponent_sc)
    min(maxabsuniteigcomponent_sc)
    
    shiftstring
    % Compute nearness of unity eig to 1 (assuming it's real-valued)
    max(unityeigs-1)
    min(unityeigs-1)
    
    max(nearunityeigs)
    min(nearunityeigs)
    
    % List state names assoc. with larger components
    statenames(maxuniteigcomponentind_sc,:)
    statenames(maxnearuniteigcomponentind_sc,:)
    %statenames(sortedabsnearuniteigcomponentind(1,:),:)
    statenames(sortedabsnearuniteigcomponentind_sc(2,:),:)
    
    % For compactness, strings to screen in matrix form:
    [statenames(maxalteigcomponentind_sc,:) statenames(maxuniteigcomponentind_sc,:) statenames(maxnearuniteigcomponentind_sc,:) statenames(sortedabsnearuniteigcomponentind_sc(2,:),:)]
    
    % Compute average component size (scaled system) over modal index
    altmodalind = mode(maxalteigcomponentind_sc)
    bclinds_altmode = bclinds(maxalteigcomponentind_sc == altmodalind)
    meanmaxabsaltmode_sc(shiftctr) = mean(maxabsalteigcomponent_sc(bclinds_altmode));
    disp('Press any key to continue. Avg. modal alternans eigenvector component is below.')
    meanmaxabsaltmode_sc(shiftctr)
    %   pause
    
    
    unitmodalind = mode(maxabsuniteigcomponent_sc)
    bclinds_unitmode = bclinds(maxabsuniteigcomponent_sc == unitmodalind)
    meanmaxabsunitmode_sc = mean(maxabsuniteigcomponent_sc(bclinds_unitmode));
    meanmaxabsunitmode_sc
    
    % time = (1:length(Y))*outstep;
    if eveccompplotflag
        figure(hv1alt)
        hold on;
        %    plot(shiftctr,statename_symbols(maxalteigcomponentind,:))
        for i=1:nbcls
            %         if exist('offset','var')
            %           if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
            %               offset_time = time(i_repol(i));
            %           else
            %               offset_time = offset;
            %           end
            %         else
            %             offset_time = 0;
            %         end
            %        legendobjs(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsalteigcomponent(i),statename_symbols(maxalteigcomponentind(i),:));
            %        legendobjs(shiftctr,i)=scatter(selected_bcls_for_fps(i),offset_time,ma,maxabsalteigcomponent(i),statename_symbols(maxalteigcomponentind(i),:));
            %        legendobjs(shiftctr,i)=scatter(selected_bcls_for_fps(i),approxshifttimes(shiftctr),ma,maxabsalteigcomponent(i),statename_symbols(maxalteigcomponentind(i),:));
            legendobjs(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsalteigcomponent(i),statename_symbols(maxalteigcomponentind(i),:));
            legendobjs(shiftctr,i).LineWidth = lw;
            allmaxalteiginds(shiftctr,i) = maxalteigcomponentind(i);
        end
        
        figure(hv1alt_sc)
        hold on;
        for i=1:nbcls
            %        legendobjs_sc(shiftctr,i)=scatter(selected_bcls_for_fps(i),approxshifttimes(shiftctr),ma,maxabsalteigcomponent_sc(i),statename_symbols(maxalteigcomponentind_sc(i),:));
            legendobjs_sc(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsalteigcomponent_sc(i),statename_symbols(maxalteigcomponentind_sc(i),:));
            legendobjs_sc(shiftctr,i).LineWidth = lw;
            allmaxalteiginds_sc(shiftctr,i) = maxalteigcomponentind_sc(i);
        end
        
        
        
        figure(hv1m1)
        hold on;
        for i=1:nbcls
            %        legendobjsm1(shiftctr,i)=scatter(selected_bcls_for_fps(i),approxshifttimes(shiftctr),ma,maxabsuniteigcomponent(i),statename_symbols(maxuniteigcomponentind(i),:));
            legendobjsm1(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsuniteigcomponent(i),statename_symbols(maxuniteigcomponentind(i),:));
            legendobjsm1(shiftctr,i).LineWidth = lw;
            allmaxuniteiginds(shiftctr,i) = maxuniteigcomponentind(i);
        end
        
        figure(hv1m1_sc)
        hold on;
        for i=1:nbcls
            %        legendobjsm1_sc(shiftctr,i)=scatter(selected_bcls_for_fps(i),approxshifttimes(shiftctr),ma,maxabsuniteigcomponent_sc(i),statename_symbols(maxuniteigcomponentind_sc(i),:));
            legendobjsm1_sc(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsuniteigcomponent_sc(i),statename_symbols(maxuniteigcomponentind_sc(i),:));
            legendobjsm1_sc(shiftctr,i).LineWidth = lw;
            allmaxuniteiginds_sc(shiftctr,i) = maxuniteigcomponentind_sc(i);
        end
        
        
        figure(hv1m2)
        hold on;
        for i=1:nbcls
            %        legendobjsm2(shiftctr,i)=scatter(selected_bcls_for_fps(i),approxshifttimes(shiftctr),ma,maxabsnearuniteigcomponent(i),statename_symbols(maxnearuniteigcomponentind(i),:));
            legendobjsm2(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsnearuniteigcomponent(i),statename_symbols(maxnearuniteigcomponentind(i),:));
            legendobjsm2(shiftctr,i).LineWidth = lw;
            allmaxnearuniteiginds(shiftctr,i) = maxnearuniteigcomponentind(i);
        end
        
        figure(hv1m2_sc)
        hold on;
        for i=1:nbcls
            %        legendobjsm2_sc(shiftctr,i)=scatter(selected_bcls_for_fps(i),approxshifttimes(shiftctr),ma,maxabsnearuniteigcomponent_sc(i),statename_symbols(maxnearuniteigcomponentind_sc(i),:));
            legendobjsm2_sc(shiftctr,i)=scatter(selected_bcls_for_fps(i),shiftctr,ma,maxabsnearuniteigcomponent_sc(i),statename_symbols(maxnearuniteigcomponentind_sc(i),:));
            legendobjsm2_sc(shiftctr,i).LineWidth = lw;
            allmaxnearuniteiginds_sc(shiftctr,i) = maxnearuniteigcomponentind_sc(i);
        end
    end
    
    if individualfigflag
        %% Plotting Magnitude of Each Eigenvalue
        figure
        hold on;
        %grid on;
        for i=1:nbcls
            for j=1:numstate
                % plot negative-real-part eigenvalues with a different color and
                % symbol (Elizabeth's suggestion)
                if real(alleigs{i}(j)) >= 0
                    p(1) = scatter(selected_bcls_for_fps(i), alleigsabs{i}(j), ma, 'bo');
                    p(1).LineWidth = lw;
                else
                    p(2) = scatter(selected_bcls_for_fps(i), alleigsabs{i}(j), ma, 'rd');
                    p(2).LineWidth = lw;
                end
            end
        end
        %set(gcf,'DefaultTextFontSize', fs)
        %title(['Eigenvalue moduli for ' param ' parameters, epsilon = 10^{' num2str(logepsln) '}']);
        title(shiftstrings{shiftctr}(2:end))
        ylabel('      $|\lambda|$');
        xlabel('$T$, ms');
        if param == 'def'
            text(750,0.9, 'default ','fontsize',fontsize)
            figname = ['eigsdef' shiftstring];
        elseif param == 'adj'
            text(750,0.9, 'adjusted ','fontsize',fontsize)
            figname = ['eigsadj' shiftstring];
        end
        axis([40 1030 0 1.25])
        leghandle=legend([p(1), p(2)], 'Re$(\lambda) \geq 0$', 'Re$(\lambda) < 0$','Location','east');
        set(leghandle,'Position',[0.6584    0.4090    0.2304    0.1152])
        set(gca,'FontSize',fontsize)
        %set(gcf,'units','normalized','outerposition',[0 0 1 1])
        set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.55    0.79])
        set(gca,'Position', [0.165    0.1618    0.7551    0.7632])
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,figname,'-dpdf')
        saveas(fig,figname)
        
        %% Plot eigenvalues on complex plane
        figure
        hold on;
        clear pe; 
        nbclforeigs = length(selected_bcls_eig_complex); 
        leglabel_e = cell(1,nbclforeigs);
        for i=1:nbclforeigs%nbcls
            sbindex = bclinds(selected_bcls_for_fps == selected_bcls_eig_complex(i));
            leglabel_e{i} = ['$T=$ ' num2str(selected_bcls_for_fps(sbindex)) ' ms']; 
            for j=1:numstate
                pe(i) = plot(real(alleigs{sbindex}(j)),imag(alleigs{sbindex}(j)),eigsymbols(i,:),'MarkerSize',mksz,'LineWidth',lw); 
            end
        end
        title(shiftstrings{shiftctr}(2:end))
        ylabel('Im($\lambda$)');
        xlabel('Re($\lambda$)');
        if param == 'def'
            text(-1.1,0.2, 'default ','fontsize',fontsize)
            figname = ['eigscomplexpldef' shiftstring];
        elseif param == 'adj'
            text(-1.1,0.2, 'adjusted ','fontsize',fontsize)
            figname = ['eigscomplexpladj' shiftstring];
        end
        % Add graph of unit circle: 
        th = 0:pi/50:2*pi;
        xunit = cos(th);
        yunit = sin(th);
        plot(xunit, yunit,'--');
        leg=legend(pe,leglabel_e);
        set(leg,'position', [0.5853    0.6047    0.3168    0.3213])

        axis equal
        axis([-1.5 1.5 -.5 .5])
        set(gca,'FontSize',fontsize)
        set(gcf,'units','normalized','outerposition',[0.1838    0.0014    0.5398    0.4993])
        set(gca,'Position', [0.165    0.1618    0.7551    0.7632])
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,figname,'-dpdf')
        saveas(fig,figname)
        
        
        %% Compute and plot largest component of largest eigenvalue, which happens to be an alternans
        % % eigenvalue in alternans range of BCLs.
        % if param == 'def'
        %     startaltbcl = 240;
        %     endaltbcl = 150;
        % else
        %     startaltbcl = 300;
        %     endaltbcl = 120;
        % end
        % bclind = 1:nbcls;
        % startaltind = bclind(selected_bcls_for_fps == startaltbcl);
        % endaltind = bclind(selected_bcls_for_fps == endaltbcl);
        % altbcls = selected_bcls_for_fps(startaltind:endaltind);
        % maxabsaltcomponent = zeros(1,length(altbcls));
        % maxaltcomponentind = zeros(1,length(altbcls));
        % for ii = startaltind:endaltind
        %     % find index of largest eigenvalue
        %     [maxeigabs, maxeigabsind] = max(alleigsabs{ii});
        %     % find largest element of eigenvector belonging to largest eigenvalue
        %     absvofmaxeig = abs(allv{ii}(:,maxeigabsind));
        %     [maxabsaltcomponent(ii-startaltind+1),maxaltcomponentind(ii-startaltind+1)] = max(absvofmaxeig);
        % end
        % statenames(maxaltcomponentind,:)
        % figure
        % plot(altbcls,maxabsaltcomponent,'x')
        % xlabel('BCL, ms')
        % ylabel('Modulus of largest component of v for |\lambda|_{max}')
        
        
        figure
        hold on;
        plot(selected_bcls_for_fps,unityeigs,'o')
        xlabel('BCL, ms')
        ylabel('$\lambda_{m1}$')
        
        figure
        hold on;
        plot(selected_bcls_for_fps,nearunityeigs,'o')
        xlabel('BCL, ms')
        ylabel('$\lambda_{m2}$')
        
        figure
        plot(selected_bcls_for_fps,maxabsuniteigcomponent_sc,'x')
        xlabel('BCL, ms')
        ylabel('Modulus of largest component of $\bar{v}$ for $|\lambda_{m1}| = 1$')
        
        figure
        hold on;
        plot(selected_bcls_for_fps,maxabsnearuniteigcomponent_sc,'x')
        %plot(selected_bcls_for_fps,sortedabsnearuniteigcomponents(1,:),'go')
        plot(selected_bcls_for_fps,sortedabsnearuniteigcomponents_sc(2,:),'rs')
        xlabel('BCL, ms')
        ylabel('Moduli of largest components of $\bar{v}$ for $|\lambda_{m2}|$')
        
        % Combine plots
        figure
        hold on;
        %p(1) = plot(altbcls(2:end-1),maxabsaltcomponent(2:end-1),'b^','markersize',mksz);
        %         p(1) = plot(selected_bcls_for_fps,maxabsalteigcomponent,'r^','markersize',mksz);
        %         p(2) = plot(selected_bcls_for_fps,maxabsuniteigcomponent,'gs','markersize',mksz);
        %         p(3) = plot(selected_bcls_for_fps,maxabsnearuniteigcomponent,'bv','markersize',mksz);
        %         p(4) = plot(selected_bcls_for_fps,sortedabsnearuniteigcomponents(2,:),'kh','markersize',mksz);
        p(1) = plot(selected_bcls_for_fps,maxabsalteigcomponent_sc,'r^','markersize',mksz);
        p(2) = plot(selected_bcls_for_fps,maxabsuniteigcomponent_sc,'gs','markersize',mksz);
        p(3) = plot(selected_bcls_for_fps,maxabsnearuniteigcomponent_sc,'bv','markersize',mksz);
        p(4) = plot(selected_bcls_for_fps,sortedabsnearuniteigcomponents_sc(2,:),'kh','markersize',mksz);
        leghandle = legend(p,'$\bar{v}_{1,alt}$', '$\bar{v}_{1,m1}$', '$\bar{v}_{1,m2}$', '$\bar{v}_{2,m2}$');
        %leghandle = legend(p,'$v_{1,alt}$', '$v_{1,m1}$', '$v_{1,m2}$');
        p(1).LineWidth = lw;
        p(2).LineWidth = lw;
        p(3).LineWidth = lw;
        p(4).LineWidth = lw;
        xlabel('$T$, ms')
        %ylabel('Moduli of largest components of eigenvectors')
        ylabel('$|\bar{v}_{i,k}|$')
        if param == 'def'
            text(750,0.856, 'default ','fontsize',fontsize)
            figname = ['vcompsdef' shiftstring];
        elseif param == 'adj'
            text(750,0.856, 'adjusted ','fontsize',fontsize)
            figname = ['vcompsadj' shiftstring];
        end
        title(shiftstrings{shiftctr}(2:end))
        axis([40 1030 0 1.05])
        % Add annotations, if all directions are the same for a given component.
        %altcomps = maxaltcomponentind(2:end-1); % Exclude endpoints since alternans eig crosses other values there
        %if max(altcomps) == min(altcomps)
        %        if max(maxalteigcomponentind) == min(maxalteigcomponentind)
        if max(maxalteigcomponentind_sc) == min(maxalteigcomponentind_sc)
            
            arrowxcoords = [0.4425    0.3658];%[220 200]/(1000-70);
            arrowycoords = [0.8918    0.8490];%[1.05-0.65 1-0.65]/(1.05-0.65);
            %annotation('textarrow',arrowxcoords,arrowycoords,'String',statenames_latex(maxalteigcomponentind(1),:),'Interpreter','latex','FontSize',fontsize*.8)
            annotation('textarrow',arrowxcoords,arrowycoords,'String',['$\bar{v}_{1,alt}$ dir. ' statenames_latex(maxalteigcomponentind_sc(1),:)],'Interpreter','latex','FontSize',fontsize*.8)
        end
        %        if max(maxuniteigcomponentind) == min(maxuniteigcomponentind)
        if max(maxuniteigcomponentind_sc) == min(maxuniteigcomponentind_sc)
            arrowxcoords = [0.6776    0.6081];
            arrowycoords = [0.7341    0.7766];
            %            annotation('textarrow',arrowxcoords,arrowycoords,'String',statenames_latex(maxuniteigcomponentind(1),:),'Interpreter','latex','FontSize',fontsize*.8)
            annotation('textarrow',arrowxcoords,arrowycoords,'String',['$\bar{v}_{1,m1}$ dir. ' statenames_latex(maxuniteigcomponentind_sc(1),:)],'Interpreter','latex','FontSize',fontsize*.8)
        end
        %        if max(maxnearuniteigcomponentind) == min(maxnearuniteigcomponentind)
        if max(maxnearuniteigcomponentind_sc) == min(maxnearuniteigcomponentind_sc)
            arrowxcoords = [0.6776    0.6081];
            arrowycoords = [.3307 .2882];
            %            annotation('textarrow',arrowxcoords,arrowycoords,'String',statenames_latex(maxnearuniteigcomponentind(1),:),'Interpreter','latex','FontSize',fontsize*.8)
            annotation('textarrow',arrowxcoords,arrowycoords,'String',['$\bar{v}_{1,m2}$ dir. ' statenames_latex(maxnearuniteigcomponentind_sc(1),:)],'Interpreter','latex','FontSize',fontsize*.8)
        end
        if max(sortedabsnearuniteigcomponentind_sc(2,:)) == min(sortedabsnearuniteigcomponentind_sc(2,:))
            arrowxcoords = [0.4425    0.3658];
            arrowycoords = [0.1911    0.2339];
            annotation('textarrow',arrowxcoords,arrowycoords,'String',['$\bar{v}_{2,m2}$ dir. ' statenames_latex(sortedabsnearuniteigcomponentind_sc(2,1),:)],'Interpreter','latex','FontSize',fontsize*.8)
        end
        
        set(leghandle,'Position',[0.2871    0.3602    0.2304    0.3458])
        set(gca,'FontSize',fontsize)
        set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.55    0.79])
        set(gca,'Position', [0.165    0.1618    0.7551    0.7632])
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,figname,'-dpdf')
        saveas(fig,figname)
        
        close(1:6)
        %% Plotting Real vs Imaginary Part for each Eigenvalue
        figure
        title(['$\lambda$ for ' param ' parameters']);
        ylabel('Im($\lambda$)');
        xlabel('Real($\lambda$)');
        grid on;
        hold on
        %for i=1:length(bcls)
        scatter(real(alleigs{33}), imag(alleigs{33}), 'o');
        %end
    end
    
end
if 0
    % alt eig plots
    [uniquealtinds,uniquealtindinds,tempinds] = unique(allmaxalteiginds);
    numuniqueind = length(uniquealtinds);
    leglabel = cell(1,numuniqueind);
    for jj=1:numuniqueind
        leglabel{jj} = [statenames_latex(uniquealtinds(jj),:) ' dir.'];
    end
    
    % alternans eig plots
    figure(hv1alt)
    colormap winter
    xlabel('BCL, ms')
    %ylabel('Shift index')
    ylabel('Approx. time (ms)')
    c = colorbar;
    c.Label.String = '|v_{1,alt}|=|v_{alt}|_{\infty}';
    leg=legend(legendobjs(uniquealtindinds),leglabel);
    set(leg,'Interpreter','latex')
    if param == 'def'
        text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['v1altdef'];
    elseif param == 'adj'
        text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['v1altadj'];
    end
    %axis([40 1030 0 1.25])
    %set(leg,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    set(gca,'Position', [.165 .1618 .57 .7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    
    
    [uniquealtinds_sc,uniquealtindinds_sc,tempinds_sc] = unique(allmaxalteiginds_sc);
    numuniqueind_sc = length(uniquealtinds_sc);
    leglabel_sc = cell(1,numuniqueind_sc);
    for jj=1:numuniqueind_sc
        leglabel_sc{jj} = [statenames_latex(uniquealtinds_sc(jj),:) ' dir.'];
    end
    
    figure(hv1alt_sc)
    colormap winter
    xlabel('BCL, ms')
    %ylabel('Shift index')
    ylabel('Approx. time (ms)')
    title('Scaled eigenvector components')
    c = colorbar;
    c.Label.String = '|v_{1,alt}|=|v_{alt}|_{\infty}';
    leg_sc=legend(legendobjs_sc(uniquealtindinds_sc),leglabel_sc);
    set(leg_sc,'Interpreter','latex')
    if param == 'def'
        text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['v1altdef_sc'];
    elseif param == 'adj'
        text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['v1altadj_sc'];
    end
    %axis([40 1030 0 1.25])
    %set(leg_sc,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    set(gca,'Position', [.165 .1618 .57 .7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    
    
    % unity eig plots
    [uniqueunitinds,uniqueunitindinds,tempinds] = unique(allmaxuniteiginds);
    numuniqueunitind = length(uniqueunitinds);
    leglabel_unit = cell(1,numuniqueunitind);
    for jj=1:numuniqueunitind
        leglabel_unit{jj} = [statenames_latex(uniqueunitinds(jj),:) ' dir.'];
    end
    
    figure(hv1m1)
    colormap winter
    xlabel('BCL, ms')
    %ylabel('Shift index')
    ylabel('Approx. time (ms)')
    c = colorbar;
    c.Label.String = '|v_{1,m1}|=|v_{m1}|_{\infty}';
    leg=legend(legendobjsm1(uniqueunitindinds),leglabel_unit);
    set(leg,'Interpreter','latex')
    if param == 'def'
        text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['v1m1def'];
    elseif param == 'adj'
        text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['v1m1adj'];
    end
    %axis([40 1030 0 1.25])
    %set(leg,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    set(gca,'Position', [.165 .1618 .57 .7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    
    
    [uniqueunitinds_sc,uniqueunitindinds_sc,tempinds_sc] = unique(allmaxuniteiginds_sc);
    numuniqueunitind_sc = length(uniqueunitinds_sc);
    leglabel_unit_sc = cell(1,numuniqueunitind_sc);
    for jj=1:numuniqueunitind_sc
        leglabel_unit_sc{jj} = [statenames_latex(uniqueunitinds_sc(jj),:) ' dir.'];
    end
    
    figure(hv1m1_sc)
    colormap winter
    xlabel('BCL, ms')
    %ylabel('Shift index')
    ylabel('Approx. time (ms)')
    title('Scaled eigenvector components')
    c = colorbar;
    c.Label.String = '|v_{1,m1}|=|v_{m1}|_{\infty}';
    leg_sc=legend(legendobjsm1_sc(uniqueunitindinds_sc),leglabel_unit_sc);
    set(leg_sc,'Interpreter','latex')
    if param == 'def'
        text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['v1m1def_sc'];
    elseif param == 'adj'
        text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['v1m1adj_sc'];
    end
    %axis([40 1030 0 1.25])
    %set(leg_sc,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    set(gca,'Position', [.165 .1618 .57 .7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    
    
    % near-unity eig plots
    [uniquenearunitinds,uniquenearunitindinds,tempinds] = unique(allmaxnearuniteiginds);
    numuniquenearunitind = length(uniquenearunitinds);
    leglabel_nearunit = cell(1,numuniquenearunitind);
    for jj=1:numuniquenearunitind
        leglabel_nearunit{jj} = [statenames_latex(uniquenearunitinds(jj),:) ' dir.'];
    end
    
    figure(hv1m2)
    colormap winter
    xlabel('BCL, ms')
    %ylabel('Shift index')
    ylabel('Approx. time (ms)')
    c = colorbar;
    c.Label.String = '|v_{1,m2}|=|v_{m2}|_{\infty}';
    leg=legend(legendobjsm2(uniquenearunitindinds),leglabel_nearunit);
    set(leg,'Interpreter','latex')
    if param == 'def'
        text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['v1m2def'];
    elseif param == 'adj'
        text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['v1m2adj'];
    end
    %axis([40 1030 0 1.25])
    %set(leg,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    set(gca,'Position', [.165 .1618 .57 .7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    
    
    [uniquenearunitinds_sc,uniquenearunitindinds_sc,tempinds_sc] = unique(allmaxnearuniteiginds_sc);
    numuniquenearunitind_sc = length(uniquenearunitinds_sc);
    leglabel_nearunit_sc = cell(1,numuniquenearunitind_sc);
    for jj=1:numuniquenearunitind_sc
        leglabel_nearunit_sc{jj} = [statenames_latex(uniquenearunitinds_sc(jj),:) ' dir.'];
    end
    
    figure(hv1m2_sc)
    colormap winter
    xlabel('BCL, ms')
    %ylabel('Shift index')
    ylabel('Approx. time (ms)')
    title('Scaled eigenvector components')
    c = colorbar;
    c.Label.String = '|v_{1,m2}|=|v_{m2}|_{\infty}';
    leg_sc=legend(legendobjsm2_sc(uniquenearunitindinds_sc),leglabel_nearunit_sc);
    set(leg_sc,'Interpreter','latex')
    if param == 'def'
        text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['v1m2def_sc'];
    elseif param == 'adj'
        text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['v1m2adj_sc'];
    end
    %axis([40 1030 0 1.25])
    %set(leg_sc,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    set(gca,'Position', [.165 .1618 .57 .7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
end

mean(mean(allobsvmags_valt_sc))
mean(mean(allobsvmags_vm1_sc))
mean(mean(allobsvmags_vm2_sc))
% approx avg during alternans (bcl range is based more on default params than
% adjusted ones);
mean(mean(allobsvmags_valt_sc(:,26:39)))

if surfplotflag
    %imagesc doesn't work with unevenly spaced data, and pcolor won't show
    %z-values correctly in data indicator, so that leaves surf.
    %hcdotv_sc = figure;
    hobsv_valt_sc = figure;
    % augment matrix since surf won't otherwise show the last row or column
    %allcdotvmags_valt_sc_aug = [allcdotvmags_valt_sc allcdotvmags_valt_sc(:,end); allcdotvmags_valt_sc(end,:) allcdotvmags_valt_sc(end,end)];
    %allcdotvmags_valt_sc_aug = [allcdotvmags_valt_sc allcdotvmags_valt_sc(:,end) allcdotvmags_valt_sc(:,end); allcdotvmags_valt_sc(end,:) allcdotvmags_valt_sc(end,end) allcdotvmags_valt_sc(end,end); allcdotvmags_valt_sc(end,:) allcdotvmags_valt_sc(end,end) allcdotvmags_valt_sc(end,end)];
    %allobsvmags_valt_sc_aug = [allobsvmags_valt_sc allobsvmags_valt_sc(:,end) allobsvmags_valt_sc(:,end); allobsvmags_valt_sc(end,:) allobsvmags_valt_sc(end,end) allobsvmags_valt_sc(end,end); allobsvmags_valt_sc(end,:) allobsvmags_valt_sc(end,end) allobsvmags_valt_sc(end,end)];
    allobsvmags_valt_sc_aug = [allobsvmags_valt_sc allobsvmags_valt_sc(:,end); allobsvmags_valt_sc(end,:) allobsvmags_valt_sc(end,end)];
    %su = surf([selected_bcls_for_fps selected_bcls_for_fps(end) selected_bcls_for_fps(end)-10],[1:nshiftstrings nshiftstrings nshiftstrings+1],allobsvmags_valt_sc_aug);
    su = surf([selected_bcls_for_fps selected_bcls_for_fps(end)-10],[shiftaxisticks shiftaxisticks(end)+surfyaxisoffset],allobsvmags_valt_sc_aug);
    view([0 90])
    su.EdgeColor = 'none';
    tempstr = strtrim(statenames_latex(kc,:));
    c = colorbar;
    caxis([0 0.52])
    %c.Label.String = '|Cv|/||v||';
    %c.Label.String = '|Cv|/||C||||v||';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    %c.Label.String = '$|\cos \phi_{k,alt}|$';
    c.Label.String = ['$|\cos \phi_{\,' tempstr(2:end-1) ',\; alt}|$'];
    %colormap winter
    %xlabel('BCL, ms')
    xlabel('$T$, ms')
    % The following modifications to x-axis labels are hardcoded and will not work properly if the list of BCLs
    % changes:
    %xlold = xticklabels;
    %xlold{1} = '70';
    %xlold{2} = '240';
    %xtold = xticks;
    %xtold(1) = 70;
    %xtold(2) = 240;
    %xticks([selected_bcls_for_fps(1:13)+50 selylabel('Normalized V')ected_bcls_for_fps(14:end)+5 selected_bcls_for_fps(end)-5]);
    xticks([xtold(1:3)-surfxaxisoffset_shortbcl xtold(4:end)-surfxaxisoffset_longbcl]);
    xticklabels(xlold);
    %ylabel('Shift Index')
    %ylabel('Normalized V')
    ylabel('$\breve{V}$')
    yticks([shiftaxisticks(1:8)+surfyaxisoffset shiftaxisticks(9:10)+surfyaxisoffset_nearrest]);
    yticklabels(shiftaxislabels);
    %title(['|Cv|/||v||, scaled system, meas. index = ' num2str(kc)])
    %title(['|Cv|/||C||||v||, scaled system, meas. index = ' num2str(kc)])
    %title(['$|\cos \phi_{' tempstr(2:end-1) ',alt}|$, scaled system'])
    if param == 'def'
        %    text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['obsvaltdef_sc_measind_' num2str(kc)];
    elseif param == 'adj'
        %    text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['obsvaltadj_sc_measind_' num2str(kc)];
    end
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    %set(gca,'Position', [.165 .1618 .57 .7632])
    set(gca,'Position', [.18 .1618 .58 .7632])
    axis tight
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %print(fig,figname,'-dpdf') % surf prints white lines to pdf but not jpg or png
    print(fig,figname,'-djpeg')
    saveas(fig,figname)
    
    
    % m1 surf plot
    hobsv_m1_sc = figure;
    % augment matrix since surf won't otherwise show the last row or column
    allobsvmags_vm1_sc_aug = [allobsvmags_vm1_sc allobsvmags_vm1_sc(:,end); allobsvmags_vm1_sc(end,:) allobsvmags_vm1_sc(end,end)];
    su = surf([selected_bcls_for_fps selected_bcls_for_fps(end)-10],[shiftaxisticks shiftaxisticks(end)+surfyaxisoffset],allobsvmags_vm1_sc_aug);
    view([0 90])
    su.EdgeColor = 'none';
    c = colorbar;
    %caxis([0 12.5e-4]) % for meas = V
    caxis([.9999 .99997]) % for meas = K
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    %c.Label.String = '$|\cos \phi_{k,m1}|$';
    c.Label.String = ['$|\cos \phi_{\,' tempstr(2:end-1) ',\; m1}|$'];
    %colormap winter
    %xlabel('BCL, ms')
    xlabel('$T$, ms')
    % The following modifications to x-axis labels are hardcoded and will not work properly if the list of BCLs
    % changes:
    %xlold = xticklabels;
    %xlold{1} = '70';
    %xlold{2} = '240';
    %xtold = xticks;
    %xtold(1) = 70;
    %xtold(2) = 240;
    %xticks([selected_bcls_for_fps(1:13)+50 selylabel('Normalized V')ected_bcls_for_fps(14:end)+5 selected_bcls_for_fps(end)-5]);
    xticks([xtold(1:3)-surfxaxisoffset_shortbcl xtold(4:end)-surfxaxisoffset_longbcl]);
    xticklabels(xlold);
    %ylabel('Shift Index')
    %ylabel('Normalized V')
    ylabel('$\breve{V}$')
    yticks([shiftaxisticks(1:8)+surfyaxisoffset shiftaxisticks(9:10)+surfyaxisoffset_nearrest]);
    yticklabels(shiftaxislabels);
    %title(['$|\cos \phi_{' tempstr(2:end-1) ',m1}|$, scaled system'])
    if param == 'def'
        %    text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['obsvm1def_sc_measind_' num2str(kc)];
    elseif param == 'adj'
        %    text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['obsvm1adj_sc_measind_' num2str(kc)];
    end
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    %set(gca,'Position', [.165 .1618 .57 .7632])
    set(gca,'Position', [.18 .1618 .58 .7632])
    axis tight
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %print(fig,figname,'-dpdf') % surf prints white lines to pdf but not jpg or png
    print(fig,figname,'-djpeg')
    saveas(fig,figname)
    
    % m2 surf plot
    hobsv_m2_sc = figure;
    % augment matrix since surf won't otherwise show the last row or column
    allobsvmags_vm2_sc_aug = [allobsvmags_vm2_sc allobsvmags_vm2_sc(:,end); allobsvmags_vm2_sc(end,:) allobsvmags_vm2_sc(end,end)];
    su = surf([selected_bcls_for_fps selected_bcls_for_fps(end)-10],[shiftaxisticks shiftaxisticks(end)+surfyaxisoffset],allobsvmags_vm2_sc_aug);
    view([0 90])
    su.EdgeColor = 'none';
    c = colorbar;
    %caxis([0 12.5e-4]) % for meas = V
    caxis([.9771 .9777]) % for meas = K
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    %c.Label.String = '$|\cos \phi_{k,m2}|$';
    c.Label.String = ['$|\cos \phi_{\,' tempstr(2:end-1) ',\; m2}|$'];
    %colormap winter
    %xlabel('BCL, ms')
    xlabel('$T$, ms')
    % The following modifications to x-axis labels are hardcoded and will not work properly if the list of BCLs
    % changes:
    %xlold = xticklabels;
    %xlold{1} = '70';
    %xlold{2} = '240';
    %xtold = xticks;
    %xtold(1) = 70;
    %xtold(2) = 240;
    %xticks([selected_bcls_for_fps(1:13)+50 selylabel('Normalized V')ected_bcls_for_fps(14:end)+5 selected_bcls_for_fps(end)-5]);
    xticks([xtold(1:3)-surfxaxisoffset_shortbcl xtold(4:end)-surfxaxisoffset_longbcl]);
    xticklabels(xlold);
    %ylabel('Shift Index')
    %ylabel('Normalized V')
    ylabel('$\breve{V}$')
    yticks([shiftaxisticks(1:8)+surfyaxisoffset shiftaxisticks(9:10)+surfyaxisoffset_nearrest]);
    yticklabels(shiftaxislabels);
    %title(['$|\cos \phi_{' tempstr(2:end-1) ',m2}|$, scaled system'])
    if param == 'def'
        %    text(750,0.9, 'default ','fontsize',fontsize)
        figname = ['obsvm2def_sc_measind_' num2str(kc)];
    elseif param == 'adj'
        %    text(750,0.9, 'adjusted ','fontsize',fontsize)
        figname = ['obsvm2adj_sc_measind_' num2str(kc)];
    end
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.75    0.79])
    %set(gca,'Position', [.165 .1618 .57 .7632])
    set(gca,'Position', [.18 .1618 .58 .7632])
    axis tight
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %print(fig,figname,'-dpdf') % surf prints white lines to pdf but not jpg or png
    print(fig,figname,'-djpeg')
    saveas(fig,figname)
end