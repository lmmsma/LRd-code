% Use this program to produce and plot several cycles of a high-resolution
% trajectory. Based on fp_compile2_shiftedtraj.m.
%
% Revised version (after first round of reviews): overlay examples of fixed
% points on trajectories. Note! Skip to section 2 to just produce plots
% instead of compiling fixed points. I don't think the code has been tested 
% (or will necessarily work) if shiftstring is other than ''. 

clear variables;
tstartfpsearch = tic;

%offset = 6.5; % offset of stimulus time, in ms
%offset = 0.75; % offset of stimulus time, in ms
%offset = 7; % this a voltage threshold in mV, not an offset in ms
%offset = -50; % this a voltage threshold in mV, not an offset in ms
shiftstring = ''; % if using default shift of data.dt
%shiftstring = '_shift0p75ms';
%shiftstring = '_shift-50mVrepol';
%shiftstring = '_shift1Vnormdepol';
%shiftstring = '_shift0p2Vnormdepol';
%shiftstring = '_shift0p4Vnormdepol';
%shiftstring = '_shift0p6Vnormdepol';
%shiftstring = '_shift0p8Vnormdepol';
%shiftstring = '_shift0p8Vnormrepol';
%shiftstring = '_shift0p6Vnormrepol';
%shiftstring = '_shift0p4Vnormrepol';
%shiftstring = '_shift0p2Vnormrepol';
%shiftstring = '_shift0p001Vnormrepol';

%Reads in files of located fixed points
%sourcefolder = 'lrddata_fromfixedpointsearch/';%This folder is the source for fixed points
sourcefolder = 'lrddata_fromfixedpointsearch_highres011719/'; %alternative location if using fp_compile2.m to produce high-resolution trajectories
%sourcefolder = 'C:/Users/laura/''Google Drive''/''REU 2017''/Data/''Fixed
%Points''/''Default Parameters- Ryan''/';%This folder is the source for
%fixed points, but these are wrong (swapped parameter set) and the double
%quote format no longer works in Matlab 2018
%fixedpointfolder = 'fixedpoints/';%Save compiled fixed points here
%fixedpointfolder = 'lrddata_fromfixedpointsearch_highres011719/'; %alternative location if using fp_compile2.m to produce high-resolution trajectories
%fixedpointfolder = 'lrddata_fromfixedpointsearch_highres_shift6p5ms/'; %alternative location if using fp_compile2.m to produce high-resolution trajectories
%fixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '/']; %store data here
fixedpointfolder = ['lrddata_fromfixedpointsearch_plothighres' shiftstring '/']; %store data here
%fixedpointfolder = 'fixedpoints_fromaandr/';%Save compiled fixed points here
modelinputfolder = 'lrddata/'; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script.
% Next part is no longer needed if threshold times are stored in a
% different file.
% if ~exist(fixedpointfolder,'dir')
%     mkdir(fixedpointfolder)
% end

%selected_bcls_for_fps = [1000:-50:500 490:-10:70];%Vector of bcls
% Change to match available FPs for adjusted data set:
selected_bcls_for_fps = [1000:-50:400 390:-10:70];%Vector of bcls for Anthony's adjusted parameters
% There are more fixed points than the ones listed in selected_bcls...
% Just choose BCLs where Jacobians are desired. They should be available
% at least every 10 ms, starting at 810 ms.
% The original list for A&R's fixed points is below:
%bcls = [1000:-50:300 290:-10:70 69:-1:50]; %revised full pacedown (more points leading up to 200ms)

L = length(selected_bcls_for_fps);
modelname = 'lrd';
nsolitol = 1e-12; % abs & rel tolerance for Newton-Krylov solver
tolfactor = 4; % for offset trajectories, accept error within tolfactor*nsolitol
% The actual stop-criterion tolerance can be larger than nsolitol, so I'm
% trying to account for that here.

ncyc = 4; % Run the model for this number of cycles per bcl setting.

%Compile all fixed points into a single matrix, fp_found
fp_found = zeros(17, L);
fp_errnorms = zeros(1,L);
%i_repol = zeros(1,L);
% Load repol and depol times:
%load([fixedpointfolder 'threshold_times.mat'],'i_depol','i_repol');
stimstartsteps = 1;

%In case program crashed partway through, check for previously saved fixed
%points:
if exist([fixedpointfolder 'compiled_fp.mat'],'file')
    eval(['load ' fixedpointfolder 'compiled_fp.mat fp_found']);
end

startindex = 35; % Choose BCL = 180ms for plot
% Find index of bcl where program crashed
%failedbcl = 130; %ms
if exist('failedbcl','var')
    bclindices = 1:L;
    startindex = bclindices(selected_bcls_for_fps == failedbcl);
end

for i = startindex%:L
    bcl = selected_bcls_for_fps(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    
    fname = [sourcefolder 'lrddata_1cell_b' num2str(bcl)];
    %  fp = load(fname, 'Y'); % I can't get this to work properly with folder
    %  names that contain spaces
    %   fp_found(:, i) = fp.Y(:, 1);
    %    eval(['load ' fname ' data Y;']);
    eval(['load ' fname ' data Y outstep;']);
    
    %    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    subdiv_per_cyc = bcl/data.dt; % sampling subdivision; high-res option
    
    % All of the following if/else statement assumes that outstep=data.dt:
    %     if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
    %         % borrow falling-edge detection from plot_restitution.m:
    %         V = Y(1,:); % membrane potential
    %         % Search for falling edge crossing of threshold
    %         i_repol(i) = find(V(1:end-1)>offset & V(2:end)<=offset);
    %         % Stimulus start time:
    %         % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
    %         % (bcl + data.dt) – (offset-data.dt). I verified that this formula
    %         % produces matched values for BCL = 1000ms.
    %         stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
    %         % Convert stimstart to time index
    %         stimstartsteps = i_repol(i);
    %     else % otherwise the offset is in ms (right now there is no condition for a rising-edge event detection)
    %         % Stimulus start time:
    %         % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
    %         % (bcl + data.dt) – (offset-data.dt). I verified that this formula
    %         % produces matched values for BCL = 1000ms.
    %         stimstart = bcl + 2*data.dt - offset;
    %         % Convert stimstart to time index
    %         stimstartsteps = round(offset/outstep);
    %     end
    stimstart = stimstartsteps*data.dt;
    
    % All of the following if/else statements assume that outstep=data.dt:
    if strfind(shiftstring,'repol') % in this case the offset is in normalized membrane potential units on the repolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
        stimstartsteps = i_repol(i);
    elseif strfind(shiftstring,'depol') % in this case the offset is in normalized membrane potential units on the depolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_depol(i)*data.dt;
        stimstartsteps = i_depol(i);
    end
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since nsoli expects our function to only
    % have one input, which is the initial condition.
    %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
    
    %    fp_found(:, i) = Y(:, 1);
    fp_found(:, i) = Y(:, stimstartsteps);
    
    %     % Compute fixed-point error:
    eval(['fperr = ' modelname '_p2pzero(fp_found(:,i));']); % compute fixed-point error
    fperrnorm = norm(fperr)
    % from nsoli.m, the stop criterion formula is stop_tol = atol + rtol*fnrm, computed
    % once based on the inital error. Hence, it doesn't really make sense to
    % put in an extra check on the error size, though it can be recorded.
    
    %     if fperrnorm > tolfactor*nsolitol
    % %        tstartfp=tic;
    %         % Find fixed point:
    %         %        [fixedpt, it_hist, ierr, x_hist] = nsoli(fp_found(:, i),[modelname '_p2pzero'],[nsolitol nsolitol])
    %         [fixedpt, it_hist, ierr, x_hist] = nsoli(fp_found(:, i),[modelname '_p2pzero'],tolfactor*[nsolitol nsolitol])
    % %        toc(tstartfp);
    %
    %         if ierr == 0
    %             fp_found(:, i) = fixedpt; % fixed point
    %             eval(['fperr = ' modelname '_p2pzero(fixedpt);']); % compute fixed-point error
    %             fperrnorm = norm(fperr)
    %             %Save the matrix to a separate file in the folder specified above
    %             eval(['save ' fixedpointfolder 'compiled_fp.mat *']);
    %         else
    %             % Exit the "for" loop, and print an error message, if ierr ~= 0
    %             disp(['Fixed-point search failed, ierr = ' num2str(ierr)])
    %             save fperrdumpfile *
    %             sendmail('lmmsma@rit.edu',['Fixed point search error, BCL = ' num2str(bcl)]);
    %             break % exit the for loop
    %         end
    %
    %
    %     end
    %     fp_errnorms(i) = fperrnorm;
    
end

if ~exist(fixedpointfolder,'dir')
    mkdir(fixedpointfolder)
end
%Save the matrix to a separate file in the folder specified above
eval(['save ' fixedpointfolder 'compiled_fp.mat *']);
%Move generated simulation data to fixed point folder
eval(['movefile ' modelinputfolder '* ' fixedpointfolder]);

% % Display minimum and maximum errors:
% max(fp_errnorms)
% min(fp_errnorms)

%Display full matrix if desired (uncomment)
%disp(fp_found);
toc(tstartfpsearch)
%%
clear variables; 
% Adapted from plot_periodic_trajectories.m. This plotting section can be
% run independently from the foregoing simulation section.
fs = 20; % fontsize
lw = 2; % linewidth
markersize = 10; 
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

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
%shiftstring = '';
%shiftstring = '_shift1Vnormdepol'; 
%shiftstring = '_shift0p2Vnormdepol'; shiftlabel = '0.2 \uparrow'; 
%shiftstring = '_shift0p4Vnormdepol'; shiftlabel = '0.4 \uparrow$'; 
%shiftstring = '_shift0p6Vnormdepol'; shiftlabel = '0.6 \uparrow$'; 
%shiftstring = '_shift0p8Vnormdepol'; shiftlabel = '0.8 \uparrow$'; 
%shiftstring = '_shift0p8Vnormrepol'; shiftlabel = '0.8 \downarrow'; 
%shiftstring = '_shift0p6Vnormrepol'; shiftlabel = '0.6 \downarrow'; 
%shiftstring = '_shift0p4Vnormrepol'; shiftlabel = '0.4 \downarrow';
shiftstring = '_shift0p2Vnormrepol'; shiftlabel = '0.2 \downarrow';
%shiftstring = '_shift0p001Vnormrepol'; shiftlabel = '0.001 \downarrow';

% Always load the unshifted trajectories first: 
fixedpointfolder = ['lrddata_fromfixedpointsearch_plothighres/']; %store data here
if ~isempty(shiftstring)
    shiftedfixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '/']; %store data here
end

% This folder is the source for simulations
% where ICs were chosen as fixed points

bcls = 180;
defaultstimstartsteps = 1;

plotvarindices = [1 8:10]; % Choose state variable indicies of variables to plot

simoutfname = ['lrddata_1cell_b' num2str(bcls(1))];
eval(['load ' fixedpointfolder simoutfname ' Y data subdiv_per_cyc ncyc outstep;']);
% Load shift info (i_repol or i_depol) 
if ~isempty(shiftstring)
    eval(['load ' shiftedfixedpointfolder 'compiled_fp.mat i_* selected_bcls*;']);
end
bclinds = 1:length(selected_bcls_for_fps); 
selectedbclind = bclinds(selected_bcls_for_fps == bcls);

time = (1:size(Y,2))*bcls(1)/subdiv_per_cyc;
stimindexoffsets = [0:(ncyc-1)]*bcls/outstep; 
shiftindex = defaultstimstartsteps; 

% Compute stimulus time (tau)
if ~isempty(shiftstring) % copied from above
    if strfind(shiftstring,'repol') % in this case the offset is in normalized membrane potential units on the repolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcls + 2*data.dt - i_repol(selectedbclind)*data.dt;
    elseif strfind(shiftstring,'depol') % in this case the offset is in normalized membrane potential units on the depolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcls + 2*data.dt - i_depol(selectedbclind)*data.dt;
    end
end

for ii=1:length(plotvarindices)
    h(ii) = figure;
    hold on; 
    plot(time,Y(plotvarindices(ii),:))
    p(1)=plot(time(defaultstimstartsteps+stimindexoffsets),Y(plotvarindices(ii),defaultstimstartsteps+stimindexoffsets),'ro','MarkerSize',markersize); 
    if ~isempty(shiftstring)
        if exist('i_depol','var') && i_depol(selectedbclind)
           p(2) = plot(time(i_depol(selectedbclind)+stimindexoffsets),Y(plotvarindices(ii),i_depol(selectedbclind)+stimindexoffsets),'ks','MarkerSize',markersize);
           shiftindex = i_depol(selectedbclind); 
        elseif exist('i_repol','var') && i_repol(selectedbclind)
           p(2) = plot(time(i_repol(selectedbclind)+stimindexoffsets),Y(plotvarindices(ii),i_repol(selectedbclind)+stimindexoffsets),'ks','MarkerSize',markersize); 
           shiftindex = i_repol(selectedbclind); 
        end
    end
    t1 = time(defaultstimstartsteps+stimindexoffsets(3)); 
    t2 = time(defaultstimstartsteps+stimindexoffsets(4));
    Y2 = Y(plotvarindices(ii),defaultstimstartsteps+stimindexoffsets(4));

    xlim([-20 740])
    if ~isempty(shiftstring) && ii==1
        t3 = time(shiftindex+stimindexoffsets(3));
        t4 = time(shiftindex+stimindexoffsets(4));
        Y3 = Y(plotvarindices(ii),shiftindex+stimindexoffsets(3));
        % For 0.8 dn
%         line([t3 t3],[Y3-40-9.35 Y3-3],'color','k','linestyle','--')
%         line([t2-2 t2-2],[Y2+40+9.35 Y2+3],'color','k','linestyle','--') % incl horiz offset to make line more visible
%         arrowxcoords = [0.6330 0.76];%[0.555 0.697];%[t1 t2]/760;
%         arrowycoords = .4961*[1 1];
        line([t3 t3],[Y3-25 Y3-3],'color','k','linestyle','--')
        line([t2-2 t2-2],[Y2+25 Y2+3],'color','k','linestyle','--') % incl horiz offset to make line more visible
        arrowxcoords = [0.68 0.76];%[0.555 0.697];%[t1 t2]/760;
        arrowycoords = .305*[1 1];
        an1=annotation('doublearrow',arrowxcoords,arrowycoords);
        set(an1,'Head2Length',4); % Make arrow format match older arrows
        set(an1,'Head1Length',4); 
        % Here I'm assuming that tau = stimstart = t2-t3: 
        if round(stimstart) ~= round(t2-t3)
            disp('Possible problem with tau calculation. Script paused.')
            pause
        end
        % For 0.8 dn
%        text(439,-29,{'$\tau = $',[num2str(round(t2-t3)) ' ms']},'Fontsize',0.6*fs)
        text(460,-60,{'$\tau = $',[num2str(round(t2-t3))], ' ms'},'Fontsize',0.6*fs)
%         arrowxcoords = [0.5092 0.5999];
%         arrowycoords = [0.3306 0.3306]; 
%         annotation('doublearrow',arrowxcoords,arrowycoords)
    end
    %title([statenames_latex(plotvarindices(ii),:) ' vs. time'],'Interpreter','latex')
    xlabel('time, ms')
    ylabel([statenames_latex(plotvarindices(ii),:) ', mmol/L'],'Interpreter','latex')
    if ~isempty(shiftstring)
        %legend(p,{'no shift',shiftstring})
%        legend(p,{'$V^*_{T=180,\breve{V}=0}$',['$V^*_{T=180,\breve{V}=' shiftlabel '}$']},'orientation','horizontal','location','best')
        if plotvarindices(ii) == 1
            legend(p,{'$V^*_{T=180,\breve{V}=0}$',['$V^*_{T=180,\breve{V}=' shiftlabel '}$']},'orientation','horizontal','location','best')
        elseif plotvarindices(ii) == 8
            legend(p,{'$Ca^*_{T=180,\breve{V}=0}$',['$Ca^*_{T=180,\breve{V}=' shiftlabel '}$']},'orientation','horizontal','location','best')
        elseif plotvarindices(ii) == 9
            legend(p,{'$Na^*_{T=180,\breve{V}=0}$',['$Na^*_{T=180,\breve{V}=' shiftlabel '}$']},'orientation','horizontal','location','south')
        elseif plotvarindices(ii) == 10
            legend(p,{'$K^*_{T=180,\breve{V}=0}$',['$K^*_{T=180,\breve{V}=' shiftlabel '}$']},'orientation','horizontal','location','best')
        end            
        %    legend boxoff;
    end
    if plotvarindices(ii) == 1
        ylabel([statenames_latex(plotvarindices(ii),:) ', mV'],'Interpreter','latex')
        figname = ['vtimeseries_b' num2str(bcls) shiftstring];
    elseif plotvarindices(ii) == 8
        figname = ['catimeseries_b' num2str(bcls) shiftstring];
    elseif plotvarindices(ii) == 9
        figname = ['natimeseries_b' num2str(bcls) shiftstring];
        yticks([19.22 19.23 19.24]);
    elseif plotvarindices(ii) == 10
        yticks([136.834 136.835 136.836 136.837]);
        figname = ['ktimeseries_b' num2str(bcls) shiftstring];
    end
%    figname = ['vtimeseries_b' num2str(bcls) '_labels']
    set(gca,'position',[.25 .18 .7 .7769])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    if ii==1
        disp('Copy-paste APD, DI, T arrows and annotations from original figure. Press any key to continue.')
        pause
    end
    print(fig,figname,'-dpdf')
    saveas(fig,figname)
    %    set(gca,'FontName', 'cmr12'); % Font changes don't work on desktop. Maybe cmr12 isn't installed? Can't tell.
end


