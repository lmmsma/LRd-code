%LMM: perform fixed-point search on LRD model. With enhancements from A.
%Guzman and R. Vogt. 
clear variables;
%% Settings for Fixed Point search
ms = 10; fs = 14; % marker size and font size

modelinputfolder = ['lrddata/']; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script. 

fixedptfolder = 'fixedpoints/'; % save data here

setpref('Internet','E_mail','lmmsma@rit.edu');
setpref('Internet','SMTP_Server','mymail.ad.rit.edu');

constantsLRd_strand % load settings and parameters for model

% The variable data.stimflag is defined in constantsLRd_strand.
% If data.stimflag is 0, apply biphasic stimuli to the cell through
% V dynamics. Otherwise, apply monophasic stimuli through K+-dynamics.
if data.stimflag
    systemselect = 'kmonovsolem12'
    %The systemselect string is used in certain filenames
else
    systemselect = 'solem12'
end

%statenames = strvcat('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
%numstate = length(statenames);

%modelname = 'karma_sim_estimtest';
modelname = 'lrd';

%numpart = 1; % number of partitions (cells)

% List of stimulus start times: Make a decision about what time, relative
% to the beginning of a BCL, where the main stimulus will be applied.
% In an older version of the code, a list of start times was contained
% in allstimstart, and we could consider revisiting this if we want to
% loop over relative start times.
allstimstart = data.dt;
%allstimstart = data.dt + 200;
%allstimstart = data.dt + [400 600 800];
stimstart = allstimstart(1);% stimulus start time, ms

%save lrdparams L numpart bcl stimstart % these will be read in by lrd_p2p.m

nsolitol = 1e-12; % abs & rel tolerance for Newton-Krylov solver
%useoldfp = 1; % set to zero to always solve for fixed point; set to 1 to re-use fixed point
maxabsfperr = NaN; % maximum absolute value of p2p error (use to determine whether OL f.p. is a valid CL f.p.)

% % initialize nsoli output variables, in case usoldfp==1
% it_hist = -1;
% ierr=NaN;
% x_hist=NaN;

%% Set up BCL's to be evaluated by prompt or manually (hard-coded)
%USE ONLY ONE bcls vector; Comment out all other lines
%---------------------Prompts--------------------------------%
% bcl_max= input('Enter max BCL\n');
% bcl_step= -input('Enter BCL step size\n');
% bcl_num= input('Enter number of BCLs to evaluate\n');
% bcls= [bcl_max:bcl_step:bcl_max+bcl_step*(bcl_num-1)]


%---------------------Manual Setup---------------------------%
% BCL steps loosely based on Koller/Riccio/Gilmour dynamic protocol:
% bcls = [600 550 500]
%bcls = [1000:-50:200 190:-10:70 69:-1:50]; %full pacedown; cycle lengths in ms
%bcls = [1000:-50:300 290:-10:70 69:-1:50]; %revised full pacedown (more points leading up to 200ms)
%bcls = [400:-50:300 290:-10:70]; %revised shortened pacedown
%bcls = [50:-1:40]; %cycle lengths in ms
% Below is A&R's bcls_tested vector from fp_compile.m. Is this for the
% default parameters? 
%bcls = [1000:-50:600 590:-10:510 505 500:-10:70];%Vector of bcls
% Reduce step size after failure at 850ms: 
%bcls = [875:-25:600 590:-10:510 505 500:-10:70];%Vector of bcls
% Reduce step size after failure at 800ms: 
%bcls = [810:-10:510 505 500:-10:70];%Vector of bcls
% Reduce step size after failure at 570ms: 
%bcls = [575:-5:70];%Vector of bcls
% Reduce step size after failure at 530ms: 
%bcls = 532:-2:70;%Vector of bcls
% Reduce step size after failure at 522ms: 
%bcls = 523:-1:70;%Vector of bcls
% Increase maxarm from 20 to 30 in nsoli.m, after ierr=2 failure at 512ms: 
%bcls = 513:-1:70;%Vector of bcls
% Increase maxit from 40 to 50 in nsoli.m, after ierr=1 failure at 507ms: 
%bcls = 507:-1:70;%Vector of bcls
% Increase maxarm from 30 to 40 in nsoli.m, after ierr=2 failure at 505ms: 
%bcls = 505:-1:70;%Vector of bcls
% Reduce step size, after ierr=1 failure at 502ms: 
%bcls = 502:-0.5:70;%Vector of bcls
% Restart search at 500 using 505 fixed point, after another failure at 502ms.
% Increase bcl step size, even though this is counterintuitive:
%bcls = 500:-5:70;%Vector of bcls
% Continue search at 450ms, after ierr=2 error:  
%bcls = 450:-10:70;%Vector of bcls
% Continue search at 410ms, after ierr=2 error. Increase maxarm from 540 to 0 in nsoli.m: at  
bcls = 410:-10:70;%Vector of bcls

%% Loading Initial Conditions for starting BCL, if they exist
if bcls(1) == 1000
    bcl = bcls(1);
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b' num2str(bcl)]; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    else % load fixed points from older study
        if strcmp(systemselect, 'solem12')
            unpertfilename = 'b1000fsolem12_fwde_shift0_newpulse';
        elseif strcmp(systemselect,'kmonovsolem12')
            unpertfilename = 'b1000fkmonovsolem12_fwde_shift0_newpulse_relpert1e-7';
        end
        % Load the unperturbed fixed point
        eval(['load ' unpertfilename ' ' systemselect])
        yinit = eval(systemselect);
    end
elseif bcls(1) == 875
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b900']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 810
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b825']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 575
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b580']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 532
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b535']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 523
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b524']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 513
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b513']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 507
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b507_old']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 505
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b505_old']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 502
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b502_old']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 500
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b505']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
elseif bcls(1) == 450
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b450_old']; %
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    end
else
    bcl = bcls(1);
    % Filenames and label settings based on fixed point type:
    fname = [fixedptfolder 'lrddata_1cell_b' num2str(bcl)]; 
    if exist([fname '.mat'],'file') % if a fixed point was found on a previous run, load it here
        eval(['load ' fname ' Y']) % load the state vectors
        yinit = Y(:,end); % load final condition of previous recording        
    else % load fixed points from older study
        if strcmp(systemselect, 'solem12')
            unpertfilename = 'b1000fsolem12_fwde_shift0_newpulse';
        elseif strcmp(systemselect,'kmonovsolem12')
            unpertfilename = 'b1000fkmonovsolem12_fwde_shift0_newpulse_relpert1e-7';
        end
        % Load the unperturbed fixed point
        eval(['load ' unpertfilename ' ' 'yin'])
        yinit = eval(systemselect);
    end
end

%% Find Fixed Points
ncyc = 1; % Run the model for this number of cycles per bcl setting.
%To identify a period-n trajectory, should set ncyc = n.

% initialize matrices
allfp = zeros(length(yinit),length(bcls)); % fixed points
allierr = NaN*ones(1,length(bcls)); % nsoli exit conditions
allfperr = NaN*ones(1,length(bcls)); % fixed-point errors

fpstarttimer = tic;
for i = 1:length(bcls)
    bcl = bcls(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
%    ncyc = ncycs(i);
    
    %    subdiv_per_cyc = bcl/data.dt;
    %    subdiv_per_cyc = bcl;
    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    %    subdiv_per_cyc = 5*bcl;
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since nsoli expects our function to only
    % have one input, which is the initial condition.
    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    
    % The following portion should eventually be revised to look for and load fixed
    % points from previous runs, if available(if useoldfp == 1), or to
    % always try to recalculate fixed points (if useoldfp == 0).
    %if ~useoldfp
    if exist('fpfolder')
        eval(['load ' fpfolder '/1 Y']);
        ivshiftindex = bcl/data.dt - stimstart/data.dt + 2;
        if ivshiftindex > bcl/data.dt
            ivshiftindex = ivshiftindex - bcl/data.dt;
        end
        yinit = Y(:,ivshiftindex);
        tstartfp=tic;
        % Find fixed point:
        [fixedpt, it_hist, ierr, x_hist] = nsoli(yinit,[modelname '_p2pzero'],[nsolitol nsolitol])
        toc(tstartfp);
    else
        tstartfp=tic;
        % Find fixed point:
       [fixedpt, it_hist, ierr, x_hist] = nsoli(yinit,[modelname '_p2pzero'],[nsolitol nsolitol])
%        [fixedpt, it_hist, ierr, x_hist] = nsoli(yinit,[modelname '_p2pzero'],[nsolitol nsolitol],[40, 40, .9, 2, 20])
        toc(tstartfp);
    end
    %else
    %    fixedpt = ???; % load fixed pt value from previous run
    %end
    
    % Store results in matrices
    allierr(i) = ierr; % exit condition
    
    if ierr == 0
        allfp(:,i) = fixedpt; % fixed point
    else 
        % Exit the "for" loop, and print an error message, if ierr ~= 0
        disp(['Fixed-point search failed, ierr = ' num2str(ierr)])
        save fperrdumpfile *
        sendmail('lmmsma@rit.edu',['Fixed point search error, BCL = ' num2str(bcl)]);
        break % exit the for loop
    end
    
    eval(['fperr = ' modelname '_p2pzero(fixedpt);']); % compute fixed-point error
    maxabsfperr = max(abs(fperr))
    if maxabsfperr > nsolitol
        disp('Fixed point error exceeds tolerance.')
        break
    end
        
    % Store fixed-point error
    allfperr(i) = maxabsfperr;
    
    % Next IC = fixed point identified at previous BCL
    yinit = fixedpt;
    
    % Save settings and fixed-point values
    eval(['save ' fixedptfolder 'fpfile *'])
end

% Time elapsed during search
fpstoptimer = toc(fpstarttimer)

