% Use this program to compile all the fixed points into a single matrix
% Modified version of Anthony and Ryan's script to include error size
% storage.
% Modified again to use to produce higher-resolution trajectories.
% Main changes: different fixedpointfolder name, rearranged contents of for
% loop to load data.dt, set subdiv_per_bcl to bcl/data.dt.
% This is a modified version of fp_compile2_highrestraj.m that
% re-identifies fixed points that don't meet the tolerance.

clear variables;

%offset = 6.5; % offset of stimulus time, in ms
%offset = 0.75; % offset of stimulus time, in ms
%offset = 7; % this a voltage threshold in mV, not an offset in ms
offset = -50; % this a voltage threshold in mV, not an offset in ms
%shiftstring = ''; % if using default shift of data.dt
%shiftstring = '_shift0p75ms';
shiftstring = '_shift-50mVrepol';

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
fixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '/']; %alternative location if using fp_compile2.m to produce high-resolution trajectories
%fixedpointfolder = 'fixedpoints_fromaandr/';%Save compiled fixed points here
modelinputfolder = 'lrddata/'; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script.
if ~exist(fixedpointfolder,'dir')
    mkdir(fixedpointfolder)
end

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

ncyc = 1; % Run the model for this number of cycles per bcl setting.

%Compile all fixed points into a single matrix, fp_found
fp_found = zeros(17, L);
fp_errnorms = zeros(1,L);
i_repol = zeros(1,L);

%In case program crashed partway through, check for previously saved fixed
%points:
if exist([fixedpointfolder 'compiled_fp.mat'],'file')
    eval(['load ' fixedpointfolder 'compiled_fp.mat fp_found']);
end

% Find index of bcl where program crashed
startindex = 1;
%failedbcl = 130; %ms
if exist('failedbcl','var')
    bclindices = 1:L;
    startindex = bclindices(selected_bcls_for_fps == failedbcl);
end

for i = startindex:L %1:L
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
    if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
        % borrow falling-edge detection from plot_restitution.m:
        V = Y(1,:); % membrane potential
        % Search for falling edge crossing of threshold
        i_repol(i) = find(V(1:end-1)>offset & V(2:end)<=offset);
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
        % Convert stimstart to time index
        stimstartsteps = i_repol(i);
    else % otherwise the offset is in ms (right now there is no condition for a rising-edge event detection)
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - offset;
        % Convert stimstart to time index
        stimstartsteps = round(offset/outstep);
    end
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since nsoli expects our function to only
    % have one input, which is the initial condition.
    %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
    
    %    fp_found(:, i) = Y(:, 1);
    fp_found(:, i) = Y(:, stimstartsteps);
    
    % Compute fixed-point error:
    eval(['fperr = ' modelname '_p2pzero(fp_found(:,i));']); % compute fixed-point error
    fperrnorm = norm(fperr)
    % from nsoli.m, the stop criterion formula is stop_tol = atol + rtol*fnrm, computed
    % once based on the inital error. Hence, it doesn't really make sense to
    % put in an extra check on the error size, though it can be recorded.
    
    if fperrnorm > tolfactor*nsolitol
        tstartfp=tic;
        % Find fixed point:
        %        [fixedpt, it_hist, ierr, x_hist] = nsoli(fp_found(:, i),[modelname '_p2pzero'],[nsolitol nsolitol])
        [fixedpt, it_hist, ierr, x_hist] = nsoli(fp_found(:, i),[modelname '_p2pzero'],tolfactor*[nsolitol nsolitol])
        toc(tstartfp);
        
        if ierr == 0
            fp_found(:, i) = fixedpt; % fixed point
            eval(['fperr = ' modelname '_p2pzero(fixedpt);']); % compute fixed-point error
            fperrnorm = norm(fperr)
            %Save the matrix to a separate file in the folder specified above
            eval(['save ' fixedpointfolder 'compiled_fp.mat *']);
        else
            % Exit the "for" loop, and print an error message, if ierr ~= 0
            disp(['Fixed-point search failed, ierr = ' num2str(ierr)])
            save fperrdumpfile *
            sendmail('lmmsma@rit.edu',['Fixed point search error, BCL = ' num2str(bcl)]);
            break % exit the for loop
        end
        
        
    end
    fp_errnorms(i) = fperrnorm;
    
end

%Save the matrix to a separate file in the folder specified above
eval(['save ' fixedpointfolder 'compiled_fp.mat *']);
%Move generated simulation data to fixed point folder
eval(['movefile ' modelinputfolder '* ' fixedpointfolder]);

% Display minimum and maximum errors:
max(fp_errnorms)
min(fp_errnorms)

%Display full matrix if desired (uncomment)
%disp(fp_found);
