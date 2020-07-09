% Script for computing Jacobians for LRd (and possibly other) models.
% This code was adapted from lrd_batch_eig.m.
% Laura Munoz, June 2017
% Input: epsln = relative perturbation size for use with diffjac_mod
% Outputs are stored in jacfolder. 
% March 2019: adapt to work with new offset timing computations from
% fp_compile2_shiftedtraj.m and plotV_allbcls.m 

function computejacobians(epsln)

%shiftstring = ''; % if using default shift of data.dt
%shiftstring = '_shift0p75ms';
%shiftstring = '_shift6p5ms';
%shiftstring = '_shift7mVrepol';
%shiftstring = '_shift-50mVrepol';
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
shiftstring = '_shift0p001Vnormrepol';

modelinputfolder = ['lrddata/']; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script.

%fpfolder = ['fixedpoints/' ]; % folder where fixed points are stored
%fpfolder = ['lrddata_fromfixedpointsearch_highres_shift6p5ms/' ]; % folder where fixed points are stored
fpsimfolder = 'lrddata_fromfixedpointsearch_highres';%This folder is the source for simulations
fpsimfolder = [fpsimfolder shiftstring '/'];
% Unshifted:
if isempty(shiftstring)
    fpfolder = 'fixedpoints/';%This folder is the source for fixed points
    % Shifted:
else
    fpfolder = fpsimfolder;%This folder is the source for fixed points
end

%jacfolder = ['jacobians/']; %folder where Jacobians will be saved.
%jacfolder = ['jacobians_shift6p5ms/']; %folder where Jacobians will be saved.
jacfolder = ['jacobians' shiftstring '/']; %folder where Jacobians will be saved.
if ~exist(jacfolder,'dir')
    mkdir(jacfolder)
end


jacfilename = ['jacfile' num2str(log10(epsln)) '.mat']; % file name should depend on epsilon

ncyc = 1; % Run the model for this number of cycles per bcl setting,
% when computing Jacobians.
%epsln = 1e-4; % relative perturbation size for use with diffjac_mod

% Here I'm assuming a MAT file exists that contains the following
% quantities:
% fp_found (matrix where each column is a fixed-point vector)
% selected_bcls_for_fps (vector of corresponding bcls; same number of columns as compiled_fp)
% data (this is the structure created by constantsLRd_strand)
% modelname (a string such as 'lrd')
% allstimstart (stimulus start times; will use later)

% Load all of the above quantities:
eval(['load ' fpfolder 'compiled_fp *'])
% Load repol and depol times:
load([fixedpointfolder 'threshold_times.mat'],'i_depol','i_repol');

% Initialize cell arrays
alljacs = cell(1,length(selected_bcls_for_fps)); % Store Jacobians here. Could instead use
% a 3D array, but I think then you have to use "squeeze" to extract things.

jacstarttimer = tic;
for i = 1:length(selected_bcls_for_fps)
    bcl = selected_bcls_for_fps(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    %    subdiv_per_cyc = bcl/data.dt;
    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since diffjac_mod (which is called by
    % jacobian_cd) does not currently allow these quantities to be inputs,
    % but diffjac_mod could be modified further.
    
    % Add code to check for stimulus offset, which (if it exists) was stored in
    % compiled_fp mat file
%     if exist('offset','var')
%         if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
%             % Stimulus start time:
%             % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
%             % (bcl + data.dt) – (offset-data.dt). I verified that this formula
%             % produces matched values for BCL = 1000ms.
%             % Indexing of i_repol must match selected_bcls_for_fps
%             stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
%         else % otherwise the offset is in ms (right now there is no condition for a rising-edge event detection)
%             stimstart = bcl + 2*data.dt - offset;
%         end
%         
%         % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
%         % Pass-by-file is being used since nsoli expects our function to only
%         % have one input, which is the initial condition.
%         %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
%         eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
%     else
%         eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
%     end
    
    % All of the following if/else statements assume that outstep=data.dt:
    if strfind(shiftstring,'repol') % in this case the offset is in normalized membrane potential units on the repolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
    elseif strfind(shiftstring,'depol') % in this case the offset is in normalized membrane potential units on the depolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_depol(i)*data.dt;
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
    else
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    end
    
    alljacs{i} = jacobian_cd(fp_found(:,i),epsln,modelname);
    
    % Save settings and Jacobians (OK to overwrite each time)
    eval(['save ' jacfolder jacfilename ' *'])
end

% Time elapsed during computation
jacstoptimer = toc(jacstarttimer)



