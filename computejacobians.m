% Script for computing Jacobians for LRd (and possibly other) models.
% This code was adapted from lrd_batch_eig.m.
% Laura Munoz, June 2017
% Input: epsln = relative perturbation size for use with diffjac_mod
% Outputs are stored in 

function computejacobians(epsln)

modelinputfolder = ['lrddata/']; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script.

fpfolder = ['fixedpoints/' ]; % folder where fixed points are stored

jacfolder = ['jacobians/']; %folder where Jacobians will be saved. 

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
    
    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    
    alljacs{i} = jacobian_cd(fp_found(:,i),epsln,modelname);
    
    % Save settings and Jacobians (OK to overwrite each time) 
    eval(['save ' jacfolder jacfilename ' *'])
end

% Time elapsed during computation
jacstoptimer = toc(jacstarttimer)



