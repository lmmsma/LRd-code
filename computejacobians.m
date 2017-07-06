% Script for computing Jacobians for LRd (and possibly other) models.
% This code was adapted from lrd_batch_eig.m.
% Laura Munoz, June 2017

clear variables;

fpfolder = ['lrddata/']; % folder where fixed points are stored

jacfolder = ['lrddata/']; %folder where Jacobians will be saved. 
% Currently this folder must match the one where lrd_p2p.m looks for 
% inputs. 

ncyc = 1; % Run the model for this number of cycles per bcl setting,
% when computing Jacobians.
epsln = 1e-7; % relative perturbation size for use with diffjac_mod
% LRd: 1e-5 may be better for biphasic V, 1e-7 for monophasic K+ stimulus

% Here I'm assuming a MAT file exists that contains the following
% quantities:
% allfp (matrix where each column is a fixed-point vector)
% bcls (vector of corresponding bcls; same number of columns as allfp)
% data (this is the structure created by constantsLRd_strand)
% modelname (a string such as 'lrd')
% allstimstart (stimulus start times; will use later)
%
% Load all of the above quantities:
eval(['load ' fpfolder 'fpfile *'])

% Initialize cell arrays
alljacs = cell(1,length(bcls)); % Store Jacobians here. Could instead use
% a 3D array, but I think then you have to use "squeeze" to extract things.

jacstarttimer = tic;
for i = 1:length(bcls)
    bcl = bcls(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    %    subdiv_per_cyc = bcl/data.dt;
    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since diffjac_mod (which is called by
    % jacobian_cd) does not currently allow these quantities to be inputs,
    % but diffjac_mod could be modified further.
    
    eval(['save ' jacfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    
    alljacs{i} = jacobian_cd(allfp(:,i),epsln,modelname);
end

% Time elapsed during computation
jacstoptimer = toc(jacstarttimer)

% Save settings and Jacobians
eval(['save ' jacfolder 'jacfile *'])


