% Script for computing Jacobians for LRd (and possibly other) models.
% This code was adapted from lrd_batch_eig.m.
% Laura Munoz, June 2017
% Input: epsln = relative perturbation size for use with diffjac_mod
% Outputs are stored in jacfolder. 
% March 2019: adapt to work with new offset timing computations from
% fp_compile2_shiftedtraj.m and plotV_allbcls.m 
% Apr 2020: adapt to compute CL Jacobians (just save in a different folder)

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
shiftstring = '_shift0p2Vnormrepol';
%shiftstring = '_shift0p001Vnormrepol';

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
jacfolder = ['jacobians_cl' shiftstring '/']; %folder where Jacobians will be saved.
if ~exist(jacfolder,'dir')
    mkdir(jacfolder)
end


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

selected_bcls_for_perturbations = 200;%200;%330;%190;%1000;%70; %[1000];%Vector of bcls
numselbcls = length(selected_bcls_for_perturbations);
bclfpindices = 1:length(selected_bcls_for_fps);

% Initialize cell arrays
alljacs = cell(1,length(selected_bcls_for_fps)); % Store Jacobians here. Could instead use
% a 3D array, but I think then you have to use "squeeze" to extract things.

% This gain is for bcl = 1000ms, shift = 0.2dn, measurementindex = 1
measurementindex = 1; 

jacstarttimer = tic;


for i = 1:numselbcls%1:length(selected_bcls_for_fps)
    bcl = selected_bcls_for_perturbations(i);
    bclind = bclfpindices(selected_bcls_for_fps == bcl);
    
    if bcl == 1000 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 1
        Llu = [0.8906   -0.0036    0.0069   -0.0015    0.0000   -0.0058    0.0049   0.0001   -0.0729   -4.9068    0.0274    0.0139    0.0017    0.0050   -0.0034    0.0011    0.0000]; % for BCL = 1000ms, 0.2dn, unscaled, measurementindex =1
    elseif bcl == 200 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 1
% high-precision format only seemed to be necessary for bcl=100 test cases,
% but using it here just in case
%        Llu = [-0.0058   -0.0058    0.0001   -0.0028    0.0000   -0.0240    0.0221    0.0006   -1.8410  -88.0942    0.1702    0.0914     0.0298    0.0189   -0.0075    0.0327    0.0001]; 
        Llu = [-0.005750641503498  -0.005763357664634   0.000144145593421  -0.002775791569397   0.000031526085499  -0.023982849157018   0.022071810256199   0.000564080635719  -1.840997251403683 -88.094177613076809   0.170188064295353   0.091422577584066   0.029808134438865   0.018880643058383  -0.007535007734080   0.032690307284951   0.000065828675278];
    elseif bcl == 200 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 10
% high-precision format only seemed to be necessary for bcl=100 test cases,
% but using it here just in case
%        Llu = 1.0e+04*[-2.1538    0.0086   -0.0149    0.0034   -0.0001    0.0121  -0.0131    0.0003    0.0133    0.0001   -0.0825    0.0019  -0.0021   -0.0119    0.0074   -0.0006   -0.0000]; 
        Llu = 1.0e+04*[-2.153775638316071   0.008642037382794  -0.014948842556174   0.003435418563187 -0.000064257363951   0.012094311408771  -0.013064809060953   0.000292254463773   0.013266330987352   0.000105514852743  -0.082504075852370   0.001857273772603  -0.002079659054639  -0.011892992039428   0.007354671108559  -0.000574379933953  -0.000039383516605]; 
    elseif bcl == 100 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 1
        Llu = 1.0e+03*[0.000389751311207  -0.000093237536550 0.000005619060195  -0.000045241960194  0.000000502975428  -0.000499549611300  0.000606910864915  -0.000012097865701 -0.028104784284541  -1.214325117276513 -0.000141203282165  -0.000582614472652  0.000410382237115   0.000330869242412  -0.000140684839695   0.000427576873080   0.000006777685206]; 
    elseif bcl == 100 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 10
        Llu = 1.0e+03*[-1.097925810265023   0.003637262000388  -0.007774919161288   0.001379787464483  -0.000030463631231  -0.001660867070130   0.003101636295633  -0.000891059164431   0.575366283569920   0.001218617080828  -0.010322908419630  -0.064933082658678  -0.004246377455183  -0.000350153571477   0.002136706192036  -0.000529096785608  -0.000014687317062]; 
    elseif bcl == 100 && strcmp(shiftstring,'_shift1Vnormdepol') && measurementindex == 1     
        Llu = 1.0e+03*[0.000369234530581   0.004460211749332  -0.000002114922454   0.018044286497205  -0.026058394996194   0.002715438654225  -0.007518137614274   0.000158664470606   0.270742001642501   4.424574882768074   0.007946070291199   0.012767879478589  -0.002879644989076  -0.007794061132781   0.001457074821435  -0.002893719441405  -0.003913750149751]; 
    end    
    
    jacfilename = ['jacfile_b' num2str(bcl) '_meas' num2str(measurementindex) '.mat']; % file name should depend on epsilon
    
    fp_for_fbk = fp_found(:,bclind); 
%    bcl = selected_bcls_for_fps(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    subdiv_per_cyc = bcl/data.dt;
    %subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    
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
        stimstart = bcl + 2*data.dt - i_repol(bclind)*data.dt;
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Llu'])
    elseif strfind(shiftstring,'depol') % in this case the offset is in normalized membrane potential units on the depolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_depol(bclind)*data.dt;
%        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Llu'])
    else
%        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
         eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Llu'])
   end
    
    alljacs{bclind} = jacobian_cd(fp_found(:,bclind),epsln,modelname);
    
    % Save settings and Jacobians (OK to overwrite each time)
    eval(['save ' jacfolder jacfilename ' *'])
    eig(alljacs{bclind})
end

% Time elapsed during computation
jacstoptimer = toc(jacstarttimer)



