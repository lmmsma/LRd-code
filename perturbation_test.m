% Apply perturbation to initial condition, run model for specified number of
% cycles, then analyze and plot results.

clear variables;

%shiftstring = '' % if using default shift of data.dt
%shiftstring = '_shift0p75ms'
%shiftstring = '_shift6p5ms'
%shiftstring = '_shift7mVrepol'
shiftstring = '_shift-50mVrepol'

%fpsimfolder = 'lrddata_fromfixedpointsearch_highres/';%This folder is the source for simulations
fpsimfolder = 'lrddata_fromfixedpointsearch_highres';%This folder is the source for simulations
fpsimfolder = [fpsimfolder shiftstring '/'];
% where ICs were chosen as fixed points
% Unshifted:
if isempty(shiftstring)
    fixedpointfolder = 'fixedpoints/';%This folder is the source for fixed points
    % Shifted:
else
    fixedpointfolder = fpsimfolder;%This folder is the source for fixed points
end
%eigfolder = 'eigenvalues/def/'; %This folder is the source for eigenvectors and eigenvalues
eigfolder = ['eigenvalues' shiftstring '/def/']; %This folder is the source for eigenvectors and eigenvalues
pertfolder = 'perturbation_tests/';%Save perturbation outputs here
modelinputfolder = 'lrddata/'; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script. It is also the default location where
% simulation output is stored.

pertflag = 0; %0 to perturb IC along an eigenvector,
% 1 to select a point along a closed (fixed-point) trajectory

selected_bcls_for_perturbations = 1000;%330;%190;%1000;%70; %[1000];%Vector of bcls
numselbcls = length(selected_bcls_for_perturbations);

modelname = 'lrd';

epsln = 10^-5; % This is the epsilon value used in Jacobian computation.

%eigvalselect = 1.0 % To perturb along an eigenvector, select the corresponding
% eigenvalue here. For now, this value needs to be a real number,
% with real part within +/- epsln of the true eigenvalue.
%eigvalselect = 0.718622401 % a 1000ms BCL eigval, for some shifts
eigvalselect = 0.71859 %0.718587654698170 this is the corresponding value for the 50mV, 1000ms case
%eigvalselect = 0.111706954019 % a 1000ms BCL eigval
%eigvalselect = 0.0195661939 % a 330ms eigval 
%eigvalselect = -1.073215029 % a 190ms BCL eigval
%eigvalselect = 0.9647198918 % a 70ms eigval, for some shifts
%eigvalselect = 0.96471 %0.964707494652390 this is the corresponding value for the 7mV, 70ms case

%pertsize = -10^6*epsln % This is the size of the perturbation that will be applied
pertsize = epsln % This is the size of the perturbation that will be applied
%pertsize = 0; %10^2*epsln % This is the size of the perturbation that will be applied
% to the default IC

% load fixed points
if isempty(shiftstring)
    eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps'])
else
    if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
        eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps data i_repol'])
    else
        eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps data'])
    end
end

matobj = matfile([fixedpointfolder 'compiled_fp.mat']);
if isprop(matobj,'offset')
    eval(['load ' fixedpointfolder 'compiled_fp offset'])
end

bclfpindices = 1:length(selected_bcls_for_fps);

% load eigenvalues and eigenvectors
%eval(['load ' eigfolder '/eigfile' num2str(log10(epsln)) ' alleigs allv'])
eval(['load ' eigfolder 'eigfile' num2str(log10(epsln)) ' alleigs allv'])

%fp_found = zeros(17, numselbcls);
%fp_errnorms = zeros(1,numselbcls);
initcond = zeros(17, numselbcls); % store initial conditions here
inputdiff = NaN*ones(17, numselbcls); % initial deviations from fixed point
outputdiff = NaN*ones(17, numselbcls); % final deviations from fixed point
outin_mag_ratio = NaN*ones(1, numselbcls); % ratio of output magnitude to input magnitude
outin_ang = NaN*ones(1, numselbcls); % angle between output and input vectors
pertnorm = cell(1,numselbcls); % store difference between system state and FP here

for i = 1:numselbcls
    bcl = selected_bcls_for_perturbations(i);
    bclind = bclfpindices(selected_bcls_for_fps == bcl);
    
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    %    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    subdiv_per_cyc = bcl/0.005;
    subdiv_per_cyc_pert = subdiv_per_cyc; %store this value
    
    ncyc = 1; % Run the model for this number of cycles per bcl setting.
    %ncyc = 5143; % For 70ms BCL, 6 min is approximately 5143 cycles
    %ncyc = 360; % For 1000ms BCL, 6 min is 360 cycles
    %    ncyc = round(6*60*1000/bcl); % convert 6 min to number of BCLs
    %    ncyc = round(10*60*1000/bcl); % convert 10 min to number of BCLs
    %    ncyc = round(5*60*1000/bcl); % convert 5 min to number of BCLs
    %ncyc = 3;
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since nsoli expects our function to only
    % have one input, which is the initial condition.
    %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    % Add code to check for stimulus offset, which (if it exists) was stored in
    % compiled_fp mat file
    if exist('offset','var')
        if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
            % Stimulus start time:
            % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
            % (bcl + data.dt) – (offset-data.dt). I verified that this formula
            % produces matched values for BCL = 1000ms.
            % Indexing of i_repol must match selected_bcls_for_fps
            %stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
            stimstart = bcl + 2*data.dt - i_repol(bclind)*data.dt;
        else % otherwise the offset is in ms (right now there is no condition for a rising-edge event detection)
            stimstart = bcl + 2*data.dt - offset;
        end
        
        % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
        % Pass-by-file is being used since nsoli expects our function to only
        % have one input, which is the initial condition.
        %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
    else
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    end
    
    % Find eigenvector, from Jacobian for this bcl, corresponding to
    % eigenvalue closest to eigvalselect
    eigind = find(real(alleigs{bclind})< eigvalselect+epsln & eigvalselect-epsln <real(alleigs{bclind}) ...
        & abs(imag(alleigs{bclind})) < epsln);
    if ~isempty(eigind)
        direction = allv{bclind}(:,eigind);
    else
        disp(['Error: Eigenvalue at ' num2str(eigvalselect) ' not found.'])
        break;
    end
    
    % filename for both unperturbed and perturbed trajectories (they will
    % be distinguished by folder name)
    simoutfname = ['lrddata_1cell_b' num2str(selected_bcls_for_perturbations(i))];
    
    if ~pertflag % perturb IC along eigenvector
        % Load fixed point matching selected bcl, and
        % add a perturbation of size epsln to the fixed point, in the direction
        % corresponding to the selected eigenvector
        initcond(:, i) = fp_found(:, bclind) + pertsize*direction;
    elseif pertflag == 1 % perturb along closed fixed-point trajectory
        % load fixed-point trajectory
        eval(['load ' fpsimfolder simoutfname ' Y subdiv_per_cyc;']);
        Yfp = Y;
        % Ask user where (along the closed fixed-point trajectory) to
        % choose the new IC
        pertsteps = input(['The closed fixed-point trajectory had a sampling interval \nof ' num2str(bcl/subdiv_per_cyc) ' ms, with ' num2str(subdiv_per_cyc*ncyc) ' total sample steps. To choose the new IC, \nenter a whole number corresponding to the number of \nsampling intervals to move along the trajectory: ']);
        initcond(:,i) = Yfp(:,pertsteps+1); % The first column contains the IC, so displace one more step to get to the next set of values
    end
    
    inputdiff(:, i) = initcond(:,i) - fp_found(:,bclind);
    
    % Compute error relative to IC:
    eval(['perterr = ' modelname '_p2pzero(initcond(:,i));']); % compute fixed-point error
    perterrnorm = norm(perterr)
    
    % Load trajectory produced by model
    eval(['load ' modelinputfolder simoutfname ' Y;']);
    
    outputdiff(:, i) = Y(:,end) - fp_found(:,bclind);
    
    % compute magnitude ratio and angle between input and output
    outin_mag_ratio(i) = norm(outputdiff(:,i))/norm(inputdiff(:,i))
    % if perturbation is along an eigenvector, the above ratio should match the
    % eigenvalue
    
    costheta = dot(outputdiff(:,i),inputdiff(:,i))/(norm(outputdiff(:,i))*norm(inputdiff(:,i)));
    outin_ang(i) = acosd(costheta)
    % if perturbation is along an eigenvector, the above angle should be close to zero
    
    % Copy model output to perturbation folder:
    eval(['movefile ' modelinputfolder simoutfname '.mat ' pertfolder simoutfname '.mat'])
    %h(i) = figure;
    
    numtimes = size(Y,2);
    pertnorm{i} = zeros(1,numtimes);
    for kk=1:numtimes
        pertnorm{i}(kk) = norm(Y(:,kk)-fp_found(:,bclind));
    end
    figure
    hold on;
    plot((1:numtimes)*(bcl/subdiv_per_cyc), pertnorm{i})
    plot(((0:ncyc-1)*subdiv_per_cyc +1)*bcl/subdiv_per_cyc,pertnorm{i}((0:ncyc-1)*subdiv_per_cyc +1),'ro')
    
end

%Save the matrix to a separate file in the folder specified above
%eval(['save ' pertfolder 'perturbation_results *']);
eval(['save ' pertfolder 'perturbation_results' shiftstring ' *']);

