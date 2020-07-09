% Modified version of perturbation_test.m to try to apply
% Luenberger gains to single-cell nonlinear system. 
% Apply perturbation to initial condition, run model for specified number of
% cycles, then analyze and plot results.

clear variables;

%shiftstring = '' % if using default shift of data.dt
%shiftstring = '_shift0p75ms'
%shiftstring = '_shift6p5ms'
%shiftstring = '_shift7mVrepol'
%shiftstring = '_shift-50mVrepol'
%shiftstring = '_shift1Vnormdepol';
shiftstring = '_shift0p2Vnormrepol'
rng(1);

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
numstate = length(statenames); 
%%fpsimfolder = 'lrddata_fromfixedpointsearch_highres/';%This folder is the source for simulations
%fpsimfolder = 'lrddata_fromfixedpointsearch_highres';%This folder is the source for simulations
%fpsimfolder = [fpsimfolder shiftstring '/'];

ncyc = 30; % Run the model for this number of cycles per bcl setting.

% where ICs were chosen as fixed points
if ncyc == 1
    if isempty(shiftstring)
        fixedpointfolder = ['lrddata_fromfixedpointsearch_highres011719/'];
    else
        fixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '/']; %alternative location if using fp_compile2.m to produce high-resolution trajectories
    end
else
    fixedpointfolder = ['lrddata_fromfixedpointsearch_highres' shiftstring '_ncyc' num2str(ncyc) '/']; %alternative location if multi-cycle simulation
end
%eigfolder = 'eigenvalues/def/'; %This folder is the source for eigenvectors and eigenvalues
eigfolder = ['eigenvalues' shiftstring '/def/']; %This folder is the source for eigenvectors and eigenvalues
pertfolder = 'nonlinear_luen_tests/';%Save perturbation outputs here
modelinputfolder = 'lrddata/'; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script. It is also the default location where
% simulation output is stored.

% load approximate "state normalization" scaling matrix
load b1000fsolem12variable_amplitudes varamp
% Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
Smat = diag(1./varamp);

pertflag = 4; %0 to perturb IC along an eigenvector,
% 1 to select a point along a closed (fixed-point) trajectory
% 2 to perturb along a state variable direction
% 3 to perturb only along Na and K in a charge-conservative manner
% 4 to apply random perturbation to all state variable directions

%pertsize = -10^6*epsln % This is the size of the perturbation that will be applied
%pertsize = epsln % This is the size of the perturbation that will be applied
pertsize = 10^-5;%10^-2; %10^2*epsln % This is the size of the perturbation that will be applied
% to the default IC

% if pertflag = 2 (perturb along state variable direction), indicate the
% direction in which the IC will be perturbed: 
sv_pert_dir = zeros(17,1); 
if pertflag == 2
    sv_pert_ind = 1;
%    sv_pert_ind = 10;    
    sv_pert_dir(sv_pert_ind) = 1; 
elseif pertflag == 3
    sv_pert_ind = 'NaK'; 
    sv_pert_dir(9) = -1; % Na dir
    sv_pert_dir(10) = 1; % K dir
elseif pertflag == 4
    sv_pert_ind = 'rand';
    sv_pert_dir = inv(Smat)*(rand(numstate,1)-0.5*ones(numstate,1))
end

% Select BCL: 
selected_bcls_for_perturbations = 500;%330;%190;%1000;%70; %[1000];%Vector of bcls

numselbcls = length(selected_bcls_for_perturbations);

measurementindex=  input('Enter Measurement Index (enter 0 for open-loop): ');
if measurementindex
    measlabel =  statenames_latex(measurementindex, :); 
else
    measlabel = 'OL'; 
end

adj_yn = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if adj_yn
    param = 'adj';
else
    param = 'def';
end

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


% load fixed points
if isempty(shiftstring)
    eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps data'])
else
%    if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
        eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps data i_repol i_depol'])
%    else
%        eval(['load ' fixedpointfolder 'compiled_fp fp_found selected_bcls_for_fps data'])
%    end
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
    
    fp_for_fbk = fp_found(:,bclind); 
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    %    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    subdiv_per_cyc = bcl/data.dt;
    subdiv_per_cyc_pert = subdiv_per_cyc; %store this value

%    ncyc = 10; 
%    ncyc = 5; 
%    ncyc = 20;
%    ncyc = 1; % Run the model for this number of cycles per bcl setting.
    %ncyc = 5143; % For 70ms BCL, 6 min is approximately 5143 cycles
    %ncyc = 360; % For 1000ms BCL, 6 min is 360 cycles
    %    ncyc = round(6*60*1000/bcl); % convert 6 min to number of BCLs
    %    ncyc = round(10*60*1000/bcl); % convert 10 min to number of BCLs
    %    ncyc = round(5*60*1000/bcl); % convert 5 min to number of BCLs
    %ncyc = 3;
    if bcl == 1000 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 1
%        Llu = [0.8906   -0.0036    0.0069   -0.0015    0.0000   -0.0058    0.0049   0.0001   -0.0729   -4.9068    0.0274    0.0139    0.0017    0.0050   -0.0034    0.0011    0.0000]; % for BCL = 1000ms, 0.2dn, unscaled, measurementindex =1
        Llu = [0.890607864346635  -0.003647503422018   0.006892381791094  -0.001516846279912   0.000042136130225  -0.005799602166817   0.004895716049793   0.000134901249450  -0.072861565844498  -4.906848956901893   0.027400078932940   0.013877025306040   0.001656342713882   0.004975173021898  -0.003377442125606   0.001111871473395   0.000000035532622];
    elseif bcl == 1000 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 10
        Llu = 1.0e+03*[1.979361164614477  -0.007638557410004   0.015297467913533  -0.003145312228526   0.000090000191887  -0.012334031007641   0.010982833813179  -0.000429103892945   0.023784473807704   0.001252464083921   0.037132691040578  -0.042652461052620   0.004020131389316   0.010865138940394  -0.006937123954572   0.001587196762620   0.000000066732027]; 
    elseif bcl == 500 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 1
        Llu = [0.783427356944253  -0.003746312094000   0.005888133087999  -0.001581789776874   0.000037317326505  -0.006926015853800   0.005450088289279   0.000361110850909  -0.230376631980297 -14.052735276475454   0.050748647591461   0.038500953659538   0.003229660195650   0.006010431447710  -0.003594727620941   0.003795619115661   0.000000233869891]; 
    elseif bcl == 500 && strcmp(shiftstring,'_shift0p2Vnormrepol') && measurementindex == 10
        Llu = 1.0e+02*[-1.460080709332431  -0.004956621492517  -0.010597932234220  -0.003079755698088   0.000001335203782  -0.082787249316491   0.079744349153248  -0.019908029644496   0.284648501245726   0.011621494863116   5.308723521345258   0.029952196674306   0.036141441911792   0.086526473826801  -0.010914674513371   0.206809600129045   0.000017590452200];
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
    elseif measurementindex == 0
        Llu = zeros(1,numstate); 
    end    

    % filename for both unperturbed and perturbed trajectories (they will
    % be distinguished by folder name)
    simoutfname = ['lrddata_1cell_b' num2str(selected_bcls_for_perturbations(i))];
    eval(['load ' fixedpointfolder simoutfname ' Y;']);
    Yfp_traj = Y; 

    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since nsoli expects our function to only
    % have one input, which is the initial condition.
    %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    % Add code to check for stimulus offset, which (if it exists) was stored in
    % compiled_fp mat file
%     if exist('offset','var')
% %        if strfind(shiftstring,'mVrepol') % in this case the offset is in mV and should be detected on the falling edge
%             % Stimulus start time:
%             % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
%             % (bcl + data.dt) – (offset-data.dt). I verified that this formula
%             % produces matched values for BCL = 1000ms.
%             % Indexing of i_repol must match selected_bcls_for_fps
%             %stimstart = bcl + 2*data.dt - i_repol(i)*data.dt;
%             stimstart = bcl + 2*data.dt - i_repol(bclind)*data.dt;
% %         else % otherwise the offset is in ms (right now there is no condition for a rising-edge event detection)
% %             stimstart = bcl + 2*data.dt - offset;
% %         end
%         
%         % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
%         % Pass-by-file is being used since nsoli expects our function to only
%         % have one input, which is the initial condition.
%         %    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
%         eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart'])
    % All of the following if/else statements assume that outstep=data.dt:
    if strfind(shiftstring,'repol') % in this case the offset is in normalized membrane potential units on the repolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_repol(bclind)*data.dt;
        stimstartsteps = i_repol(bclind);
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Yfp_traj Llu'])
%        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Llu'])
    elseif strfind(shiftstring,'depol') % in this case the offset is in normalized membrane potential units on the depolarization edge
        % Stimulus start time:
        % This is (simulation duration) – (timing between original stimulus time (or simulation start time) and new IC time).
        % (bcl + data.dt) – (offset-data.dt). I verified that this formula
        % produces matched values for BCL = 1000ms.
        stimstart = bcl + 2*data.dt - i_depol(bclind)*data.dt;
        stimstartsteps = i_depol(bclind);
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Yfp_traj Llu'])
%        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc stimstart measurementindex fp_for_fbk Llu'])
    else
        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc measurementindex fp_for_fbk Yfp_traj Llu'])
%        eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc measurementindex fp_for_fbk Llu'])
    end
    
    if pertflag == 1
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
    end
    

    if ~pertflag % perturb IC along eigenvector
        % Load fixed point matching selected bcl, and
        % add a perturbation of size epsln to the fixed point, in the direction
        % corresponding to the selected eigenvector
        initcond(:, i) = fp_found(:, bclind) + pertsize*direction;
    elseif pertflag == 1 % perturb along closed fixed-point trajectory
        % load fixed-point trajectory
%        eval(['load ' fpsimfolder simoutfname ' Y subdiv_per_cyc;']);
        eval(['load ' fixedpointfolder simoutfname 'subdiv_per_cyc;']);
%        Yfp_traj = Y;
        % Ask user where (along the closed fixed-point trajectory) to
        % choose the new IC
        pertsteps = input(['The closed fixed-point trajectory had a sampling interval \nof ' num2str(bcl/subdiv_per_cyc) ' ms, with ' num2str(subdiv_per_cyc*ncyc) ' total sample steps. To choose the new IC, \nenter a whole number corresponding to the number of \nsampling intervals to move along the trajectory: ']);
        initcond(:,i) = Yfp_traj(:,pertsteps+1); % The first column contains the IC, so displace one more step to get to the next set of values
    else % if pertflag==2 perturb along a selected state variable direction
        % or pertflag == 3 or 4
        initcond(:,i) = fp_found(:,bclind) + pertsize*sv_pert_dir; 
    end
    
    inputdiff(:, i) = initcond(:,i) - fp_found(:,bclind);
    
    % Compute error relative to IC:
    eval(['perterr = ' modelname '_p2pzero(initcond(:,i));']); % compute fixed-point error
    perterrnorm = norm(perterr)
    
    % Load trajectory produced by model
    eval(['load ' modelinputfolder simoutfname ' Y;']);
    
%    outputdiff(:, i) = Y(:,end) - fp_found(:,bclind);
    trajectorydiff = Y-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1));
%    trajectorydiffrel = (Y-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1)))./Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1));
    % To avoid div by zero (or close to it), try relative percent difference (not really a percent?):
    trajectorydiffrpd = 2*(Y-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1)))./(abs(Y)+abs(Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1))));
%    trajectorydiffinf = (Y-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1)))./max(abs(Y),abs(Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1))));

% compute magnitude ratio and angle between input and output
%    outin_mag_ratio(i) = norm(outputdiff(:,i))/norm(inputdiff(:,i))
    % if perturbation is along an eigenvector, the above ratio should match the
    % eigenvalue
    
%     costheta = dot(outputdiff(:,i),inputdiff(:,i))/(norm(outputdiff(:,i))*norm(inputdiff(:,i)));
%     outin_ang(i) = acosd(costheta)
    % if perturbation is along an eigenvector, the above angle should be close to zero
    
    % Copy model output to perturbation folder:
    eval(['movefile ' modelinputfolder simoutfname '.mat ' pertfolder simoutfname '.mat'])
    %h(i) = figure;
    
     numtimes = size(Y,2);
%     pertnorm{i} = zeros(1,numtimes);
%     for kk=1:numtimes
%         pertnorm{i}(kk) = norm(Y(:,kk)-fp_found(:,bclind));
%     end
%     figure
%     hold on;
%     plot((1:numtimes)*(bcl/subdiv_per_cyc), pertnorm{i})
%     plot(((0:ncyc-1)*subdiv_per_cyc +1)*bcl/subdiv_per_cyc,pertnorm{i}((0:ncyc-1)*subdiv_per_cyc +1),'ro')
    
    figure
    hold on;
    plot((1:numtimes)*(bcl/subdiv_per_cyc), Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1)))
    plot((1:numtimes)*(bcl/subdiv_per_cyc), Y,'--')
    xlabel('time, ms')
    ylabel('Reference and perturbed trajectories, original units')
    title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 

%     figure
%     hold on;
%     plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiff)
%     xlabel('time, ms')
%     ylabel('estimation errors, original units')

%     figure
%     hold on;
%     plot((1:numtimes)*(bcl/subdiv_per_cyc), Smat*trajectorydiff)
%     xlabel('time, ms')
%     ylabel('estimation errors, nondimensionalized')
%     title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 
    
    figure
    hold on;
    plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiffrpd)
    xlabel('time, ms')
    ylabel('RPD estimation errors')
    title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 

    errorendncyc = ncyc; % ending cycle index over which to compute error norm
    errorstartindex = round(bcl*(errorendncyc-5)/data.dt)+1; %1; % Enter a starting index for error norm computation
    errorendindex = round(bcl*errorendncyc/data.dt)+1;
    error_nondim_norm = norm(Smat*trajectorydiff(:,errorstartindex:errorendindex))
    error_nondim_mean = mean(mean(abs(Smat*trajectorydiff(:,errorstartindex:errorendindex))))
%    errorrel_nondim_norm = norm(trajectorydiffrel(:,errorstartindex:errorendindex)) 
    errorrpd_norm = norm(trajectorydiffrpd(:,errorstartindex:errorendindex)) 
    errorrpd_mean = mean(mean(abs(trajectorydiffrpd(:,errorstartindex:errorendindex)))) 
end

%Save the matrix to a separate file in the folder specified above
%eval(['save ' pertfolder 'perturbation_results *']);
if pertflag == 2
    % To save space, save everything except for Yfp 
%    eval(['save ' pertfolder 'nonl_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(ncyc) '_pfl' num2str(pertflag) '_pdir' num2str(sv_pert_ind) '_psz' num2str(round(log10(pertsize))) ' *']);
    eval(['save ' pertfolder 'nonl_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(ncyc) '_pfl' num2str(pertflag) '_pdir' num2str(sv_pert_ind) '_psz' num2str(round(log10(pertsize))) ' -regexp ^(?!(Yfp_traj)$).']);
elseif (pertflag == 3) || (pertflag == 4)
    eval(['save ' pertfolder 'nonl_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(ncyc) '_pfl' num2str(pertflag) '_pdir' sv_pert_ind '_psz' num2str(round(log10(pertsize))) ' -regexp ^(?!(Yfp_traj)$).']);
else
    disp('You should save the results manually.') 
end

if 0
    clear variables
    pertfolder = 'nonlinear_luen_tests/';%Save perturbation outputs here
%    loadfilename = 'nonl_luen_results_b500_m0_shift0p2Vnormrepol_ncyc20_pfl4_pdirrand_psz-5.mat';
%    loadfilename = 'nonl_luen_results_b500_m1_shift0p2Vnormrepol_ncyc20_pfl4_pdirrand_psz-5.mat';
%    loadfilename = 'nonl_luen_results_b500_m10_shift0p2Vnormrepol_ncyc20_pfl4_pdirrand_psz-5.mat';
    loadfilename = 'nonl_luen_results_b200_m1_shift0p2Vnormrepol_ncyc30_pfl4_pdirrand_psz-5.mat'; 
%    loadfilename = 'nonl_luen_results_b200_m0_shift0p2Vnormrepol_ncyc30_pfl4_pdirrand_psz-5.mat'; 
    eval(['load ' pertfolder loadfilename]);
 
%     figure
%     hold on;
%     plot((1:numtimes)*(bcl/subdiv_per_cyc), Smat*trajectorydiff)
%     xlabel('time, ms')
%     ylabel('estimation errors, nondimensionalized')
%     title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 
    if ~exist('trajectorydiffrpd','var')
        Yold = Y; 
        eval(['load ' fixedpointfolder simoutfname ' Y;']);
        Yfp_traj = Y; 
        trajectorydiffrpd = 2*(Yold-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1)))./(abs(Yold)+abs(Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1))));
    end
    if ~exist('trajectorydiff','var')
        Yold = Y; 
        eval(['load ' fixedpointfolder simoutfname ' Y;']);
        Yfp_traj = Y; 
        trajectorydiff = Yold-Yfp_traj(:,1:(round(ncyc*bcl/data.dt)+1));
    end

    figure
    hold on;
    plot((1:numtimes)*(bcl/subdiv_per_cyc), trajectorydiffrpd)
    xlabel('time, ms')
    ylabel('RPD estimation errors')
    title(['$T = $' num2str(bcl) 'ms, Meas = ' measlabel ', ' shiftstring(7:end) ', ' param],'Interpreter','latex'); 

    errorendncyc = ncyc; % ending cycle index over which to compute error norm
    errorstartindex = round(bcl*(errorendncyc-5)/data.dt)+1; %1; % Enter a starting index for error norm computation
    errorendindex = round(bcl*errorendncyc/data.dt)+1;
    error_nondim_norm = norm(Smat*trajectorydiff(:,errorstartindex:errorendindex))
    error_nondim_mean = mean(mean(abs(Smat*trajectorydiff(:,errorstartindex:errorendindex))))
    errorrpd_norm = norm(trajectorydiffrpd(:,errorstartindex:errorendindex)) 
    errorrpd_mean = mean(mean(abs(trajectorydiffrpd(:,errorstartindex:errorendindex)))) 
end

