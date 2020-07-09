% Test out linear Kalman filter (KF) designs for selected
% conditions, using LRd model Jacobians. Code based on lrd_linear_sim.m
% and compute_obsv_ctrb.m.
% Laura Munoz, July 2017

clear variables;

%shiftstring = '' % if using default shift of data.dt
%shiftstring = '_shift0p75ms'
%shiftstring = '_shift6p5ms'
%shiftstring = '_shift7mVrepol'
%shiftstring = '_shift-50mVrepol'
%shiftstrings = {''};
shiftstrings = {'_shift0p2Vnormrepol'}
%shiftstrings = {'_shift1Vnormdepol'}
%shiftstrings = {'_shift0p001Vnormrepol'}
%shiftstrings = {'','_shift0p75ms','_shift6p5ms','_shift7mVrepol','_shift-50mVrepol'};
%shiftstrings = {'','_shift1Vnormdepol','_shift0p2Vnormdepol','_shift0p4Vnormdepol','_shift0p6Vnormdepol','_shift0p8Vnormdepol','_shift0p2Vnormrepol','_shift0p4Vnormrepol','_shift0p6Vnormrepol','_shift0p8Vnormrepol','_shift0p001Vnormrepol'};


rng(1, 'twister'); % reset random number generator

markersize = 80; fontsize = 14; linewidth = 2; % marker size, font size, line width
shorttitleflag = 1; % set to 1 to reduce title length for presentation figures

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
% plot symbols
symbols = char('ko-','bs-','rp-','m*-','g^-','cv-','yh-');

selected_logepsln = -5; % Base 10 log of step size used in Jacobian computation

kf_folder = 'kalman/';
%OCvalues = 'OCvalues/'; %folder where obsv and ctrb values will be saved

adj_yn = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if adj_yn
    param = 'adj';
else
    param = 'def';
end

% Not all of the following statements need to be inside the loop.
for shiftctr = 1:length(shiftstrings)
    shiftstring = shiftstrings{shiftctr}
    
    %jacfolder = ['jacobians/' param '/']; % folder where jacobians are stored
    jacfolder = ['jacobians' shiftstring '/' param '/']; %This folder is where jacobians are stored
    
    eval(['load ' jacfolder 'jacfile' num2str(selected_logepsln) ' *']); %Load Jacobians
    
    numstate = size(alljacs{1},1); % number of state variables
    
    % number of bcls
    nbcls = length(selected_bcls_for_fps);
    
    % Simulate the system for some fixed amount of time, and figure out the
    % corresponding number of BCLs later.
    % What is a reasonable amount of time? How quickly do we want
    % the noise-free errors to converge?
    %simtime = 60000; % ms
    %simtime = 20000; % ms
    simtime = 10000000; % ms
    
    noisecycles = round(simtime/60); % number of cycles for which noise will be generated.
    % Pick a number greater than the number of cycles implied by simtime.
    % The above should work as long as BCL >= 60ms.
    
    % input matrix for all possible individual inputs
    B = eye(numstate); % This isn't necessarily right, just a placeholder. In general, B isn't square and its nonzero elements aren't necessarily 1.
    
    % output matrix for all possible individual measurements
    C = eye(numstate); % This is also a placeholder.
    
    % load approximate "state normalization" scaling matrix
    load b1000fsolem12variable_amplitudes varamp
    % Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
    Smat = diag(1./varamp);
    Smatinv = inv(Smat); % only need to compute once
    umax = diag(ones(1,size(B,2))); % This is also incorrect and just a placeholder.
    % For inputs, also need to know approximate maximum values of each element
    % of input vector. Note that u is a deviational quantity that may or
    % may not attain the size of the stimulus.
    % Each umax diagonal element should probably be the size of
    % deviational max for integrated system.
    Bs = Smat*B*umax; % scaled B matrix
    Cs = C*Smatinv; % scaled C matrix (could use Smat*C*Smatinv instead, if the output should also be normalized?)
    
    measurementindex=  input('Enter Measurement Index: ');
    disp(['Measurement: ' statenames(measurementindex, :)])
    % measurementindex = 1; % Select an integer from 1 to numstate.
    % This is the index of the state variable that is measured.
    
    % Reconstruction indices: these are the indices of the variables that
    % are judged to be the most important to reconstruct. These variables will
    % be highlighted in the plots.
    %reconstructionindices = [8 11 12]; % Calcium concentration variables are 8, 11, 12
    reconstructionindices = [8];
    
    controlindex = 1;  % Select an integer from 1 to numstate.
    % This is the index of the state variable that is controlled.
    % It should not affect KF design, unless process noise matrix
    % Bw is redesigned to make use of this index.
    
    % Initial condition
    %x0 = 0.01*ones(numstate,1);
    %Warning!!! Sizes probably don't make sense
    %(variables have different ranges, also xhat0 = 0.01*Smatinv*randn(numstate,1); % Try to make proportional to state variable amplitudes
    
    % perturbations of size 1 prob. not in linear regime)
    % This should be replaced with a random IC (e.g., xhat system should have
    % random IC) after the xhat system is separated out.
    % Choose a random initial condition for "ground truth" (=gt) system;
    % according to standard KF theory this should be a random variable.
    %xgtlin0 = 0.001*randn(numstate,1); % Random IC.
    %xgtlin0 = 0.001*Smatinv*randn(numstate,1); % Random IC, but try to make proportional to state variable amplitudes.
    xgtlin0 = 0.01*Smatinv*randn(numstate,1); % Random IC, but try to make proportional to state variable amplitudes.
    
    % Initialize state estimate (doesn't have to be a random variable)
    %xhat0 = xgtlin0; % Can match ICs to investigate the effect of noise-based
    %error in the absence of IC error
    %xhat0 = 0.01*ones(numstate,1);
    %xhat0 = 0.01*Smatinv*ones(numstate,1); % Try to make proportional to state variable amplitudes
    %xhat0 = 0.005*Smatinv*randn(numstate,1); % Try to make proportional to state variable amplitudes
    %xhat0 = 0.01*Smatinv*randn(numstate,1); % Try to make proportional to state variable amplitudes
    xhat0 = zeros(numstate,1); % If we don't know anything about the true IC, just initialize estiamtor at zero
    
    % Pre-generate unit-variance noise signals so that different Q/R
    % ratios can be tested with the same shape of w and v trajectories
    
    %wu = randn(1,noisecycles);
    %wu = randn(numstate,noisecycles);
    %vu = randn(1,noisecycles);
    % Later, I attempted to choose different sizes for the noise signals based
    % on the amplitudes of the corresponding variables. The noise size is reduced later via
    % Qn and Rn, when w and v are computed.
    wu = Smatinv*randn(numstate,noisecycles); % Try to make input noise proportional to state variable amplitudes
    vu = Smatinv(measurementindex,measurementindex)*randn(1,noisecycles);
    
    % Units of noise signals will depend on which variables are
    % measured and/or receive process noise.
    
    %Qnscalar = 0.1; % process noise factor
    %Qnscalar = 0.01; % process noise factor
    %Qnscalar = 0.000001; % process noise factor
    %Qnscalar = 100; % process noise factor
    %Qnscalar = 0.00000001; % process noise factor
%    Qnscalar = 0.0000001; % process noise factor
    %Qnscalar = 0.00001; % process noise factor
    Qnscalar = 1; % process noise factor
    
    %Rn = 100*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)
%    Rn = 0.001*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)
    Rn = 1e-14*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)
    
 %   Qn = Qnscalar*eye(numstate); % process noise covariance
    Qn = Qnscalar; % process noise covariance
    
    % For the 2013 study, I added the same noise signal to every channel, which
    % isn't particularly realistic. A better test could be to set Bw = eye and
    % add separate random process noise signals to each channel.
    %Bw = ones(numstate,1); % Process noise matrix
    Bw = eye(numstate); % Process noise matrix
    % For the unscaled system, may want to scale the noise signals (i.e.,
    % larger-sized variable should receive larger-sized noise signal)
    
    % Adjust noise signals to reflect covariances chosen above
    w=sqrt(Qn)*wu; % Process noise
    v=sqrt(Rn)*vu; % Measurement noise
    
    % Compute empirical noise variances
    covw_empirical = cov(w');
    covv_empirical = cov(v');
    
    % Select a BCL of interest. Trying to keep the number of BCLs low for now,
    % to reduce the number of plots
    bclindices = 1:nbcls;
    bclselect = input('Select BCL: '); % ms
    
    %%% Code after this line could be put inside a BCL loop, if more than one
    %%% BCL were selected
    
    % Find Jacobian index matching the BCL you chose
    bclselectindex = bclindices(selected_bcls_for_fps==bclselect);
    
    % Exit if you didn't find the right index
    if bclselect ~= selected_bcls_for_fps(bclselectindex)
        disp('Error: BCL index mismatch')
        return;
    end
    
    bcl = selected_bcls_for_fps(bclselectindex);
    % print current BCL to screen
    disp(['BCL = ' num2str(selected_bcls_for_fps(bclselectindex)) ' ms'])
    
    nsteps = round(simtime./bclselect); % number of cycles for simtime
    
    % Load Jacobian for selected BCL
    jaccd = alljacs{bclselectindex};
    % For debugging, produce a random (real-valued) Jacobian
    % tempmat = rand(17,17);
    % jaccd = tempmat + tempmat';
    % jaccd = jaccd/20;
    jaccd_scaled = Smat*jaccd*Smatinv;
    
    % Set up a discrete-time state-space system for use with "kalman" function
    %sys = ss(jaccd, [B(:,controlindex) Bw], C(measurementindex,:), [0 zeros(1,size(Bw,2))], -1);
    % Need to define Ts properly for use with freqsep:
    sys = ss(jaccd, [B(:,controlindex) Bw], C(measurementindex,:), [0 zeros(1,size(Bw,2))], bcl/1000);
    
    % Should the KFs be designed for scaled or unscaled systems? KF results for
    % scaled systems should be more consistent with scaled observabilty
    % predictions, but if we eventually develop any nonlinear observers, any
    % normalization step may need to wait until after the design and
    % implementation stage (e.g., only for the purpose of plotting or averaging
    % the estimation errors).
    
    % Set up a discrete-time state-space system for use with "kalman" function
    % Use scaled matrices this time
    %sys_scaled = ss(Smat*jaccd*Smatinv, [Bs(:,controlindex) Smat*Bw], Cs(measurementindex,:), [0 zeros(1,size(Bw,2))], -1);
    % Need to define Ts properly for use with freqsep:
    sys_scaled = ss(jaccd_scaled, [Bs(:,controlindex) Smat*Bw], Cs(measurementindex,:), [0 zeros(1,size(Bw,2))], bcl/1000);
    % I'm not doing anything with the scaled system (yet), but it may be useful
    % for debugging purposes.
    
    % Compute KF gain (= Lkf) just for unscaled system for now
    %[kest,Lkf,P] = kalman(sys,Qn,Rn,0);
    %[kest,Lkf,P] = kalman(sys,Qn,Rn,zeros(size(w,1),1),'delayed'); % excluding the �delayed� seemed to have no impact on Lkf in the one test I ran on 11/19.
%    [kest,Lkf,P,M,Z] = kalman(sys,Qn,Rn,zeros(size(w,1),1));
    [kest,Lkf,P,M,Z] = kalman(sys,Qn,Rn,zeros(size(Qn,1),1));
    
    % % For debugging purposes, can check computation of Lkf.
    % % Should get Lkf = Ltemp1' = Ltemp2.
    % [Ms,cleig1,Ltemp1] = dare(jaccd',C(measurementindex,:)',Qn,Rn);
    % Ltemp2 = jaccd*Ms*C(measurementindex,:)'/(C(measurementindex,:)*Ms*C(measurementindex,:)' + Rn);
    
%    [kestscaled,Lkfscaled,Pscaled,Mscaled,Zscaled] = kalman(sys_scaled,Qn,Rn,zeros(size(w,1),1));
    [kestscaled,Lkfscaled,Pscaled,Mscaled,Zscaled] = kalman(sys_scaled,Qn,Rn,zeros(size(Qn,1),1));
        
    % Open-loop eigenvalues
    oleig = eig(jaccd);
    
    % Closed-loop eigenvalues
    cleig = eig(jaccd-Lkf*C(measurementindex,:));
    
    % Should be the same as above
    cleigsc = eig(jaccd_scaled-Lkfscaled*Cs(measurementindex,:));
    
    % Display open and closed-loop eigenvalues
    [oleig sort(cleig,'descend','ComparisonMethod','abs') sort(cleigsc,'descend','ComparisonMethod','abs')]
    
    % Smat*L*Cs should be dimensionless:
    SmatLkfCsnorm = norm(Smat*Lkf*Cs(measurementindex,:))
    % This should be the same number as the previous one:
    LkfscCsnorm = norm(Lkfscaled*Cs(measurementindex,:))

    % Compute quantities for reduced-order scaled system
    %cutoff_freq = 5; % rad/s
    cutoff_freq_disc = .01;%.1; % eigenvalue modulus threshold
    cutoff_freq_cts = abs(log(cutoff_freq_disc))/(bcl/1000)
    %cutoff_freq_cts = abs(log(cutoff_freq_disc))/(2*pi*bcl/1000)
    [Gs_scaled, Gf_scaled] = freqsep(sys_scaled, cutoff_freq_cts);
    sort(abs(eig(Gs_scaled.A)))
    Gs_scaled.Ts
    Gs_scaleddim = length(Gs_scaled.C);
    [kestscaled_lf,Lkfscaled_lf,Pscaled_lf,Mscaled_lf,Zscaled_lf] = kalman(Gs_scaled,Qn(1:Gs_scaleddim,1:Gs_scaleddim),Rn);
    
    minsv_obsv_scaled_lf = min(svd(obsv(Gs_scaled.A,Gs_scaled.C)))
    Pscaledlfnorm = norm(Pscaled_lf)
    Zscaledlfnorm=norm(Zscaled_lf)
    
    
    figure;
    t= {[param ': Eigenvalue mags for OL and CL'], ['BCL =' num2str(bcl) 'ms, Meas= ' statenames_latex(measurementindex, :)]};
    plot(1:17, [abs(oleig) abs(cleig)])
    legend('Open Loop', 'Closed Loop');
    grid on
    title(t,'Interpreter','latex');
    %set(gca,'XTick',1:17,'XTickLabel',statenames)
    %xlabel('Measurement');
    xlabel('Eigenvalue index');
    ylabel('\lambda Magnitude');
    saveas(gcf,[kf_folder param '/kfeigs_' param '_' strtrim(statenames(measurementindex,:)) '_' num2str(bcl)])
    %saveas(gcf,[kf_folder 'kfeigs_' param '_' statenames(measurementindex, :) '_' num2str(bcl) '.jpeg'])
    
    
    % % open-loop system without noise (error dynamics only)
    % xall = zeros(numstate,nsteps);
    % xall(:,1) = x0;
    % for ii = 1:nsteps-1
    %     xall(:,ii+1) = jaccd*xall(:,ii);
    % end
    %
    % % open-loop system with noise (error dynamics only)
    % xnall = zeros(numstate,nsteps);
    % xnall(:,1) = x0;
    % for ii = 1:nsteps-1
    %     xnall(:,ii+1) = jaccd*xnall(:,ii) + Bw*w(:,ii);
    % end
    %
    % % pseudo-KF (no noise) error dynamics
    % xkfall = zeros(numstate,nsteps);
    % xkfall(:,1) = x0;
    % for ii = 1:nsteps-1
    %     xkfall(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*xkfall(:,ii);
    % end
    %
    % % KF with noise (error dynamics only)
    % xkfnall = zeros(numstate,nsteps);
    % xkfnall(:,1) = x0;
    % for ii = 1:nsteps-1
    %     xkfnall(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*xkfnall(:,ii) + Bw*w(:,ii) - Lkf*v(ii);
    % end
    
    % % Compute "ground truth" linear system state, without noise
    % xgtlinall_noiseless = zeros(numstate,nsteps);
    % xgtlinall_noiseless(:,1) = xgtlin0;
    % for ii = 1:nsteps-1
    %     xgtlinall_noiseless(:,ii+1) = jaccd*xgtlinall_noiseless(:,ii);
    % end
    % ygtlin_noiseless = C(measurementindex,:)*xgtlinall_noiseless; % measurement without noise
    
    % % Compute "ground truth" linear system state, with noise
    % xgtlinall = zeros(numstate,nsteps);
    % xgtlinall(:,1) = xgtlin0;
    % for ii = 1:nsteps-1
    %     xgtlinall(:,ii+1) = jaccd*xgtlinall(:,ii)+ Bw*w(:,ii);
    % end
    %%ygtlin = C(measurementindex,:)*xgtlinall + v(ii); % measurement with noise
    %%(incorrect indexing of v)
    %ygtlin = C(measurementindex,:)*xgtlinall + v(1:nsteps); % measurement with noise
    
    % % Compute pseudo-KF (no noise) state estimate, open loop
    % % This is almost the same as xgtlin_noiseless, but with the option of a non-random IC
    % xhatolall_noiseless = zeros(numstate,nsteps);
    % xhatolall_noiseless(:,1) = xhat0;
    % for ii = 1:nsteps-1
    %     xhatolall_noiseless(:,ii+1) = jaccd*xhatolall_noiseless(:,ii);
    % end
    
    % Compute "ground truth" linear system state, without noise
    xgtlinall_noiseless = zeros(numstate,nsteps);
    xgtlinall_noiseless(:,1) = xgtlin0;
    ygtlin_noiseless = zeros(length(measurementindex),nsteps);
    ygtlin_noiseless(1) = C(measurementindex,:)*xgtlinall_noiseless(:,1); % measurement without noise
    
    % Compute "ground truth" linear system state, with noise
    xgtlinall = zeros(numstate,nsteps);
    xgtlinall(:,1) = xgtlin0;
    ygtlin = zeros(length(measurementindex),nsteps);
    ygtlin(1) = C(measurementindex,:)*xgtlinall(:,1); % measurement with noise
    
    % Compute pseudo-KF (no noise) state estimate, closed loop
    xkfall_noiseless = zeros(numstate,nsteps);
    xkfall_noiseless(:,1) = xhat0;
    xkfapostall_noiseless = zeros(numstate,nsteps);
    xkfapostall_noiseless(:,1) = xkfall_noiseless(:,1) + M*(ygtlin_noiseless(1)-C(measurementindex,:)*xkfall_noiseless(:,1));
    
    % Compute pseudo-KF (no noise) state estimate, open loop
    % This is almost the same as xgtlin_noiseless, but with the option of a non-random IC
    xhatolall_noiseless = zeros(numstate,nsteps);
    xhatolall_noiseless(:,1) = xhat0;
    
    % Compute KF state estimate, closed loop. xkfall(:, ii+1) = xhat(ii+1|ii),
    % in other words, it is the a-priori predictor-corrector error.
    xkfall = zeros(numstate,nsteps);
    xkfall(:,1) = xhat0;
    xkfapostall = zeros(numstate,nsteps);
    xkfapostall(:,1) = xkfall(:,1) + M*(ygtlin(1)-C(measurementindex,:)*xkfall(:,1));
    
    err_apri_directsim_noiseless = zeros(numstate,nsteps);
    err_apri_directsim_noiseless(:,1) = xgtlinall_noiseless(:,1)-xkfall_noiseless(:,1);
    err_apri_directsim = zeros(numstate,nsteps);
    err_apri_directsim(:,1) = xgtlinall(:,1)-xkfall(:,1);
    
    
    for ii = 1:nsteps-1
        err_apri_directsim_noiseless(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*err_apri_directsim_noiseless(:,ii);
        err_apri_directsim(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*err_apri_directsim(:,ii) + Bw*w(:,ii) -Lkf*v(ii);
        
        xgtlinall_noiseless(:,ii+1) = jaccd*xgtlinall_noiseless(:,ii);
        
        ygtlin_noiseless(ii+1) = C(measurementindex,:)*xgtlinall_noiseless(:,ii+1); % measurement without noise
        
        
        xgtlinall(:,ii+1) = jaccd*xgtlinall(:,ii)+ Bw*w(:,ii);
        
        ygtlin(ii+1) = C(measurementindex,:)*xgtlinall(:,ii+1) + v(ii+1); % measurement with noise
        
        
        xhatolall_noiseless(:,ii+1) = jaccd*xhatolall_noiseless(:,ii);
        
        %    xkfall_noiseless(:,ii+1) = jaccd*xkfall_noiseless(:,ii) + Lkf*(ygtlin_noiseless(ii)-C(measurementindex,:)*xkfall_noiseless(:,ii));
        xkfall_noiseless(:,ii+1) = jaccd*xkfapostall_noiseless(:,ii);
        
        xkfapostall_noiseless(:,ii+1) = xkfall_noiseless(:,ii+1) + M*(ygtlin_noiseless(ii+1)-C(measurementindex,:)*xkfall_noiseless(:,ii+1));
        
        
        %    xkfall(:,ii+1) = jaccd*xkfall(:,ii) + Lkf*(ygtlin(ii)-C(measurementindex,:)*xkfall(:,ii));
        xkfall(:,ii+1) = jaccd*xkfapostall(:,ii);
        
        xkfapostall(:,ii+1) = xkfall(:,ii+1) + M*(ygtlin(ii+1)-C(measurementindex,:)*xkfall(:,ii+1));
        
    end
    %ygtlin = C(measurementindex,:)*xgtlinall + v(1:nsteps); % measurement with noise
    
    % Compute estimation errors
    eol_noiseless = xgtlinall_noiseless - xhatolall_noiseless;
    eol = xgtlinall - xhatolall_noiseless;
    ekf_noiseless = xgtlinall_noiseless - xkfall_noiseless;
    ekfapost_noiseless = xgtlinall_noiseless - xkfapostall_noiseless;
    ekf = xgtlinall - xkfall;
    ekfapost = xgtlinall - xkfapostall;
    
    
    Lkfnorm = norm(Lkf) % size of gain is related to strength of feedback
    Pnorm = norm(P) % steady-state a-priori cov norm
    Znorm = norm(Z) % steady-state a-posteriori cov norm
    
    Lkfscalednorm = norm(Lkfscaled) % size of gain is related to strength of feedback
    Pscalednorm = norm(Pscaled) % steady-state a-priori cov norm
    Zscalednorm = norm(Zscaled) % steady-state a-posteriori cov norm
    
    maxabscleig = max(abs(cleig)) % size of least stable eigenvalue
    
    % Compute error covariance norms
    ekfcov = cov(ekf');
    size(ekfcov) % Make sure this is nxn, where n is the number of states.
    % The presence or absence of the transpose operator matters. Cov assumes
    % that columns are variables and rows are observations.
    eolcovnorm_noiseless = norm(cov(eol_noiseless'));
    ekfcovnorm_noiseless = norm(cov(ekf_noiseless'));
    ekfapostcovnorm_noiseless = norm(cov(ekfapost_noiseless'));
    eolcovnorm = norm(cov(eol'));
    ekfcovnorm = norm(ekfcov);
    ekfapostcovnorm = norm(cov(ekfapost'));
    
    % Check a value against manual computation
    ekfvarcheck = sum((ekf(reconstructionindices(1),:)-mean(ekf(reconstructionindices(1),:))).^2)/(nsteps-1)
    ekfcov(reconstructionindices(1),reconstructionindices(1))
    
    % Exclude transient period from covariances?
    % No obvious fixed time when this should be; maybe first ten cycles?
    %postransindex = 11;
    % Revise this to use last 5% of simulation?
    % eolcovnorm_noiseless_ss = norm(cov(eol_noiseless(:,postransindex:end)'));
    % ekfcovnorm_noiseless_ss = norm(cov(ekf_noiseless(:,postransindex:end)'));
    % ekfapostcovnorm_noiseless_ss = norm(cov(ekfapost_noiseless(:,postransindex:end)'));
    % eolcovnorm_ss = norm(cov(eol(:,postransindex:end)'));
    % ekfcovnorm_ss = norm(cov(ekf(:,postransindex:end)'));
    % ekfapostcovnorm_ss = norm(cov(ekfapost(:,postransindex:end)'));
    index_last5pct = round(0.05*nsteps);
    index_last_n = 10;
    
    eolcovnorm_noiseless_ss = norm(cov(eol_noiseless(:,end-index_last_n:end)'));
    ekfcovnorm_noiseless_ss = norm(cov(ekf_noiseless(:,end-index_last_n:end)'));
    ekfapostcovnorm_noiseless_ss = norm(cov(ekfapost_noiseless(:,end-index_last_n:end)'));
    eolcovnorm_ss = norm(cov(eol(:,end-index_last_n:end)'));
    ekfcovnorm_ss = norm(cov(ekf(:,end-index_last_n:end)'));
    ekfapostcovnorm_ss = norm(cov(ekfapost(:,end-index_last_n:end)'));
    
    
    % Compute scaled error covariance norms
    eolcovnorm_noiseless_scaled = norm(cov((Smat*eol_noiseless)'));
    ekfcovnorm_noiseless_scaled = norm(cov((Smat*ekf_noiseless)'));
    ekfapostcovnorm_noiseless_scaled = norm(cov((Smat*ekfapost_noiseless)'));
    eolcovnorm_scaled = norm(cov((Smat*eol)'));
    ekfcovnorm_scaled = norm(cov((Smat*ekf)'));
    ekfapostcovnorm_scaled = norm(cov((Smat*ekfapost)'));
    
    eolcovnorm_noiseless_ss_scaled = norm(cov((Smat*eol_noiseless(:,end-index_last_n:end))'));
    ekfcovnorm_noiseless_ss_scaled = norm(cov((Smat*ekf_noiseless(:,end-index_last_n:end))'));
    ekfapostcovnorm_noiseless_ss_scaled = norm(cov((Smat*ekfapost_noiseless(:,end-index_last_n:end))'));
    eolcovnorm_ss_scaled = norm(cov((Smat*eol(:,end-index_last_n:end))'));
    ekfcovnorm_ss_scaled = norm(cov((Smat*ekf(:,end-index_last_n:end))'));
    ekfapostcovnorm_ss_scaled = norm(cov((Smat*ekfapost(:,end-index_last_n:end))'));
    
    
    eolcovnorm
    ekfcovnorm
    ekfapostcovnorm
    ekfcovnorm_noiseless
    ekfapostcovnorm_noiseless
    ekfcovnorm_noiseless_ss
    ekfapostcovnorm_noiseless_ss
    ekfcovnorm_ss
    ekfapostcovnorm_ss
    
    eolcovnorm_scaled
    ekfcovnorm_scaled
    ekfapostcovnorm_scaled
    eolcovnorm_ss_scaled
    ekfcovnorm_ss_scaled
    ekfapostcovnorm_ss_scaled
    
    
    % Compare errors that were computed in different ways
    figure
    hold on;
    plot(1:nsteps, err_apri_directsim(1,:));
    plot(1:nsteps, ekf(1,:),'m:x');
    ylabel('V est error')
    xlabel('time index')
    legend('e(j+1) = (A-LC)e(j) + Bw*w(j) - Lv(j)', 'e(j) = xtrue(j) - xhat(j|j-1)')
    
    figure
    hold on;
    plot(1:nsteps, err_apri_directsim_noiseless(1,:));
    plot(1:nsteps, ekf_noiseless(1,:),'m:x');
    ylabel('V est error')
    xlabel('time index')
    legend('e(j+1) = (A-LC)e(j)', 'e(j) = xtrue noiseless(j) - xhat noiseless(j|j-1)')
    
    % Compute scaled errors (may be better for plotting, to prevent larger variables
    %from dominating the plot?)
    % xsnall = Smat*xnall;
    % xskfall = Smat*xkfall;
    % xskfnall = Smat*xkfnall;
    
    % Maximum absolute errors, per state variable
    eolmax = max(abs(eol_noiseless),[],2); % max open-loop error
    eolmaxn = max(abs(eol),[],2); % max open-loop error, with noise
    ekfmax = max(abs(ekf_noiseless),[],2); % max closed-loop error
    ekfmaxn = max(abs(ekf),[],2); % max closed-loop error, with noise
    ekfapostmax = max(abs(ekfapost_noiseless),[],2); % max closed-loop a-post error
    ekfapostmaxn = max(abs(ekfapost),[],2); % max a-post closed-loop error, with noise
    
    % Choose a y-axis limit for error plots
    if max(abs(oleig)) <= 1 % Open-loop system is (arguably) at least marginally stable
        ylim = max(eolmax); % choose max open-loop error
        ylimn = max(eolmaxn); % choose max open-loop error, with noise
    else % OL system is unstable
        ylim = max(ekfmax); % choose max closed-loop error
        ylimn = max(ekfmaxn); % choose max closed-loop error, with noise
    end
    
    % Plot absolute noiseless errors
    figure
    subplot(2,1,1)
    hold
    plot(bclselect*(1:nsteps)/1000,eol_noiseless,'.-')
    grid
    ylabel('Error w/o fbk')
    axis([0 simtime/1000 -ylim ylim])
    title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]},'Interpreter', 'latex');
    subplot(2,1,2)
    hold
    plot(bclselect*(1:nsteps)/1000,ekf_noiseless,'.-')
    plot(bclselect*(1:nsteps)/1000,ekfapost_noiseless,'m--')
    grid
    legend('a-priori','a-posteriori')
    xlabel('time (sec)')
    ylabel('Error w/ fbk')
    axis([0 simtime/1000 -ylim ylim])
    
    saveas(gcf,[kf_folder param '/kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    %saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
    
    % Plot absolute errors with noise
    figure
    subplot(2,1,1)
    hold
    plot(bclselect*(1:nsteps)/1000,eol,'.-')
    grid
    ylabel('Error w/o fbk')
    title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex') ;
    axis([0 simtime/1000 -ylimn ylimn])
    subplot(2,1,2)
    hold
    plot(bclselect*(1:nsteps)/1000,ekf,'.-')
    plot(bclselect*(1:nsteps)/1000,ekfapost,'m--')
    legend('a-priori','a-posteriori')
    grid
    ylabel('Error w/ fbk')
    xlabel('time (sec)')
    axis([0 simtime/1000 -ylimn ylimn])
    
    saveas(gcf,[kf_folder param '/kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    %saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
    
    % Plot normalized noiseless errors
    figure
    subplot(2,1,1)
    hold
    plot(bclselect*(1:nsteps)/1000,eol_noiseless./eolmax,'.-')
    grid
    ylabel('Normalized err. w/o fbk')
    axis([0 simtime/1000 -1.1 1.1])
    title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
    subplot(2,1,2)
    hold
    plot(bclselect*(1:nsteps)/1000,ekf_noiseless./ekfmax,'.-')
    grid
    xlabel('time (sec)')
    ylabel('Normalized err. w/ fbk')
    axis([0 simtime/1000 -1.1 1.1])
    
    saveas(gcf,[kf_folder param '/kfplot_wonoise_normalized_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    %saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
    
    % Plot normalized errors, with noise
    figure
    subplot(2,1,1)
    hold
    plot(bclselect*(1:nsteps)/1000,eol./eolmaxn,'.-')
    grid
    ylabel('Normalized err. w/o fbk')
    title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
    axis([0 simtime/1000 -1.1 1.1])
    subplot(2,1,2)
    hold
    plot(bclselect*(1:nsteps)/1000,ekf./ekfmaxn,'.-')
    grid
    ylabel('Normalized err. w/ fbk')
    xlabel('time (sec)')
    axis([0 simtime/1000 -1.1 1.1])
    
    saveas(gcf,[kf_folder param '/kfplot_noise_normalized_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    %saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
    
    % Plot absolute errors of variables targeted for "reconstruction"
    figure
    subplot(2,1,1)
    hold
    for k = 1:length(reconstructionindices)
        plot(bclselect*(1:nsteps)/1000,eol(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
    end
    grid
    %ylabel('Error w/o fbk')
    ylabel('Error w/o fbk, mmol / L') % Use if concentrations were chosen
    title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
    leg = legend(statenames_latex(reconstructionindices,:));
    set(leg,'Interpreter', 'latex')
    %axis([0 simtime/1000 -ylimn ylimn])
    subplot(2,1,2)
    hold
    for k = 1:length(reconstructionindices)
        plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
        plot(bclselect*(1:nsteps)/1000,ekfapost(reconstructionindices(k),:),[symbols(k,1) ':'], 'Linewidth', linewidth/2)
    end
    grid
    %ylabel('Error w/ fbk')
    ylabel('A-pri (-) and a-post (--) errors w/ fbk, mmol / L') % Use if concentrations were chosen
    xlabel('time (sec)')
    leg = legend(statenames_latex(reconstructionindices,:));
    set(leg,'Interpreter', 'latex')
    %axis([0 simtime/1000 -ylimn ylimn])
    saveas(gcf,[kf_folder param '/kfplot_targetedvariables_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    
    
    % Plot absolute CL errors of variables targeted for "reconstruction"
    figure
    if shorttitleflag
        title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
    else
        title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
    end
    hold
    for k = 1:length(reconstructionindices)
        plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
        plot(bclselect*(1:nsteps)/1000,ekfapost(reconstructionindices(k),:),[symbols(k,1) ':'], 'Linewidth', linewidth/2)
    end
    grid
    %ylabel('Error w/ fbk')
    %ylabel('Estimation error, mmol / L') % Use if concentrations were chosen
    ylabel('A-pri (-) and a-post (--) errors w/ fbk, mmol / L') % Use if concentrations were chosen
    xlabel('time (sec)')
    leg = legend(statenames_latex(reconstructionindices,:));
    set(leg,'Interpreter', 'latex')
    set(gca,'fontsize',fontsize)
    set(gca,'linewidth',linewidth)
    %axis([0 simtime/1000 -ylimn ylimn])
    %saveas(gcf,[kf_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    saveas(gcf,[kf_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    
    % Plot normalized CL errors of variables targeted for "reconstruction"
    figure
    if shorttitleflag
        title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
    else
        title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
    end
    hold
    for k = 1:length(reconstructionindices)
        plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:)./ekfmaxn(reconstructionindices(k)),symbols(k,:), 'Linewidth', linewidth)
    end
    grid
    ylabel('Normalized est. err.')
    xlabel('time (sec)')
    leg = legend(statenames_latex(reconstructionindices,:));
    set(leg,'Interpreter', 'latex')
    set(gca,'fontsize',fontsize)
    set(gca,'linewidth',linewidth)
    %axis([0 simtime/1000 -ylimn ylimn])
    %saveas(gcf,[kf_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    saveas(gcf,[kf_folder param '/kfplot_targetedvariables_normederr_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    
    % Plot normalized errors, with noise, and highlight first targeted variable
    figure
    if shorttitleflag
        title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
    else
        title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
    end
    hold
    plot(bclselect*(1:nsteps)/1000,ekf./ekfmaxn,'.-')
    grid
    ylabel('Normalized est. err.')
    xlabel('time (sec)')
    axis([0 simtime/1000 -1.1 1.1])
    for k = 1:1
        p=plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:)./ekfmaxn(reconstructionindices(k)),symbols(k,:), 'Linewidth', 2*linewidth);
        leg = legend(p,statenames_latex(reconstructionindices(k),:));
    end
    grid
    set(leg,'Interpreter', 'latex')
    set(gca,'fontsize',fontsize)
    set(gca,'linewidth',linewidth)
    saveas(gcf,[kf_folder param '/kfplot_normederr_overlay_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    
    % Can use next portion if no Lkf was found. Compute eol first, then compute
    % associated OL norms and save command.
    if 0
        % Plot absolute CL errors of variables targeted for "reconstruction"
        figure
        if shorttitleflag
            title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
        else
            title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
        end
        
        hold
        for k = 1:length(reconstructionindices)
            plot(bclselect*(1:nsteps)/1000,eol(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
        end
        grid
        ylabel('Estimation error, mmol / L') % Use if concentrations were chosen
        xlabel('time (sec)')
        leg = legend(statenames_latex(reconstructionindices,:));
        set(leg,'Interpreter', 'latex')
        set(gca,'fontsize',fontsize)
        set(gca,'linewidth',linewidth)
        saveas(gcf,[kf_folder param '/kfplot_targetedvariables_ol_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    end
    
    %save kfwv1to10 *
    eval(['save ' kf_folder param '/kalmanfile_measind_' num2str(measurementindex) '_bcl_' num2str(bcl) ' *'])
    
end
