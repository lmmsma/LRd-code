% Test out linear Kalman filter (KF) designs for selected 
% conditions, using LRd model Jacobians. Code based on lrd_linear_sim.m
% and compute_obsv_ctrb.m.
% Laura Munoz, July 2017

clear variables;

rng('default'); % reset random number generator

ms = 10; fs = 14; % marker size and font size

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
selected_logepsln = -5; % Base 10 log of step size used in Jacobian computation

kf_folder = 'kalman/';
%OCvalues = 'OCvalues/'; %folder where obsv and ctrb values will be saved

adj_yn = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if adj_yn
    param = 'adj';
else
    param = 'def';
end

jacfolder = ['jacobians/' param '/']; % folder where jacobians are stored

eval(['load ' jacfolder 'jacfile' num2str(selected_logepsln) ' *']); %Load Jacobians

numstate = size(alljacs{1},1); % number of state variables

% number of bcls
nbcls = length(selected_bcls_for_fps);

% Simulate the system for some fixed amount of time, and figure out the
% corresponding number of BCLs later. 
% What is a reasonable amount of time? How quickly do we want 
% the noise-free errors to converge? 
%simtime = 60000; % ms
simtime = 30000; % ms

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
display(['Measurement: ' statenames(measurementindex, :)])
% measurementindex = 1; % Select an integer from 1 to numstate. 
% This is the index of the state variable that is measured. 

controlindex = 1;  % Select an integer from 1 to numstate.  
% This is the index of the state variable that is controlled. 
% It should not affect KF design, unless process noise matrix 
% Bw is redesigned to make use of this index. 

% Initial condition
%x0 = 0.01*ones(numstate,1); 
% This should be replaced with a random IC (e.g., xhat system should have 
% random IC) after the xhat system is separated out. 
% Choose a random initial condition for "ground truth" system; 
% according to standard KF theory this should be a random variable. 
%xgtlin0 = 0.001*randn(numstate,1); % Random IC. 
xgtlin0 = 0.001*Smatinv*randn(numstate,1); % Random IC, but try to make proportional to state variable amplitudes.  
%Warning!!! Sizes probably don't make sense 
%(variables have different ranges, also 
% perturbations of size 1 prob. not in linear regime)

% Initialize state estimate (doesn't have to be a random variable) 
%xhat0 = 0.01*ones(numstate,1); 
%xhat0 = 0.01*Smatinv*ones(numstate,1); % Try to make proportional to state variable amplitudes
%xhat0 = 0.01*Smatinv*randn(numstate,1); % Try to make proportional to state variable amplitudes
xhat0 = 0.005*Smatinv*randn(numstate,1); % Try to make proportional to state variable amplitudes


% Pre-generate unit-variance noise signals so that different Q/R 
% ratios can betested with the same shape of w and v trajectories

%rng('default'); % Set random number seed to default for now, to make it 
% easier to compare across different runs? This didn't work the way I
% expected. 

%wu = randn(1,noisecycles); 
%wu = randn(numstate,noisecycles); 
wu = Smatinv*randn(numstate,noisecycles); % Try to make input noise proportional to state variable amplitudes
vu = randn(1,noisecycles); 

% Units of noise signals will depend on which variables are 
% measured and/or receive process noise. 

%Qnscalar = 0.1; % process noise factor 
%Qnscalar = 0.01; % process noise factor 
%Qnscalar = 0.000001; % process noise factor 
Qnscalar = 0.0000001; % process noise factor 
%Qnscalar = 100; % process noise factor 

Rn = 0.001*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)
%Rn = .01*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)

Qn = Qnscalar*eye(numstate); % process noise covariance 

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

% Set up a discrete-time state-space system for use with "kalman" function
sys = ss(jaccd, [B(:,controlindex) Bw], C(measurementindex,:), [0 zeros(1,size(Bw,2))], -1);

% Should the KFs be designed for scaled or unscaled systems? KF results for
% scaled systems should be more consistent with scaled observabilty
% predictions, but if we eventually develop any nonlinear observers, any
% normalization step may need to wait until after the design and 
% implementation stage (e.g., only for the purpose of plotting or averaging
% the estimation errors). 

% Set up a discrete-time state-space system for use with "kalman" function
% Use scaled matrices this time
%sys_scaled = ss(Smat*jaccd*Smatinv, [Bs(:,controlindex) Smat*Bw], Cs(measurementindex,:), [0 0], -1);
% I'm not doing anything with the scaled system (yet), but it may be useful
% for debugging purposes. 

% Compute KF gain (= Lkf) just for unscaled system for now
%[kest,Lkf,P] = kalman(sys,Qn,Rn,0);
[kest,Lkf,P] = kalman(sys,Qn,Rn,zeros(size(w,1),1));
%[kest,Lkf,P] = kalman(sys,Qn,Rn,zeros(size(w,1),1),'delayed'); % excluding the �delayed� seemed to have no impact on Lkf in the one test I ran on 11/19. 

% % For debugging purposes, can check computation of Lkf. 
% % Should get Lkf = Ltemp1' = Ltemp2. 
% [Ms,cleig1,Ltemp1] = dare(jaccd',C(measurementindex,:)',Qn,Rn);  
% Ltemp2 = jaccd*Ms*C(measurementindex,:)'/(C(measurementindex,:)*Ms*C(measurementindex,:)' + Rn); 

% Open-loop eigenvalues
oleig = eig(jaccd); 

% Closed-loop eigenvalues
cleig = eig(jaccd-Lkf*C(measurementindex,:)); 

% Display open and closed-loop eigenvalues
[oleig cleig];
figure;
t= {[param ': Eigenvalue mags for OL and CL'], ['BCL =' num2str(bcl) 'ms, Meas= ' statenames(measurementindex, :)]};
plot(1:17, [abs(oleig) abs(cleig)])
legend('Open Loop', 'Closed Loop');
grid on
title(t);
%set(gca,'XTick',1:17,'XTickLabel',statenames)
%xlabel('Measurement');
xlabel('Eigenvalue index');
ylabel('\lambda Magnitude');
%saveas(gcf,[kf_folder 'kfeigs_' param '_' statenames(measurementindex,:) '_' num2str(bcl)])
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

% Compute "ground truth" linear system state, without noise
xgtlinall_noiseless = zeros(numstate,nsteps);
xgtlinall_noiseless(:,1) = xgtlin0; 
for ii = 1:nsteps-1
    xgtlinall_noiseless(:,ii+1) = jaccd*xgtlinall_noiseless(:,ii); 
end
ygtlin_noiseless = C(measurementindex,:)*xgtlinall_noiseless; % measurement without noise

% Compute "ground truth" linear system state, without noise
xgtlinall = zeros(numstate,nsteps);
xgtlinall(:,1) = xgtlin0; 
for ii = 1:nsteps-1
    xgtlinall(:,ii+1) = jaccd*xgtlinall(:,ii)+ Bw*w(:,ii); 
end
ygtlin = C(measurementindex,:)*xgtlinall + v(ii); % measurement with noise

% Compute pseudo-KF (no noise) state estimate, open loop
% This is almost the same as xgtlin_noiseless, but with a non-random IC
xhatolall_noiseless = zeros(numstate,nsteps);
xhatolall_noiseless(:,1) = xhat0; 
for ii = 1:nsteps-1
    xhatolall_noiseless(:,ii+1) = jaccd*xhatolall_noiseless(:,ii); 
end

% Compute pseudo-KF (no noise) state estimate, closed loop
xkfall_noiseless = zeros(numstate,nsteps);
xkfall_noiseless(:,1) = xhat0; 
for ii = 1:nsteps-1
    xkfall_noiseless(:,ii+1) = jaccd*xkfall_noiseless(:,ii) + Lkf*(ygtlin_noiseless(ii)-C(measurementindex,:)*xkfall_noiseless(:,ii)); 
end

% Compute KF state estimate, closed loop
xkfall = zeros(numstate,nsteps);
xkfall(:,1) = xhat0; 
for ii = 1:nsteps-1
    xkfall(:,ii+1) = jaccd*xkfall(:,ii) + Lkf*(ygtlin(ii)-C(measurementindex,:)*xkfall(:,ii));  
end

% Compute estimation errors
eol_noiseless = xgtlinall_noiseless - xhatolall_noiseless; 
eol = xgtlinall - xhatolall_noiseless; 
ekf_noiseless = xgtlinall_noiseless - xkfall_noiseless; 
ekf = xgtlinall - xkfall; 

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
title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]});
subplot(2,1,2)
hold 
plot(bclselect*(1:nsteps)/1000,ekf_noiseless,'.-')
grid 
xlabel('time (sec)') 
ylabel('Error w/ fbk')
axis([0 simtime/1000 -ylim ylim])

%saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
%saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])

% Plot absolute errors with noise
figure
subplot(2,1,1)
hold
plot(bclselect*(1:nsteps)/1000,eol,'.-')
grid 
ylabel('Error w/o fbk')
title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]});
axis([0 simtime/1000 -ylimn ylimn])
subplot(2,1,2)
hold
plot(bclselect*(1:nsteps)/1000,ekf,'.-')
grid 
ylabel('Error w/ fbk')
xlabel('time (sec)') 
axis([0 simtime/1000 -ylimn ylimn])

%saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
%saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])

% Plot normalized noiseless errors 
figure
subplot(2,1,1)
hold
plot(bclselect*(1:nsteps)/1000,eol_noiseless./eolmax,'.-')
grid 
ylabel('Normalized err. w/o fbk')
axis([0 simtime/1000 -1.1 1.1])
title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]});
subplot(2,1,2)
hold 
plot(bclselect*(1:nsteps)/1000,ekf_noiseless./ekfmax,'.-')
grid 
xlabel('time (sec)') 
ylabel('Normalized err. w/ fbk')
axis([0 simtime/1000 -1.1 1.1])

%saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
%saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])

% Plot normalized errors, with noise
figure
subplot(2,1,1)
hold
plot(bclselect*(1:nsteps)/1000,eol./eolmaxn,'.-')
grid 
ylabel('Normalized err. w/o fbk')
title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]});
axis([0 simtime/1000 -1.1 1.1])
subplot(2,1,2)
hold
plot(bclselect*(1:nsteps)/1000,ekf./ekfmaxn,'.-')
grid 
ylabel('Normalized err. w/ fbk')
xlabel('time (sec)') 
axis([0 simtime/1000 -1.1 1.1])

%saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
%saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])

% Compute error covariance
ekfcovnorm = norm(cov(ekf))

%save save kfwv1to10 *


