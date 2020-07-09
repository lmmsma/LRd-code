% Test out LQR designs for selected
% conditions, using LRd model Jacobians. Code based on test_linear_kf.m.
% Laura Munoz, August 2018

clear variables;

rng(1, 'twister'); % reset random number generator

markersize = 80; fontsize = 14; linewidth = 2; % marker size, font size, line width
shorttitleflag = 1; % set to 1 to reduce title length for presentation figures

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
% plot symbols
symbols = char('ko-','bs-','rp-','m*-','g^-','cv-','yh-');

selected_logepsln = -5; % Base 10 log of step size used in Jacobian computation

%kf_folder = 'kalman/';

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
simtime = 20000; % ms
%simtime = 10000000; % ms

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

controlinputindex=  input('Enter Control Input Index: ');
display(['Input applied to: ' statenames(controlinputindex, :)])
% This is the index of the dynamical equation to which the input will be
% applied. 

%measurementindex = 1; % This shouldn't affect LQR design, but is included to produce a SISO system

% Initial condition
x0 = 0.01*Smatinv*randn(numstate,1); % Random IC, but try to make proportional to state variable amplitudes.

% Qscalar = 100; % weighting factor for Q matrix
% 
% R = 0.01*Qscalar;% penalty matrix (or scalar) for input
% 
% Q = Qscalar*eye(numstate); % penalty matrix for state vector

Q = diag(1./varamp.^2); % penalty matrix for state vector
R = 0.01*norm(Q);% penalty matrix (or scalar) for input. 
% R should be more like 1/||umax,i||^2, where i is the index of the input
% channel

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

% Set up a discrete-time state-space system for use with "dlqr" function
%sys = ss(jaccd, B(:,controlindex), C(measurementindex,:), 0, -1);

% Set up a discrete-time state-space system for use with "dlqr" function
% Use scaled matrices this time
%sys_scaled = ss(Smat*jaccd*Smatinv, Bs(:,controlindex), Cs(measurementindex,:), 0, -1);

% Compute LQR gain (= Klq) for unscaled system 
[Klq,Sinf,cleiglq] = dlqr(jaccd,B(:,controlinputindex),Q,R);

% % For debugging purposes, can check computation of Klq (use dare).
% % Should get Klq = Ktemp1.
% [???,cleig1,Ktemp1] = dare(jaccd,B(:,controlindex),Q,R);

[Klqscaled,Sinfscaled,cleiglqscaled] = dlqr(Smat*jaccd*Smatinv, Bs(:,controlinputindex),Q,R); 

% Open-loop eigenvalues
oleig = eig(jaccd);

% Closed-loop eigenvalues
cleig = eig(jaccd - B(:,controlinputindex)*Klq);

% Display open and closed-loop eigenvalues
[oleig cleig];
figure;
t= {[param ': Eigenvalue mags for OL and CL'], ['BCL =' num2str(bcl) 'ms, Input appl. to ' statenames_latex(controlinputindex, :)]};
plot(1:17, [abs(oleig) abs(cleig)])
legend('Open Loop', 'Closed Loop');
grid on
title(t,'Interpreter','latex');
%set(gca,'XTick',1:17,'XTickLabel',statenames)
%xlabel('Measurement');
xlabel('Eigenvalue index');
ylabel('\lambda Magnitude');
%saveas(gcf,[kf_folder param '/kfeigs_' param '_' strtrim(statenames(measurementindex,:)) '_' num2str(bcl)])


% open-loop system (error dynamics only)
err_ol_directsim = zeros(numstate,nsteps);
err_ol_directsim (:,1) = x0;
Jsimol = 0; % "empirical" cost function
for ii = 1:nsteps-1
    err_ol_directsim(:,ii+1) = jaccd*err_ol_directsim(:,ii);
    Jsimol = Jsimol + err_ol_directsim(:,ii)'*Q*err_ol_directsim(:,ii); 
end

% closed-loop system (error dynamics only) 
err_cl_directsim = zeros(numstate,nsteps);
err_cl_directsim (:,1) = x0;
ulq = zeros(1,nsteps); 
Jsim = 0; % "empirical" cost function
for ii = 1:nsteps-1
    ulq(ii) = -Klq*err_cl_directsim(:,ii); 
    err_cl_directsim(:,ii+1) = jaccd*err_cl_directsim(:,ii)+ B(:,controlinputindex)*ulq(ii);
    Jsim = Jsim + err_cl_directsim(:,ii)'*Q*err_cl_directsim(:,ii) + ulq(ii)'*R*ulq(ii); 
end


Klqnorm = norm(Klq) % size of gain is related to strength of feedback

Klqscalednorm = norm(Klqscaled) % size of gain is related to strength of feedback

maxabscleig = max(abs(cleig)) % size of least stable eigenvalue

Sinfnorm = norm(Sinf)

Sinfscalednorm = norm(Sinfscaled)

Jsimol

Jsim

% % Compare errors
figure
hold on;
p1=plot(bclselect*(1:nsteps)/1000, err_ol_directsim, 'g--');
p2=plot(bclselect*(1:nsteps)/1000, err_cl_directsim,'b-');
ylabel('V error')
xlabel('time (sec)')
legend([p1(1) p2(1)], 'OL', 'CL')

figure
subplot(2,1,1)
hold on;
%plot(bclselect*(1:nsteps), err_ol_directsim(1,:));
%plot(bclselect*(1:nsteps), err_cl_directsim(1,:),'g:x');
plot(bclselect*(1:nsteps)/1000, err_ol_directsim(1,:), 'g--');
plot(bclselect*(1:nsteps)/1000, err_cl_directsim(1,:),'b-');
ylabel('V-V_0 (mV)')
leg = legend('u \equiv 0', 'u = -Kx','Interpreter','latex');
title({['Simulated x values, input appl. to ' strtrim(statenames_latex(controlinputindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(norm(Q)), ', R = ' num2str(R)]},'Interpreter', 'latex');
subplot(2,1,2)
hold on;
%plot(bclselect*(1:nsteps), err_ol_directsim(1,:));
%plot(bclselect*(1:nsteps), err_cl_directsim(1,:),'g:x');
plot(bclselect*(1:nsteps)/1000, err_ol_directsim(8,:), 'g--');
plot(bclselect*(1:nsteps)/1000, err_cl_directsim(8,:),'b-');
ylabel('[Ca^{2+}]_{i,t}-[Ca^{2+}]_{i,t,0} (mmol/L)')
xlabel('time (sec)')
% 
% figure
% hold on;
% plot(1:nsteps, err_apri_directsim_noiseless(1,:));
% plot(1:nsteps, ekf_noiseless(1,:),'m:x');
% ylabel('V est error')
% xlabel('time index')
% legend('e(j+1) = (A-LC)e(j)', 'e(j) = xtrue noiseless(j) - xhat noiseless(j|j-1)')
% 
% % Compute scaled errors (may be better for plotting, to prevent larger variables
% %from dominating the plot?)
% % xsnall = Smat*xnall;
% % xskfall = Smat*xkfall;
% % xskfnall = Smat*xkfnall;
% 
% % Maximum absolute errors, per state variable
% eolmax = max(abs(eol_noiseless),[],2); % max open-loop error
% eolmaxn = max(abs(eol),[],2); % max open-loop error, with noise
% ekfmax = max(abs(ekf_noiseless),[],2); % max closed-loop error
% ekfmaxn = max(abs(ekf),[],2); % max closed-loop error, with noise
% ekfapostmax = max(abs(ekfapost_noiseless),[],2); % max closed-loop a-post error
% ekfapostmaxn = max(abs(ekfapost),[],2); % max a-post closed-loop error, with noise
% 
% % Choose a y-axis limit for error plots
% if max(abs(oleig)) <= 1 % Open-loop system is (arguably) at least marginally stable
%     ylim = max(eolmax); % choose max open-loop error
%     ylimn = max(eolmaxn); % choose max open-loop error, with noise
% else % OL system is unstable
%     ylim = max(ekfmax); % choose max closed-loop error
%     ylimn = max(ekfmaxn); % choose max closed-loop error, with noise
% end
% 
% % Plot absolute noiseless errors
% figure
% subplot(2,1,1)
% hold
% plot(bclselect*(1:nsteps)/1000,eol_noiseless,'.-')
% grid
% ylabel('Error w/o fbk')
% axis([0 simtime/1000 -ylim ylim])
% title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]},'Interpreter', 'latex');
% subplot(2,1,2)
% hold
% plot(bclselect*(1:nsteps)/1000,ekf_noiseless,'.-')
% plot(bclselect*(1:nsteps)/1000,ekfapost_noiseless,'m--')
% grid
% legend('a-priori','a-posteriori')
% xlabel('time (sec)')
% ylabel('Error w/ fbk')
% axis([0 simtime/1000 -ylim ylim])
% 
% saveas(gcf,[kf_folder param '/kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% %saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
% 
% % Plot absolute errors with noise
% figure
% subplot(2,1,1)
% hold
% plot(bclselect*(1:nsteps)/1000,eol,'.-')
% grid
% ylabel('Error w/o fbk')
% title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex') ;
% axis([0 simtime/1000 -ylimn ylimn])
% subplot(2,1,2)
% hold
% plot(bclselect*(1:nsteps)/1000,ekf,'.-')
% plot(bclselect*(1:nsteps)/1000,ekfapost,'m--')
% legend('a-priori','a-posteriori')
% grid
% ylabel('Error w/ fbk')
% xlabel('time (sec)')
% axis([0 simtime/1000 -ylimn ylimn])
% 
% saveas(gcf,[kf_folder param '/kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% %saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
% 
% % Plot normalized noiseless errors
% figure
% subplot(2,1,1)
% hold
% plot(bclselect*(1:nsteps)/1000,eol_noiseless./eolmax,'.-')
% grid
% ylabel('Normalized err. w/o fbk')
% axis([0 simtime/1000 -1.1 1.1])
% title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
% subplot(2,1,2)
% hold
% plot(bclselect*(1:nsteps)/1000,ekf_noiseless./ekfmax,'.-')
% grid
% xlabel('time (sec)')
% ylabel('Normalized err. w/ fbk')
% axis([0 simtime/1000 -1.1 1.1])
% 
% saveas(gcf,[kf_folder param '/kfplot_wonoise_normalized_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% %saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
% 
% % Plot normalized errors, with noise
% figure
% subplot(2,1,1)
% hold
% plot(bclselect*(1:nsteps)/1000,eol./eolmaxn,'.-')
% grid
% ylabel('Normalized err. w/o fbk')
% title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
% axis([0 simtime/1000 -1.1 1.1])
% subplot(2,1,2)
% hold
% plot(bclselect*(1:nsteps)/1000,ekf./ekfmaxn,'.-')
% grid
% ylabel('Normalized err. w/ fbk')
% xlabel('time (sec)')
% axis([0 simtime/1000 -1.1 1.1])
% 
% saveas(gcf,[kf_folder param '/kfplot_noise_normalized_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% %saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])
% 
% % Plot absolute errors of variables targeted for "reconstruction"
% figure
% subplot(2,1,1)
% hold
% for k = 1:length(reconstructionindices)
%     plot(bclselect*(1:nsteps)/1000,eol(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
% end
% grid
% %ylabel('Error w/o fbk')
% ylabel('Error w/o fbk, mmol / L') % Use if concentrations were chosen
% title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
% leg = legend(statenames_latex(reconstructionindices,:));
% set(leg,'Interpreter', 'latex')
% %axis([0 simtime/1000 -ylimn ylimn])
% subplot(2,1,2)
% hold
% for k = 1:length(reconstructionindices)
%     plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
%     plot(bclselect*(1:nsteps)/1000,ekfapost(reconstructionindices(k),:),[symbols(k,1) ':'], 'Linewidth', linewidth/2)
% end
% grid
% %ylabel('Error w/ fbk')
% ylabel('A-pri (-) and a-post (--) errors w/ fbk, mmol / L') % Use if concentrations were chosen
% xlabel('time (sec)')
% leg = legend(statenames_latex(reconstructionindices,:));
% set(leg,'Interpreter', 'latex')
% %axis([0 simtime/1000 -ylimn ylimn])
% saveas(gcf,[kf_folder param '/kfplot_targetedvariables_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% 
% 
% % Plot absolute CL errors of variables targeted for "reconstruction"
% figure
% if shorttitleflag
%     title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
% else
%     title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
% end
% hold
% for k = 1:length(reconstructionindices)
%     plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
%     plot(bclselect*(1:nsteps)/1000,ekfapost(reconstructionindices(k),:),[symbols(k,1) ':'], 'Linewidth', linewidth/2)
% end
% grid
% %ylabel('Error w/ fbk')
% %ylabel('Estimation error, mmol / L') % Use if concentrations were chosen
% ylabel('A-pri (-) and a-post (--) errors w/ fbk, mmol / L') % Use if concentrations were chosen
% xlabel('time (sec)')
% leg = legend(statenames_latex(reconstructionindices,:));
% set(leg,'Interpreter', 'latex')
% set(gca,'fontsize',fontsize)
% set(gca,'linewidth',linewidth)
% %axis([0 simtime/1000 -ylimn ylimn])
% %saveas(gcf,[kf_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% saveas(gcf,[kf_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% 
% % Plot normalized CL errors of variables targeted for "reconstruction"
% figure
% if shorttitleflag
%     title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
% else
%     title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
% end
% hold
% for k = 1:length(reconstructionindices)
%     plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:)./ekfmaxn(reconstructionindices(k)),symbols(k,:), 'Linewidth', linewidth)
% end
% grid
% ylabel('Normalized est. err.')
% xlabel('time (sec)')
% leg = legend(statenames_latex(reconstructionindices,:));
% set(leg,'Interpreter', 'latex')
% set(gca,'fontsize',fontsize)
% set(gca,'linewidth',linewidth)
% %axis([0 simtime/1000 -ylimn ylimn])
% %saveas(gcf,[kf_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% saveas(gcf,[kf_folder param '/kfplot_targetedvariables_normederr_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% 
% % Plot normalized errors, with noise, and highlight first targeted variable
% figure
% if shorttitleflag
%     title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
% else
%     title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
% end
% hold
% plot(bclselect*(1:nsteps)/1000,ekf./ekfmaxn,'.-')
% grid
% ylabel('Normalized est. err.')
% xlabel('time (sec)')
% axis([0 simtime/1000 -1.1 1.1])
% for k = 1:1
%     p=plot(bclselect*(1:nsteps)/1000,ekf(reconstructionindices(k),:)./ekfmaxn(reconstructionindices(k)),symbols(k,:), 'Linewidth', 2*linewidth);
%     leg = legend(p,statenames_latex(reconstructionindices(k),:));
% end
% grid
% set(leg,'Interpreter', 'latex')
% set(gca,'fontsize',fontsize)
% set(gca,'linewidth',linewidth)
% saveas(gcf,[kf_folder param '/kfplot_normederr_overlay_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% 
% % Can use next portion if no Lkf was found. Compute eol first, then compute
% % associated OL norms and save command.
% if 0
%     % Plot absolute CL errors of variables targeted for "reconstruction"
%     figure
%     if shorttitleflag
%         title({['Sim. estim. errors, w/ noise, meas. = ' strtrim(statenames_latex(measurementindex,:)) ', BCL = ' num2str(bcl) 'ms'], ' '}, 'Interpreter', 'latex');
%     else
%         title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames_latex(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qnscalar), ', R = ' num2str(Rn)]}, 'Interpreter', 'latex');
%     end
%     
%     hold
%     for k = 1:length(reconstructionindices)
%         plot(bclselect*(1:nsteps)/1000,eol(reconstructionindices(k),:),symbols(k,:), 'Linewidth', linewidth)
%     end
%     grid
%     ylabel('Estimation error, mmol / L') % Use if concentrations were chosen
%     xlabel('time (sec)')
%     leg = legend(statenames_latex(reconstructionindices,:));
%     set(leg,'Interpreter', 'latex')
%     set(gca,'fontsize',fontsize)
%     set(gca,'linewidth',linewidth)
%     saveas(gcf,[kf_folder param '/kfplot_targetedvariables_ol_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
% end
% 
% %save kfwv1to10 *
% eval(['save ' kf_folder param '/kalmanfile_measind_' num2str(measurementindex) '_bcl_' num2str(bcl) ' *'])
% 
% 
