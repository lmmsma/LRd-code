% Design closed-loop controllers and observers for linearized LRd single-cell model.
% Generate open-loop responses for comparison. 
% Closed-loop tests: Luenberger observer (pole placement), Kalman filter,
% state feedback control (pole placement), Linear Quadratic Regulator
% (LQR). An older(?) version is lrd_linear_sim_cinc13_abstract.m. 
%
% Assume we can measure and perturb membrane potential (mp).
% Question: If you switch to stimulating via potassium current, does Bmp
% need to change? 

clear variables;
load lrddata_1cell_b1000_solem12_eig_relpert1e-005_start1_cdiff_Kp0
deltat = .005; % ms
numstate = size(jaccd,1);
fs = 18; %fontsize

%Bmp = deltat*[1; zeros(numstate-1,1)]; % Warning!!! This scaling is probably nonsensical. Why not leave out deltat for periodic control? 
Bmp = [1; zeros(numstate-1,1)]; 
Cmp = [1 zeros(1,numstate-1)];

%x0 = 0.1*rand(numstate,1); % Random IC. Warning!!! Sizes probably don't make sense (variables have different ranges, also perturbations of size 1 prob. not in linear regime)
x0 = 0.1*randn(numstate,1); % Random IC. Warning!!! Sizes probably don't make sense (variables have different ranges, also perturbations of size 1 prob. not in linear regime)

% Open-loop system
ncyc = 400000; % Simulate more than one cycle length
xall = zeros(numstate,ncyc);
xall(:,1) = x0; 
for ii = 1:ncyc-1
    xall(:,ii+1) = jaccd*xall(:,ii); 
end

figure
plot(1:ncyc,xall,'.')
xlabel('number of cycles') 
ylabel('Linearized LRd state variables')
axis([0 100 -0.1 2])

% For comparison, use 'initial' function to check results of earlier open-loop
% simulation 
syslrd = ss(jaccd, Bmp, eye(numstate), 0, deltat);
figure
initial(syslrd,rand(numstate,1))

% Check observability: 
% Use a larger rankcutoff value for reduced-order observer 
% to try to avoid using weakly observable modes in feedback 
rankcutoff = 1e-5; % What value should I choose? 
[Ab,Bb,Cb,T,kb] = obsvf(jaccd,Bmp,Cmp,rankcutoff);

% Observable part of state matrix in canonical form
Abo = Ab((numstate-sum(kb)+1):end,(numstate-sum(kb)+1):end); 
% Non-observable part of state matrix in canonical form
Abno = Ab(1:(numstate-sum(kb)),1:(numstate-sum(kb)));

% Observable part of output matrix in canonical form
Cbo = Cb(:,(numstate-sum(kb)+1):end); 

% "Strongly" observable eigenvalues (as defined by rankcutoff) 
po = eig(Abo,'nobalance');

% Desired eigenvalues
pod = po;

% Move largest eigenvalues closer to origin
%pod(1) = .96;
%pod(2) = .98; 
pod(1) = .8;
pod(2) = .9; 

% Reduced-order observer gain in canonical coordinates
L = place(Abo',Cbo',pod).'

% Check: are closed-loop eigenvalues in the correct locations? 
eig(Abo-L*Cbo,'nobalance')

% Full-order observer gain in canonical coordinates (pad with zeros to
% avoid feeding back weakly observable modes) 
Lb = [zeros(numstate-sum(kb),1);L]; 

% Transform gain back into original coordinates
Lfull = T'*Lb;
 
% Check: are closed-loop eigenvalues in the correct locations? 
eig(jaccd-Lfull*Cmp,'nobalance')

% Warning!!! You could use the predictor-corrector A-LCA format instead, and I'm not really sure why you
% aren't. 

% Closed-loop reduced-order Luenberger observer simulation 
% (linear system only) 
xcall = zeros(numstate,ncyc);
xcall(:,1) = x0; 
for ii = 1:ncyc-1
    xcall(:,ii+1) = (jaccd-Lfull*Cmp)*xcall(:,ii); 
end

figure
plot(1:ncyc,xcall)
xlabel('number of cycles') 
ylabel('Linearized LRd state variables')


% Kalman filter
% Pre-generate unscaled noise signals so that different ratios can be
% tested with the same shape of w and v trajectories
% set(0,'DefaultAxesLineStyleOrder','x-|s--|o-.|^:|-|s--|o-.|^:|')
  
wu = randn(ncyc,1);
vu = randn(ncyc,1); 

%Qn = .01;%0.1; % process noise covariance (units?) 
%Rn = 0.001;%0.01; % meas. noise covariance (mV^2?)
Qn = .01;%0.1; % process noise covariance (units?) 
%Rn = Qn/10;% meas. noise covariance (mV^2?)
Rn = 10*Qn;% meas. noise covariance (mV^2?)
Bw = [1;ones(numstate-1,1)]; % Input noise matrix, noise only added to first channel
%Bw = ones(numstate,1); % Input noise matrix
% For the abstract, I added the same noise signal to every channel, which
% isn't particularly realistic. A better test could be to set Bw = eye and
% add separate random process noise signals to each channel. For now, I've
% used the same structure as B (input to first channel only). 

%w=sqrt(Qn)*randn(ncyc,1);
%v=sqrt(Rn)*randn(ncyc,1);
w=sqrt(Qn)*wu;
v=sqrt(Rn)*vu;

figure
plot(v)
title('meas noise')

syslrdkf = ss(jaccd, [Bmp Bw], Cmp, [0 0], deltat);

% open-loop system with noise (error dynamics only)
xnall = zeros(numstate,ncyc);
xnall(:,1) = x0; 
for ii = 1:ncyc-1
    xnall(:,ii+1) = jaccd*xnall(:,ii) + Bw*w(ii); 
end

figure
hold
plot(1:ncyc,xnall,'.-')
xlabel('number of cycles') 
%ylabel('Linearized LRd state variables')
ylabel('open-loop estimation error')
axis([0 300 -10 10])

figure
hold
p1=plot(1:ncyc,xnall(1,:),'b+-','linewidth',2);
p2=plot(1:ncyc,xnall(5,:),'g--','linewidth',2);
p3=plot(1:ncyc,xnall(8,:),'r-.','linewidth',2);
p4=plot(1:ncyc,xnall(9,:),'k:','linewidth',2);
p5=plot(1:ncyc,xnall(10,:),'m-','linewidth',2);
legend([p1 p2 p3 p4 p5],'V','d','[Ca^{2+}]_{i,t}','[Na^+]_i','[K^+]_i')
xlabel('number of cycles','fontsize',fs) 
ylabel(['open-loop estimation error'],'fontsize',fs) 
set(gca,'fontsize',fs)
axis([0 300 -5 10])
set(gcf,'position',[403 100 560 550])
YTL = get(gca,'yticklabel');
set(gca,'yticklabel',[YTL,repmat(' ',size(YTL,1),1)])


[kest,Lkf,P] = kalman(syslrdkf,Qn,Rn,0);

% pseudo-KF (no noise) error dynamics
xkfall = zeros(numstate,ncyc);
xkfall(:,1) = x0; 
for ii = 1:ncyc-1
    xkfall(:,ii+1) = (jaccd-Lkf*Cmp)*xkfall(:,ii); 
end

figure
plot(1:ncyc,xkfall,'.-')
xlabel('number of cycles') 
%ylabel('Linearized LRd state variables')
ylabel('closed-loop estimation error, no input or output noise')
axis([0 300 -10 10])

% KF with noise (error dynamics only)
xkfnall = zeros(numstate,ncyc);
xkfnall(:,1) = x0; 
for ii = 1:ncyc-1
    xkfnall(:,ii+1) = (jaccd-Lkf*Cmp)*xkfnall(:,ii) + Bw*w(ii) - Lkf*v(ii); 
end

figure
plot(1:ncyc,xkfnall,'.')
xlabel('number of cycles') 
%ylabel('Linearized LRd state variables')
ylabel('Kalman filter estimation error')
axis([0 300 -10 10])

figure
hold
p1=plot(1:ncyc,xkfnall(1,:),'b+-','linewidth',2);
p2=plot(1:ncyc,xkfnall(5,:),'g--','linewidth',2);
p3=plot(1:ncyc,xkfnall(8,:),'r-.','linewidth',2);
p4=plot(1:ncyc,xkfnall(9,:),'k:','linewidth',2);
p5=plot(1:ncyc,xkfnall(10,:),'m-','linewidth',2);
legend([p1 p2 p3 p4 p5],'V','d','[Ca^{2+}]_{i,t}','[Na^+]_i','[K^+]_i')
xlabel('number of cycles','fontsize',fs) 
ylabel(['Kalman filter estimation error, W / V = ' num2str(Qn/Rn)],'fontsize',fs)
set(gca,'fontsize',fs)
axis([0 300 -5 10])
set(gcf,'position',[403 100 560 550])
YTL = get(gca,'yticklabel');
set(gca,'yticklabel',[YTL,repmat(' ',size(YTL,1),1)])


% Compare error covariances over some number of cycles
col = cov(xnall(:,1:300)');
norm(col)
ckf = cov(xkfnall(:,1:300)');
norm(ckf)
%norm(col)/norm(ckf)
norm(ckf)/norm(col)

%save save kfwv1to10 *



% Check controllability: 
[Abc,Bbc,Cbc,Tc,kbc] = ctrbf(jaccd,Bmp,Cmp,rankcutoff);

% Controllable part of state matrix in canonical form
Abcc = Abc((numstate-sum(kbc)+1):end,(numstate-sum(kbc)+1):end); 
% Non-controllable part of state matrix in canonical form
Abcnc = Abc(1:(numstate-sum(kbc)),1:(numstate-sum(kbc)));

% Controllable part of input matrix in canonical form
Bbcc = Bbc((numstate-sum(kbc)+1):end,:); 

% "Strongly" controllable eigenvalues (as defined by rankcutoff) 
pc = eig(Abcc,'nobalance');

% Are these the same eigenvalues? Probably.  
[sort(pc) sort(po)]
sort(diag(d))

% Desired eigenvalues
pcd = pc;
pcd(2) = pod(1);
pcd(3) = pod(2);


% Reduced-order feedback gain in canonical coordinates
K = place(Abcc,Bbcc,pcd)

% Check: are closed-loop eigenvalues in the correct locations? 
eig(Abcc-Bbcc*K,'nobalance')

% Full-order controller gain in canonical coordinates (pad with zeros to
% avoid feeding back weakly controllable modes) 
Kb = [zeros(1,numstate-sum(kbc)) K]; 

% Transform gain back into original coordinates
Kfull = Kb*Tc;
 
% Check: are closed-loop eigenvalues in the correct locations? 
eig(jaccd-Bmp*Kfull,'nobalance')


% Closed-loop reduced-order pole-placement controller simulation 
% (linear system only) 
xcpall = zeros(numstate,ncyc);
xcpall(:,1) = x0; 
for ii = 1:ncyc-1
    xcpall(:,ii+1) = (jaccd-Bmp*Kfull)*xcpall(:,ii); 
end

figure
plot(1:ncyc,xcpall)
xlabel('number of cycles') 
ylabel('Linearized LRd state variables')

% Try LQR (given that (jaccd,Bmp) is stabilizable)
Q = inv(Smat.^2);
R = 0.1; 
[Klqr,S,e] = dlqr(jaccd,Bmp,Q,R,0);

eig(jaccd-Bmp*Klqr,'nobalance')

% Closed-loop LQR simulation 
% (linear system only) 
xlqall = zeros(numstate,ncyc);
xlqall(:,1) = x0; 
for ii = 1:ncyc-1
    xlqall(:,ii+1) = (jaccd-Bmp*Klqr)*xlqall(:,ii); 
end

figure
plot(1:ncyc,xlqall)
xlabel('number of cycles') 
ylabel('Linearized LRd state variables')


%save cinc13_abs_results *

