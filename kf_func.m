%Runs the test_linear_kf script but asks for the measurement index as an
%input. Allows the kf to be calculated/plotted in a for loop very easily
function kf_func( meas )

for adj_yn = 0:1
    for bcl_in = [200]
        ms = 10; fs = 14; % marker size and font size

        statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');

        Jacobians = 'Jacobians/'; % folder where jacobians are stored

        kf_folder = 'Kalman/';

        if adj_yn
            param= 'adj';
        else
            param = 'def';
        end

        eval(['load ' Jacobians 'jacfile_' param ' *']) %Load Jacobians

        numstate = size(alljacs{1},1); % number of state variables

        nbcls = length(bcls); % number of bcls

        % Simulate the system for some fixed amount of time, and figure out the
        % corresponding number of BCLs later. 
        % What is a reasonable amount of time? How quickly do we want 
        % the noise-free errors to converge? 
        simtime = 60000; % ms

        noisecycles = round(simtime/60); % number of cycles for which noise will be generated. 
        % Pick a number greater than the number of cycles implied by simtime. 
        % The above should work as long as BCL >= 60ms.

        % input matrix for all possible individual inputs
        B = data.dt*eye(numstate); % This isn't necessarily right, just a placeholder. In general, B isn't square and its nonzero elements aren't 1.

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

        measurementindex=  meas;
        % measurementindex = 1; % Select an integer from 1 to numstate. 
        % This is the index of the state variable that is measured. 

        controlindex = 1;  % Select an integer from 1 to numstate.  
        % This is the index of the state variable that is controlled. 
        % It should not affect KF design, unless process noise matrix 
        % Bw is redesigned to make use of this index. 

        % Initial condition
        x0 = 0.01*ones(numstate,1); 
        % This should be replaced with a random IC (e.g., xhat system should have 
        % random IC) after the xhat system is separated out. 
        % Choose a random initial condition 
        %x0 = 0.01*randn(numstate,1); % Random IC. 
        %Warning!!! Sizes probably don't make sense 
        %(variables have different ranges, also 
        % perturbations of size 1 prob. not in linear regime)

        % Pre-generate unit-variance noise signals so that different Q/R 
        % ratios can betested with the same shape of w and v trajectories

        %rng('default'); % Set random number seed to default for now, to make it 
        % easier to compare across different runs? This didn't work the way I
        % expected. 

        wu = randn(1,noisecycles); 
        vu = randn(1,noisecycles); 

        % Units of noise signals will depend on which variables are 
        % measured and/or receive process noise. 

        Qn = 0.1; % process noise covariance 
        Rn = 0.001*Qn;% meas. noise covariance (units of mV^2 if V is the measurement?)

        Bw = data.dt*ones(numstate,1); % Process noise matrix
        %Bw = data.dt*eye(numstate); % Process noise matrix
        % For the 2013 study, I added the same noise signal to every channel, which
        % isn't particularly realistic. A better test could be to set Bw = eye and
        % add separate random process noise signals to each channel. 

        % Adjust noise signals to reflect covariances chosen above
        w=sqrt(Qn)*wu; % Process noise 
        v=sqrt(Rn)*vu; % Measurement noise

        % Select a BCL of interest. Trying to keep the number of BCLs low for now, 
        % to reduce the number of plots 
        bclindices = 1:nbcls; 
        bclselect = bcl_in; % ms

        %%% Code after this line could be put inside a BCL loop, if more than one
        %%% BCL were selected 

        % Find Jacobian index matching the BCL you chose
        bclselectindex = bclindices(bcls==bclselect); 

        % Exit if you didn't find the right index
        if bclselect ~= bcls(bclselectindex)
            disp('Error: BCL index mismatch')
            return;
        end

        bcl = bcls(bclselectindex);
        % print current BCL to screen
        disp(['BCL = ' num2str(bcls(bclselectindex)) ' ms'])

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
        [kest,Lkf,P] = kalman(sys,Qn,Rn,0);

        % Open-loop eigenvalues
        oleig = eig(jaccd); 

        % Closed-loop eigenvalues
        cleig = eig((jaccd-Lkf*C(measurementindex,:))); 

        % Display open and closed-loop eigenvalues
        [oleig cleig];

        % open-loop system without noise (error dynamics only)
        xall = zeros(numstate,nsteps);
        xall(:,1) = x0; 
        for ii = 1:nsteps-1
            xall(:,ii+1) = jaccd*xall(:,ii); 
        end

        % open-loop system with noise (error dynamics only)
        xnall = zeros(numstate,nsteps);
        xnall(:,1) = x0; 
        for ii = 1:nsteps-1
            xnall(:,ii+1) = jaccd*xnall(:,ii) + Bw*w(:,ii); 
        end

        % pseudo-KF (no noise) error dynamics
        xkfall = zeros(numstate,nsteps);
        xkfall(:,1) = x0; 
        for ii = 1:nsteps-1
            xkfall(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*xkfall(:,ii); 
        end

        % KF with noise (error dynamics only)
        xkfnall = zeros(numstate,nsteps);
        xkfnall(:,1) = x0; 
        for ii = 1:nsteps-1
            xkfnall(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*xkfnall(:,ii) + Bw*w(:,ii) - Lkf*v(ii); 
        end

        % Compute scaled errors (may be better for plotting, to prevent larger variables
        %from dominating the plot?) 
        % xsnall = Smat*xnall; 
        % xskfall = Smat*xkfall;
        % xskfnall = Smat*xkfnall; 

        % Choose a y-axis limit for error plots
        if max(abs(oleig)) <= 1 % Open-loop system is (arguably) at least marginally stable
            ylim = max(max(abs(xall))); % choose max open-loop error
            ylimn = max(max(abs(xnall))); % choose max open-loop error, with noise
        else % OL system is unstable
            ylim = max(max(abs(xkfall))); % choose max closed-loop error
            ylimn = max(max(abs(xkfnall))); % choose max closed-loop error, with noise
        end

        figure
        subplot(2,1,1)
        hold;
        plot(bclselect*(1:nsteps)/1000,xall,'.-')
        grid 
        ylabel('Error without feedback')
        axis([0 simtime/1000 -ylim ylim])
        title({[param ': Simulated estimation errors, with no noise, meas. = ' strtrim(statenames(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qn), ', R = ' num2str(Rn)]});
        subplot(2,1,2)
        hold 
        plot(bclselect*(1:nsteps)/1000,xkfall,'.-')
        grid 
        xlabel('time (sec)') 
        ylabel('Error with feedback')
        axis([0 simtime/1000 -ylim ylim])

        saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
        saveas(gcf,[kf_folder 'kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl) '.jpeg'])

        figure
        subplot(2,1,1)
        hold;
        plot(bclselect*(1:nsteps)/1000,xnall,'.-')
        grid 
        ylabel('Error without feedback')
        title({[param ': Simulated estimation errors, with noise, meas. = ' strtrim(statenames(measurementindex,:))], ['BCL = ' num2str(bcl) 'ms, Q = ' num2str(Qn), ', R = ' num2str(Rn)]});
        axis([0 simtime/1000 -ylimn ylimn])
        subplot(2,1,2)
        hold
        plot(bclselect*(1:nsteps)/1000,xkfnall,'.-')
        grid 
        ylabel('Error with feedback')
        xlabel('time (sec)') 
        axis([0 simtime/1000 -ylimn ylimn])

        saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex, :) '_' num2str(bcl)])
        saveas(gcf,[kf_folder 'kfplot_noise_' param '_' statenames(measurementindex, :) '_' num2str(bcl) '.jpeg'])
    end
end