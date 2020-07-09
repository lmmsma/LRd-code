% Test out linear Luenberger observer (pole placement) designs for selected
% conditions, using LRd model Jacobians. Code based on test_linear_kf.m
% Laura Munoz, Sep 2018

clear variables;
%shiftstring = '' % if using default shift of data.dt
%shiftstring = '_shift0p75ms'
%shiftstring = '_shift6p5ms'
%shiftstring = '_shift7mVrepol'
%shiftstring = '_shift-50mVrepol'
%shiftstrings = {''};
%shiftstrings = {'_shift0p2Vnormrepol'}
%shiftstrings = {'_shift0p001Vnormrepol'}
%shiftstrings = {'','_shift0p75ms','_shift6p5ms','_shift7mVrepol','_shift-50mVrepol'};
%shiftstrings = {'','_shift1Vnormdepol','_shift0p2Vnormdepol','_shift0p4Vnormdepol','_shift0p6Vnormdepol','_shift0p8Vnormdepol','_shift0p2Vnormrepol','_shift0p4Vnormrepol','_shift0p6Vnormrepol','_shift0p8Vnormrepol','_shift0p001Vnormrepol'};
%allshiftstrings = {'_shift1Vnormdepol'}
allshiftstrings = {'_shift0p2Vnormrepol'}
%allshiftstrings = {''}; 

%rng(1, 'twister'); % reset random number generator
rng(1); % reset random number generator

markersize = 80; fontsize = 14; linewidth = 2; % marker size, font size, line width
shorttitleflag = 1; % set to 1 to reduce title length for presentation figures

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{I}$','$[Na^+]_I$','$[K^+]_I$','$[Ca^{2+}]_{J}$','$[Ca^{2+}]_N$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');
% plot symbols
symbols = char('ko-','bs-','rp-','m*-','g^-','cv-','yh-');
numstate = length(statenames);

selected_logepsln = -5; % Base 10 log of step size used in Jacobian computation

%lu_folder = 'luenberger/';
lu_folder = 'linear_luen_tests/';
%OCvalues = 'OCvalues/'; %folder where obsv and ctrb values will be saved

adj_yn = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if adj_yn
    param = 'adj';
else
    param = 'def';
end
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

Qnscalar = 1; % process noise factor

%Rn = 100*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)
Rn = 0.001*Qnscalar;% meas. noise covariance (units of mV^2 if V is the measurement?)

Qn = Qnscalar*eye(numstate); % process noise covariance

% For the 2013 study, I added the same noise signal to every channel, which
% isn't particularly realistic. A better test could be to set Bw = eye and
% add separate random process noise signals to each channel.
%Bw = ones(numstate,1); % Process noise matrix
Bw = eye(numstate); % Process noise matrix

cleigfactor_moveonepole = 0.9; % When moving one pole, multiply original pole by this factor
cleigfactor = 0.9; % Multiply KF CL poles by this factor, before using in pole placement
eigthresh = 0.1; % OL eigenvalues below this size won't be moved
eigtol = 1e-6; % use to search for eigenvalues near a reference value

for shiftctr = 1:length(allshiftstrings)
    shiftstring = allshiftstrings{shiftctr}
    
    %jacfolder = ['jacobians/' param '/']; % folder where jacobians are stored
    jacfolder = ['jacobians' shiftstring '/' param '/']; %This folder is where jacobians are stored
    
    
    eval(['load ' jacfolder 'jacfile' num2str(selected_logepsln) ' *']); %Load Jacobians
    
    
    % number of bcls
    nbcls = length(selected_bcls_for_fps);
    
    
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
    Cscaled = C(measurementindex,:)*Smatinv;
    
    
    
    % Need to define Ts properly for use with freqsep:
    sys = ss(jaccd, B(:,controlindex), C(measurementindex,:), 0, bcl/1000);
    
    % Use scaled matrices this time
    sys_scaled = ss(jaccd_scaled, Bs(:,controlindex), Cscaled, 0, bcl/1000);
    
    
    sys_scaledforkf = ss(jaccd_scaled, [Bs(:,controlindex) Smat*Bw], Cscaled, [0 zeros(1,size(Bw,2))], bcl/1000);
    
    [kestscaled,Lkfscaled,Pscaled,Mscaled,Zscaled] = kalman(sys_scaledforkf,Qn,Rn,zeros(size(Qn,1),1));
    
    
    %oleig = eig(jaccd);
    oleig = sort(eig(jaccd),'descend','ComparisonMethod','abs');
    
    cleigsckf = sort(eig(jaccd_scaled-Lkfscaled*Cscaled),'descend','ComparisonMethod','abs');
    
    % Open-loop eigenvalues
    
    % Move one pole:
    deseig_moveonepole = oleig;
    [uniteigrow, uniteigcol] = find(abs(oleig) < 1.0+eigtol & abs(oleig) > 1.0-eigtol);
    deseig_moveonepole(uniteigrow) = cleigfactor_moveonepole*deseig_moveonepole(uniteigrow);
    deseig_moveonepole
    largesteigindices = find(abs(oleig) >= eigthresh);
    %deseig(largesteigindices) = 0.9*deseig(largesteigindices);
    %deseig(1) = 0.95;
    %    deseig(1) =0.7*oleig(1);
    %     % Try to move top 3 eigenvalues inward
    %     if real(oleig(1) < 0) % alternans eigenvalue is largest eigenvalue
    %         if abs(oleig(1))> 0.9
    %             deseig(1) = 0.9*oleig(1);
    %         end
    %         deseig(2) = 0.9+0.1i;
    %         deseig(3) = 0.9-0.1i;
    %     else
    %         deseig(1) = 0.9+0.1i;
    %         deseig(2) = 0.9-0.1i;
    %         if abs(oleig(3))> 0.9
    %             deseig(3) = 0.9*oleig(3);
    %         end
    %     end
    
    % Move multiple poles:
    %    deseig = 0.9*cleigsckf; % This is kludgey, but multiply KF CL eigs by a scaling factor to achieve a desired max modulus
    deseig = oleig;
    % Check for complex conjugates in CL KF eigs: 
    if largesteigindices < numstate
        % If the next eig on the list is within eigtol of previous eig,
        % they are probably conjugates, and both should be included
        if (abs(cleigsckf(largesteigindices(end)+1)) < abs(cleigsckf(largesteigindices(end)))+eigtol) && (abs(cleigsckf(largesteigindices(end)+1)) > abs(cleigsckf(largesteigindices(end)))-eigtol)
            largesteigindices = [largesteigindices; largesteigindices(end) + 1]; 
        end
    end
    deseig(largesteigindices) = cleigfactor*cleigsckf(largesteigindices);
    
    % Display observability values
    minsv_obsv = min(svd(obsv(jaccd,C(measurementindex,:))))
    minsv_obsvsc = min(svd(obsv(jaccd_scaled,Cscaled)))
    
    % Compare with modal values for the given condition:
    ocfolder = ['ocvalues' shiftstring '/']; %folder where obsv and ctrb values will be saved
    eval(['load ' ocfolder param '/ocfile sortedeigsc eigof obsv_avgd_over_modes_sc']);
    
    disp('Eigenvalues and scaled obsv, no averaging: ')
    [sortedeigsc(bclselectindex,:)' log10(eigof.obsvmag_sc{bclselectindex}(measurementindex,:))']
    
    disp('Scaled obsv avgd over modes but not BCLs or shifts: ')
    log10(obsv_avgd_over_modes_sc(measurementindex,bclselectindex))
    
    if 0
        % Recompute modal observability measure for largest mode(these should match the values from
        % compute_obsv_ctrb.m)
        [v,d,w]=eig(jaccd,'nobalance');
        [sortedeig,eigsortind] = sort(abs(diag(d)),'descend');
        sortedv = v(:,eigsortind);
        sortedw = w(:,eigsortind);
        
        % Compute eigenvectors for scaled Jacobian:
        [vsc,dsc,wsc]=eig(jaccd_scaled,'nobalance');
        [sortedeigsc,eigsortindsc] = sort(abs(diag(dsc)),'descend');
        sortedvsc = vsc(:,eigsortindsc);
        sortedwsc = wsc(:,eigsortindsc);
        obsvmag = abs(C(measurementindex,:)*sortedv(:,1))/(norm(C(measurementindex,:))*norm(sortedv(:,1)))
        obsvmagscnumerator = abs(Cscaled*sortedvsc(:,1))
        obsvmagsc = abs(Cscaled*sortedvsc(:,1))/(norm(Cscaled)*norm(sortedvsc(:,1)))
    end
    % Put a check (for debugging purposes) on results from compute_obsv_ctrb.m
    if 0
        eigthresh = 0.9
        controlindexdebug = 1;
        Bscaled = Smat*B(:,controlindexdebug);
        % Compute individual (per mode, per BCL) obsv and ctrb measures
        startbcl = 110 %ms, this should be the longer BCL
        endbcl = 100 %ms, this should be the shorter BCL
        allbclindices = bclindices;
        
        startbclindex = allbclindices(selected_bcls_for_fps == startbcl)
        endbclindex = allbclindices(selected_bcls_for_fps == endbcl)
        
        jaccdstart = alljacs{startbclindex};
        jaccdend = alljacs{endbclindex};
        
        jaccdstart_scaled = Smat*jaccdstart*Smatinv;
        jaccdend_scaled = Smat*jaccdend*Smatinv;
        
        % Compute obsv and ctrb magnitudes
        
        [vs,ds,ws]=eig(jaccdstart,'nobalance');
        [sortedeigs,eigsortinds] = sort(abs(diag(ds)),'descend');
        sortedvs = vs(:,eigsortinds);
        sortedws = ws(:,eigsortinds);
        
        [ve,de,we]=eig(jaccdend,'nobalance');
        [sortedeige,eigsortinde] = sort(abs(diag(de)),'descend');
        sortedve = ve(:,eigsortinde);
        sortedwe = we(:,eigsortinde);
        
        [vssc,dssc,wssc]=eig(jaccdstart_scaled,'nobalance');
        [sortedeigssc,eigsortindssc] = sort(abs(diag(dssc)),'descend');
        sortedvssc = vssc(:,eigsortindssc);
        sortedwssc = wssc(:,eigsortindssc);
        
        [vesc,desc,wesc]=eig(jaccdend_scaled,'nobalance');
        [sortedeigesc,eigsortindesc] = sort(abs(diag(desc)),'descend');
        sortedvesc = vesc(:,eigsortindesc);
        sortedwesc = wesc(:,eigsortindesc);
        
        % Since this is a "manual" check, I'm trying to hardcode certain items.
        % Examination of eigenvalues shows that there are 3 eigenvalues with
        % moduli > 0.9 for BCLs of 110 and 100ms.
        nmodesabovethresh = 3;
        % Compute magnitudes for top modes:
        obsvmags = zeros(1,nmodesabovethresh);
        obsvmage = zeros(1,nmodesabovethresh);
        obsvmagssc = zeros(1,nmodesabovethresh);
        obsvmagesc = zeros(1,nmodesabovethresh);
        ctrbmags = zeros(1,nmodesabovethresh);
        ctrbmage = zeros(1,nmodesabovethresh);
        ctrbmagssc = zeros(1,nmodesabovethresh);
        ctrbmagesc = zeros(1,nmodesabovethresh);
        
        for ii = 1:nmodesabovethresh
            obsvmags(ii) = abs(C(measurementindex,:)*sortedvs(:,ii))/(norm(C(measurementindex,:))*norm(sortedvs(:,ii)));
            obsvmagssc(ii) = abs(Cscaled*sortedvssc(:,ii)); % Leave out normalization for scaled values
            obsvmage(ii) = abs(C(measurementindex,:)*sortedve(:,ii))/(norm(C(measurementindex,:))*norm(sortedve(:,ii)));
            obsvmagesc(ii) = abs(Cscaled*sortedvesc(:,ii)); % Leave out normalization for scaled values
            ctrbmags(ii) = abs(sortedws(:,ii)'*B(:,controlindexdebug))/(norm(sortedws(:,ii))*norm(B(:,controlindexdebug)));
            ctrbmagssc(ii) = abs(sortedwssc(:,ii)'*Bscaled); % Leave out normalization for scaled values
            ctrbmage(ii) = abs(sortedwe(:,ii)'*B(:,controlindexdebug))/(norm(sortedwe(:,ii))*norm(B(:,controlindexdebug)));
            ctrbmagesc(ii) = abs(sortedwesc(:,ii)'*Bscaled); % Leave out normalization for scaled values
        end
        obsvmags
        obsvmagssc
        obsvmage
        obsvmagesc
        ctrbmags
        ctrbmagssc
        ctrbmage
        ctrbmagesc
        
        % Compute obsv and ctrb averages for same BCL over top modes
        avgobsv_overmodes_s = mean(obsvmags)
        avgobsv_overmodes_ssc = mean(obsvmagssc)
        avgobsv_overmodes_e = mean(obsvmage)
        avgobsv_overmodes_esc = mean(obsvmagesc)
        
        avgctrb_overmodes_s = mean(ctrbmags)
        avgctrb_overmodes_ssc = mean(ctrbmagssc)
        avgctrb_overmodes_e = mean(ctrbmage)
        avgctrb_overmodes_esc = mean(ctrbmagesc)
        
        % Compute obsv and ctrb averages for top 2 BCLs over top mode
        avgobsv_overbcls_topmode = mean([obsvmags(1) obsvmage(1)])
        avgobsv_overbcls_topmode_sc = mean([obsvmagssc(1) obsvmagesc(1)])
        avgctrb_overbcls_topmode = mean([ctrbmags(1) ctrbmage(1)])
        avgctrb_overbcls_topmode_sc = mean([ctrbmagssc(1) ctrbmagesc(1)])
        
        % Compute obsv and ctrb averages for top 2 BCLs over top modes
        avgobsv = mean([avgobsv_overmodes_s avgobsv_overmodes_e])
        avgobsvsc = mean([avgobsv_overmodes_ssc avgobsv_overmodes_esc])
        avgctrb = mean([avgctrb_overmodes_s avgctrb_overmodes_e])
        avgctrbsc = mean([avgctrb_overmodes_ssc avgctrb_overmodes_esc])
        
        log10(avgobsv)
        log10(avgobsvsc)
        log10(avgctrb)
        log10(avgctrbsc)
    end
    
    [Llu,prec,message] = place(jaccd',C(measurementindex,:)',deseig);
    [Lluscaled_moveonepole,precscaled_moveonepole,messagescaled_moveonepole] = place(jaccd_scaled',Cscaled',deseig_moveonepole);
    [Lluscaled,precscaled,messagescaled] = place(jaccd_scaled',Cscaled',deseig);
    
    disp('Scaled normalized gain product (Ls*Cs), move one pole only:')
    LluscCsnorm_moveonepole = norm(Lluscaled_moveonepole'*Cscaled)
    disp('Log10 Scaled normalized gain product (Ls*Cs), move one pole only:')
    log10(LluscCsnorm_moveonepole)
    
    % Smat*L*Cs should be dimensionless:
    SmatLluCsnorm = norm(Smat*Llu'*Cscaled)
    disp('Scaled normalized gain product (Ls*Cs):')
    % This should be the same number as the previous one:
    LluscCsnorm = norm(Lluscaled'*Cscaled)
    disp('Log 10 Scaled normalized gain product (Ls*Cs):')
    log10(LluscCsnorm)
    %
    %     % Compare with Ackermann's method
    %     La = acker(jaccd',C(measurementindex,:)',deseig);
    %     Las = acker(jaccd_scaled',Cscaled',deseig);
    %
    %     Llunorm = norm(Llu) % size of gain is related to strength of feedback
    %     Lluscalednorm = norm(Lluscaled) % size of gain is related to strength of feedback
    %     %Llusc_divbyCsc = norm(Lluscalednorm)/norm(Cscaled)
    %
    %
    %     Lanorm = norm(La);
    %     Lasnorm = norm(Las);
    %
    
    % Closed-loop eigenvalues
    %cleig = eig(jaccd-Llu'*C(measurementindex,:));
    cleig = sort(eig(jaccd-Llu'*C(measurementindex,:)),'descend','ComparisonMethod','abs');
    disp('Max CL eig, unscaled sys check:')
    cleig(1)
    cleigscaled = sort(eig(jaccd_scaled-Lluscaled'*Cscaled),'descend','ComparisonMethod','abs'); % These should match the unscaled case
    disp('Max CL eig difference norm, unscaled vs. scaled (should be small):')
    max(abs(sort(cleig,'descend','ComparisonMethod','abs')-sort(cleigscaled,'descend','ComparisonMethod','abs')))
    cleigscaled_moveonepole = sort(eig(jaccd_scaled-Lluscaled_moveonepole'*Cscaled),'descend','ComparisonMethod','abs'); % These should match the unscaled case
    
    [oleig cleigscaled_moveonepole cleigscaled]
    
    disp('Max CL eig modulus, move one pole: ')
    max(abs(cleigscaled_moveonepole))
    disp('Max CL eig modulus, move all poles: ')
    max(abs(cleigscaled))
    

    % Simulate linear closed-loop scaled system
    numsimsteps = 30; 
    err_scaled = zeros(numstate,numsimsteps); % estimation error for scaled system
%    err_scaled(:,1) = 0.001*randn(numstate,1); % IC
    pertsize = 10^-5;%10^-2; %10^2*epsln % This is the size of the perturbation that will be applied
    sv_pert_ind = 'rand';
%    sv_pert_dir = inv(Smat)*(rand(numstate,1)-0.5*ones(numstate,1))
    % Don't need factor of inv(Smat) for the scaled system)
    sv_pert_dir = (rand(numstate,1)-0.5*ones(numstate,1)); 
    err_scaled(:,1) = pertsize*sv_pert_dir; % IC

    for ii = 1:numsimsteps-1
        err_scaled(:,ii+1) = (jaccd_scaled-Lluscaled'*Cscaled)*err_scaled(:,ii);    
    end
    figure
    hold on;
    plot(1:numsimsteps, err_scaled);
    ylabel('estimation error')
    xlabel('beat index')
    title(['$T = $' num2str(bcl) 'ms, Meas = ' statenames_latex(measurementindex, :) ', ' shiftstring(7:end) ', ' param],'Interpreter','Latex'); 
    
    errorendindex = numsimsteps; % ending cycle index over which to compute error norm
    errorstartindex = errorendindex-4; %1; % Enter a starting index for error norm computation
    error_lin_nondim_mean = mean(mean(abs(err_scaled(:,errorstartindex:errorendindex))))

    % Rel pct diff just appears to double the rel error here, since the
    % reference trajectory should just be zero? I think RPD doesn't make as
    % much sense if you're just working directly with an error quantity. 
    % Breaking err into xhat and x still doesn't seem like it would help 
    % since x = X-X* will still be zero if x is taken from the fixed point
    % trajectory. The following will usually just be -2 or 2, not taking
    % any intermediate values: 
    % err_scaled_rpd = 2*(err_scaled-0)./(abs(err_scaled)+abs(0));
    % A workaround is to rebuild everything: X=X*+x, then compute 
    % 2*(X-X*)/(|X|+|X*|)
    
    savefilename = ['lin_luen_results_b' num2str(bcl) '_m' num2str(measurementindex) shiftstring '_ncyc' num2str(numsimsteps) '_pdirrand_psz' num2str(round(log10(pertsize)))];
    eval(['save ' lu_folder savefilename ' *']);

    
    %% Try ROO design
    rankcutoff = 1e-5; % below this level, singular values don't contribute to the rank
    [Abars,Bbars,Cbars,Ts,ks] = obsvf(jaccd_scaled,Bs(:,controlindex),Cscaled,rankcutoff);
    rankofsc = sum(ks)
    Abarso = Abars((numstate-rankofsc+1):end,(numstate-rankofsc+1):end);
    Cbarso = Cbars(:,(numstate-rankofsc+1):end);
    eigofsc = eig(Abarso,'nobalance') % eigenvalues of Ao
    %This part is incomplete. Using default rank cutoff can sometimes just
    %yield entire system (no states removed).
    
    
    %% Compute quantities for reduced-order scaled system
    disp('Problem! Freqsep does not capture negative eigenvalues, even if they are large. \n Press any key to continue')
    pause
    cutoff_freq_disc = .1; % eigenvalue modulus threshold
    %cutoff_freq_disc = .01; % eigenvalue modulus threshold
    cutoff_freq_cts = abs(log(cutoff_freq_disc))/(bcl/1000)
    
    [Gs, Gf] = freqsep(sys, cutoff_freq_cts);
    [Gs_scaled, Gf_scaled] = freqsep(sys_scaled, cutoff_freq_cts);
    
    %minsv_obsv_scaled_lf = min(svd(obsv(Gs_scaled.A,Gs_scaled.C)))
    
    oleigGs = sort(eig(Gs.A),'descend','ComparisonMethod','abs');
    [vGs,dGs]=eig(Gs.A);
    % Sort eigenvalues and eigenvectors in descending order of size
    [sortedeigGs,eigGssortind] = sort(abs(diag(dGs)),'descend');
    sortedvGs = vGs(:,eigGssortind);
    % This is only for the first eigenvalue:
    cdotvGs = Gs.C*sortedvGs(:,1)
    obsvmagGs = abs(cdotvGs)/(norm(Gs.C)*norm(sortedvGs(:,1)))
    
    [vGsscaled,dGsscaled]=eig(Gs_scaled.A);
    % Sort eigenvalues and eigenvectors in descending order of size
    [sortedeigGsscaled,eigGssortindscaled] = sort(abs(diag(dGsscaled)),'descend');
    sortedvGsscaled = vGs(:,eigGssortindscaled);
    % This is only for the first eigenvalue:
    cdotvGsscaled = Gs_scaled.C*sortedvGsscaled(:,1)
    obsvmagGsscaled = abs(cdotvGsscaled)/(norm(Gs_scaled.C)*norm(sortedvGsscaled(:,1)))
    
    % Pole placement
    deseigGs = oleigGs;
    %deseigGs(1) = 0.95;
    %deseigGs(1) = 0.9*oleigGs(1);
    deseigGs(1) = 0.9*oleigGs(1);
    deseigGs(2) = 0.85*oleigGs(2);
    deseigGs(3) = 0.8*oleigGs(3);
    
    [LluGs,precGs,messageGs] = place(Gs.A',Gs.C',deseigGs);
    [LluGsscaled,precGsscaled,messageGsscaled] = place(Gs_scaled.A',Gs_scaled.C',deseigGs);
    
    %  LluGsnorm = norm(LluGs)
    % LluGsscalednorm = norm(LluGsscaled)
    
    % Are reduced-order quantities still dimensionless?
    LluGsscaledCsnorm = norm(LluGsscaled'*Gs_scaled.C)
    
    % Closed-loop eigenvalues
    %cleig = eig(jaccd-Llu'*C(measurementindex,:));
    cleigGs = sort(eig(Gs.A-LluGs'*Gs.C),'descend','ComparisonMethod','abs');
    cleigGsscaled = sort(eig(Gs_scaled.A-LluGsscaled'*Gs_scaled.C),'descend','ComparisonMethod','abs');
    max(cleigGs-sort(deseigGs,'descend','ComparisonMethod','abs'))
    max(cleigGsscaled-sort(deseigGs,'descend','ComparisonMethod','abs'))
    [cleigGs cleigGsscaled]
    
    % Display open and closed-loop eigenvalues
    [oleig cleig];
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
    %saveas(gcf,[lu_folder param '/kfeigs_' param '_' strtrim(statenames(measurementindex,:)) '_' num2str(bcl)])
    
    
    
    
    % % Compute "ground truth" linear system state, without noise
    % xgtlinall_noiseless = zeros(numstate,nsteps);
    % xgtlinall_noiseless(:,1) = xgtlin0;
    % ygtlin_noiseless = zeros(length(measurementindex),nsteps);
    % ygtlin_noiseless(1) = C(measurementindex,:)*xgtlinall_noiseless(:,1); % measurement without noise
    %
    % % Compute "ground truth" linear system state, with noise
    % xgtlinall = zeros(numstate,nsteps);
    % xgtlinall(:,1) = xgtlin0;
    % ygtlin = zeros(length(measurementindex),nsteps);
    % ygtlin(1) = C(measurementindex,:)*xgtlinall(:,1); % measurement with noise
    %
    % % Compute pseudo-KF (no noise) state estimate, closed loop
    % xkfall_noiseless = zeros(numstate,nsteps);
    % xkfall_noiseless(:,1) = xhat0;
    % xkfapostall_noiseless = zeros(numstate,nsteps);
    % xkfapostall_noiseless(:,1) = xkfall_noiseless(:,1) + M*(ygtlin_noiseless(1)-C(measurementindex,:)*xkfall_noiseless(:,1));
    %
    % % Compute pseudo-KF (no noise) state estimate, open loop
    % % This is almost the same as xgtlin_noiseless, but with the option of a non-random IC
    % xhatolall_noiseless = zeros(numstate,nsteps);
    % xhatolall_noiseless(:,1) = xhat0;
    %
    % % Compute KF state estimate, closed loop. xkfall(:, ii+1) = xhat(ii+1|ii),
    % % in other words, it is the a-priori predictor-corrector error.
    % xkfall = zeros(numstate,nsteps);
    % xkfall(:,1) = xhat0;
    % xkfapostall = zeros(numstate,nsteps);
    % xkfapostall(:,1) = xkfall(:,1) + M*(ygtlin(1)-C(measurementindex,:)*xkfall(:,1));
    %
    % err_apri_directsim_noiseless = zeros(numstate,nsteps);
    % err_apri_directsim_noiseless(:,1) = xgtlinall_noiseless(:,1)-xkfall_noiseless(:,1);
    % err_apri_directsim = zeros(numstate,nsteps);
    % err_apri_directsim(:,1) = xgtlinall(:,1)-xkfall(:,1);
    %
    %
    % for ii = 1:nsteps-1
    %     err_apri_directsim_noiseless(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*err_apri_directsim_noiseless(:,ii);
    %     err_apri_directsim(:,ii+1) = (jaccd-Lkf*C(measurementindex,:))*err_apri_directsim(:,ii) + Bw*w(:,ii) -Lkf*v(ii);
    %
    %     xgtlinall_noiseless(:,ii+1) = jaccd*xgtlinall_noiseless(:,ii);
    %
    %     ygtlin_noiseless(ii+1) = C(measurementindex,:)*xgtlinall_noiseless(:,ii+1); % measurement without noise
    %
    %
    %     xgtlinall(:,ii+1) = jaccd*xgtlinall(:,ii)+ Bw*w(:,ii);
    %
    %     ygtlin(ii+1) = C(measurementindex,:)*xgtlinall(:,ii+1) + v(ii+1); % measurement with noise
    %
    %
    %     xhatolall_noiseless(:,ii+1) = jaccd*xhatolall_noiseless(:,ii);
    %
    % %    xkfall_noiseless(:,ii+1) = jaccd*xkfall_noiseless(:,ii) + Lkf*(ygtlin_noiseless(ii)-C(measurementindex,:)*xkfall_noiseless(:,ii));
    %     xkfall_noiseless(:,ii+1) = jaccd*xkfapostall_noiseless(:,ii);
    %
    %     xkfapostall_noiseless(:,ii+1) = xkfall_noiseless(:,ii+1) + M*(ygtlin_noiseless(ii+1)-C(measurementindex,:)*xkfall_noiseless(:,ii+1));
    %
    %
    % %    xkfall(:,ii+1) = jaccd*xkfall(:,ii) + Lkf*(ygtlin(ii)-C(measurementindex,:)*xkfall(:,ii));
    %     xkfall(:,ii+1) = jaccd*xkfapostall(:,ii);
    %
    %     xkfapostall(:,ii+1) = xkfall(:,ii+1) + M*(ygtlin(ii+1)-C(measurementindex,:)*xkfall(:,ii+1));
    %
    % end
    % %ygtlin = C(measurementindex,:)*xgtlinall + v(1:nsteps); % measurement with noise
    %
    % % Compute estimation errors
    % eol_noiseless = xgtlinall_noiseless - xhatolall_noiseless;
    % eol = xgtlinall - xhatolall_noiseless;
    % ekf_noiseless = xgtlinall_noiseless - xkfall_noiseless;
    % ekfapost_noiseless = xgtlinall_noiseless - xkfapostall_noiseless;
    % ekf = xgtlinall - xkfall;
    % ekfapost = xgtlinall - xkfapostall;
    %
    %
    % Llunorm = norm(Llu) % size of gain is related to strength of feedback
    % Pnorm = norm(P) % steady-state a-priori cov norm
    % Znorm = norm(Z) % steady-state a-posteriori cov norm
    %
    % Lluscalednorm = norm(Lluscaled) % size of gain is related to strength of feedback
    % Pscalednorm = norm(Pscaled) % steady-state a-priori cov norm
    % Zscalednorm = norm(Zscaled) % steady-state a-posteriori cov norm
    %
    % maxabscleig = max(abs(cleig)) % size of least stable eigenvalue
    %
    % % Compute error covariance norms
    % ekfcov = cov(ekf');
    % size(ekfcov) % Make sure this is nxn, where n is the number of states.
    % % The presence or absence of the transpose operator matters. Cov assumes
    % % that columns are variables and rows are observations.
    % eolcovnorm_noiseless = norm(cov(eol_noiseless'));
    % ekfcovnorm_noiseless = norm(cov(ekf_noiseless'));
    % ekfapostcovnorm_noiseless = norm(cov(ekfapost_noiseless'));
    % eolcovnorm = norm(cov(eol'));
    % ekfcovnorm = norm(ekfcov);
    % ekfapostcovnorm = norm(cov(ekfapost'));
    %
    % % Check a value against manual computation
    % ekfvarcheck = sum((ekf(reconstructionindices(1),:)-mean(ekf(reconstructionindices(1),:))).^2)/(nsteps-1)
    % ekfcov(reconstructionindices(1),reconstructionindices(1))
    %
    % % Exclude transient period from covariances?
    % % No obvious fixed time when this should be; maybe first ten cycles?
    % %postransindex = 11;
    % % Revise this to use last 5% of simulation?
    % % eolcovnorm_noiseless_ss = norm(cov(eol_noiseless(:,postransindex:end)'));
    % % ekfcovnorm_noiseless_ss = norm(cov(ekf_noiseless(:,postransindex:end)'));
    % % ekfapostcovnorm_noiseless_ss = norm(cov(ekfapost_noiseless(:,postransindex:end)'));
    % % eolcovnorm_ss = norm(cov(eol(:,postransindex:end)'));
    % % ekfcovnorm_ss = norm(cov(ekf(:,postransindex:end)'));
    % % ekfapostcovnorm_ss = norm(cov(ekfapost(:,postransindex:end)'));
    % index_last5pct = round(0.05*nsteps);
    % index_last_n = 10;
    %
    % eolcovnorm_noiseless_ss = norm(cov(eol_noiseless(:,end-index_last_n:end)'));
    % ekfcovnorm_noiseless_ss = norm(cov(ekf_noiseless(:,end-index_last_n:end)'));
    % ekfapostcovnorm_noiseless_ss = norm(cov(ekfapost_noiseless(:,end-index_last_n:end)'));
    % eolcovnorm_ss = norm(cov(eol(:,end-index_last_n:end)'));
    % ekfcovnorm_ss = norm(cov(ekf(:,end-index_last_n:end)'));
    % ekfapostcovnorm_ss = norm(cov(ekfapost(:,end-index_last_n:end)'));
    %
    %
    % % Compute scaled error covariance norms
    % eolcovnorm_noiseless_scaled = norm(cov((Smat*eol_noiseless)'));
    % ekfcovnorm_noiseless_scaled = norm(cov((Smat*ekf_noiseless)'));
    % ekfapostcovnorm_noiseless_scaled = norm(cov((Smat*ekfapost_noiseless)'));
    % eolcovnorm_scaled = norm(cov((Smat*eol)'));
    % ekfcovnorm_scaled = norm(cov((Smat*ekf)'));
    % ekfapostcovnorm_scaled = norm(cov((Smat*ekfapost)'));
    %
    % eolcovnorm_noiseless_ss_scaled = norm(cov((Smat*eol_noiseless(:,end-index_last_n:end))'));
    % ekfcovnorm_noiseless_ss_scaled = norm(cov((Smat*ekf_noiseless(:,end-index_last_n:end))'));
    % ekfapostcovnorm_noiseless_ss_scaled = norm(cov((Smat*ekfapost_noiseless(:,end-index_last_n:end))'));
    % eolcovnorm_ss_scaled = norm(cov((Smat*eol(:,end-index_last_n:end))'));
    % ekfcovnorm_ss_scaled = norm(cov((Smat*ekf(:,end-index_last_n:end))'));
    % ekfapostcovnorm_ss_scaled = norm(cov((Smat*ekfapost(:,end-index_last_n:end))'));
    %
    %
    % eolcovnorm
    % ekfcovnorm
    % ekfapostcovnorm
    % ekfcovnorm_noiseless
    % ekfapostcovnorm_noiseless
    % ekfcovnorm_noiseless_ss
    % ekfapostcovnorm_noiseless_ss
    % ekfcovnorm_ss
    % ekfapostcovnorm_ss
    %
    % eolcovnorm_scaled
    % ekfcovnorm_scaled
    % ekfapostcovnorm_scaled
    % eolcovnorm_ss_scaled
    % ekfcovnorm_ss_scaled
    % ekfapostcovnorm_ss_scaled
    %
    %
    % % Compare errors that were computed in different ways
    % figure
    % hold on;
    % plot(1:nsteps, err_apri_directsim(1,:));
    % plot(1:nsteps, ekf(1,:),'m:x');
    % ylabel('V est error')
    % xlabel('time index')
    % legend('e(j+1) = (A-LC)e(j) + Bw*w(j) - Lv(j)', 'e(j) = xtrue(j) - xhat(j|j-1)')
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
    % saveas(gcf,[lu_folder param '/kfplot_wonoise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_noise_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_wonoise_normalized_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_noise_normalized_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_targetedvariables_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_targetedvariables_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_targetedvariables_normederr_cl_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    % saveas(gcf,[lu_folder param '/kfplot_normederr_overlay_' param '_' statenames(measurementindex) '_' num2str(bcl)])
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
    %     saveas(gcf,[lu_folder param '/kfplot_targetedvariables_ol_only_' param '_' statenames(measurementindex) '_' num2str(bcl)])
    % end
    %
    % %save kfwv1to10 *
    % eval(['save ' lu_folder param '/kalmanfile_measind_' num2str(measurementindex) '_bcl_' num2str(bcl) ' *'])
    %
    %
end