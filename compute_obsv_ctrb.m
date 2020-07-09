% Compute and plot observable eigenvalues for LRd Jacobians, for
% different choices of measured variables (i.e., outputs).
% In addition, compute controllable eigenvalues for different choices
% of inputs.
% Code adapted from lrd_batch_eig.m and computeeig.m.
%
%Laura's changes:
%   Edited to use 'nobalance' option in eigenvalue computation.
%   Changed certain variable and folder names to be more consistent with revisions in
%   other files, and/or to be more descriptive.
%   Removed scaling factor from B matrix.
%   For now, commented out portion about checking for repeated eigenvalues,
%   since we haven't yet settled on an appropriate threshold.
%   Recomputed eigenvalues internally due to concerns about mismatches.
%   Added portions to compute and store modal observability and
%   controllability magnitudes.

clear variables;

% I can't find a way to change LaTeX-interpreted labels back to the default
% font (Helvetica), so I can try changing everything else to the LaTeX
% font, cmr12. Update: this doesn't work either since the Matlab pdf
% printer doesn't interpret the font correctly. Set defaults globally
% instead.
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


%shiftstring = '' % if using default shift of data.dt
%shiftstring = '_shift0p75ms'
%shiftstring = '_shift6p5ms'
%shiftstring = '_shift7mVrepol'
%shiftstring = '_shift-50mVrepol'
%shiftstrings = {''};
%shiftstrings = {'_shift0p8Vnormdepol'}
%shiftstrings = {'_shift1Vnormdepol'}
shiftstrings = {'_shift0p001Vnormrepol'}
%shiftstrings = {'','_shift0p75ms','_shift6p5ms','_shift7mVrepol','_shift-50mVrepol'};
%shiftstrings = {'','_shift1Vnormdepol','_shift0p2Vnormdepol','_shift0p4Vnormdepol','_shift0p6Vnormdepol','_shift0p8Vnormdepol','_shift0p2Vnormrepol','_shift0p4Vnormrepol','_shift0p6Vnormrepol','_shift0p8Vnormrepol','_shift0p001Vnormrepol'};

numberoftopmodes = 8;%5; % plot n largest eigenvalues, where n = numberoftopmodes
%numberoftopmodes = 17; % plot n largest eigenvalues, where n = numberoftopmodes
eigthresh = 0.9 % eigenvalues with magnitudes greater than or equal to this threshold will
%eigthresh = 0.5; % eigenvalues with magnitudes greater than or equal to this threshold will
%eigthresh = 0; % eigenvalues with magnitudes greater than or equal to this threshold will
% count toward averaged observability and controllability measures
selected_logepsln = -5; % Base 10 log of step size used in Jacobian computation
debugflag = 0 % if nonzero, overwrite actual Jacobians with a simple
%structure that makes it easier to tell if obsv/ctrb values are right or
%wrong. Echo value to screen.

paramflag = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if paramflag
    param= 'adj';
else
    param = 'def';
end

ma = 80; fontsize = 28; linewidth = 2; % marker size, font size, line width

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');

%--------------------------------------------------------------------------%
%thresh = .0001; %Threshold for determining whether eigenvalues are identical
rankcutoff = 1e-14; % below this level, singular values don't contribute to the rank
%--------------------------------------------------------------------------%

% load approximate "state normalization" scaling matrix
load b1000fsolem12variable_amplitudes varamp
% Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
Smat = diag(1./varamp); % scaling matrix
Smatinv = inv(Smat); % only need to compute once

numstate = size(Smat,1);
% input matrix for all possible individual inputs
B = eye(numstate);

% output matrix for all possible individual measurements
C = eye(numstate);

umax = diag(ones(1,size(B,2))); % This is just a placeholder.
% For inputs, also need to know approximate maximum values of each element
% of input vector. Note that u is a deviational quantity that may or
% may not attain the size of the stimulus.
% Each umax diagonal element should probably be the size of
% deviational max for integrated system.
Bs = Smat*B*umax; % scaled B matrix
Cs = C*Smatinv; % scaled C matrix

% %halto = figure;
% halto_sc = figure;
% hold on;
%
% hm1o_sc = figure;
% hold on;
%
% hm2o_sc = figure;
% hold on;
%
% haltc_sc = figure;
% hold on;
%
% hm1c_sc = figure;
% hold on;
%
% hm2c_sc = figure;
% hold on;

for shiftctr = 1:length(shiftstrings)
    shiftstring = shiftstrings{shiftctr}
    
    %jacfolder = ['jacobians/' param '/']; % folder where jacobians are stored
    jacfolder = ['jacobians' shiftstring '/' param '/']; %This folder is where jacobians are stored
    
    %eigfolder = ['eigenvalues/' param '/']; %folder where eigenvalues are stored
    
    %ocfolder = 'ocvalues/'; %folder where obsv and ctrb values will be saved
    ocfolder = ['ocvalues' shiftstring '/']; %folder where obsv and ctrb values will be saved
    if ~exist(ocfolder,'dir')
        mkdir([ocfolder param '/'])
    end
    
    eval(['load ' jacfolder 'jacfile' num2str(selected_logepsln) ' *']); %Load Jacobians
    %eval(['load ' eigenvalues 'eigfile_' param ' allv alleigs']);
    
    %numstate = size(alljacs{1},1); % number of state variables
    
    %selected_bcls_for_oc = 200;
    % number of bcls
    nbcls = length(selected_bcls_for_fps);
    rankcheck = zeros(1,nbcls); % store right eigenvector matrix ranks here
    
    % Specify range of BCL indices for averaged values
    % startbcl = 110; %ms, this should be the longer BCL
    % endbcl = 100; %ms, this should be the shorter BCL
    % allbclindices = 1:nbcls;
    % startbclindex = allbclindices(selected_bcls_for_fps == startbcl);
    % endbclindex = allbclindices(selected_bcls_for_fps == endbcl);
    % bclindrangeforavgs = [startbclindex endbclindex] % startbclindex:endbclindex % choose spacing between indices here
    bclindrangeforavgs = 1:nbcls; % startbclindex:endbclindex % choose spacing between indices here
%    bclindrangeforavgs = [32 35]; % for BCLs of 210 and 180ms
    nbclsforavgs = length(bclindrangeforavgs);
    disp('Check BCL range, then press any key to continue')
    pause
    
    % Intitalize matrices
    eigcf.cval = cell(nbcls,numstate); % controllable eigenvalues
    eigof.oval = cell(nbcls,numstate); % observable eigenvalues
    eigcf.ncval = cell(nbcls,numstate); % uncontrollable eigenvalues
    eigof.noval = cell(nbcls,numstate); % unobservable eigenvalues
    eigcf.rank = zeros(nbcls,numstate); % controllability matrix rank
    eigof.rank = zeros(nbcls,numstate); % observability matrix rank
    eigcf.cvec = cell(nbcls, numstate); % controllable eigenvectors
    eigcf.ncvec = cell(nbcls, numstate); % uncontrollable eigenvectors
    eigof.ovec = cell(nbcls, numstate); % observable eigenvectors
    eigof.novec = cell(nbcls, numstate); % unobservable eigenvectors
    
    eigof.cdotv = cell(nbcls,1);
    eigof.cmagvmag = cell(nbcls,1);
    eigof.obsvmag = cell(nbcls,1);
    eigof.cdotv_sc = cell(nbcls,1);
    eigof.cmagvmag_sc = cell(nbcls,1);
    eigof.obsvmag_sc = cell(nbcls,1);
    
    
    % Initialize scaled-system matrices:
    eigcf.cval_sc = cell(nbcls,numstate);
    eigof.oval_sc = cell(nbcls,numstate);
    eigcf.ncval_sc = cell(nbcls,numstate);
    eigof.noval_sc = cell(nbcls,numstate);
    eigcf.rank_sc = zeros(nbcls,numstate);
    eigof.rank_sc = zeros(nbcls,numstate);
    eigcf.cvec_sc = cell(nbcls, numstate);
    eigcf.ncvec_sc = cell(nbcls, numstate);
    eigof.ovec_sc = cell(nbcls, numstate);
    eigof.novec_sc = cell(nbcls, numstate);
    
    eigcf.bdotw = cell(nbcls,1);
    eigcf.bmagwmag = cell(nbcls,1);
    eigcf.ctrbmag = cell(nbcls,1);
    eigcf.bdotw_sc = cell(nbcls,1);
    eigcf.bmagwmag_sc = cell(nbcls,1);
    eigcf.ctrbmag_sc = cell(nbcls,1);
    
    sortedeig = zeros(nbcls,numstate);
    sortedv = cell(nbcls,1);
    sortedw = cell(nbcls,1);
    
    sortedeigsc = zeros(nbcls,numstate);
    sortedvsc = cell(nbcls,1);
    sortedwsc = cell(nbcls,1);
    
    for ii = 1:nbcls
        %% Controllability
        bcl = selected_bcls_for_fps(ii);
        % print current BCL to screen
        disp(['BCL = ' num2str(bcl) ' ms'])
        
        % For debugging:
        if debugflag
            jaccd = eye(17);
            jaccd(17,:) = ones(1,17);
        else
            jaccd = alljacs{ii};
        end
        jaccd_sc = Smat*jaccd*Smatinv; % scaled Jacobian
        
        % The eigenvalues could be loaded instead of recomputed here,
        % but recomputing doesn't take much time, whereas reloading
        % increases risk of a mismatch between eigenvalues and Jacobians.
        [v,d,w]=eig(jaccd,'nobalance');
        rankcheck(ii) = rank(v); 
        % Sort eigenvalues and eigenvectors in descending order of size
        %    [sortedeig(ii,:),eigsortind] = sort(abs(diag(d)),'descend');
        [sortedeig(ii,:),eigsortind] = sort(diag(d),'descend','ComparisonMethod','abs');
        sortedv{ii} = v(:,eigsortind);
        sortedw{ii} = w(:,eigsortind);
        
        % Compute eigenvectors for scaled Jacobian:
        [vsc,dsc,wsc]=eig(jaccd_sc,'nobalance');
        %    [sortedeigsc(ii,:),eigsortindsc] = sort(abs(diag(dsc)),'descend');
        [sortedeigsc(ii,:),eigsortindsc] = sort(diag(dsc),'descend','ComparisonMethod','abs');
        sortedvsc{ii} = vsc(:,eigsortindsc);
        sortedwsc{ii} = wsc(:,eigsortindsc);
        
        %    rankcutoffdefault = max(size(jaccd))*eps(norm(jaccd));% This is Matlab's default rank cutoff. Can't easily control cutoffs in orth /null, so to interpret these results, would need to know default
        
        % Cycle through possible inputs
        for kb = 1:size(B,2)
            % Compute controllable and uncontrollable subspaces, keeping C
            % fixed
            [Abar,Bbar,Cbar,T,k] = ctrbf(jaccd,B(:,kb),C(1,:),rankcutoff);
            eigcf.rank(ii,kb) = sum(k); % rank of controllable subspace
            % eigenvalues that are controllable from this input
            eigcf.cval{ii,kb} = eig(Abar((numstate-eigcf.rank(ii,kb)+1):end,(numstate-eigcf.rank(ii,kb)+1):end),'nobalance'); % eigenvalues of Ac
            % eigenvalues that are not controllable from this input
            eigcf.ncval{ii,kb} = eig(Abar(1:(numstate-eigcf.rank(ii,kb)),1:(numstate-eigcf.rank(ii,kb))),'nobalance'); % eigenvalues of Anc
            % repeat calculations for scaled system
            [Abars,Bbars,Cbars,Ts,ks] = ctrbf(jaccd_sc,Bs(:,kb),Cs(1,:),rankcutoff);
            eigcf.rank_sc(ii,kb) = sum(ks);
            eigcf.cval_sc{ii,kb} = eig(Abars((numstate-eigcf.rank_sc(ii,kb)+1):end,(numstate-eigcf.rank_sc(ii,kb)+1):end),'nobalance'); % eigenvalues of Ac
            eigcf.ncval_sc{ii,kb} = eig(Abars(1:(numstate-eigcf.rank_sc(ii,kb)),1:(numstate-eigcf.rank_sc(ii,kb))),'nobalance'); % eigenvalues of Anc
            
            % Compute modal controllability values.
            % Warning: computational method should be altered for repeated
            % eigenvalues.
            for kk = 1:numstate
                eigcf.bdotw{ii}(kb,kk) = sortedw{ii}(:,kk)'*B(:,kb);
                eigcf.bmagwmag{ii}(kb,kk) = norm(B(:,kb))*norm(sortedw{ii}(:,kk));
                if eigcf.bmagwmag{ii}(kb,kk) > 0
                    eigcf.ctrbmag{ii}(kb,kk) = abs(eigcf.bdotw{ii}(kb,kk))/eigcf.bmagwmag{ii}(kb,kk);
                else
                    eigcf.ctrbmag{ii}(kb,kk) = NaN;
                end
                
                eigcf.bdotw_sc{ii}(kb,kk) = sortedwsc{ii}(:,kk)'*Bs(:,kb);
                eigcf.bmagwmag_sc{ii}(kb,kk) = norm(Bs(:,kb))*norm(sortedwsc{ii}(:,kk));
                if eigcf.bmagwmag_sc{ii}(kb,kk) > 0
                    eigcf.ctrbmag_sc{ii}(kb,kk) = abs(eigcf.bdotw_sc{ii}(kb,kk))/eigcf.bmagwmag_sc{ii}(kb,kk);
                else
                    eigcf.ctrbmag_sc{ii}(kb,kk) = NaN;
                end
            end
        end
        
        
        %% Observability
        % Cycle through possible outputs
        for kc = 1:size(C,1)
            % Compute observable and unobservable subspaces, keeping B
            % fixed
            [Abar,Bbar,Cbar,T,k] = obsvf(jaccd,B(:,1),C(kc,:),rankcutoff);
            eigof.rank(ii,kc) = sum(k); % rank of observable subspace
            % eigenvalues that are observable from this output
            eigof.oval{ii,kc} = eig(Abar((numstate-eigof.rank(ii,kc)+1):end,(numstate-eigof.rank(ii,kc)+1):end),'nobalance'); % eigenvalues of Ao
            %         eigof.ovec{ii ,kc}= zeros(17, eigof.rank(ii, kc));
            %         for j=1:eigof.rank(ii, kc)
            %             f= allv{1, ii}(:,abs(real(alleigs{1, ii})-eigof.oval{ii, kc}(j)) < abs(real(eigof.oval{ii, kc}(j)*thresh)));
            %             if ~isempty(f)
            %                 eigof.ovec{ii,kc}(:, j)= f;
            %             end
            %         end
            % eigenvalues that are not observable from this output
            eigof.noval{ii,kc} = eig(Abar(1:(numstate-eigof.rank(ii,kc)),1:(numstate-eigof.rank(ii,kc))),'nobalance'); % eigenvalues of Ano
            %         eigof.novec{ii ,kc}= zeros(17, numstate-eigof.rank(ii, kc));
            %         for j=1:(numstate- eigof.rank(ii, kc))
            %             f= allv{1, ii}(:,abs(real(alleigs{1, ii})-eigof.noval{ii, kc}(j)) < abs(real(eigof.noval{ii, kc}(j)*thresh)));
            %             if ~isempty(f)
            %                 eigof.novec{ii,kc}(:, j)=  f;
            %             end
            %         end
            % repeat calculations for scaled system
            [Abars,Bbars,Cbars,Ts,ks] = obsvf(jaccd_sc,Bs(:,1),Cs(kc,:),rankcutoff);
            eigof.rank_sc(ii,kc) = sum(ks);
            eigof.oval_sc{ii,kc} = eig(Abars((numstate-eigof.rank_sc(ii,kc)+1):end,(numstate-eigof.rank_sc(ii,kc)+1):end),'nobalance'); % eigenvalues of Ao
            %        eigof.ovec_sc{ii ,kc}= zeros(17, eigof.rank_sc(ii, kc));
            %         for j=1:eigof.rank_sc(ii, kc)
            %             f= allv{1, ii}(:, abs(real(alleigs{1, ii})-eigof.oval_sc{ii, kc}(j)) < abs(real(eigof.oval_sc{ii, kc}(j)*thresh)));
            %             if ~isempty(f)
            %                 eigof.ovec_sc{ii,kc}(:, j)= f;
            %             end
            %         end
            eigof.noval_sc{ii,kc} = eig(Abars(1:(numstate-eigof.rank_sc(ii,kc)),1:(numstate-eigof.rank_sc(ii,kc))),'nobalance'); % eigenvalues of Ano
            %         eigof.novec_sc{ii ,kc}= zeros(17, numstate-eigof.rank_sc(ii, kc));
            %         for j=1:(numstate-eigof.rank_sc(ii, kc))
            %             f= allv{1, ii}(:, abs(real(alleigs{1, ii})-eigof.noval_sc{ii, kc}(j)) < abs(real(eigof.noval_sc{ii, kc}(j)*thresh)));
            %             if ~isempty(f)
            %                 eigof.novec_sc{ii,kc}(:, j)= f;
            %             end
            %         end
            
            % Compute modal observability values.
            % Warning: computational method should be altered for repeated
            % eigenvalues.
            for kk = 1:numstate
                eigof.cdotv{ii}(kc,kk) = C(kc,:)*sortedv{ii}(:,kk);
                eigof.cmagvmag{ii}(kc,kk) = norm(C(kc,:))*norm(sortedv{ii}(:,kk));
                if eigof.cmagvmag{ii}(kc,kk) > 0
                    eigof.obsvmag{ii}(kc,kk) = abs(eigof.cdotv{ii}(kc,kk))/eigof.cmagvmag{ii}(kc,kk);
                else
                    eigof.obsvmag{ii}(kc,kk) = NaN;
                end
                
                eigof.cdotv_sc{ii}(kc,kk) = Cs(kc,:)*sortedvsc{ii}(:,kk);
                eigof.cmagvmag_sc{ii}(kc,kk) = norm(Cs(kc,:))*norm(sortedvsc{ii}(:,kk));
                if eigof.cmagvmag_sc{ii}(kc,kk) > 0
                    eigof.obsvmag_sc{ii}(kc,kk) = abs(eigof.cdotv_sc{ii}(kc,kk))/eigof.cmagvmag_sc{ii}(kc,kk);
                else
                    eigof.obsvmag_sc{ii}(kc,kk) = NaN;
                end
            end
        end
    end
    
    % % Compute averages over BCLs and modes
    largesteigindices = cell(nbcls,1); % store indices of largest eigenvalues here
    numlargesteigs = zeros(nbcls,1); % store number of largest eigenvalues here
    largesteigindices_sc = cell(nbcls,1); % store indices of largest eigenvalues here
    numlargesteigs_sc = zeros(nbcls,1); % store number of largest eigenvalues here
    % largesteigindices = cell(nrbi,1); % store indices of largest eigenvalues here
    % numlargesteigs = zeros(nrbi,1); % store number of largest eigenvalues here
    % largesteigindices_sc = cell(nrbi,1); % store indices of largest eigenvalues here
    % numlargesteigs_sc = zeros(nrbi,1); % store number of largest eigenvalues here
    
    % We're averaging over cell indices, so can't use mean function
    % directly. There's probably a better way, but for now, initalize
    % zero matrices, add up matrices, and divide by the number of bcls.
    tempoa = zeros(size(C,1),numstate);
    tempca = zeros(size(B,2),numstate);
    tempoasc = zeros(size(C,1),numstate);
    tempcasc = zeros(size(B,2),numstate);
    
    % set up zero matrices for averages over BCLs and largest eigenvalues
    obsv_avgd_over_modes = zeros(size(C,1),nbcls);
    ctrb_avgd_over_modes = zeros(size(B,2),nbcls);
    obsv_avgd_over_modes_sc = zeros(size(C,1),nbcls);
    ctrb_avgd_over_modes_sc = zeros(size(B,2),nbcls);
    tempoam = zeros(size(C,1),1);
    tempcam = zeros(size(B,2),1);
    tempoamsc = zeros(size(C,1),1);
    tempcamsc = zeros(size(B,2),1);
    
    %for ii = 1:nbcls
    for ii = bclindrangeforavgs %startbclindex:endbclindex
        tempoa = tempoa + eigof.obsvmag{ii};
        tempca = tempca + eigcf.ctrbmag{ii};
        tempoasc = tempoasc + eigof.obsvmag_sc{ii};
        tempcasc = tempcasc + eigcf.ctrbmag_sc{ii};
        % Based on test_linear_luen.m results, for scaled systems,
        % the modal ctrb/obsv numerators are more consistent with placement
        % gain sizes than the normalized dot products. For the unscaled
        % systems, B and C have unit size (and eigenvectors are normalized), but this isn't true for the scaled
        % systems, so for unscaled systems the normalization step shouldn't matter.
        %tempoasc = tempoasc + abs(eigof.cdotv_sc{ii});
        %tempcasc = tempcasc + abs(eigcf.bdotw_sc{ii});
        
        largesteigindices{ii} = find(abs(sortedeig(ii,:)) >= eigthresh);
        largesteigindices_sc{ii} = find(abs(sortedeigsc(ii,:)) >= eigthresh);
        numlargesteigs(ii) = length(largesteigindices{ii});
        numlargesteigs_sc(ii) = length(largesteigindices_sc{ii});
        
        obsv_avgd_over_modes(:,ii) = mean(eigof.obsvmag{ii}(:,largesteigindices{ii}),2);
        ctrb_avgd_over_modes(:,ii) = mean(eigcf.ctrbmag{ii}(:,largesteigindices{ii}),2);
        obsv_avgd_over_modes_sc(:,ii) = mean(eigof.obsvmag_sc{ii}(:,largesteigindices_sc{ii}),2);
        ctrb_avgd_over_modes_sc(:,ii) = mean(eigcf.ctrbmag_sc{ii}(:,largesteigindices_sc{ii}),2);
        %         % Based on test_linear_luen.m results, for scaled systems,
        %         % the modal ctrb/obsv numerators are more consistent with placement
        %         % gain sizes than the normalized dot products. For the unscaled
        %         % systems, B and C have unit size (and eigenvectors are normalized), but this isn't true for the scaled
        %         % systems, so for unscaled systems the normalization step shouldn't matter.
        %         obsv_avgd_over_modes_sc(:,ii) = mean(abs(eigof.cdotv_sc{ii}(:,largesteigindices_sc{ii})),2);
        %         ctrb_avgd_over_modes_sc(:,ii) = mean(abs(eigcf.bdotw_sc{ii}(:,largesteigindices_sc{ii})),2);
        tempoam = tempoam + obsv_avgd_over_modes(:,ii);
        tempcam = tempcam + ctrb_avgd_over_modes(:,ii);
        tempoamsc = tempoamsc + obsv_avgd_over_modes_sc(:,ii);
        tempcamsc = tempcamsc + ctrb_avgd_over_modes_sc(:,ii);
    end
    
    % obsv_avgd_over_bcls = tempoa/nbcls;
    % ctrb_avgd_over_bcls = tempca/nbcls;
    % obsv_avgd_over_bcls_sc = tempoasc/nbcls;
    % ctrb_avgd_over_bcls_sc = tempcasc/nbcls;
    %
    % obsv_avgd_over_bcls_and_modes = tempoam/nbcls;
    % ctrb_avgd_over_bcls_and_modes = tempcam/nbcls;
    % obsv_avgd_over_bcls_and_modes_sc = tempoamsc/nbcls;
    % ctrb_avgd_over_bcls_and_modes_sc = tempcamsc/nbcls;
    
    obsv_avgd_over_bcls = tempoa/nbclsforavgs;
    ctrb_avgd_over_bcls = tempca/nbclsforavgs;
    obsv_avgd_over_bcls_sc = tempoasc/nbclsforavgs;
    ctrb_avgd_over_bcls_sc = tempcasc/nbclsforavgs;
    
    obsv_avgd_over_bcls_and_modes = tempoam/nbclsforavgs;
    ctrb_avgd_over_bcls_and_modes = tempcam/nbclsforavgs;
    obsv_avgd_over_bcls_and_modes_sc = tempoamsc/nbclsforavgs;
    ctrb_avgd_over_bcls_and_modes_sc = tempcamsc/nbclsforavgs;
    
    
    % obsv_avgd_over_bcls_and_modes = mean(obsv_avgd_over_bcls(:,1:numberoftopmodes),2);
    % ctrb_avgd_over_bcls_and_modes = mean(ctrb_avgd_over_bcls(:,1:numberoftopmodes),2);
    % obsv_avgd_over_bcls_and_modes_sc = mean(obsv_avgd_over_bcls_sc(:,1:numberoftopmodes),2);
    % ctrb_avgd_over_bcls_and_modes_sc = mean(ctrb_avgd_over_bcls_sc(:,1:numberoftopmodes),2);

    if 0
    disp('Unscaled observability for top mode, averaged over bcls:')
    [sortedobsv,obsvsortindex] = sort(abs(obsv_avgd_over_bcls(:,1)),'descend');
    statenames_latex(obsvsortindex,:)
    %obsv_avgd_over_bcls(obsvsortindex) % typo since there should be 2 indices
    %obsv_avgd_over_bcls(obsvsortindex,1) % this works for top mode
    sortedobsv %should be same as above, since average is over magnitudes
    
    disp('Unscaled observability, averaged over bcls and modes above threshold:')
    [sortedobsvam,obsvsortindexam] = sort(abs(obsv_avgd_over_bcls_and_modes),'descend');
    statenames_latex(obsvsortindexam,:)
    obsv_avgd_over_bcls_and_modes(obsvsortindexam)
    end
    
    disp('Scaled observability for top mode, averaged over bcls:')
    [sortedobsv_sc,obsvsortindex_sc] = sort(abs(obsv_avgd_over_bcls_sc(:,1)),'descend');
    statenames_latex(obsvsortindex_sc,:)
    %obsv_avgd_over_bcls_sc(obsvsortindex_sc)
    obsv_avgd_over_bcls_sc(obsvsortindex_sc,1)
    
    disp('Scaled observability, averaged over bcls and modes above threshold:')
    [sortedobsvam_sc,obsvsortindexam_sc] = sort(abs(obsv_avgd_over_bcls_and_modes_sc),'descend');
    statenames_latex(obsvsortindexam_sc,:)
    obsv_avgd_over_bcls_and_modes_sc(obsvsortindexam_sc)
    
    if 0 % bug, why is this labeled as a modal average?
        disp('Unscaled observability, averaged over bcls and modes above threshold:')
        % Export LaTeX tables
        obsv_avgd_over_bcls_table = cell(numstate,2);
        obsv_avgd_over_bcls_table(:,1) = cellstr(statenames_latex(obsvsortindex,:));
        obsv_avgd_over_bcls_table(:,2) = cellstr(num2str(log10(obsv_avgd_over_bcls(obsvsortindex)),'%f'));
        latextableinput.data = obsv_avgd_over_bcls_table;
        latexTable(latextableinput)
        
        disp('Scaled observability, averaged over bcls and modes above threshold:')
        obsv_avgd_over_bcls_sc_table = cell(numstate,2);
        obsv_avgd_over_bcls_sc_table(:,1) = cellstr(statenames_latex(obsvsortindex_sc,:));
        %obsv_avgd_over_bcls_sc_table(:,2) = cellstr(num2str(obsv_avgd_over_bcls_sc(obsvsortindex_sc),'%0.2e'));
        obsv_avgd_over_bcls_sc_table(:,2) = cellstr(num2str(log10(obsv_avgd_over_bcls_sc(obsvsortindex_sc)),'%f'));
        %latex(cell2sym(obsv_avgd_over_bcls_sc_table))
        % Wasn't able to get native latex function to format numbers the way I
        % wanted it to, so I'm using a custom script, latexTable.m, instead:
        latextableinput.data = obsv_avgd_over_bcls_sc_table;
        latexTable(latextableinput)
    end
    
%     disp('Unscaled observability, averaged over bcls and modes above threshold:')
%     % Export LaTeX tables
%     obsv_avgd_over_bcls_and_modes_table = cell(numstate,2);
%     obsv_avgd_over_bcls_and_modes_table(:,1) = cellstr(statenames_latex(obsvsortindexam,:));
%     obsv_avgd_over_bcls_and_modes_table(:,2) = cellstr(num2str(log10(obsv_avgd_over_bcls_and_modes(obsvsortindexam)),'%f'));
%     latextableinput.data = obsv_avgd_over_bcls_and_modes_table;
%     latexTable(latextableinput)
    
    disp('Scaled observability, averaged over bcls and modes above threshold:')
    obsv_avgd_over_bcls_and_modes_sc_table = cell(numstate,2);
    obsv_avgd_over_bcls_and_modes_sc_table(:,1) = cellstr(statenames_latex(obsvsortindexam_sc,:));
    %obsv_avgd_over_bcls_sc_table(:,2) = cellstr(num2str(obsv_avgd_over_bcls_sc(obsvsortindex_sc),'%0.2e'));
    obsv_avgd_over_bcls_and_modes_sc_table(:,2) = cellstr(num2str(log10(obsv_avgd_over_bcls_and_modes_sc(obsvsortindexam_sc)),'%f'));
    %latex(cell2sym(obsv_avgd_over_bcls_sc_table))
    % Wasn't able to get native latex function to format numbers the way I
    % wanted it to, so I'm using a custom script, latexTable.m, instead:
    latextableinput.data = obsv_avgd_over_bcls_and_modes_sc_table;
    latexTable(latextableinput)
    
%     % Save table (read into Excel later) 
%     Update: writecell only works with 2019 version of Matlab and later. 
%     obsv_avgd_over_bcls_and_modes_sc_filedata = cell(numstate,3);
%     ofilename = [ocfolder param '/obsv_avgd_over_bcls_and_modes_sc_table_b' num2str(selectedbclsforfps(bclindrangeforavgs(end))) 'to' num2str(selectedbclsforfps(bclindrangeforavgs(1))) ' .txt']
%     writecell(obsv_avgd_over_bcls_and_modes_sc_table, ofilename);
    eigthreshstring = num2str(eigthresh,'%1.2f');
    xlsfname = [ocfolder param '/octable_sc_' param '_b' num2str(bclindrangeforavgs(1)) 'to' num2str(bclindrangeforavgs(end)) '_eigthresh' eigthreshstring(1) 'p' eigthreshstring(3:end) shiftstring '.xlsx'];
    % Need to embed string column headers in cells to get more than the
    % first character to appear in the Excel cell. 
    xlswrite(xlsfname,{['BCLs (ms), ' shiftstring ', ' param ', ' num2str(eigthresh)]},'obsv','A1:A1');
    xlswrite(xlsfname,{'var'},'obsv','B1:B1'); 
    xlswrite(xlsfname,{'log10 avg obsv'},'obsv','C1:C1'); 
    xlswrite(xlsfname,{'avg obsv'},'obsv','D1:D1'); 
    xlswrite(xlsfname, selected_bcls_for_fps(bclindrangeforavgs)','obsv','A2');
    xlswrite(xlsfname,obsv_avgd_over_bcls_and_modes_sc_table(:,1),'obsv','B2'); 
    xlswrite(xlsfname,obsv_avgd_over_bcls_and_modes_sc_table(:,2),'obsv','C2'); 
    xlswrite(xlsfname,obsv_avgd_over_bcls_and_modes_sc(obsvsortindexam_sc),'obsv','D2'); 
    
    % Controllability
    if 0
    disp('Unscaled controllability for top mode, averaged over bcls:')
    [sortedctrb,ctrbsortindex] = sort(abs(ctrb_avgd_over_bcls(:,1)),'descend');
    statenames_latex(ctrbsortindex,:)
    %ctrb_avgd_over_bcls(ctrbsortindex)
    ctrb_avgd_over_bcls(ctrbsortindex,1)
    
    disp('Unscaled controllability, averaged over bcls and modes above threshold:')
    [sortedctrbam,ctrbsortindexam] = sort(abs(ctrb_avgd_over_bcls_and_modes),'descend');
    statenames_latex(ctrbsortindexam,:)
    ctrb_avgd_over_bcls_and_modes(ctrbsortindexam)
    end
    
    disp('Scaled controllability for top mode, averaged over bcls:')
    [sortedctrb_sc,ctrbsortindex_sc] = sort(abs(ctrb_avgd_over_bcls_sc(:,1)),'descend');
    statenames_latex(ctrbsortindex_sc,:)
    %ctrb_avgd_over_bcls_sc(ctrbsortindex_sc)
    ctrb_avgd_over_bcls_sc(ctrbsortindex_sc,1)
    
    disp('Scaled controllability, averaged over bcls and modes above threshold:')
    [sortedctrbam_sc,ctrbsortindexam_sc] = sort(abs(ctrb_avgd_over_bcls_and_modes_sc),'descend');
    statenames_latex(ctrbsortindexam_sc,:)
    ctrb_avgd_over_bcls_and_modes_sc(ctrbsortindexam_sc)
    
    if 0 % bug, should refer to modal averages
        disp('Unscaled controllabilty, averaged over bcls and modes above threshold:')
        % Export LaTeX tables
        ctrb_avgd_over_bcls_table = cell(numstate,2);
        ctrb_avgd_over_bcls_table(:,1) = cellstr(statenames_latex(ctrbsortindex,:));
        ctrb_avgd_over_bcls_table(:,2) = cellstr(num2str(log10(ctrb_avgd_over_bcls(ctrbsortindex)),'%f'));
        latextableinput.data = ctrb_avgd_over_bcls_table;
        latexTable(latextableinput)
        
        disp('Scaled controllabilty, averaged over bcls and modes above threshold:')
        ctrb_avgd_over_bcls_sc_table = cell(numstate,2);
        ctrb_avgd_over_bcls_sc_table(:,1) = cellstr(statenames_latex(ctrbsortindex_sc,:));
        ctrb_avgd_over_bcls_sc_table(:,2) = cellstr(num2str(log10(ctrb_avgd_over_bcls_sc(ctrbsortindex_sc)),'%f'));
        latextableinput.data = ctrb_avgd_over_bcls_sc_table;
        latexTable(latextableinput)
    end
%     disp('Unscaled controllabilty, averaged over bcls and modes above threshold:')
%     % Export LaTeX tables
%     ctrb_avgd_over_bcls_and_modes_table = cell(numstate,2);
%     ctrb_avgd_over_bcls_and_modes_table(:,1) = cellstr(statenames_latex(ctrbsortindexam,:));
%     ctrb_avgd_over_bcls_and_modes_table(:,2) = cellstr(num2str(log10(ctrb_avgd_over_bcls_and_modes(ctrbsortindexam)),'%f'));
%     latextableinput.data = ctrb_avgd_over_bcls_and_modes_table;
%     latexTable(latextableinput)
    
    disp('Scaled controllabilty, averaged over bcls and modes above threshold:')
    ctrb_avgd_over_bcls_and_modes_sc_table = cell(numstate,2);
    ctrb_avgd_over_bcls_and_modes_sc_table(:,1) = cellstr(statenames_latex(ctrbsortindexam_sc,:));
    ctrb_avgd_over_bcls_and_modes_sc_table(:,2) = cellstr(num2str(log10(ctrb_avgd_over_bcls_and_modes_sc(ctrbsortindexam_sc)),'%f'));
    latextableinput.data = ctrb_avgd_over_bcls_and_modes_sc_table;
    latexTable(latextableinput)

    % Should use writecell starting with 2019 Matlab
    %    writecell(ctrb_avgd_over_bcls_and_modes_sc_table, [ocfolder param '/ctrb_avgd_over_bcls_and_modes_sc_table.txt']);
    % Need to embed string column headers in cells to get more than the
    % first character to appear in the Excel cell. 
    xlswrite(xlsfname,{['BCLs (ms), ' shiftstring ', ' param ', ' num2str(eigthresh)]},'ctrb','A1:A1');
    xlswrite(xlsfname,{'var'},'ctrb','B1:B1'); 
    xlswrite(xlsfname,{'log10 avg ctrb'},'ctrb','C1:C1'); 
    xlswrite(xlsfname,{'avg ctrb'},'ctrb','D1:D1'); 
    xlswrite(xlsfname, selected_bcls_for_fps(bclindrangeforavgs)','ctrb','A2');
    xlswrite(xlsfname,ctrb_avgd_over_bcls_and_modes_sc_table(:,1),'ctrb','B2'); 
    xlswrite(xlsfname,ctrb_avgd_over_bcls_and_modes_sc_table(:,2),'ctrb','C2'); 
    xlswrite(xlsfname,ctrb_avgd_over_bcls_and_modes_sc(ctrbsortindexam_sc),'ctrb','D2'); 
   
    
    %clear paramflag bcl B C Cs eigfolder f fs i j jacfolder k kb kc ks ms Smat Smatinv umax varamp;
    if ~debugflag
%        eval(['save ' ocfolder param '/ocfile *'])
        eval(['save ' ocfolder param '/ocfile_b' num2str(bclindrangeforavgs(1)) 'to' num2str(bclindrangeforavgs(end)) '_eigthresh' eigthreshstring(1) 'p' eigthreshstring(3:end) ' *'])
    end
    %Use plot_obs_cont to plot these results

    %kc = 10
    %kb = 1;
    kc = input('Enter a measurement index from 1 to 17: ');
    kb = input('Enter a control input index from 1 to 17: ');

    % plot symbols
    symbols = char('bo-','rs-','gp-','m*-','k^-','cv-','yh-');
    %
        
    if 1
        % Plot scaled values
        figure
        plot(selected_bcls_for_fps,obsv_avgd_over_modes_sc(kc,:),'linewidth',linewidth);
        xlabel('$T$, ms')
        ylabel(['$|cos(\phi)|$ averaged over $|\lambda| > $' num2str(eigthresh)])
        title(['Avg. obsv. magnitude, meas. ' strtrim(statenames_latex(kc,:)) ', shift = ' shiftstring])
        set(gca,'fontsize',fontsize)
        set(gcf,'Position',[433 243 938 615])
        saveas(gcf,[ocfolder param '/obsvavgplot_' param '_measindex_' num2str(kc) '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)])
        
        figure
        plot(selected_bcls_for_fps,mean(obsv_avgd_over_modes_sc),'linewidth',linewidth);
        xlabel('$T$, ms')
        ylabel(['$|cos(\phi)| averaged over |\lambda| > $' num2str(eigthresh)])
        title(['Avg. obsv. magnitude over all measurements, ' 'shift = ' shiftstring])
        set(gca,'fontsize',fontsize)
        set(gcf,'Position',[433 243 938 615])
        saveas(gcf,[ocfolder param '/obsvavgplot_' param '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)])
        
        % Plot scaled values (ctrb)
        figure
        plot(selected_bcls_for_fps,ctrb_avgd_over_modes_sc(kb,:),'linewidth',linewidth);
        xlabel('$T$, ms')
        ylabel(['$|cos(\theta)| averaged over |\lambda| > $' num2str(eigthresh)])
        title(['Avg. ctrb. magnitude, control index = ' num2str(kb) ' shift = ' shiftstring])
        set(gca,'fontsize',fontsize)
        set(gcf,'Position',[433 243 938 615])
        saveas(gcf,[ocfolder param '/ctrbavgplot_' param '_ctrlindex_' num2str(kb) '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)])
        
        figure
        plot(selected_bcls_for_fps,mean(ctrb_avgd_over_modes_sc),'linewidth',linewidth);
        xlabel('$T$, ms')
        ylabel(['$|cos(\theta)| averaged over |\lambda| > $' num2str(eigthresh)])
        title(['Avg. ctrb. magnitude over all control inputs, ' 'shift = ' shiftstring])
        set(gca,'fontsize',fontsize)
        set(gcf,'Position',[433 243 938 615])
        saveas(gcf,[ocfolder param '/ctrbavgplot_' param '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)])
    end

        % hor = figure;
    % hold on;
    %
    % hoi = figure;
    % hold on;
    %
    hoa = figure;
    hold on;
    %
    % hcr = figure;
    % hold on;
    %
    % hci = figure;
    % hold on;
    %
    hca = figure;
    hold on;

    for ii = 1:nbcls
        if 0
            for jj = 1:numberoftopmodes
                figure(hor)
                plot(selected_bcls_for_fps(ii),real(eigof.cdotv_sc{ii}(kc,jj))/eigof.cmagvmag_sc{ii}(kc,jj),symbols(jj,:));
                figure(hoi)
                plot(selected_bcls_for_fps(ii),imag(eigof.cdotv_sc{ii}(kc,jj))/eigof.cmagvmag_sc{ii}(kc,jj),symbols(jj,:));
                
                
                figure(hcr)
                plot(selected_bcls_for_fps(ii),real(eigcf.bdotw_sc{ii}(kb,jj))/eigcf.bmagwmag_sc{ii}(kb,jj),symbols(jj,:));
                figure(hci)
                plot(selected_bcls_for_fps(ii),imag(eigcf.bdotw_sc{ii}(kb,jj))/eigcf.bmagwmag_sc{ii}(kb,jj),symbols(jj,:));
            end
        end
        for jj = 1:numberoftopmodes%numstate
            figure(hoa)
            so=scatter(selected_bcls_for_fps(ii),abs(sortedeigsc(ii,jj)),ma,eigof.obsvmag_sc{ii}(kc,jj));
            % Replace with numerators:
            %           so=scatter(selected_bcls_for_fps(ii),abs(sortedeigsc(ii,jj)),ma,abs(eigof.cdotv_sc{ii}(kc,jj)));
            so.LineWidth = linewidth;
            % plot negative-real-part eigenvalues with a different
            % symbol (Elizabeth's suggestion)
            if real(sortedeigsc(ii,jj)) < 0
                so.Marker = 'd';
            end
            
            figure(hca)
            sc = scatter(selected_bcls_for_fps(ii),abs(sortedeigsc(ii,jj)),ma,eigcf.ctrbmag_sc{ii}(kb,jj));
            % Replace with numerators
            %            sc = scatter(selected_bcls_for_fps(ii),abs(sortedeigsc(ii,jj)),ma, abs(eigcf.bdotw_sc{ii}(kb,kk)));
            set(sc,'linewidth',linewidth)
            sc.LineWidth = linewidth;
            % plot negative-real-part eigenvalues with a different
            % symbol (Elizabeth's suggestion)
            if real(sortedeigsc(ii,jj)) < 0
                sc.Marker = 'd';
            end
        end
        
    end
    
    % figure(hor)
    % xlabel('BCL, ms')
    % ylabel('Re C*v')
    %
    % figure(hoi)
    % xlabel('BCL, ms')
    % ylabel('Im C*v')

    figure(hoa)
    colormap winter
    %xlabel('BCL, ms')
    %ylabel('Eigenvalue modulus, |\lambda|')
    c = colorbar;
    %c.Label.String = '|Cv|';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
%    c.Label.String = '$|\cos \phi_{ki}|$';
    tempstr = strtrim(statenames_latex(kc,:)); 
    c.Label.String = ['$|\cos \phi_{' tempstr(2:end-1) ',\; k}|$'];
    %    title({['Observability magnitude, meas. index = ' num2str(kc)], [' shift = ' shiftstring]})
%    title({['Observability magnitude, meas. index = ' num2str(kc) ' shift = ' shiftstring]})
    title({['Observability magnitude, meas. ' tempstr ', shift = ' shiftstring(2:end)]}, 'Interpreter','latex')
    %set(gca,'fontsize',fontsize)
    %set(gcf,'Position',[550 243 821 615])
    %saveas(gcf,[ocfolder param '/obsveigplot_' param '_measindex_' num2str(kc)])
    ylabel('      $|\lambda|$');
    xlabel('$T$, ms');
    % if param == 'def'
    %     text(750,0.9, 'default ','fontsize',fontsize)
    %     figname = 'obsveigdef';
    % elseif param == 'adj'
    %     text(750,0.9, 'adjusted ','fontsize',fontsize)
    %     figname = 'obsveigadj';
    % end
    axis([40 1030 0 1.25])
    %leghandle=legend([p(1), p(2)], 'Re(\lambda) \geq 0', 'Re(\lambda) < 0','Location','east');
    %set(leghandle,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.55    0.79])
    set(gca,'Position', [0.165    0.1618    0.6    0.7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %print(fig,figname,'-dpdf')
    %saveas(fig,figname)
    saveas(fig,[ocfolder param '/obsveigplot_' param '_measindex_' num2str(kc) '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)])
    print(fig,[ocfolder param '/obsveigplot_' param '_measindex_' num2str(kc) '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)],'-dpdf')
    
    
    % figure(hcr)
    % xlabel('BCL, ms')
    % ylabel('Re w^T*B')
    %
    % figure(hci)
    % xlabel('BCL, ms')
    % ylabel('Im w^T*B')
    
    figure(hca)
    colormap winter
    %xlabel('BCL, ms')
    %ylabel('Eigenvalue modulus, |\lambda|')
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    tempstr = strtrim(statenames_latex(kb,:)); 
    %c.Label.String = '|w^TB|';
%    c.Label.String = '|cos \theta_{ki}|';
    c.Label.String = ['$|\cos \theta_{' tempstr(2:end-1) ',\; k}|$'];
    title(['Controllability magnitude, input index = ' num2str(kb) ' shift = ' shiftstring(2:end)],'Interpreter','latex')
    %set(gca,'fontsize',fontsize)
    %set(gcf,'Position',[550 243 821 615])
    %saveas(gcf,[ocfolder param '/ctrbeigplot_' param '_ctrlindex_' num2str(kb)])
    ylabel('      $|\lambda|$');
    xlabel('$T$, ms');
    % if param == 'def'
    %     text(750,0.9, 'default ','fontsize',fontsize)
    %     figname = 'obsveigdef';
    % elseif param == 'adj'
    %     text(750,0.9, 'adjusted ','fontsize',fontsize)
    %     figname = 'obsveigadj';
    % end
    axis([40 1030 0 1.25])
    %leghandle=legend([p(1), p(2)], 'Re(\lambda) \geq 0', 'Re(\lambda) < 0','Location','east');
    %set(leghandle,'Position',[0.6584    0.4090    0.2304    0.1152])
    set(gca,'FontSize',fontsize)
    set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.55    0.79])
    set(gca,'Position', [0.165    0.1618    0.6    0.7632])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %print(fig,figname,'-dpdf')
    %saveas(fig,figname)
    saveas(fig,[ocfolder param '/ctrbeigplot_' param '_ctrlindex_' num2str(kb) '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)])
    print(fig,[ocfolder param '/ctrbeigplot_' param '_ctrlindex_' num2str(kb) '_eigthr_'  eigthreshstring(1) 'p' eigthreshstring(3:end)],'-dpdf')
    
    %     %%%% Plot unscaled ctrb values
    %     hcaunsc = figure;
    %     hold on;
    %
    %     for ii = 1:nbcls
    %         for jj = 1:numberoftopmodes%numstate
    %             figure(hcaunsc)
    %             sc = scatter(selected_bcls_for_fps(ii),abs(sortedeig(ii,jj)),ma,eigcf.ctrbmag{ii}(kb,jj));
    %             set(sc,'linewidth',linewidth)
    %             sc.LineWidth = linewidth;
    %         end
    %     end
    %
    %
    %     figure(hcaunsc)
    %     colormap winter
    %     xlabel('BCL, ms')
    %     ylabel('Eigenvalue modulus, |\lambda|')
    %     c = colorbar;
    %     c.Label.String = '|cos \theta_{ki}|';
    %     title(['Controllability magnitude, input index = ' num2str(kb) ', unscaled' ', shift = ' shiftstring])
    %     set(gca,'fontsize',fontsize)
    %     set(gcf,'Position',[550 243 821 615])
    %     saveas(gcf,[ocfolder param '/ctrbeigplot_' param '_ctrlindex_' num2str(kb) 'unscaled'])
    %
    %     % Plot unscaled values (ctrb)
    %     figure
    %     plot(selected_bcls_for_fps,ctrb_avgd_over_modes(kb,:),'linewidth',linewidth);
    %     xlabel('BCL, ms')
    %     ylabel(['|cos(\theta)| averaged over |\lambda| > ' num2str(eigthresh)])
    %     title(['Avg. ctrb. magnitude, control index = ' num2str(kb) ', unscaled' ', shift = ' shiftstring])
    %     set(gca,'fontsize',fontsize)
    %     set(gcf,'Position',[433 243 938 615])
    %     saveas(gcf,[ocfolder param '/ctrbavgplot_' param '_ctrlindex_' num2str(kb) '_unscaled'])
    %
    %     figure
    %     plot(selected_bcls_for_fps,mean(ctrb_avgd_over_modes),'linewidth',linewidth);
    %     xlabel('BCL, ms')
    %     ylabel(['|cos(\theta)| averaged over |\lambda| > ' num2str(eigthresh) ', unscaled'])
    %     title(['Avg. ctrb. magnitude over all control inputs,' ' shift = ' shiftstring])
    %     set(gca,'fontsize',fontsize)
    %     set(gcf,'Position',[433 243 938 615])
    %     saveas(gcf,[ocfolder param '/ctrbavgplot_' param '_unscaled'])
    
    
end
