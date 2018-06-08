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

numberoftopmodes = 5; % plot n largest eigenvalues, where n = numberoftopmodes
eigthresh = 0.9; % eigenvalues with magnitudes greater than or equal to this threshold will 
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

markersize = 80; fontsize = 20; linewidth = 2; % marker size, font size, line width

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
statenames_latex = char('$V$','$h$','$m$','$j$','$d$','$f$','$x_r$','$[Ca^{2+}]_{i,t}$','$[Na^+]_i$','$[K^+]_i$','$[Ca^{2+}]_{j,t}$','$[Ca^{2+}]_n$','$x_{s1}$','$b$','$g$','$x_{s2}$','$I_{rel}$');

jacfolder = ['jacobians/' param '/']; % folder where jacobians are stored

%eigfolder = ['eigenvalues/' param '/']; %folder where eigenvalues are stored

ocfolder = 'ocvalues/'; %folder where obsv and ctrb values will be saved

eval(['load ' jacfolder 'jacfile' num2str(selected_logepsln) ' *']); %Load Jacobians
%eval(['load ' eigenvalues 'eigfile_' param ' allv alleigs']);

numstate = size(alljacs{1},1); % number of state variables

% number of bcls
nbcls = length(selected_bcls_for_fps);
%--------------------------------------------------------------------------%
%thresh = .0001; %Threshold for determining whether eigenvalues are identical
rankcutoff = 1e-14; % below this level, singular values don't contribute to the rank
%--------------------------------------------------------------------------%
% input matrix for all possible individual inputs
B = eye(numstate);

% output matrix for all possible individual measurements
C = eye(numstate);

% load approximate "state normalization" scaling matrix
load b1000fsolem12variable_amplitudes varamp
% Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
Smat = diag(1./varamp); % scaling matrix
Smatinv = inv(Smat); % only need to compute once
umax = diag(ones(1,size(B,2))); % This is just a placeholder.
% For inputs, also need to know approximate maximum values of each element
% of input vector. Note that u is a deviational quantity that may or
% may not attain the size of the stimulus.
% Each umax diagonal element should probably be the size of
% deviational max for integrated system.
Bs = Smat*B*umax; % scaled B matrix
Cs = C*Smatinv; % scaled C matrix

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

for i = 1:nbcls
    %% Controllability
    bcl = selected_bcls_for_fps(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    
    % For debugging:
    if debugflag
        jaccd = eye(17);
        jaccd(17,:) = ones(1,17);
    else
        jaccd = alljacs{i};
    end
    jaccd_sc = Smat*jaccd*Smatinv; % scaled Jacobian
    
    % The eigenvalues could be loaded instead of recomputed here,
    % but recomputing doesn't take much time, whereas reloading
    % increases risk of a mismatch between eigenvalues and Jacobians.
    [v,d,w]=eig(jaccd,'nobalance');
    % Sort eigenvalues and eigenvectors in descending order of size
    [sortedeig(i,:),eigsortind] = sort(abs(diag(d)),'descend');
    sortedv{i} = v(:,eigsortind);
    sortedw{i} = w(:,eigsortind);
    
    % Compute eigenvectors for scaled Jacobian:
    [vsc,dsc,wsc]=eig(jaccd_sc,'nobalance');
    [sortedeigsc(i,:),eigsortindsc] = sort(abs(diag(dsc)),'descend');
    sortedvsc{i} = vsc(:,eigsortindsc);
    sortedwsc{i} = wsc(:,eigsortindsc);
    
    %    rankcutoffdefault = max(size(jaccd))*eps(norm(jaccd));% This is Matlab's default rank cutoff. Can't easily control cutoffs in orth /null, so to interpret these results, would need to know default
    
    % Cycle through possible inputs
    for kb = 1:size(B,2)
        % Compute controllable and uncontrollable subspaces, keeping C
        % fixed
        [Abar,Bbar,Cbar,T,k] = ctrbf(jaccd,B(:,kb),C(1,:),rankcutoff);
        eigcf.rank(i,kb) = sum(k); % rank of controllable subspace
        % eigenvalues that are controllable from this input
        eigcf.cval{i,kb} = eig(Abar((numstate-eigcf.rank(i,kb)+1):end,(numstate-eigcf.rank(i,kb)+1):end),'nobalance'); % eigenvalues of Ac
        % eigenvalues that are not controllable from this input
        eigcf.ncval{i,kb} = eig(Abar(1:(numstate-eigcf.rank(i,kb)),1:(numstate-eigcf.rank(i,kb))),'nobalance'); % eigenvalues of Anc
        % repeat calculations for scaled system
        [Abars,Bbars,Cbars,Ts,ks] = ctrbf(jaccd_sc,Bs(:,kb),Cs(1,:),rankcutoff);
        eigcf.rank_sc(i,kb) = sum(ks);
        eigcf.cval_sc{i,kb} = eig(Abars((numstate-eigcf.rank_sc(i,kb)+1):end,(numstate-eigcf.rank_sc(i,kb)+1):end),'nobalance'); % eigenvalues of Ac
        eigcf.ncval_sc{i,kb} = eig(Abars(1:(numstate-eigcf.rank_sc(i,kb)),1:(numstate-eigcf.rank_sc(i,kb))),'nobalance'); % eigenvalues of Anc
        
        % Compute modal controllability values.
        % Warning: computational method should be altered for repeated
        % eigenvalues.
        for kk = 1:numstate
            eigcf.bdotw{i}(kb,kk) = sortedw{i}(:,kk)'*B(:,kb);
            eigcf.bmagwmag{i}(kb,kk) = norm(B(:,kb))*norm(sortedw{i}(:,kk));
            if eigcf.bmagwmag{i}(kb,kk) > 0
                eigcf.ctrbmag{i}(kb,kk) = abs(eigcf.bdotw{i}(kb,kk))/eigcf.bmagwmag{i}(kb,kk);
            else
                eigcf.ctrbmag{i}(kb,kk) = NaN;
            end
            
            eigcf.bdotw_sc{i}(kb,kk) = sortedwsc{i}(:,kk)'*Bs(:,kb);
            eigcf.bmagwmag_sc{i}(kb,kk) = norm(Bs(:,kb))*norm(sortedwsc{i}(:,kk));
            if eigcf.bmagwmag_sc{i}(kb,kk) > 0
                eigcf.ctrbmag_sc{i}(kb,kk) = abs(eigcf.bdotw_sc{i}(kb,kk))/eigcf.bmagwmag_sc{i}(kb,kk);
            else
                eigcf.ctrbmag_sc{i}(kb,kk) = NaN;
            end
        end
    end
    
    
    %% Observability
    % Cycle through possible outputs
    for kc = 1:size(C,1)
        % Compute observable and unobservable subspaces, keeping B
        % fixed
        [Abar,Bbar,Cbar,T,k] = obsvf(jaccd,B(:,1),C(kc,:),rankcutoff);
        eigof.rank(i,kc) = sum(k); % rank of observable subspace
        % eigenvalues that are observable from this output
        eigof.oval{i,kc} = eig(Abar((numstate-eigof.rank(i,kc)+1):end,(numstate-eigof.rank(i,kc)+1):end),'nobalance'); % eigenvalues of Ao
        %         eigof.ovec{i ,kc}= zeros(17, eigof.rank(i, kc));
        %         for j=1:eigof.rank(i, kc)
        %             f= allv{1, i}(:,abs(real(alleigs{1, i})-eigof.oval{i, kc}(j)) < abs(real(eigof.oval{i, kc}(j)*thresh)));
        %             if ~isempty(f)
        %                 eigof.ovec{i,kc}(:, j)= f;
        %             end
        %         end
        % eigenvalues that are not observable from this output
        eigof.noval{i,kc} = eig(Abar(1:(numstate-eigof.rank(i,kc)),1:(numstate-eigof.rank(i,kc))),'nobalance'); % eigenvalues of Ano
        %         eigof.novec{i ,kc}= zeros(17, numstate-eigof.rank(i, kc));
        %         for j=1:(numstate- eigof.rank(i, kc))
        %             f= allv{1, i}(:,abs(real(alleigs{1, i})-eigof.noval{i, kc}(j)) < abs(real(eigof.noval{i, kc}(j)*thresh)));
        %             if ~isempty(f)
        %                 eigof.novec{i,kc}(:, j)=  f;
        %             end
        %         end
        % repeat calculations for scaled system
        [Abars,Bbars,Cbars,Ts,ks] = obsvf(jaccd_sc,Bs(:,1),Cs(kc,:),rankcutoff);
        eigof.rank_sc(i,kc) = sum(ks);
        eigof.oval_sc{i,kc} = eig(Abars((numstate-eigof.rank_sc(i,kc)+1):end,(numstate-eigof.rank_sc(i,kc)+1):end),'nobalance'); % eigenvalues of Ao
        %        eigof.ovec_sc{i ,kc}= zeros(17, eigof.rank_sc(i, kc));
        %         for j=1:eigof.rank_sc(i, kc)
        %             f= allv{1, i}(:, abs(real(alleigs{1, i})-eigof.oval_sc{i, kc}(j)) < abs(real(eigof.oval_sc{i, kc}(j)*thresh)));
        %             if ~isempty(f)
        %                 eigof.ovec_sc{i,kc}(:, j)= f;
        %             end
        %         end
        eigof.noval_sc{i,kc} = eig(Abars(1:(numstate-eigof.rank_sc(i,kc)),1:(numstate-eigof.rank_sc(i,kc))),'nobalance'); % eigenvalues of Ano
        %         eigof.novec_sc{i ,kc}= zeros(17, numstate-eigof.rank_sc(i, kc));
        %         for j=1:(numstate-eigof.rank_sc(i, kc))
        %             f= allv{1, i}(:, abs(real(alleigs{1, i})-eigof.noval_sc{i, kc}(j)) < abs(real(eigof.noval_sc{i, kc}(j)*thresh)));
        %             if ~isempty(f)
        %                 eigof.novec_sc{i,kc}(:, j)= f;
        %             end
        %         end
        
        % Compute modal observability values.
        % Warning: computational method should be altered for repeated
        % eigenvalues.
        for kk = 1:numstate
            eigof.cdotv{i}(kc,kk) = C(kc,:)*sortedv{i}(:,kk);
            eigof.cmagvmag{i}(kc,kk) = norm(C(kc,:))*norm(sortedv{i}(:,kk));
            if eigof.cmagvmag{i}(kc,kk) > 0
                eigof.obsvmag{i}(kc,kk) = abs(eigof.cdotv{i}(kc,kk))/eigof.cmagvmag{i}(kc,kk);
            else
                eigof.obsvmag{i}(kc,kk) = NaN;
            end
            
            eigof.cdotv_sc{i}(kc,kk) = Cs(kc,:)*sortedvsc{i}(:,kk);
            eigof.cmagvmag_sc{i}(kc,kk) = norm(Cs(kc,:))*norm(sortedvsc{i}(:,kk));
            if eigof.cmagvmag_sc{i}(kc,kk) > 0
                eigof.obsvmag_sc{i}(kc,kk) = abs(eigof.cdotv_sc{i}(kc,kk))/eigof.cmagvmag_sc{i}(kc,kk);
            else
                eigof.obsvmag_sc{i}(kc,kk) = NaN;
            end
        end
    end
end

% Compute averages over BCLs and modes
largesteigindices = cell(nbcls,1); % store indices of largest eigenvalues here
numlargesteigs = zeros(nbcls,1); % store number of largest eigenvalues here
largesteigindices_sc = cell(nbcls,1); % store indices of largest eigenvalues here
numlargesteigs_sc = zeros(nbcls,1); % store number of largest eigenvalues here

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

for ii = 1:nbcls
    tempoa = tempoa + eigof.obsvmag{ii};
    tempca = tempca + eigcf.ctrbmag{ii};
    tempoasc = tempoasc + eigof.obsvmag_sc{ii};
    tempcasc = tempcasc + eigcf.ctrbmag_sc{ii};
    
    largesteigindices{ii} = find(abs(sortedeig(ii,:)) >= eigthresh); 
    largesteigindices_sc{ii} = find(abs(sortedeigsc(ii,:)) >= eigthresh);
    numlargesteigs(ii) = length(largesteigindices{ii}); 
    numlargesteigs_sc(ii) = length(largesteigindices_sc{ii}); 
    
    obsv_avgd_over_modes(:,ii) = mean(eigof.obsvmag{ii}(:,largesteigindices{ii}),2);
    ctrb_avgd_over_modes(:,ii) = mean(eigcf.ctrbmag{ii}(:,largesteigindices{ii}),2);
    obsv_avgd_over_modes_sc(:,ii) = mean(eigof.obsvmag_sc{ii}(:,largesteigindices_sc{ii}),2);
    ctrb_avgd_over_modes_sc(:,ii) = mean(eigcf.ctrbmag_sc{ii}(:,largesteigindices_sc{ii}),2);

    tempoam = tempoam + obsv_avgd_over_modes(:,ii);
    tempcam = tempcam + ctrb_avgd_over_modes(:,ii);
    tempoamsc = tempoamsc + obsv_avgd_over_modes_sc(:,ii);
    tempcamsc = tempcamsc + ctrb_avgd_over_modes_sc(:,ii);
end

obsv_avgd_over_bcls = tempoa/nbcls;
ctrb_avgd_over_bcls = tempca/nbcls;
obsv_avgd_over_bcls_sc = tempoasc/nbcls;
ctrb_avgd_over_bcls_sc = tempcasc/nbcls;

obsv_avgd_over_bcls_and_modes = tempoam/nbcls; 
ctrb_avgd_over_bcls_and_modes = tempcam/nbcls; 
obsv_avgd_over_bcls_and_modes_sc = tempoamsc/nbcls; 
ctrb_avgd_over_bcls_and_modes_sc = tempcamsc/nbcls; 


% obsv_avgd_over_bcls_and_modes = mean(obsv_avgd_over_bcls(:,1:numberoftopmodes),2);
% ctrb_avgd_over_bcls_and_modes = mean(ctrb_avgd_over_bcls(:,1:numberoftopmodes),2);
% obsv_avgd_over_bcls_and_modes_sc = mean(obsv_avgd_over_bcls_sc(:,1:numberoftopmodes),2);
% ctrb_avgd_over_bcls_and_modes_sc = mean(ctrb_avgd_over_bcls_sc(:,1:numberoftopmodes),2);

disp('Unscaled observability for top mode, averaged over bcls:')
[sortedobsv,obsvsortindex] = sort(abs(obsv_avgd_over_bcls(:,1)),'descend');
statenames_latex(obsvsortindex,:)
obsv_avgd_over_bcls(obsvsortindex)

disp('Unscaled observability, averaged over bcls and modes above threshold:')
[sortedobsvam,obsvsortindexam] = sort(abs(obsv_avgd_over_bcls_and_modes),'descend');
statenames_latex(obsvsortindexam,:)
obsv_avgd_over_bcls_and_modes(obsvsortindexam)

disp('Scaled observability for top mode, averaged over bcls and modes above threshold:')
[sortedobsv_sc,obsvsortindex_sc] = sort(abs(obsv_avgd_over_bcls_sc(:,1)),'descend');
statenames_latex(obsvsortindex_sc,:)
obsv_avgd_over_bcls_sc(obsvsortindex_sc)

disp('Scaled observability, averaged over bcls and modes above threshold:')
[sortedobsvam_sc,obsvsortindexam_sc] = sort(abs(obsv_avgd_over_bcls_and_modes_sc),'descend');
statenames_latex(obsvsortindexam_sc,:)
obsv_avgd_over_bcls_and_modes_sc(obsvsortindexam_sc)

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

% Controllability
disp('Unscaled controllability for top mode, averaged over bcls:')
[sortedctrb,ctrbsortindex] = sort(abs(ctrb_avgd_over_bcls(:,1)),'descend');
statenames_latex(ctrbsortindex,:)
ctrb_avgd_over_bcls(ctrbsortindex)

disp('Unscaled controllability, averaged over bcls and modes above threshold:')
[sortedctrbam,ctrbsortindexam] = sort(abs(ctrb_avgd_over_bcls_and_modes),'descend');
statenames_latex(ctrbsortindexam,:)
ctrb_avgd_over_bcls_and_modes(ctrbsortindexam)

disp('Scaled controllability for top mode, averaged over bcls and modes above threshold:')
[sortedctrb_sc,ctrbsortindex_sc] = sort(abs(ctrb_avgd_over_bcls_sc(:,1)),'descend');
statenames_latex(ctrbsortindex_sc,:)
ctrb_avgd_over_bcls_sc(ctrbsortindex_sc)

disp('Scaled controllability, averaged over bcls and modes above threshold:')
[sortedctrbam_sc,ctrbsortindexam_sc] = sort(abs(ctrb_avgd_over_bcls_and_modes_sc),'descend');
statenames_latex(ctrbsortindexam_sc,:)
ctrb_avgd_over_bcls_and_modes_sc(ctrbsortindexam_sc)

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



%clear paramflag bcl B C Cs eigfolder f fs i j jacfolder k kb kc ks ms Smat Smatinv umax varamp;
if ~debugflag
    eval(['save ' ocfolder param '/ocfile *'])
end
%Use plot_obs_cont to plot these results


% plot symbols
symbols = char('bo-','rs-','gp-','m*-','k^-','cv-','yh-');

hor = figure;
hold on;

hoi = figure;
hold on;

hoa = figure;
hold on;

hcr = figure;
hold on;

hci = figure;
hold on;

hca = figure;
hold on;

%kc = 10
%kb = 1;
kc = input('Enter a measurement index from 1 to 17: ');
kb = input('Enter a control input index from 1 to 17: ');

for ii = 1:nbcls
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
    for jj = 1:numberoftopmodes%numstate
        figure(hoa)
        so=scatter(selected_bcls_for_fps(ii),abs(sortedeigsc(ii,jj)),markersize,eigof.obsvmag_sc{ii}(kc,jj));
        so.LineWidth = linewidth; 
        figure(hca)
        sc = scatter(selected_bcls_for_fps(ii),abs(sortedeigsc(ii,jj)),markersize,eigcf.ctrbmag_sc{ii}(kb,jj));
        set(sc,'linewidth',linewidth)
        sc.LineWidth = linewidth; 
    end
end

figure(hor)
xlabel('BCL, ms')
ylabel('Re C*v')

figure(hoi)
xlabel('BCL, ms')
ylabel('Im C*v')

figure(hoa)
colormap winter
xlabel('BCL, ms')
ylabel('Eigenvalue modulus, |\lambda|')
c = colorbar;
%c.Label.String = '|Cv|';
c.Label.Interpreter = 'latex'; 
c.Label.String = '$|\cos \phi_{ki}|$';
title(['Observability magnitude, meas. index = ' num2str(kc)])
set(gca,'fontsize',fontsize)
set(gcf,'Position',[550 243 821 615])
saveas(gcf,[ocfolder param '/obsveigplot_' param '_measindex_' num2str(kc)])

figure(hcr)
xlabel('BCL, ms')
ylabel('Re w^T*B')

figure(hci)
xlabel('BCL, ms')
ylabel('Im w^T*B')

figure(hca)
colormap winter
xlabel('BCL, ms')
ylabel('Eigenvalue modulus, |\lambda|')
c = colorbar;
c.Label.String = '|w^TB|';
title(['Controllability magnitude, input index = ' num2str(kb)])
set(gca,'fontsize',fontsize)
set(gcf,'Position',[550 243 821 615])
saveas(gcf,[ocfolder param '/ctrbeigplot_' param '_ctrlindex_' num2str(kb)])

% Plot scaled values
figure
plot(selected_bcls_for_fps,obsv_avgd_over_modes_sc(kc,:),'linewidth',linewidth);
xlabel('BCL, ms')
ylabel('|cos(\phi)| averaged over |\lambda| > 0.9')
title(['Avg. obsv. magnitude, meas. index = ' num2str(kc)])
set(gca,'fontsize',fontsize)
set(gcf,'Position',[433 243 938 615])
saveas(gcf,[ocfolder param '/obsvavgplot_' param '_measindex_' num2str(kc)])

figure
plot(selected_bcls_for_fps,mean(obsv_avgd_over_modes_sc),'linewidth',linewidth);
xlabel('BCL, ms')
ylabel('|cos(\phi)| averaged over |\lambda| > 0.9')
title(['Avg. obsv. magnitude over all measurements'])
set(gca,'fontsize',fontsize)
set(gcf,'Position',[433 243 938 615])
saveas(gcf,[ocfolder param '/obsvavgplot_' param])

% Plot scaled values (ctrb) 
figure
plot(selected_bcls_for_fps,ctrb_avgd_over_modes_sc(kb,:),'linewidth',linewidth);
xlabel('BCL, ms')
ylabel('|cos(\theta)| averaged over |\lambda| > 0.9')
title(['Avg. ctrb. magnitude, control index = ' num2str(kb)])
set(gca,'fontsize',fontsize)
set(gcf,'Position',[433 243 938 615])
saveas(gcf,[ocfolder param '/ctrbavgplot_' param '_ctrlindex_' num2str(kb)])

figure
plot(selected_bcls_for_fps,mean(ctrb_avgd_over_modes_sc),'linewidth',linewidth);
xlabel('BCL, ms')
ylabel('|cos(\theta)| averaged over |\lambda| > 0.9')
title(['Avg. ctrb. magnitude over all control inputs'])
set(gca,'fontsize',fontsize)
set(gcf,'Position',[433 243 938 615])
saveas(gcf,[ocfolder param '/ctrbavgplot_' param])

