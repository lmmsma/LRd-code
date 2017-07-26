% Compute and plot observable eigenvalues for LRd Jacobians, for
% different choices of measured variables (i.e., outputs).
% In addition, compute controllable eigenvalues for different choices
% of inputs.
% Code adapted from lrd_batch_eig.m and computeeig.m.
% Laura Munoz, July 2017

clear variables;

adj_yn = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted) ') == 1;
if adj_yn
    param= 'adj';
else
    param = 'def';
end
ms = 10; fs = 14; % marker size and font size

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');

Jacobians = 'Jacobians/'; % folder where jacobians are stored

Eigenvalues = 'Eigenvalues/'; %folder where eigenvalues are stored

OCvalues = 'OCvalues/'; %folder where obsv and ctrb values will be saved

eval(['load ' Jacobians 'jacfile_' param ' data allfp modelname bcls alljacs']); %Load Jacobians
eval(['load ' Eigenvalues 'eigfile_' param ' allv alleigs']);

numstate = size(alljacs{1},1); % number of state variables

nbcls = length(bcls); % number of bcls
%--------------------------------------------------------------------------%
thresh = .0001; %Threshold for determining whether eigenvalues are identical
rankcutoff = 1e-14; % below this level, singular values don't contribute to the rank
%--------------------------------------------------------------------------%

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
Cs = C*Smatinv; % scaled C matrix (could use Smat*C*Smatinv instead, if the output should also be normalized?

% Intitalize matrices
eigcf.cval = cell(nbcls,numstate);
eigof.oval = cell(nbcls,numstate);
eigcf.ncval = cell(nbcls,numstate);
eigof.noval = cell(nbcls,numstate);
eigcf.rank = zeros(nbcls,numstate);
eigof.rank = zeros(nbcls,numstate);
eigcf.rank_sc = zeros(nbcls,numstate);
eigof.rank_sc = zeros(nbcls,numstate);
eigcf.cval_sc = cell(nbcls,numstate);
eigof.oval_sc = cell(nbcls,numstate);
eigcf.ncval_sc = cell(nbcls,numstate);
eigof.noval_sc = cell(nbcls,numstate);
eigcf.cvec = cell(nbcls, numstate);
eigcf.ncvec = cell(nbcls, numstate);
eigof.ovec = cell(nbcls, numstate);
eigof.novec = cell(nbcls, numstate);
eigcf.cvec_sc = cell(nbcls, numstate);
eigcf.ncvec_sc = cell(nbcls, numstate);
eigof.ovec_sc = cell(nbcls, numstate);
eigof.novec_sc = cell(nbcls, numstate);

for i = 1:nbcls
    %% Controllability
    bcl = bcls(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    
    jaccd = alljacs{i};
    
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
        [Abars,Bbars,Cbars,Ts,ks] = ctrbf(Smat*jaccd*Smatinv,Bs(:,kb),Cs(1,:),rankcutoff);
        eigcf.rank_sc(i,kb) = sum(ks);
        eigcf.cval_sc{i,kb} = eig(Abars((numstate-eigcf.rank_sc(i,kb)+1):end,(numstate-eigcf.rank_sc(i,kb)+1):end),'nobalance'); % eigenvalues of Ac
        eigcf.ncval_sc{i,kb} = eig(Abars(1:(numstate-eigcf.rank_sc(i,kb)),1:(numstate-eigcf.rank_sc(i,kb))),'nobalance'); % eigenvalues of Anc
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
        eigof.ovec{i ,kc}= zeros(17, eigof.rank(i, kc));
        for j=1:eigof.rank(i, kc)
            f= allv{1, i}(:,abs(real(alleigs{1, i})-eigof.oval{i, kc}(j)) < abs(real(eigof.oval{i, kc}(j)*thresh)));
            if ~isempty(f)
                eigof.ovec{i,kc}(:, j)= f;
            end
        end
        % eigenvalues that are not observable from this output
        eigof.noval{i,kc} = eig(Abar(1:(numstate-eigof.rank(i,kc)),1:(numstate-eigof.rank(i,kc))),'nobalance'); % eigenvalues of Ano
        eigof.novec{i ,kc}= zeros(17, numstate-eigof.rank(i, kc));
        for j=1:(numstate- eigof.rank(i, kc))
            f= allv{1, i}(:,abs(real(alleigs{1, i})-eigof.noval{i, kc}(j)) < abs(real(eigof.noval{i, kc}(j)*thresh)));
            if ~isempty(f)
                eigof.novec{i,kc}(:, j)=  f;
            end
        end
        % repeat calculations for scaled system
        [Abars,Bbars,Cbars,Ts,ks] = obsvf(Smat*jaccd*Smatinv,Bs(:,1),Cs(kc,:),rankcutoff);
        eigof.rank_sc(i,kc) = sum(ks);
        eigof.oval_sc{i,kc} = eig(Abars((numstate-eigof.rank_sc(i,kc)+1):end,(numstate-eigof.rank_sc(i,kc)+1):end),'nobalance'); % eigenvalues of Ao
        eigof.ovec_sc{i ,kc}= zeros(17, eigof.rank_sc(i, kc));
        for j=1:eigof.rank_sc(i, kc)
            f= allv{1, i}(:, abs(real(alleigs{1, i})-eigof.oval_sc{i, kc}(j)) < abs(real(eigof.oval_sc{i, kc}(j)*thresh)));
            if ~isempty(f)
                eigof.ovec_sc{i,kc}(:, j)= f;
            end
        end
        eigof.noval_sc{i,kc} = eig(Abars(1:(numstate-eigof.rank_sc(i,kc)),1:(numstate-eigof.rank_sc(i,kc))),'nobalance'); % eigenvalues of Ano
        eigof.novec_sc{i ,kc}= zeros(17, numstate-eigof.rank_sc(i, kc));
        for j=1:(numstate-eigof.rank_sc(i, kc))
            f= allv{1, i}(:, abs(real(alleigs{1, i})-eigof.noval_sc{i, kc}(j)) < abs(real(eigof.noval_sc{i, kc}(j)*thresh)));
            if ~isempty(f)
                eigof.novec_sc{i,kc}(:, j)= f;
            end
        end
    end
end
clear adj_yn bcl B C Cs Eigenvalues f fs i j Jacobians k kb kc ks ms Smat Smatinv umax varamp;
eval(['save ' OCvalues param '/ocfile *'])
%Use plot_obs_cont to plot these results
