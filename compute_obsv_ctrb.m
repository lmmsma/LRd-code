% Compute and plot observable eigenvalues for LRd Jacobians, for
% different choices of measured variables (i.e., outputs).
% In addition, compute controllable eigenvalues for different choices
% of inputs.
% Code adapted from lrd_batch_eig.m and computeeig.m.
% Laura Munoz, July 2017

clear variables;

ms = 10; fs = 14; % marker size and font size

statenames = char('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');

Jacobians = 'Jacobians/'; % folder where jacobians are stored

%Eigenvalues = 'Eigenvalues/'; %folder where eigenvalues are stored

OCvalues = 'OCvalues/'; %folder where obsv and ctrb values will be saved

eval(['load ' Jacobians 'jacfile *']) %Load Jacobians

numstate = size(alljacs{1},1); % number of state variables

nbcls = length(bcls); % number of bcls

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

rankcutoff = 1e-14; % below this level, singular values don't contribute to the rank

% Intitalize matrices
evcf = cell(nbcls,numstate);
evof = cell(nbcls,numstate);
evncf = cell(nbcls,numstate);
evnof = cell(nbcls,numstate);
rankcf = zeros(nbcls,numstate);
rankof = zeros(nbcls,numstate);
rankcf_scaled = zeros(nbcls,numstate);
rankof_scaled = zeros(nbcls,numstate);
evcf_scaled = cell(nbcls,numstate);
evof_scaled = cell(nbcls,numstate);
evncf_scaled = cell(nbcls,numstate);
evnof_scaled = cell(nbcls,numstate);


for i = 1:nbcls
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
        rankcf(i,kb) = sum(k); % rank of controllable subspace
        % eigenvalues that are controllable from this input
        evcf{i,kb} = eig(Abar((numstate-rankcf(i,kb)+1):end,(numstate-rankcf(i,kb)+1):end),'nobalance'); % eigenvalues of Ac
        % eigenvalues that are not controllable from this input
        evncf{i,kb} = eig(Abar(1:(numstate-rankcf(i,kb)),1:(numstate-rankcf(i,kb))),'nobalance'); % eigenvalues of Anc
        % repeat calculations for scaled system
        [Abars,Bbars,Cbars,Ts,ks] = ctrbf(Smat*jaccd*Smatinv,Bs(:,kb),Cs(1,:),rankcutoff);
        rankcf_scaled(i,kb) = sum(ks);
        evcf_scaled{i,kb} = eig(Abars((numstate-rankcf_scaled(i,kb)+1):end,(numstate-rankcf_scaled(i,kb)+1):end),'nobalance'); % eigenvalues of Ac
        evncf_scaled{i,kb} = eig(Abars(1:(numstate-rankcf_scaled(i,kb)),1:(numstate-rankcf_scaled(i,kb))),'nobalance'); % eigenvalues of Anc
    end
    % Cycle through possible outputs
    for kc = 1:size(C,1)
        % Compute observable and unobservable subspaces, keeping B
        % fixed
        [Abar,Bbar,Cbar,T,k] = obsvf(jaccd,B(:,1),C(kc,:),rankcutoff);
        rankof(i,kc) = sum(k); % rank of observable subspace
        % eigenvalues that are observable from this output
        evof{i,kc} = eig(Abar((numstate-rankof(i,kc)+1):end,(numstate-rankof(i,kc)+1):end),'nobalance'); % eigenvalues of Ao
        % eigenvalues that are not observable from this output
        evnof{i,kc} = eig(Abar(1:(numstate-rankof(i,kc)),1:(numstate-rankof(i,kc))),'nobalance'); % eigenvalues of Ano
        % repeat calculations for scaled system
        [Abars,Bbars,Cbars,Ts,ks] = obsvf(Smat*jaccd*Smatinv,Bs(:,1),Cs(kc,:),rankcutoff);
        rankof_scaled(i,kc) = sum(ks);
        evof_scaled{i,kc} = eig(Abars((numstate-rankof_scaled(i,kc)+1):end,(numstate-rankof_scaled(i,kc)+1):end),'nobalance'); % eigenvalues of Ao
        evnof_scaled{i,kc} = eig(Abars(1:(numstate-rankof_scaled(i,kc)),1:(numstate-rankof_scaled(i,kc))),'nobalance'); % eigenvalues of Ano
    end
    
    eval(['save ' OCvalues 'ocfile *'])
end

% Plot results and save figures

% Cycle through possible outputs (measurements), unscaled system  
for kk = 1:numstate
    figure
    title(['Observable eigenvalues, measurement = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]);
    ylabel('Eigenvalue magnitude');
    xlabel('BCL (ms)');
    grid on;
    hold on;
    for i=1:nbcls
        if ~isempty(evof{i,kk})
            p1 = plot(bcls(i)*ones(size(evof{i,kk})),abs(evof{i,kk}),'b*');
        end
        if ~isempty(evnof{kk})
            p2 = plot(bcls(i)*ones(size(evnof{i,kk})),abs(evnof{i,kk})','r*');
        end
        legend([p1 p2],'observable','unobservable')
        saveas(gcf,[OCvalues 'obsvf_eig_meas' num2str(kk)])
    end
end


% Cycle through possible outputs (measurements), scaled system  
for kk = 1:numstate
    figure
    title(['Scaled system, observable eig., meas. = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]);
    ylabel('Eigenvalue magnitude');
    xlabel('BCL (ms)');
    grid on;
    hold on;
    for i=1:nbcls
        if ~isempty(evof_scaled{i,kk})
            p1 = plot(bcls(i)*ones(size(evof_scaled{i,kk})),abs(evof_scaled{i,kk}),'b*');
        end
        if ~isempty(evnof{kk})
            p2 = plot(bcls(i)*ones(size(evnof_scaled{i,kk})),abs(evnof_scaled{i,kk})','r*');
        end
        legend([p1 p2],'observable','unobservable')
        saveas(gcf,[OCvalues 'obsvf_scaled_eig_meas' num2str(kk)])
    end
end

% Having 17x4 = 68 figures open at once may cause problems, so close some
% of them, if desired
cfflag = input('Before generating controllability figures, enter 1 to close observability figures, or enter 0 to keep the figures: ')
if cfflag
    close all; 
end

% Cycle through possible control inputs, unscaled system  
for kk = 1:numstate
    figure
    title(['Controllable eigenvalues, input = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]);
    ylabel('Eigenvalue magnitude');
    xlabel('BCL (ms)');
    grid on;
    hold on;
    for i=1:nbcls
        if ~isempty(evcf{i,kk})
            p1 = plot(bcls(i)*ones(size(evcf{i,kk})),abs(evcf{i,kk}),'b*');
        end
        if ~isempty(evncf{kk})
            p2 = plot(bcls(i)*ones(size(evncf{i,kk})),abs(evncf{i,kk})','r*');
        end
        legend([p1 p2],'controllable','uncontrollable')
        saveas(gcf,[OCvalues 'ctrbf_eig_meas' num2str(kk)])
    end
end

% Cycle through possible control inputs, scaled system
for kk = 1:numstate
    figure
    title(['Scaled system, controllable eig., input = ' strtrim(statenames(kk,:)), ', rankcutoff = ' num2str(rankcutoff)]);
    ylabel('Eigenvalue magnitude');
    xlabel('BCL (ms)');
    grid on;
    hold on;
    for i=1:nbcls
        if ~isempty(evcf_scaled{i,kk})
            p1 = plot(bcls(i)*ones(size(evcf_scaled{i,kk})),abs(evcf_scaled{i,kk}),'b*');
        end
        if ~isempty(evncf_scaled{kk})
            p2 = plot(bcls(i)*ones(size(evncf_scaled{i,kk})),abs(evncf_scaled{i,kk})','r*');
        end
        legend([p1 p2],'controllable','uncontrollable')
        saveas(gcf,[OCvalues 'ctrbf_scaled_eig_meas' num2str(kk)])
    end
end
    