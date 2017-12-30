% Script developed by Anthony and Ryan. Modified by Laura to include
% epsilon values used in Jacobian computation, use 'nobalance' option, 
% and to compute and store left eigenvectors. 

clear variables;

logepsln = -3; % This is the log10 of the epsilon value used in Jacobian computation. 

jacfolder = 'jacobians/'; % folder where jacobians are stored

eigfolder = 'eigenvalues/'; %folder where eigenvalues will be saved. 

eval(['load ' jacfolder 'jacfile' num2str(logepsln) ' *']) %Load data from jacobians

alleigs = cell(1,length(selected_bcls_for_fps)); % Store eigenvalues here.
alleigsabs = cell(1,length(selected_bcls_for_fps)); % Store eigenvalue magnitudes here.
allv = cell(1,length(selected_bcls_for_fps)); % Store right eigenvectors here.
allw = cell(1,length(selected_bcls_for_fps)); % Store left eigenvectors here.

for i = 1:length(selected_bcls_for_fps)
    bcl = selected_bcls_for_fps(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    
    %eval(['save ' eigfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
    
    [allv{i}, eigv, allw{i}] = eig(alljacs{i},'nobalance');
    alleigs{i}= diag(eigv);
    alleigsabs{i}= abs(alleigs{i});
    eval(['save ' eigfolder 'eigfile' num2str(logepsln) ' *'])
end

% ploteigs