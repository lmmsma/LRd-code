clear variables;

Jacobians = 'Jacobians/'; % folder where jacobians are stored

Eigenvalues = 'Eigenvalues/'; %folder where eiganvalues will be saved. 

eval(['load ' Jacobians 'jacfile *']) %Load data from jacobians

alleigs = cell(1,length(bcls)); % Store eigenvalues here.
alleigsabs = cell(1,length(bcls)); % Store eigenvalue magnitude here
allv= cell(1,length(bcls)); % Store eigenvalues here.
for i = 1:length(bcls)
    bcl = bcls(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    
    eval(['save ' Eigenvalues 'lrdinputs bcl ncyc subdiv_per_cyc'])
    
    [allv{i}, eigv] = eig(alljacs{i});
    alleigs{i}= diag(eigv);
    alleigsabs{i}= abs(alleigs{i});
    eval(['save ' Eigenvalues 'eigfile *'])
end

ploteigs