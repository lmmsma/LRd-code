clear variables;

Jacobians = ['Jacobians/']; % folder where jacobians are stored

Eiganvalues = ['Eiganvalues/']; %folder where eiganvalues will be saved. 

eval(['load ' Jacobians 'jacfile *']) %Load data from jacobians

alleigs = cell(1,length(bcls)); % Store eiganvalues here.

for i = 1:length(bcls)
    bcl = bcls(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    
    eval(['save ' Eiganvalues 'lrdinputs bcl ncyc subdiv_per_cyc'])
    
    alleigs{i} = eig(alljacs{i});
    eval(['save ' Eiganvalues 'eigfile *'])
end
