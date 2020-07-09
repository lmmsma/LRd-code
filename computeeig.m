% Script developed by Anthony and Ryan. Modified by Laura to include
% epsilon values used in Jacobian computation, use 'nobalance' option, 
% and to compute and store left eigenvectors. 

clear variables;

logepsln = -5; % This is the log10 of the epsilon value used in Jacobian computation. 

%shiftstring = ''; % if using default shift of data.dt
%shiftstring = '_shift0p75ms';
%shiftstring = '_shift6p5ms';
%shiftstring = '_shift7mVrepol';
%shiftstring = '_shift-50mVrepol';
%shiftstring = '_shift1Vnormdepol';
%shiftstring = '_shift0p2Vnormdepol';
%shiftstring = '_shift0p4Vnormdepol';
%shiftstring = '_shift0p6Vnormdepol';
%shiftstring = '_shift0p8Vnormdepol';
%shiftstring = '_shift0p8Vnormrepol';
%shiftstring = '_shift0p6Vnormrepol';
%shiftstring = '_shift0p4Vnormrepol';
%shiftstring = '_shift0p2Vnormrepol';
shiftstring = '_shift0p001Vnormrepol';

paramflag = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if paramflag
    param = 'adj';
else
    param = 'def';
end

%jacfolder = 'jacobians/def/'; % folder where jacobians are stored
jacfolder = ['jacobians' shiftstring '/' param '/']; % folder where jacobians are stored

%eigfolder = 'eigenvalues/'; %folder where eigenvalues will be saved. 
eigfolder = ['eigenvalues' shiftstring  '/' param '/']; %folder where eigenvalues will be saved. 
if ~exist(eigfolder,'dir')
    mkdir(eigfolder)
end
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