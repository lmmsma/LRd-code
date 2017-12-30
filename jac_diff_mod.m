% Modified version of jac_diff, created to use a more conventional norm for
% comparisons. Find epsilon (epsln) value that minimizes norm of difference
% in Jacobians computed for successive epsilons. 
%
% Suspending development on this one for now, since the file contents 
% (for different epsln values)
% are too inconsistent to use a looping structure to read in data. 
% There is also a key error in one file, where the recorded epsln does not 
% match the value (10^-5) in the name of the folder. 
% Laura Munoz, Aug. 2017

clear variables; 

% Can't get exist to return the right result with full path name, but
% somehow it works when I cd into the folder first. 

% Home folder: 
homefolder = 'C:\Users\laura\''Google Drive''\LRd-code-master_mod\'; 

% Folder where Jacobians are stored: 
sourcefolder = 'C:\Users\laura\''Google Drive''\''REU 2017''\Data\''Jacobian Data''\''Default Parameters''\'; 

% List logs of epsilon that were tested, e.g. -5 stands for
% epsln = 10^-5 : 

exps = -5; 
%exps = -8:-3; 

% How are Jacobians stored within each folder?     
% Looks like files with "more" data could either be stored in
    % jacfile_70.mat, or just jacfile.mat, if the former doesn't exist? 

alljac1000 = cell(length(exps),1); 
alljac200 = cell(length(exps),1); 
alljac70 = cell(length(exps),1); 

for ii = 1:length(exps)
    expfolder = [sourcefolder 'eps=10_' num2str(exps(ii)) '\']; 
    eval(['cd ' expfolder]); 

    if exist('jacfile_70.mat', 'file')
        load jacfile_70.mat; 
        disp('Found jacfile_70.mat') 
    else
        load jacfile.mat; 
%        eval(['load ' expfolder 'jacfile.mat']); 
    end
    eval(['cd ' homefolder]); 
    
    disp(['epsln = ' num2str(epsln)])
    if log10(epsln) ~= exps(ii)
        disp(['Error: epsilon mismatch, selected index = ' num2str(exps(ii))])
%        return;
    end
   % disp(['BCL = ' num2str(bcls(bcls==1000))])
    alljac1000{ii} = alljacs{bcls==1000}; 
    
    
    disp(['BCL = ' num2str(bcls(bcls==200))])
    %alljac200{ii} = alljacs{bcls==200}; 
    

    disp(['BCL = ' num2str(bcls(bcls==70))])
    %alljac70{ii} = alljacs{bcls==70}; 
  
end

