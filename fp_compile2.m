% Use this program to compile all the fixed points into a single matrix
% Modified version of Anthony and Ryan's script to include error size
% storage. 

clear variables; 

%Reads in files of located fixed points
%sourcefolder = 'lrddata_fromfixedpointsearch/';%This folder is the source for fixed points
sourcefolder = 'C:/Users/laura/''Google Drive''/''REU 2017''/Data/''Fixed Points''/''Default Parameters- Ryan''/';%This folder is the source for fixed points
%fixedpointfolder = 'fixedpoints/';%Save compiled fixed points here
fixedpointfolder = 'fixedpoints_fromaandr/';%Save compiled fixed points here
modelinputfolder = 'lrddata/'; % this folder name is currently hardcoded into
% lrd_p2p.m. The M-file will look in this folder for the lrdinputs files
% specified later in this script. 

selected_bcls_for_fps = [1000:-50:500 490:-10:70];%Vector of bcls
% There are more fixed points than the ones listed in selected_bcls...
% Just choose BCLs where Jacobians are desired. They should be available
% at least every 10 ms, starting at 810 ms.  
% The original list for A&R's fixed points is below: 
%bcls = [1000:-50:300 290:-10:70 69:-1:50]; %revised full pacedown (more points leading up to 200ms)

L = length(selected_bcls_for_fps);
modelname = 'lrd'; 
nsolitol = 1e-12; % abs & rel tolerance for Newton-Krylov solver
ncyc = 1; % Run the model for this number of cycles per bcl setting.

%Compile all fixed points into a single matrix, fp_found
fp_found = zeros(17, L);
fp_errnorms = zeros(1,L);

for i = 1:L
    bcl = selected_bcls_for_fps(i);
    % print current BCL to screen
    disp(['BCL = ' num2str(bcl) ' ms'])
    subdiv_per_cyc = 2*bcl; % sampling subdivision; record state vector this many times per cycle
    
    % Pass inputs (bcl, ncyc, and subdiv_per_cyc) to lrd_p2p by file.
    % Pass-by-file is being used since nsoli expects our function to only
    % have one input, which is the initial condition.
    eval(['save ' modelinputfolder 'lrdinputs bcl ncyc subdiv_per_cyc'])
   
   fname = [sourcefolder 'lrddata_1cell_b' num2str(selected_bcls_for_fps(i))];
 %  fp = load(fname, 'Y'); % I can't get this to work properly with folder
 %  names that contain spaces
%   fp_found(:, i) = fp.Y(:, 1);
   eval(['load ' fname ' Y;']);
   fp_found(:, i) = Y(:, 1);
   
   % Compute fixed-point error: 
   eval(['fperr = ' modelname '_p2pzero(fp_found(:,i));']); % compute fixed-point error
   fperrnorm = norm(fperr)
   % from nsoli.m, the stop criterion formula is stop_tol = atol + rtol*fnrm, computed
   % once based on the inital error. Hence, it doesn't really make sense to
   % put in an extra check on the error size, though it can be recorded. 
   fp_errnorms(i) = fperrnorm; 
end

%Save the matrix to a separate file in the folder specified above
eval(['save ' fixedpointfolder 'compiled_fp.mat *']);

% Display minimum and maximum errors:
max(fp_errnorms)
min(fp_errnorms)

%Display full matrix if desired (uncomment)
%disp(fp_found);
