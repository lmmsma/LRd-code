%This script is to be run after computeeig.m
%It takes in the calculated eigenvalues and plots either their a)Magnitudes
%vs. BCL or b)their Real vs. Imaginary parts
% Script developed by Anthony and Ryan. Modified by Laura to include
% epsilon values used in Jacobian computation. 

clear variables; 

logepsln = -5; % This is the log10 of the epsilon value used in Jacobian computation. 

paramflag = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted) ') == 1;
if paramflag
    param= 'adj';
else
    param = 'def';
end

eigfolder = ['eigenvalues/' param];
eval(['load ' eigfolder '/eigfile' num2str(logepsln) ' *'])
%% Plotting Magnitude of Each Eigenvalue
figure
title(['Eigenvalue moduli for ' param ' parameters, epsilon = 10^{' num2str(logepsln) '}']);
ylabel('|\lambda|');
xlabel('BCL (ms)');
grid on;
hold on
for i=1:length(selected_bcls_for_fps)
    for j=1:17
        scatter(selected_bcls_for_fps(i), alleigsabs{i}(j), 'b*');
    end
end
hold off;

%% Plotting Real vs Imaginary Part for each Eigenvalue
figure
title(['\lambda for ' param ' parameters']);
ylabel('Im(\lambda)');
xlabel('Real(\lambda)');
grid on;
hold on
%for i=1:length(bcls)
    scatter(real(alleigs{33}), imag(alleigs{33}), 'o');
%end