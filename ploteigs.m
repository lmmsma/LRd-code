%This script is to be run after computeeig.m
%It takes in the calculated eigenvalues and plots either their a)Magnitudes
%vs. BCL or b)their Real vs. Imaginary parts

eigfolder = 'Eigenvalues';
eval(['load ' eigfolder '/eigfile bcls'])
%% Plotting Magnitude of Each Eigenvalue
%figure
title('\lambda Magnitudes for default parameters');
ylabel('Eigenvalue magnitude');
xlabel('BCL (ms)');
grid on;
hold on
for i=1:length(bcls)
    for j=1:17
        scatter(bcls(i), alleigsabs{i}(j), '*');
    end
end
hold off;

%% Plotting Real vs Imaginary Part for each Eigenvalue
figure
title('\lambda for default parameters');
ylabel('Im(\lambda)');
xlabel('Real(\lambda)');
grid on;
hold on
%for i=1:length(bcls)
    scatter(real(alleigs{33}), imag(alleigs{33}), 'o');
%end