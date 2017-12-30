% computejacobians_batch.m: Compute Jacobians in batch mode; loop over
% different epsilon (relative perturbation) values. 
% Laura Munoz, August 2017

clear variables;
 
%epsilons = [10^-6 10^-5 10^-4]; % These values seemed to be the 
% best ones based on preliminary tests of 2-norms of Jacobian and eigenvalue
% differences, for A&R's vlaues. 
epsilons = [10^-7 10^-3]; %  

for j=1:length(epsilons)
    computejacobians(epsilons(j)); 
end

