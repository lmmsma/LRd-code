function jacobian = jacobian_cd(fixedpt,epsln,modelname)
% Compute an approximate Jacobian for a given model about a provided 
% fixed point, using modified versions of C.T. Kelley's 
% directional derivative (dirder.m) and Jacobian computation (diffjac.m) 
% functions. A central differencing step was included below to reduce 
% the impact of fixed-point approximation error on the Jacobian. 
% Inputs: 
% fixedpt = fixed point, in vector form 
% epsln = relative perturbation size for use with diffjac_mod
% modelname = a string containing the model name, such as 'lrd'. This
% assumes that there is a [modelname '_p2p'] function available, but  
% jacobian_cd could be edited to accommodate a more general name. 
% Output: jacobian = approxiate Jacobian matrix. 
% 
% LRd: epsln = 1e-5 may be better for biphasic V, 1e-7 for monophasic K+ stimulus
% This function is based on code originally used in lrd_batch_eig.m. 
% Laura Munoz, June 2017

% Forward-difference approximation to Jacobian
jacfwd = diffjac_mod(fixedpt,[modelname '_p2p'],feval([modelname '_p2p'],fixedpt),epsln); % Solve for empirical bcl-to-bcl Jacobian. 

% Backward-difference approximation to Jacobian
jacback = diffjac_mod(fixedpt,[modelname '_p2p'],feval([modelname '_p2p'],fixedpt),-epsln); % Solve for empirical bcl-to-bcl Jacobian using perturbation of opposite sign

% Central-difference approximation to Jacobian
jacobian = (jacfwd+jacback)/2;
