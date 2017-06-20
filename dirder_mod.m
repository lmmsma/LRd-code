%function z = dirder(x,w,f,f0)
function z = dirder(x,w,f,f0,epsilon)
% Finite difference directional derivative
% Approximate f'(x) w
% Modified to change difference increment to values other than 1d-7.
% Modified again to take epsilon as an argument, and use relative perturbations scaled to each individual element of the 
% fixed-point vector, rather than scaling epsilon by the entire norm of the
% fixed point as in the original code. This only works 
% in cases where the selected fixed point element is sufficiently different
% from zero. To give an idea of the norm(x0)-based scaling,  norm(solem12)
% = 1.690052685089206e+002, so in the original code, epsnew=1e-7
% corresponds to a perturbation size of 1.69e-5 applied to all
% variables. I'm not sure what the default behavior should be when the
% fixed-point element is near zero; here it seems plausible to treat
% epsnew as an absolute perturbation, or as a perturbation relative to the
% norm of the fixed point vector, or as a perturbation relative to some
% "small" threshold value. The middle option doesn't make much sense to me,
% since the reasoning behind scaling the perturbations for each state element
% is that smaller fixed-point values should be perturbed by smaller amounts, 
% and using epsnew*norm(x) would revert to using the largest perturbation. 

% 
% C. T. Kelley, November 25, 1993
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder

%
% Hardwired difference increment.
%epsnew=1.d-7;
epsnew = epsilon; 
%epsnew=1.d-12;
%
n=length(x);
%
% scale the step
%
if norm(w) == 0
    z=zeros(n,1);
return
end
epsnew = epsnew/norm(w);
%if norm(x) > 0
%    epsnew=epsnew*norm(x);
%end
x(w==1)
if abs(x(w==1))>1e-10 % if the element size is above some threshold value
    epsnew = epsnew*abs(x(w==1))
else
    disp(['Element size is below threshold, using absolute difference of ' num2str(epsnew)])
end
    
%
% del and f1 could share the same space if storage
% is more important than clarity
%
del=x+epsnew*w;
f1=feval(f,del);
z = (f1 - f0)/epsnew;
