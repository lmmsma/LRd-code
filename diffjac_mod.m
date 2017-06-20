
%function [jac] =diffjac(x, f, f0)
function [jac] =diffjac(x, f, f0, epsilon)
% compute a forward difference Jacobian f'(x), return lu factors
% uses dirder.m to compute the columns
% Modified to return f'(x) instead of lu factors, and to take relative
% perturbation size epsilon as an argument
%
% C. T. Kelley, November 25, 1993
%
% This code comes with no guarantee or warranty of any kind.
%
%
% inputs:
%         x, f = point and function
%		  f0   = f(x), preevaluated
%
n=length(x);
for j=1:n
    zz=zeros(n,1);
    zz(j)=1;
%    jac(:,j)=dirder_mod(x,zz,f,f0);
    jac(:,j)=dirder_mod(x,zz,f,f0,epsilon);
end
%save jac jac; 
save jac jac epsilon; 
%[l, u] = lu(jac);
%atv = jac*v; 
