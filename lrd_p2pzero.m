function out = lrd_p2pzero(x)
% Subtract off initial condition (since we want to
% find an x at which lrd_p2p_rev(x)-x = 0). 

%xp1 = main_LRd_1cell(x);
xp1 = lrd_p2p(x);
out=xp1-x;
