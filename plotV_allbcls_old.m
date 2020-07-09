% Using unshifted fixed points as ICs, plot all APs on the same plot

clear variables;

%folder = ['lrddata_fromfixedpointsearch121718/'];
folder = ['lrddata_fromfixedpointsearch_highres011719/'];
%folder = ['lrddata/'];
%folder = ['lrddata_fromfixedpointsearch_highres_shift7mVrepol/'];

bcls = [1000:-50:400 390:-10:70];%Vector of bcls in ms

nbcls = length(bcls);

% Put this here since I initially didn't save all values in
% fp_compile2_shiftedtraj.m.
i_repol = zeros(1,length(bcls));
%offset = 7; % this a voltage threshold in mV, not an offset in ms
% Now add i_depol: 
i_depol = zeros(1,length(bcls));

Vthresh = 1; % threshold value for normalized voltages
i_thresh = zeros(1,length(bcls));

h=figure
hold on;

hnorm = figure
hold on;

for i = 1:nbcls
    bcl = bcls(i);
    
    fname = [folder 'lrddata_1cell_b' num2str(bcl)];
    % load simulation data for this bcl
    load(fname);
    V = Y(1,:); % membrane potential
    time = (1:length(V))*outstep;
    
    % borrow falling-edge detection from plot_restitution.m:
    % Search for falling edge crossing of threshold
    %    i_repol(i) = find(V(1:end-1)>offset & V(2:end)<=offset);
    
    ampl = abs(max(V)-min(V));
    Vmin = min(V);
    
    Vnorm = (V-Vmin)/ampl; % normalized membrane potential
    
    
    % borrow falling-edge detection from plot_restitution.m:
    % Search for falling and rising edge crossings of threshold
    idtemp = find(Vnorm(1:end-1)<Vthresh & Vnorm(2:end)>=Vthresh);
    irtemp = find(Vnorm(1:end-1)>Vthresh & Vnorm(2:end)<=Vthresh);
    if ~isempty(idtemp)
        i_depol(i) = idtemp;
    end
    if ~isempty(irtemp)
        i_repol(i) = irtemp;
    end
    
    figure(h)
    plot(time, V)
    if i_depol(i)
        plot(time(i_depol(i)),V(i_depol(i)),'rx')
    elseif i_repol(i)
        plot(time(i_repol(i)),V(i_repol(i)),'rx')
    end
    
    figure(hnorm)
    plot(time/bcl, (V-Vmin)/ampl)
    if i_depol(i)
        plot(time(i_depol(i))/bcl,Vnorm(i_depol(i)),'rx')
    elseif i_repol(i)
        plot(time(i_repol(i))/bcl,Vnorm(i_repol(i)),'rx')
    end
    
end
figure(h)
xlabel('time, ms')
ylabel('V, mV')

figure(hnorm)
xlabel('time/BCL, ms')
ylabel('(V-minV)/(maxV-minV), dimensionless')

% Be careful with this one. Just put it here since I forgot to compute it
% as an array earlier:
%save('lrddata_fromfixedpointsearch_highres_shift7mVrepol/compiled_fp.mat','i_repol','-append')