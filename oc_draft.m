% For given C and B matrices, this script computes eigenvectors using
% several different approaches, and plots the resulting modal observability 
% and controllabilty measures. The main purpose of the script is to provide
% an assessment of which computational method is "best" in terms of 
% yielding the closest agreement with the eigenvector-eigenvalue equation, 
% and to help determine whether different computational methods produce 
% similar trends in obsv or ctrb with respect to BCL. 
% The tested methods considered are as follows: 
% Right eigenvectors: 
%   Default eig method with and without balancing, default eigs method
% Left eigenvectors: 
%   Default eig method with and without balancing, inversion without
%   balancing, transpose with and without balancing, and default eigs
%   method with transpose. 
% For the tested C and B matrices (only tried one pair so far), the default
% method without balancing seems best. 

clear variables;

scalingflag = 1;
numberofmodes = 5; % focus on n largest eigenvalues, where n = numberofmodes
selected_logepsln = -5; %
jacfolder = 'jacobians/'; % folder where Jacobians are stored
%ocfolder = 'ocvalues/'; % folder where obsv & ctrb values are stored
markersize = 18; % marker size for scatter plots

C=zeros(1,17);
B = C';
C(1) = 1;
B(1) = 1;

eval(['load ' jacfolder 'jacfile' num2str(selected_logepsln) ' *']) %Load data from jacobians
nbcls_sfp = length(selected_bcls_for_fps);

load('b1000fsolem12variable_amplitudes.mat'); % load data used in scaling
Smat = diag(1./varamp); % scaling matrix
Smatinv = inv(Smat); % only need to compute once
umax = diag(ones(1,size(B,2))); % This is also incorrect and just a placeholder.
% For inputs, also need to know approximate maximum values of each element
% of input vector. Note that u is a deviational quantity that may or
% may not attain the size of the stimulus.
% Each umax diagonal element should probably be the size of
% deviational max for integrated system.
Bs = Smat*B*umax; % scaled B matrix
Cs = C*Smatinv; % scaled C matrix 

cdotv = NaN*ones(nbcls_sfp, numberofmodes);
cmagvmag = NaN*ones(nbcls_sfp, numberofmodes);
obsvmag = NaN*ones(nbcls_sfp, numberofmodes);

cdotveigs = NaN*ones(nbcls_sfp, numberofmodes);
cmagvmageigs = NaN*ones(nbcls_sfp, numberofmodes);
obsvmageigs = NaN*ones(nbcls_sfp, numberofmodes);

cdotvnb = NaN*ones(nbcls_sfp, numberofmodes);
cmagvmagnb = NaN*ones(nbcls_sfp, numberofmodes);
obsvmagnb = NaN*ones(nbcls_sfp, numberofmodes);

bdotw = NaN*ones(nbcls_sfp, numberofmodes);
bmagwmag = NaN*ones(nbcls_sfp, numberofmodes);
ctrbmag = NaN*ones(nbcls_sfp, numberofmodes);

bdotwnb = NaN*ones(nbcls_sfp, numberofmodes);
bmagwmagnb = NaN*ones(nbcls_sfp, numberofmodes);
ctrbmagnb = NaN*ones(nbcls_sfp, numberofmodes);

bdotwnbfi = NaN*ones(nbcls_sfp, numberofmodes);
bmagwmagnbfi = NaN*ones(nbcls_sfp, numberofmodes);
ctrbmagnbfi = NaN*ones(nbcls_sfp, numberofmodes);

bdotwtr = NaN*ones(nbcls_sfp, numberofmodes);
bmagwmagtr = NaN*ones(nbcls_sfp, numberofmodes);
ctrbmagtr = NaN*ones(nbcls_sfp, numberofmodes);

bdotwnbtr = NaN*ones(nbcls_sfp, numberofmodes);
bmagwmagnbtr = NaN*ones(nbcls_sfp, numberofmodes);
ctrbmagnbtr = NaN*ones(nbcls_sfp, numberofmodes);

bdotwtreigs = NaN*ones(nbcls_sfp, numberofmodes);
bmagwmagtreigs = NaN*ones(nbcls_sfp, numberofmodes);
ctrbmagtreigs = NaN*ones(nbcls_sfp, numberofmodes);


righterror = NaN*ones(1,nbcls_sfp);
righterror_nb = NaN*ones(1,nbcls_sfp);
righterror_eigs = NaN*ones(1,nbcls_sfp);

lefterror = NaN*ones(1,nbcls_sfp);
lefterror_nb = NaN*ones(1,nbcls_sfp);
lefterror_finb = NaN*ones(1,nbcls_sfp);
lefterror_tr = NaN*ones(1,nbcls_sfp);
lefterror_nbtr = NaN*ones(1,nbcls_sfp);
lefterror_treigs = NaN*ones(1,nbcls_sfp);

flageigs = NaN*ones(nbcls_sfp); 
flageigsw = NaN*ones(nbcls_sfp); 

% plot symbols
symbols = char('bo-','rs-','gp-','m*-','k^-','cv-','yh-');

hor = figure;
hold on;

hoi = figure;
hold on;

hoa = figure;
hold on;

hcr = figure;
hold on;

hci = figure;
hold on;

hca = figure;
hold on;


for ii = 1:nbcls_sfp
    if scalingflag
        jactemp = Smat*alljacs{ii}*Smatinv;
        C = Cs;
        B = Bs;
    else
        jactemp = alljacs{ii};
    end
    
    % Compute eigenvalues and eigenvectors (use several methods) 
    [v,d,w]= eig(jactemp);
    [v_nb,d_nb,w_nb]=eig(jactemp,'nobalance');
    [ve,de,flageigs(ii)]= eigs(jactemp,numberofmodes);

    [w_tr,d_tr]=eig(jactemp');
    [w_nbtr,d_nbtr]=eig(jactemp','nobalance');
    [w_tre,d_tre,flageigsw(ii)]= eigs(jactemp',numberofmodes);
    
    % Try an alternate left eigenvector computation method, since 
    % default method does not yield W'*V = I 
%    wtemp = inv(v);
    % For LRd system, computing left eigenvectors from inverse doesn't
    % appear to be very accurate, unless nobalance option is used first
%    w_frominv = wtemp'; % complex conjugate transpose, so that columns are left eigenvectors

    wtempnb = inv(v_nb); % complex conjugate transpose, so that columns are left eigenvectors
    w_frominvnb = wtempnb';
    
    % Compute right and left eigenvector-eigenvalue errors, to determine
    % which methods are more accurate
    righterror(ii) = norm(jactemp*v - v*d);
    righterror_nb(ii) = norm(jactemp*v_nb-v_nb*d_nb);
    righterror_eigs(ii) = norm(jactemp*ve-ve*de);

    lefterror(ii) = norm(w'*jactemp - d*w');
    lefterror_nb(ii) = norm(w_nb'*jactemp - d_nb*w_nb');
%    lefterror_fi = norm(w_frominv'*jactemp - d*w_frominv');
    lefterror_finb(ii) = norm(w_frominvnb'*jactemp - d_nb*w_frominvnb');
    
    lefterror_tr(ii) = norm(w_tr'*jactemp - d_tr'*w_tr');
    lefterror_nbtr(ii) = norm(w_nbtr'*jactemp - d_nbtr'*w_nbtr');
    lefterror_treigs(ii) = norm(w_tre'*jactemp - d_tre'*w_tre');

    % Sort eigenvalues and eigenvectors in order of descending eigenvalue
    % magnitude
    [sortedeig,eigsortind] = sort(abs(diag(d)),'descend');
    sortedv = v(:,eigsortind);
    sortedw = w(:,eigsortind);
%    sortedw_frominv = w_frominv(:,eigsortind);
    [sortedeige,eigsortinde] = sort(abs(diag(de)),'descend');
    sortedve = ve(:,eigsortinde); 

    [sortedeignb,eigsortindnb] = sort(abs(diag(d_nb)),'descend');
    sortedvnb = v_nb(:,eigsortindnb);
    sortedwnb = w_nb(:,eigsortindnb);
    sortedwnb_frominv = w_frominvnb(:,eigsortindnb);

    [sortedeigtr,eigsortindtr] = sort(abs(diag(d_tr)),'descend');
    sortedwtr = w_tr(:,eigsortindtr);
    
    [sortedeignbtr,eigsortindnbtr] = sort(abs(diag(d_nbtr)),'descend');
    sortedwnbtr = w_nbtr(:,eigsortindnbtr);

    [sortedeigtre,eigsortindtre] = sort(abs(diag(d_tre)),'descend');
    sortedwtre = w_tre(:,eigsortindtre);
    

    for jj = 1:numberofmodes
        
        cdotv(ii,jj) = C*sortedv(:,jj); 
        cmagvmag(ii,jj) = norm(C)*norm(sortedv(:,jj)); 
        obsvmag(ii,jj) = abs(cdotv(ii,jj))/cmagvmag(ii,jj); 

        cdotvnb(ii,jj) = C*sortedvnb(:,jj); 
        cmagvmagnb(ii,jj) = norm(C)*norm(sortedvnb(:,jj)); 
        obsvmagnb(ii,jj) = abs(cdotvnb(ii,jj))/cmagvmagnb(ii,jj); 

        cdotveigs(ii,jj) = C*sortedve(:,jj); 
        cmagvmageigs(ii,jj) = norm(C)*norm(sortedve(:,jj)); 
        obsvmageigs(ii,jj) = abs(cdotveigs(ii,jj))/cmagvmageigs(ii,jj); 
        
        bdotw(ii,jj) = sortedw(:,jj)'*B; 
        bmagwmag(ii,jj) = norm(B)*norm(sortedw(:,jj));
        ctrbmag(ii,jj) = abs(bdotw(ii,jj))/bmagwmag(ii,jj); 
        
        bdotwnb(ii,jj) = sortedwnb(:,jj)'*B; 
        bmagwmagnb(ii,jj) = norm(B)*norm(sortedwnb(:,jj));
        ctrbmagnb(ii,jj) = abs(bdotwnb(ii,jj))/bmagwmagnb(ii,jj); 

        bdotwnbfi(ii,jj) = sortedwnb_frominv(:,jj)'*B; 
        bmagwmagnbfi(ii,jj) = norm(B)*norm(sortedwnb_frominv(:,jj));
        ctrbmagnbfi(ii,jj) = abs(bdotwnbfi(ii,jj))/bmagwmagnbfi(ii,jj); 

        bdotwtr(ii,jj) = sortedwtr(:,jj)'*B; 
        bmagwmagtr(ii,jj) = norm(B)*norm(sortedwtr(:,jj));
        ctrbmagtr(ii,jj) = abs(bdotwtr(ii,jj))/bmagwmagtr(ii,jj); 
        
        bdotwnbtr(ii,jj) = sortedwnbtr(:,jj)'*B; 
        bmagwmagnbtr(ii,jj) = norm(B)*norm(sortedwnbtr(:,jj));
        ctrbmagnbtr(ii,jj) = abs(bdotwnbtr(ii,jj))/bmagwmagnbtr(ii,jj);
        
        bdotwtreigs(ii,jj) = sortedwtre(:,jj)'*B; 
        bmagwmagtreigs(ii,jj) = norm(B)*norm(sortedwtre(:,jj));
        ctrbmagtreigs(ii,jj) = abs(bdotwtreigs(ii,jj))/bmagwmagtreigs(ii,jj); 
        
        cdotvchosen = cdotveigs; 
        cmagvmagchosen = cmagvmageigs;
        obsvmagchosen = obsvmageigs; 
        sortedeigchosenv = sortedeige; 

        bdotwchosen = bdotwtreigs; 
        bmagwmagchosen = bmagwmagtreigs;
        ctrbmagchosen = ctrbmagtreigs; 
        sortedeigchosenw = sortedeigtre; 

        figure(hor)
        plot(selected_bcls_for_fps(ii),real(cdotvchosen(ii,jj))/cmagvmagchosen(ii,jj),symbols(jj,:));
        figure(hoi)
        plot(selected_bcls_for_fps(ii),imag(cdotvchosen(ii,jj))/cmagvmagchosen(ii,jj),symbols(jj,:));
        figure(hoa)
        %        plot(selected_bcls_for_fps(ii),abs(C*sortedv(:,jj)),symbols(jj,:));
        scatter(selected_bcls_for_fps(ii),abs(sortedeigchosenv(jj)),markersize,obsvmagchosen(ii,jj));
        
        figure(hcr)
        plot(selected_bcls_for_fps(ii),real(bdotwchosen(ii,jj))/bmagwmagchosen(ii,jj),symbols(jj,:));
        figure(hci)
        plot(selected_bcls_for_fps(ii),imag(bdotwchosen(ii,jj))/bmagwmagchosen(ii,jj),symbols(jj,:));
        figure(hca)
        %        plot(selected_bcls_for_fps(ii),abs(sortedw(:,jj)'*B),symbols(jj,:));
        scatter(selected_bcls_for_fps(ii),abs(sortedeigchosenw(jj)),markersize,ctrbmagchosen(ii,jj));
    end
end

figure(hor)
xlabel('BCL, ms')
ylabel('Re C*v')
if scalingflag
    title('with scaling')
end

figure(hoi)
xlabel('BCL, ms')
ylabel('Im C*v')
if scalingflag
    title('with scaling')
end

figure(hoa)
xlabel('BCL, ms')
ylabel('Eigenvalue modulus, |\lambda|')
if scalingflag
    title('with scaling')
end
c = colorbar;
c.Label.String = '|Cv|';

figure(hcr)
xlabel('BCL, ms')
ylabel('Re w^T*B')
if scalingflag
    title('with scaling')
end

figure(hci)
xlabel('BCL, ms')
ylabel('Im w^T*B')
if scalingflag
    title('with scaling')
end

figure(hca)
xlabel('BCL, ms')
ylabel('Eigenvalue modulus, |\lambda|')
if scalingflag
    title('with scaling')
end
c = colorbar;
c.Label.String = '|w^TB|';

figure
hold on;
plot(selected_bcls_for_fps,righterror)
plot(selected_bcls_for_fps,righterror_nb,'go--')
plot(selected_bcls_for_fps,righterror_eigs,'r.--')
xlabel('BCL, ms')
ylabel('right eigenvector error norm') 
if scalingflag
    title('with scaling')
end
legend('with balancing', 'no balancing', 'eigs')

figure
hold on;
plot(selected_bcls_for_fps,lefterror)
plot(selected_bcls_for_fps,lefterror_nb,'go--')
plot(selected_bcls_for_fps,lefterror_finb,'r.--')
plot(selected_bcls_for_fps,lefterror_tr,'ks--')
plot(selected_bcls_for_fps,lefterror_nbtr,'cp--')
plot(selected_bcls_for_fps,lefterror_treigs,'mh:')
xlabel('BCL, ms')
ylabel('left eigenvector error norm') 
if scalingflag
    title('with scaling')
end
legend('with balancing', 'no balancing','inverse, no balancing', 'transpose','transpose, no balancing','transpose, eigs')

figure
hold on;
plot(selected_bcls_for_fps,lefterror_nb,'go--')
plot(selected_bcls_for_fps,lefterror_finb,'r.--')
plot(selected_bcls_for_fps,lefterror_nbtr,'cp--')
plot(selected_bcls_for_fps,lefterror_treigs,'mh:')
xlabel('BCL, ms')
ylabel('left eigenvector error norm') 
if scalingflag
    title('with scaling')
end
legend('no balancing','inverse, no balancing','transpose, no balancing', 'transpose, eigs')

figure
hold on;
plot(selected_bcls_for_fps,mean(obsvmag,2))
plot(selected_bcls_for_fps,mean(obsvmagnb,2),'go--')
plot(selected_bcls_for_fps,mean(obsvmageigs,2),'r.--')
xlabel('BCL, ms')
ylabel('average observability magnitude over largest modes') 
if scalingflag
    title('with scaling')
end
legend('with balancing', 'no balancing', 'eigs')

figure
hold on;
plot(selected_bcls_for_fps,mean(ctrbmag,2))
plot(selected_bcls_for_fps,mean(ctrbmagnb,2),'go--')
plot(selected_bcls_for_fps,mean(ctrbmagnbfi,2),'r.--')
plot(selected_bcls_for_fps,mean(ctrbmagtr,2),'ks--')
plot(selected_bcls_for_fps,mean(ctrbmagnbtr,2),'cp--')
plot(selected_bcls_for_fps,mean(ctrbmagtreigs,2),'mh:')
xlabel('BCL, ms')
ylabel('average controllability magnitude over largest modes') 
if scalingflag
    title('with scaling')
end
legend('with balancing', 'no balancing','inverse, no balancing', 'transpose','transpose, no balancing', 'transpose, eigs')
