%LMM: This is  another version of a script to compute variations between
%Jacobians as epsilon is varied. Predecessors are jac_diff.m and
%jac_diff2.m. While I could have used version control instead of numbering
%the scripts, they don't have that much overlapping functionality,
%due to changes in notation.
%
%This script should be run after computejacobians.m. 
%If provided with Jacobians, the script produces plots of the norms of the
%Jacobian and eigenvalue differences, for successive values of epsilon. Our
%assumption is that ranges of epsilon values that produce less variation in
%these norms yield better estimates of the true Jacobians.

clear variables;

scalingflag = 1; % 0 to leave off pseudo-normalization, 1 to include
altplotflag = 0; % this will automatically be changed to 1 if the number of 
% curves exceeds the number of elements in "symbols" 

jacfolder = 'jacobians/def/'; % folder where Jacobians are stored

logepslns = -7:-3; % These are the log10's of the epsilon values.
% These are limited to values that correspond to files stored in the
% Jacobians folder.
numexps = length(logepslns); % number of exponents

load('b1000fsolem12variable_amplitudes.mat'); % load data used in scaling 
Smat = diag(1./varamp); % scaling matrix
Smatinv = inv(Smat); % only need to compute once

% plot symbols
symbols = char('bo-','rs-','gp-','m*-','k^-','cv-','yh-');

%selectbcls = [1000 200 70]; % select BCLs for comparison
eval(['load ' jacfolder 'jacfile' num2str(logepslns(1)) ' selected_bcls_for_fps']) %Load data from jacobians    
% the above line assumes all bcl lists are the same for every epsilon 
selectbcls = selected_bcls_for_fps; % select BCLs for comparison
numselectbcls = length(selectbcls);

if numselectbcls > size(symbols,1)
    altplotflag = 1;
    symbols = repmat('bx-',numselectbcls,1); 
end

bcllabels = ['BCL = ' num2str(selectbcls(1))];
for kk = 2:numselectbcls
    bcllabels = char(bcllabels,['BCL = ' num2str(selectbcls(kk))]);
end


% Jacobian figure
h=figure;
hold on;

% Eigenvalue figure
he=figure;
hold on;

% Jacobian element figure (similar to Mingyi Li & Niels' Fig 3 in their ABME
% 2003 paper) 
hj=figure;
hold on;

jacsatbcls = cell(numexps,numselectbcls);
jacdiffs = cell(numexps-1,numselectbcls);
jacdiffnorms = NaN*ones(numexps-1,numselectbcls);
eigsatbcls = cell(numexps,numselectbcls);
eigdiffs = cell(numexps-1,numselectbcls);
eigdiffnorms = NaN*ones(numexps-1,numselectbcls);
minjacdiff = NaN*ones(1,numselectbcls); 
mineigdiff = NaN*ones(1,numselectbcls); 
minjacdiffind = NaN*ones(1,numselectbcls); 
mineigdiffind = NaN*ones(1,numselectbcls); 

for ii = 1:numexps
    disp(['logepsln = ' num2str(logepslns(ii))])

    eval(['load ' jacfolder 'jacfile' num2str(logepslns(ii)) ' *']) %Load data from jacobians    
    
    % Check the epsln stored in the Jacobian file 
    if log10(epsln) ~= logepslns(ii)
        disp(['Error: epsilon mismatch, selected index = ' num2str(logepslns(ii))])
        %        return;
    end
         
    bclsfps = selected_bcls_for_fps; % Shorten this name for better readability
    
    for kk = 1:numselectbcls
        disp(['BCL = ' num2str(bclsfps(bclsfps==selectbcls(kk)))])
        if scalingflag && ~isempty(alljacs{bclsfps==selectbcls(kk)})
            jacsatbcls{ii,kk} = Smat*alljacs{bclsfps==selectbcls(kk)}*Smatinv;
        else
            jacsatbcls{ii,kk} = alljacs{bclsfps==selectbcls(kk)};
        end
        eigstemp = eig(jacsatbcls{ii,kk}, 'nobalance'); % alternatively, the eigenvalues
        % could be loaded from their folder instead of being recomputed
        % here, though the loading method requires more extracting and
        % checking of values 
        [sorttemp,eigsortind] = sort(abs(eigstemp));
        % Warning: sort function works differently depending on whether elements
        % are complex or not, and sometimes real numbers will have a
        % 0.0000i term tacked on.
        eigsatbcls{ii,kk} = eigstemp(eigsortind);
    end
    
end
% Compute differences and norms
for jj = 1:numexps-1
    for kk = 1:numselectbcls
        % Some of the matrices might be empty, so compare sizes
        if size(jacsatbcls{jj,kk}) == size(jacsatbcls{jj+1,kk})
            jacdiffs{jj,kk} = jacsatbcls{jj,kk} - jacsatbcls{jj+1,kk};
            jacdiffnorms(jj,kk) = norm(jacdiffs{jj,kk});
            eigdiffs{jj,kk} = eigsatbcls{jj,kk} - eigsatbcls{jj+1,kk};
            eigdiffnorms(jj,kk) = norm(eigdiffs{jj,kk});
        end
    end
end

for kk = 1:numselectbcls
    [minjacdiff(kk),minjacdiffind(kk)] = min(jacdiffnorms(:,kk));
    [mineigdiff(kk),mineigdiffind(kk)] = min(eigdiffnorms(:,kk));
end

meanminjacdiffind = mean(minjacdiffind); 
meanmineigdiffind = mean(mineigdiffind); 

logepslns(round(meanminjacdiffind))
logepslns(round(meanmineigdiffind))

% plot norms of differences
for kk = 1:numselectbcls
    figure(h)
    plot(logepslns(1:end-1)+0.5, jacdiffnorms(:,kk),symbols(kk,:));
    figure(he)
    plot(logepslns(1:end-1)+0.5, eigdiffnorms(:,kk),symbols(kk,:));
    figure(hj)
    for ii=1:numexps
        if ii==1
            p(kk) = plot(logepslns(ii), jacsatbcls{ii,kk}(1,1),symbols(kk,:));
        else
            plot(logepslns(ii), jacsatbcls{ii,kk}(1,1),symbols(kk,:));
        end
        
    end
end
figure(h)
xlabel('log_{10}(epsilon)')
set(gca,'YScale','log');
ylabel('norm of difference')
if scalingflag
    title('Jacobian differences: Revised method, with scaling')
else
    title('Jacobian differences: Revised method')
end
if ~altplotflag
    legend(bcllabels)
end
grid on;
%axis([-8 -3 10^-4 10^2])

figure(he)
xlabel('log_{10}(epsilon)')
set(gca,'YScale','log');
ylabel('norm of difference')
if scalingflag
    title('Eigenvalue differences, with scaling')
else
    title('Eigenvalue differences')
end
if ~altplotflag
    legend(bcllabels)
end
grid on;

figure(hj)
xlabel('log_{10}(epsilon)')
ylabel('J_{1,1}')
if scalingflag
    title('1,1 Jacobian element, with scaling')
else
    title('1,1 Jacobian element')
end
if ~altplotflag
    legend(bcllabels)
end
grid on;
