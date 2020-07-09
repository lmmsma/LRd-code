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

%shiftstring = ''; % if using default shift of data.dt
%shiftstring = '_shift0p75ms';
%shiftstring = '_shift6p5ms';
%shiftstring = '_shift7mVrepol';
%shiftstring = '_shift-50mVrepol';
%shiftstring = '_shift0p2Vnormdepol';
%shiftstring = '_shift0p4Vnormdepol';
%shiftstring = '_shift0p6Vnormdepol';
%shiftstring = '_shift0p8Vnormdepol';
%shiftstring = '_shift1Vnormdepol';
%shiftstring = '_shift0p8Vnormrepol';
%shiftstring = '_shift0p6Vnormrepol';
%shiftstring = '_shift0p4Vnormrepol';
%shiftstring = '_shift0p2Vnormrepol';
shiftstring = '_shift0p001Vnormrepol';

scalingflag = 0; % 0 to leave off pseudo-normalization, 1 to include
altplotflag = 0; % this will automatically be changed to 1 if the number of 
% curves exceeds the number of elements in "symbols" 
fontsize = 28; 

paramflag = input('Default or adjusted parameters (enter 0 if default, 1 if adjusted): ') == 1;
if paramflag
    param = 'adj';
else
    param = 'def';
end

%jacfolder = ['jacobians/' param '/']; % folder where Jacobians are stored
%jacfolder = ['jacobians_shift6p5ms/' param '/']; % folder where Jacobians are stored
jacfolder = ['jacobians' shiftstring '/' param '/']; % folder where Jacobians are stored

logepslns = -7:-3; % These are the log10's of the epsilon values.
% These are limited to values that correspond to files stored in the
% Jacobians folder.
numexps = length(logepslns); % number of exponents

load('b1000fsolem12variable_amplitudes.mat'); % load data used in scaling 
Smat = diag(1./varamp); % scaling matrix
Smatinv = inv(Smat); % only need to compute once

% plot symbols
symbols = char('bo-','rs-','gp-','m*-','k^-','cv-','yh-');

highlightbcls = [1000 200 70]; % Emphasize certan BCLs on the plot
eval(['load ' jacfolder 'jacfile' num2str(logepslns(1)) ' selected_bcls_for_fps']) %Load data from jacobians    
% the above line assumes all bcl lists are the same for every epsilon 
selectbcls = selected_bcls_for_fps; % select BCLs for comparison
%selectbcls = 1000:-100:100; % select BCLs for comparison
numselectbcls = length(selectbcls);
bclindices = 1:numselectbcls; 


if numselectbcls > size(symbols,1)
    altplotflag = 1;
    highlightsymbols = symbols; % Highlight certain BCLs instead
%    symbols = repmat('bx-',numselectbcls,1);
    symbols = repmat('k:',numselectbcls,1);
    bcllabels = ['BCL = ' num2str(highlightbcls(1)) ' ms'];
    for kk = 2:length(highlightbcls)
        bcllabels = char(bcllabels,['BCL = ' num2str(highlightbcls(kk)) ' ms']);
    end
else % if list of BCLs is short
    bcllabels = ['BCL = ' num2str(selectbcls(1)) ' ms'];
    for kk = 2:numselectbcls
        bcllabels = char(bcllabels,['BCL = ' num2str(selectbcls(kk)) ' ms']);
    end
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

%    eval(['load ' jacfolder 'jacfile' num2str(logepslns(ii)) ' *']) %Load data from jacobians 
% loading everything may overwrite jacfolder
    eval(['load ' jacfolder 'jacfile' num2str(logepslns(ii)) ' epsln selected_bcls_for_fps alljacs']) %Load data from jacobians    
    
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

logepsmeanminjacdiffind = logepslns(round(meanminjacdiffind)); 
logepsmeanmineigdiffind = logepslns(round(meanmineigdiffind));

logepsmeanminjacdiffind 
logepsmeanmineigdiffind 

% Compute max, mean, and min (over all BCLs) differences at the
% difference minimizing index
maxeigdiffnorms = max(eigdiffnorms(round(meanmineigdiffind),:))
meaneigdiffnorms = mean(eigdiffnorms(round(meanmineigdiffind),:))
mineigdiffnorms = min(eigdiffnorms(round(meanmineigdiffind),:))

% plot norms of differences
for kk = 1:numselectbcls
    figure(h)
    plot(logepslns(1:end-1)+0.5, jacdiffnorms(:,kk),symbols(kk,:),'linewidth',0.5);
    figure(he)
    plot(logepslns(1:end-1)+0.5, eigdiffnorms(:,kk),symbols(kk,:),'linewidth',0.5);
    figure(hj)
    for ii=1:numexps
        if ii==1
            p(kk) = plot(logepslns(ii), jacsatbcls{ii,kk}(1,1),[symbols(kk,:) 'x']);
        else
            plot(logepslns(ii), jacsatbcls{ii,kk}(1,1),[symbols(kk,:) 'x']);
        end
        
    end
end


% plot norms of differences for highlighted bcls
for kk = 1:length(highlightbcls)
    figure(h)
    pj(kk) = plot(logepslns(1:end-1)+0.5, jacdiffnorms(:,bclindices(selectbcls == highlightbcls(kk))),highlightsymbols(kk,:),'Linewidth',2);
    figure(he)
    pe(kk) = plot(logepslns(1:end-1)+0.5, eigdiffnorms(:,bclindices(selectbcls == highlightbcls(kk))),highlightsymbols(kk,:),'Linewidth',2);
end


figure(h)
xlabel('log_{10} \epsilon')
set(gca,'YScale','log');
ylabel('norm of difference')
if scalingflag
    title('Jacobian differences: Revised method, with scaling')
else
    title('Jacobian differences: Revised method')
end
%if ~altplotflag
    legend(pj, bcllabels)
%end
grid on;
%axis([-8 -3 10^-4 10^2])

figure(he)
xlabel('log_{10} \epsilon')
set(gca,'YScale','log');
%ylabel('norm of difference')
ylabel('||\Delta L||')
% if scalingflag
%     title('Eigenvalue differences, with scaling')
% else
%     title('Eigenvalue differences')
% end
%if ~altplotflag
    leghandle = legend(pe,bcllabels); 
%end
%grid on;
%set(leghandle,'Position',[0.5390    0.17    0.28   0.1415])
set(leghandle,'Location','northwest')
set(gca,'FontSize',fontsize)
set(gcf,'units','normalized','outerposition',[0.1787    0.0958    0.55    0.79])

if param == 'def'
    text(-5,10^-7.5, 'default ','fontsize',fontsize)
    figname = ['eigdiffdef' shiftstring];
elseif param == 'adj'
    text(-5,10^-7.5, 'adjusted ','fontsize',fontsize)
    figname = ['eigdiffadj' shiftstring];
end
axis([-6.5 -3.5 10^-8 10^-2])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,figname,'-dpdf')
saveas(fig,figname)

figure(hj)
xlabel('log_{10}\epsilon')
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
