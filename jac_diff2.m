%LMM: modified version of Anthony and Ryan's jac_diff.m. Some code was
%borrowed from the now defunct jac_diff_mod.m.

clear variables;

methodflag = 1; % 0 for A&R's orignal method, 1 for revised method.
scalingflag = 0; % 0 to leave off pseudo-normalization, 1 to include

% Home folder:
homefolder = 'C:\Users\laura\''Google Drive''\LRd-code-master_mod\';

% Folder where Jacobians are stored:
sourcefolder = 'C:\Users\laura\''Google Drive''\''REU 2017''\matlab\''Final Matlab Files''\Stored\';

load('b1000fsolem12variable_amplitudes.mat');
Smat = diag(1./varamp);
Smatinv = inv(Smat); % only need to compute once

if methodflag
    h=figure;
    hold on;

    hh=figure;
    hold on;

    symbols = char('bo-','rs-','gp-','m*-','k^-','cv-','yh-');
    
    % List logs of epsilon that were tested, e.g. -5 stands for
    % epsln = 10^-5 :
    
    %exps = -5;
    exps = -8:-3;
    numexps = length(exps); 
    
    selectbcls = [1000 200 150 130 100 90 70]; % select BCLs for comparison
    bcllabels = ['BCL = ' num2str(selectbcls(1))];
    numselectbcls = length(selectbcls); 
    
    for kk = 2:numselectbcls
        bcllabels = char(bcllabels,['BCL = ' num2str(selectbcls(kk))]);
    end
    
    %     alljac1000 = cell(numexps,1);
    %     alljac200 = cell(numexps,1);
    %     alljac70 = cell(numexps,1);
    jacsatbcls = cell(numexps,numselectbcls);
    jacdiffs = cell(numexps-1,numselectbcls);
    jacdiffnorms = NaN*ones(numexps-1,numselectbcls);
    eigsatbcls = cell(numexps,numselectbcls);
    eigdiffs = cell(numexps-1,numselectbcls);
    eigdiffnorms = NaN*ones(numexps-1,numselectbcls);
    
    for ii = 1:numexps
        expfolder = [sourcefolder 'eps=10_' num2str(exps(ii)) '\'];
        eval(['cd ' expfolder]);
        
        if exist('jacfile.mat', 'file')
            load jacfile.mat;
            disp('Found jacfile.mat')
        end
        eval(['cd ' homefolder]);
        
        disp(['epsln = ' num2str(epsln)])
        if log10(epsln) ~= exps(ii)
            disp(['Error: epsilon mismatch, selected index = ' num2str(exps(ii))])
            %        return;
        end
        
        for kk = 1:numselectbcls
            disp(['BCL = ' num2str(bcls(bcls==selectbcls(kk)))])
            %            eval(['alljac' num2str(bcls(bcls==selectbcls(kk))) '{ii} = alljacs{bcls==selectbcls(kk)}']);
            if scalingflag && ~isempty(alljacs{bcls==selectbcls(kk)})
                jacsatbcls{ii,kk} = Smat*alljacs{bcls==selectbcls(kk)}*Smatinv;
            else
                jacsatbcls{ii,kk} = alljacs{bcls==selectbcls(kk)};
            end
            % Warning: sort works differently depending on whether elements
            % are complex or not, and sometimes real numbers will have a
            % 0.0000i term tacked on. 
            eigstemp = eig(jacsatbcls{ii,kk}); 
            [sorttemp,eigsortind] = sort(abs(eigstemp));
            eigsatbcls{ii,kk} = eigstemp(eigsortind);
        end
        
        %         disp(['BCL = ' num2str(bcls(bcls==200))])
        %         alljac200{ii} = alljacs{bcls==200};
        %
        %         disp(['BCL = ' num2str(bcls(bcls==70))])
        %         alljac70{ii} = alljacs{bcls==70};
        
    end
    % Compute differences and norms
    for jj = 1:numexps-1
        for kk = 1:numselectbcls
            % Some of the matrices might be empty, so compare sizes
            if size(jacsatbcls{jj,kk}) == size(jacsatbcls{jj+1,kk})
                jacdiffs{jj,kk} = jacsatbcls{jj,kk} - jacsatbcls{jj+1,kk};
                jacdiffnorms(jj,kk) = norm(jacdiffs{jj,kk});
                jacdiffnorms(jj,kk) = norm(jacdiffs{jj,kk});
                eigdiffs{jj,kk} = eigsatbcls{jj,kk} - eigsatbcls{jj+1,kk};
                eigdiffnorms(jj,kk) = norm(eigdiffs{jj,kk});
            end
        end
    end
    
    for kk = 1:numselectbcls
        figure(h)
        p(i) = plot(exps(1:end-1)+0.5, jacdiffnorms(:,kk),symbols(kk,:));
        figure(hh) 
        plot(exps(1:end-1)+0.5, eigdiffnorms(:,kk),symbols(kk,:));
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
    legend(bcllabels)
    grid on;
    axis([-8 -3 10^-4 10^2])
    
    figure(hh)
    xlabel('log_{10}(epsilon)')
    set(gca,'YScale','log');
    ylabel('norm of difference')
    title('Eigenvalue differences')
    legend(bcllabels)
    grid on;

    
    %     jacA_1000 = alljac1000{2};
    %     jacB_1000 = alljac1000{3};
    %     jacC_1000 = alljac1000{4};
    %     jacD_1000 = alljac1000{5};
    %     jacE_1000 = alljac1000{6};
    
else
    eval(['load ' sourcefolder 'eps=10_-7/jacfile alljacs']);
    jacA_1000= alljacs{46};
    eval(['load ' sourcefolder 'eps=10_-6/jacfile alljacs']);
    jacB_1000= alljacs{46};
    eval(['load ' sourcefolder 'eps=10_-5/jacfile alljacs']);
    jacC_1000= alljacs{46};
    eval(['load ' sourcefolder 'eps=10_-4/jacfile alljacs']);
    jacD_1000= alljacs{46};
    eval(['load ' sourcefolder 'eps=10_-3/jacfile alljacs']);
    jacE_1000= alljacs{43};
    
    %end % use this end statement instead, to try to test out original
    %computation method with new data reading method
    
    
    diff1=jacB_1000-jacA_1000;
    diff2=jacC_1000-jacB_1000;
    diff3=jacD_1000-jacC_1000;
    diff4=jacE_1000-jacD_1000;
    
    df1= zeros(17,1); df2= zeros(17,1); df3= zeros(17,1); df4= zeros(17,1);
    
    for i= 1:17
        df1(i)= sum(diff1(i , :));
        df2(i)= sum(diff2(i , :));
        df3(i)= sum(diff3(i , :));
        df4(i)= sum(diff4(i , :));
    end
    
    norm_df1= Smat*df1;
    norm_df2= Smat*df2;
    norm_df3= Smat*df3;
    norm_df4= Smat*df4;
    n1= norm(norm_df1)
    n2= norm(norm_df2)
    n3= norm(norm_df3)
    n4= norm(norm_df4)
    
    figure
    semilogy(-6.5:-3.5, [n1 n2 n3 n4],'o')
    xlabel('log_{10}(epsilon)')
    ylabel('norm of difference')
    if methodflag
        title('Revised method')
    else
        title('A&R''s method')
    end
end