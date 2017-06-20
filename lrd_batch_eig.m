% Compute fixed point, eigenvalues, Jacobian, and ctrb/obsv of time-integrated LRd model
clear variables

ms = 10; fs = 14; % marker size and font size

statenames = strvcat('V','H','m','J','d','f','xr','ca_T','na_i','k_i','jsr_T','nsr','xs','B','G','xs2','Rel');
% Rough descriptions of the significant components of each eigenvalue
%eiglabels = strvcat('jsr_T', 'k_i', 'k_i-na_i', 'jsr_T-xs2', 'jsr_T, alt', 'V', 'V, alt', 'G', 'V', 'V', 'V-G, alt', 'G-J-H-V, alt', 'G-J-H-V, alt', 'G-m, alt', 'G-m, alt', 'G-m', 'm-G');
eiglabels_vbiph_em12 = strvcat('jsr_T', 'k_i', 'k_i-na_i', 'jsr_T-xs2', 'jsr_T, alt', 'V', 'V, alt', 'G', 'V, alt', 'V, alt', 'V, alt', 'G-J-V-H, alt', 'm-J-H-V-G', 'm-J-H-V-G', 'm-H-J-G', 'm-V-H-J, alt', 'H-J-G-V');
eiglabels_kmono_kem12 = strvcat('jsr_T', 'k_i', 'k_i-na_i', 'jsr_T-xs2', 'jsr_T, alt', 'V', 'V, alt', 'G', 'V', 'V, alt', 'V, alt', 'V, alt', 'B, alt', 'B, alt', 'B-m', 'm-xr, alt', 'm-H');

% Name of unperturbed fixed point
systemselect = 'solem12'
%systemselect = 'kmonovsolem12'

% Filenames and label settings based on fixed point type
if strcmp(systemselect, 'solem12')
    unpertfilename = 'b1000fsolem12_fwde_shift0_newpulse';
    el = eiglabels_vbiph_em12;
elseif strcmp(systemselect,'kmonovsolem12')
    unpertfilename = 'b1000fkmonovsolem12_fwde_shift0_newpulse_relpert1e-7';
    el = eiglabels_kmono_vem12;
end

% Load the unperturbed fixed point
eval(['load ' unpertfilename ' ' systemselect])
initialvalue = eval(systemselect);

numstate = length(statenames);
%modelname = 'karma_sim_estimtest';
modelname = 'lrd';

%L = squeeze(allgains(:,:,1));
numpart = 1; % number of cells
bcl = 1000; %bcl, ms
L = zeros(numstate*numpart,numpart); % observer feedback gain

deltat = .005; % ms
%stimstart = deltat; % stimulus start time, ms
allstimstart = deltat;
%allstimstart = deltat + 200;
%allstimstart = deltat + [400 600 800];
stimstart = allstimstart(1);

save lrdparams L numpart bcl stimstart % these will be read in by lrd_p2p.m

epsln = 1e-5; % relative perturbation size for use with diffjac_mod % LRd: 1e-5 may be better for biphasic V, 1e-7 for monophasic K+ stimulus
nsolitol = 1e-12; % abs & rel tolerance for Newton-Krylov solver
useoldfp = 1; % set to zero to always solve for fixed point; set to 1 to re-use fixed point
maxabsfperr = NaN; % maximum absolute value of p2p error (use to determine whether OL f.p. is a valid CL f.p.)
is2to1 = 0; % set to 1 if 2:1 block pattern
computejac = 1; % set to 1 to compute Jacobian explicitly
if numpart == 1
    computejac = 1;
end

% initialize nsoli output variables, in case usoldfp==1
it_hist = -1;
ierr=NaN;
x_hist=NaN;

fpfolder = 'lrddata_1cell_b1000_solem12_start1_Kp0'; % folder containing trajectory that passes through fixed point

% If computing Jacobians of OL system, can assess ctrb/obsv.
% input matrix
B = deltat*eye(numstate*numpart); % This isn't right, just a placeholder. In general, B isn't square and its nonzero elements aren't 1.
% output matrix
C = eye(numstate*numpart); % This is also a placeholder. There is a dimensional inconsistency with L above. If L is numstate*numpart by numpart, it is sized to fit the assumption that only V's can be measured. C should technically have the transposed dimensions of L.
%Vtmax = 120; %mV; approximate span of an AP
%ntmax = 1.1; %dimensionless; approximate max of n variable
%Smat = diag([(1/Vtmax)*ones(1,numpart) (1/ntmax)*ones(1,numpart)]); % Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
load b1000fsolem12variable_amplitudes varamp
Smat = diag(1./varamp);
Smatinv = inv(Smat); % only need to compute once
%u1max = 97; %muA/cm^2; max for current injections (based on estimate of max usf value for "smallest" stabilizing gain)
%u2max = 0.065; %dimensionless; approximate max of n-perturbations (based on estimate of max n-pert value for "smallest" stabilizing gain)
umax = diag(ones(1,size(B,2))); % This is also incorrect and just a placeholder. For inputs, also need to know approximate maximum values of each element of input vector. Note that u is a deviational quantity that may or may not attain the size of the stimulus.
%umax(1,1) = u1max;
%umax(2,2) = u2max;
% Each umax diagonal element should be size of deviational max for integrated system.
Bs = Smat*B*umax; % scaled B matrix
Cs = Smat*C*Smatinv; % scaled B matrix

rankcutoff = 1e-14; % below this level, singular values don't contribute to the rank

% Intitalize svd matrices
svdctrb = zeros(numstate*numpart);
svdctrb_scaled = zeros(numstate*numpart);
svdobsv = zeros(numstate*numpart);
svdobsv_scaled = zeros(numstate*numpart);

for i=1:length(allstimstart)%1%(size(allgains,3)-1)
    clear jacback jaccd
    % Intitalize svd matrices
    svdctrb = zeros(numstate);
    svdctrb_scaled = zeros(numstate);
    svdobsv = zeros(numstate);
    svdobsv_scaled = zeros(numstate);
    evcf = cell(numstate);
    evof = cell(numstate);
    evncf = cell(numstate);
    evnof = cell(numstate);
    rankc = zeros(1,numstate);
    ranko = zeros(1,numstate);
    rankcf = zeros(1,numstate);
    rankof = zeros(1,numstate);
    evcf_scaled = cell(numstate);
    evof_scaled = cell(numstate);
    evncf_scaled = cell(numstate);
    evnof_scaled = cell(numstate);
    rankc_scaled = zeros(1,numstate);
    ranko_scaled = zeros(1,numstate);
    rankcf_scaled = zeros(1,numstate);
    rankof_scaled = zeros(1,numstate);
    
    if numpart == 1 % Set to 1 to turn on automatic filename generation (useful for single-cell tests)
        % Select file name for saving results
        % Replace decimal with "p", otherwise it can cause problems loading files:
        for ii = 1:size(L,1)
            for jj = 1:size(L,2)
                temp = num2str(L(ii,jj)); % this can probably be vectorized but I'm not sure how (direct num2str will pad with spaces)
                kk = strfind(temp,'.'); % find decimal
                if kk
                    temp(kk) = 'p';
                end
                gainstr{ii}{jj} = temp;
            end
        end
        diffstr = '';
        % single-cell:
        fname1 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) '_' systemselect '_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_Kp' gainstr{1}{1}];
        fname2 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) '_' systemselect '_eig_relpert' num2str(-epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_Kp' gainstr{1}{1}];
        fname3 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) '_' systemselect '_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_cdiff_Kp' gainstr{1}{1}];
        % This naming convention will only work for single cells. Not sure what to do to summarize multicell feedback
    else  % Provide a filename suffix
        fnamesuffix = input('Enter filename suffix in single quotes, or the number zero to leave blank. The default suffix is ''Kp0'': ')
        diffstr = '_d0007'; % indicate diffusion coeff
        if computejac
            fname1 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_' systemselect '_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_' fnamesuffix];
            fname2 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_' systemselect '_eig_relpert' num2str(-epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_' fnamesuffix];
            fname3 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_' systemselect '_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_cdiff_' fnamesuffix];
            if is2to1
                fname1 = [fname1 '_2bcl'];
                fname2 = [fname2 '_2bcl'];
                fname3 = [fname3 '_2bcl'];
            end
        else
            fname4 = ['lrddata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_' systemselect '_eigs_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_' fnamesuffix]; % for use with eigs
            if is2to1
                fname4 = [fname4 '_2bcl'];
            end
        end
    end
    
    
    
    if ~useoldfp
        if exist('fpfolder')
            eval(['load ' fpfolder '/1 Y']);
            ivshiftindex = bcl/deltat - stimstart/deltat + 2;
            if ivshiftindex > bcl/deltat
                ivshiftindex = ivshiftindex - bcl/deltat;
            end
            initialvalue = Y(:,ivshiftindex);
            % Find fixed point:
            [sol, it_hist, ierr, x_hist] = nsoli(initialvalue,[modelname '_p2pzero'],[nsolitol nsolitol])
        else
            tstartfp=tic;
            % Find fixed point:
            [sol, it_hist, ierr, x_hist] = nsoli(initialvalue,[modelname '_p2pzero'],[nsolitol nsolitol])
            toc(tstartfp);
        end
    else
        %        eval(['load lrddata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_cdiff_Kp0 sol ierr x_hist it_hist']); % can use this version later depending on how we change the naming convention
        sol = initialvalue;
    end
    
    eval(['fperr = ' modelname '_p2pzero(sol);']); % compute fixed-point error
    maxabsfperr = max(abs(fperr))
    if maxabsfperr > nsolitol
        disp('Fixed point error exceeds tolerance.')
        break
    end
    
    if computejac
        
        tstartjacfw=tic;
        diffjac_mod(sol,[modelname '_p2p'],feval([modelname '_p2p'],sol),epsln); % Solve for empirical bcl-to-bcl Jacobian. 1e-5 (searching by decades) seemed to give the closest answer to the 1-cell analytical result.
        tendjacfw=toc(tstartjacfw)
        load jac
        [v,d] = eig(jac);
        eval(['save ' fname1 ' L bcl d epsilon epsln ierr it_hist jac numpart sol v x_hist useoldfp fperr maxabsfperr stimstart'])
        %save edata_1cell_b230_d005_tauV0p7_tauN170_Vstar4_eig_relpert1e-5_Kp-0p8 *
        
        tstartjacbk=tic;
        diffjac_mod(sol,[modelname '_p2p'],feval([modelname '_p2p'],sol),-epsln); % Solve for empirical bcl-to-bcl Jacobian using perturbation of opposite sign
        tendjacbk=toc(tstartjacbk)
        load jac
        [v,d] = eig(jac);
        eval(['save ' fname2 ' L bcl d epsilon epsln ierr it_hist jac numpart sol v x_hist useoldfp fperr maxabsfperr stimstart'])
        %save edata_1cell_b230_d005_tauV0p7_tauN170_Vstar4_eig_relpert-1e-5_Kp-0p8 *
        jacback=jac
        
        eval(['load ' fname1 '.mat jac'])
        jaccd = (jac+jacback)/2;
        [v,d] = eig(jaccd); % eigenvalues of "central difference" Jacobian
        norm(diag(d))
        eval(['save ' fname3 ' L bcl d epsilon epsln ierr it_hist jac numpart sol v x_hist jacback jaccd useoldfp fperr maxabsfperr stimstart'])
        %save edata_1cell_b230_d005_tauV0p7_tauN170_Vstar4_eig_relpert1e-5_cdiff_Kp-0p8 *
        if ~L
            % compute singular values of controllability and observability
            % matrices for specified B and C for both scaled and unscaled
            % systems
            rankcutoffdefault = max(size(jaccd))*eps(norm(jaccd));% This is Matlab's default rank cutoff. Can't easily control cutoffs in orth /null, so to interpret these results, would need to know default
            for kb = 1:size(B,2)
                P = ctrb(jaccd,B(:,kb));
                Ps = ctrb(Smat*jaccd*Smatinv,Bs(:,kb));
                %                 M = [orth(P) null(P)];
                %                 Abc = inv(M)*jaccd*M;
                %                 if ~sum(sum(isnan(Abc)))
                %                 evc{kb} = eig(Abc(1:size(orth(P),2),1:size(orth(P),2)),'nobalance');
                %                 evnc{kb} = eig(Abc((size(orth(P),2)+1):end,(size(orth(P),2)+1):end),'nobalance');
                %                 end
                svdctrb(:,kb) = svd(P);
                svdctrb_scaled(:,kb) = svd(Ps);
                rankc(kb) = sum(svdctrb(:,kb)>=rankcutoff);
                rankc_scaled(kb) = sum(svdctrb_scaled(:,kb)>=rankcutoff);
                [Abar,Bbar,Cbar,T,k] = ctrbf(jaccd,B(:,kb),C(1,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                %                rankcf(kb) = sum(k>0);
                rankcf(kb) = sum(k);
                evcf{kb} = eig(Abar((numstate-rankcf(kb)+1):end,(numstate-rankcf(kb)+1):end),'nobalance'); % eigenvalues of Ac
                evncf{kb} = eig(Abar(1:(numstate-rankcf(kb)),1:(numstate-rankcf(kb))),'nobalance'); % eigenvalues of Anc
                [Abars,Bbars,Cbars,Ts,ks] = ctrbf(Smat*jaccd*Smatinv,Bs(:,kb),Cs(1,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                %                rankcf_scaled(kb) = sum(ks>0);
                rankcf_scaled(kb) = sum(ks);
                evcf_scaled{kb} = eig(Abars((numstate-rankcf_scaled(kb)+1):end,(numstate-rankcf_scaled(kb)+1):end),'nobalance'); % eigenvalues of Ac
                evncf_scaled{kb} = eig(Abars(1:(numstate-rankcf_scaled(kb)),1:(numstate-rankcf_scaled(kb))),'nobalance'); % eigenvalues of Anc
            end
            for kc = 1:size(C,1)
                Q = obsv(jaccd,C(kc,:));
                Qs = obsv(Smat*jaccd*Smatinv,Cs(kc,:));
                %                 Q = ctrb(jaccd',C(kc,:)'); %I'm using the transposed version since I wasn't sure how to get orth/null into proper row form (transposing them doesn't work)
                %                 O = [orth(Q) null(Q)];
                %                 Abo = inv(O)*jaccd'*O;
                %                 if ~sum(sum(isnan(Abo)))
                %                 evo{kc} = eig(Abo(1:size(orth(Q),2),1:size(orth(Q),2)),'nobalance');
                %                 evno{kc} = eig(Abo((size(orth(Q),2)+1):end,(size(orth(Q),2)+1):end),'nobalance');
                %                 end
                svdobsv(:,kc) = svd(Q);
                svdobsv_scaled(:,kc) = svd(Qs);
                ranko(kc) = sum(svdobsv(:,kc)>=rankcutoff);
                ranko_scaled(kc) = sum(svdobsv_scaled(:,kc)>=rankcutoff);
                [Abar,Bbar,Cbar,T,k] = obsvf(jaccd,B(:,1),C(kc,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                % The rank of the observability matrix is supposed to be sum(k)? The result
                % doesn't seem to agree with that of the standard rank computation
                %                rankof(kc) = sum(k>0);
                rankof(kc) = sum(k);
                evof{kc} = eig(Abar((numstate-rankof(kc)+1):end,(numstate-rankof(kc)+1):end),'nobalance'); % eigenvalues of Ao
                evnof{kc} = eig(Abar(1:(numstate-rankof(kc)),1:(numstate-rankof(kc))),'nobalance'); % eigenvalues of Ano
                [Abars,Bbars,Cbars,Ts,ks] = obsvf(Smat*jaccd*Smatinv,Bs(:,1),Cs(kc,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                %                rankof_scaled(kc) = sum(ks>0);
                rankof_scaled(kc) = sum(ks);
                evof_scaled{kc} = eig(Abars((numstate-rankof_scaled(kc)+1):end,(numstate-rankof_scaled(kc)+1):end),'nobalance'); % eigenvalues of Ao
                evnof_scaled{kc} = eig(Abars(1:(numstate-rankof_scaled(kc)),1:(numstate-rankof_scaled(kc))),'nobalance'); % eigenvalues of Ano
            end
            eval(['save ' fname3 ' rank* svd* ev* B* C* umax Smat tend* -append'])
        end
    else % use eigs
        %    tic
        teigsstart = tic;
        [veigs,deigs,flageigs]=eigs(@(xpert) dirder_mod(sol,xpert,[modelname '_p2p'],feval([modelname '_p2p'],sol),epsln),numstate*numpart,6); % compute eigenvalues
        %    deigs=eigs(@(xpert) dirder_mod(sol,xpert,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln),numstate*numpart,1);
        %    toc % can't use simple tic/toc since there are others embedded in code
        toc(teigsstart)
        
        if ~useoldfp
            eval(['save ' fname4 ' L bcl epsln ierr it_hist numpart sol x_hist deigs veigs flageigs useoldfp maxabsfperr stimstart svd*'])
        else
            eval(['save ' fname4 ' L bcl epsln numpart sol deigs veigs flageigs useoldfp maxabsfperr stimstart svd*'])
        end
        
    end
    
    %L = squeeze(allgains(:,:,i+1));
    %save kseparams L -append % these will be read in by karma_sim_estimtest_p2p.m

    if i < length(allstimstart)
        stimstart=allstimstart(i+1);
        save lrdparams stimstart -append % these will be read in by karma_sim_estimtest_p2p.m
    end
    
if 1%~L
    %svdctrb
    %svdctrb_scaled
    %svdobsv
    %svdobsv_scaled
    %clerr
    
    % for kk=1:numstate
    % maxeigncf(kk) = max(abs(evcf{kk}));
    % maxeignof(kk) = max(abs(evof{kk}));
    % mineigncf(kk) = min(abs(evcf{kk}));
    % mineignof(kk) = min(abs(evof{kk}));
    % end
    
    
    figure
    plot(1:numstate, flipud(sort(abs(eig(jaccd,'nobalance')))),'b-*')
    ylabel('eigenvalue norm')
    xlabel('eigenvalue index')
    %set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
    
    figure
    plot(real(d),imag(d),'k*','markersize',ms,'linewidth',2)
    axis equal
    axis([-.3 1.3 -.2 .2])
    grid
    xlabel('real(\lambda)','fontsize',fs)
    ylabel('imag(\lambda)','fontsize',fs)
    set(gca,'fontsize',fs)
    title(['BCL = ' num2str(bcl) ' ms, start at ' num2str(stimstart) ' ms'],'fontsize',fs)

    figure
    plot(1:numstate,rankc)
    xlabel('state vector element')
    ylabel('controllability matrix rank')
    %axis([0 numstate 0 20])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_ctrb_rank'])
    
    figure
    plot(1:numstate,ranko)
    xlabel('state vector element')
    ylabel('observability matrix rank')
    %axis([0 numstate 0 100])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_obsv_rank'])
    
    figure
    plot(1:numstate,rankcf)
    xlabel('state vector element')
    ylabel('controllability-f matrix rank')
    %axis([0 numstate 0 20])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_ctrbf_rank'])

    figure
    plot(1:numstate,rankof)
    xlabel('state vector element')
    ylabel('observability-f matrix rank')
    %axis([0 numstate 0 100])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_obsvf_rank'])

    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evcf{kk})
            p1 = plot(kk,abs(evcf{kk}),'b*');
        end
        if ~isempty(evncf{kk})
            p2 = plot(kk,abs(evncf{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (ctrbf)')
    xlabel('state vector element')
    legend([p1(1) p2(1)],'controllable','uncontrollable')
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_ctrbf_eig'])
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evof{kk})
            p3 = plot(kk,abs(evof{kk}),'b*');
        end
        if ~isempty(evnof{kk})
            p4 = plot(kk,abs(evnof{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (obsvf)')
    xlabel('state vector element')
    legend([p3(1) p4(1)],'observable','unobservable')
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_obsvf_eig'])
    
    %%%%%%%%%%%%%%%%
    
    figure
    plot(1:numstate,rankc_scaled)
    xlabel('state vector element')
    ylabel('scaled controllability matrix rank')
    %axis([0 numstate 0 20])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_ctrb_rank_scaled_umax_eye'])
    
    figure
    plot(1:numstate,ranko_scaled)
    xlabel('state vector element')
    ylabel('scaled observability matrix rank')
    %axis([0 numstate 0 100])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_obsv_rank_scaled'])
    
    figure
    plot(1:numstate,rankcf_scaled)
    xlabel('state vector element')
    ylabel('scaled controllability-f matrix rank')
    %axis([0 numstate 0 20])
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_ctrbf_rank_scaled_umax_eye'])
    
    figure
    plot(1:numstate,rankof_scaled)
    xlabel('state vector element')
    ylabel('scaled observability-f matrix rank')
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_obsvf_rank_scaled'])
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evcf_scaled{kk})
            p1 = plot(kk,abs(evcf_scaled{kk}),'b*');
        end
        if ~isempty(evncf_scaled{kk})
            p2 = plot(kk,abs(evncf_scaled{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (scaled ctrbf)')
    xlabel('state vector element')
    legend([p1(1) p2(1)],'controllable','uncontrollable')
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_ctrbf_eig_scaled_umax_eye'])
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evof_scaled{kk})
            p3 = plot(kk,abs(evof_scaled{kk}),'b*');
        end
        if ~isempty(evnof_scaled{kk})
            p4 = plot(kk,abs(evnof_scaled{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (scaled obsvf)')
    xlabel('state vector element')
    legend([p3(1) p4(1)],'observable','unobservable')
    saveas(gcf,['b' num2str(bcl) '_start' num2str(round(stimstart/deltat)) '_obsvf_eig_scaled'])

end
end
