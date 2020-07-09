% LMM: made various changes to main_LRd_1cell: only run for one cycle length (BCL), load lrdparams, and derive BCL and other parameters from this file, to be
% more consistent with Karma model eigenvalue solver. In addition, the
% stimulus current (data.In) formula has been altered to make sure the
% current is conservative over one cycle.
% April 2017: Took version of code I initially sent to Kalyan, and
% updated to work in single-cycle (for fixed-point search) and
% multi-cycle (for pacedown) mode. Input order made to match Claire's code
% where possible.
% Function  inputs: 
% yin = initial condition
% bcl = basic cycle length in ms
% ncyc = number of cycles (specifies duration of simulation) 
% subdiv_per_cyc = record this number of intervals per cycle. This quantity 
% must be a positive integer. 
% Output: Changed yend to last simualted step, not last sampled step. 
% June 2017: Modified to read in any missing arguments after yin (bcl, ncyc,
% subdiv_per_cyc) from a file, if the file is available, before reverting
% to default values. Changed 'lrdparams' filename to 'lrdinputs', since the 
% former was used for different purposes than are needed here. Folder for
% saving data is now defined at the beginning.

%function [tmp,ca]= main();
% /***************************************************************************
%  *   Copyright (C) 2006 by Leonid Livshitz and Yoram Rudy  *
%  *   Email livshitz@wustl.edu   *
%  *                                                                         *
%  *   This program is free software; you can redistribute it and/or modify  *
%  *   it under the terms of the GNU General Public License as published by  *
%  *   the Free Software Foundation; either version 2 of the License, or     *
%  *   (at your option) any later version.                                   *
%  *                                                                         *
%  *   This program is distributed in the hope that it will be useful,       *
%  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%  *   GNU General Public License for more details.                          *
%  *                                                                         *
%  *   You should have received a copy of the GNU General Public License     *
%  *   along with this program; if not, write to the                         *
%  *   Free Software Foundation, Inc.,                                       *
%  *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%  ***************************************************************************/
%

function yend = main_LRd_1cell(yin, bcl, ncyc, subdiv_per_cyc)
%clear variables
%close all

constantsLRd_strand %% constants

%%data.Ist=120;% stimuli
%data.Ist=0;% stimuli; doesn't appear to do anything? At least not when pdex1pde isn't called

% The next two parameters are only needed when simulating 1D fiber
%data.cell_STIM=0.05;% fraction of the cell stimulated
%data.Diff_X=0.0007;%% Effective diffusion coefficient

folder = 'lrddata/'; % simulation output and settings are stored here 

if nargin > 4 % too many arguments
    disp('too many input arguments')
    return
end

if nargin < 4 % number of subdivisions per cycle has not been specified
    if exist([folder 'lrdinputs.mat'],'file') % this file is being used to bypass single-input restriction of nsoli
        eval(['load ' folder 'lrdinputs subdiv_per_cyc']) %LMM;
    else
        subdiv_per_cyc = bcl/data.dt; % for recordings, subdivide cycle into this many subintervals
        disp(['Using default number of sampling intervals, ' num2str(subdiv_per_cyc) ' per cycle'])
    end
    if nargin < 3 % number of cycles has not been specified
        if exist([folder 'lrdinputs.mat'],'file') % this file is being used to bypass single-input restriction of nsoli
            eval(['load ' folder 'lrdinputs ncyc']) %LMM;
        else
            ncyc = 1;
            disp(['Using default number of cycle lengths, ' num2str(ncyc)])
        end
        if nargin < 2 % bcl has not been specified
           if exist([folder 'lrdinputs.mat'],'file') % this file is being used to bypass single-input restriction of nsoli
                eval(['load ' folder 'lrdinputs bcl']) %LMM;
            else
                bcl = 1000;
                disp(['Using default BCL = ' num2str(bcl) ' ms'])
            end
        end
    end
end

N = ncyc; % number of cycles
Time = bcl; 
%    outstep = 10;%deltat;% %ms, must be a factor of Time
outstep = round((Time/subdiv_per_cyc)/data.dt)*data.dt; %ms, length of sampling interval
% adjust outstep to be an integer multiple of data.dt
yfptrajflag = 0; % Default is no fixed-point trajectory is loaded

%Tmesh=400; %% ms
data.bcl=Time;
if exist('lrdparams','file') % this file is generated by lrd_batch_eig.m
    load lrdparams stimstart %LMM; this file contains variables: bcl stimstart L numpart
else
    stimstart = data.dt;
end
% Add a new block to work with code from post-REU project: 
if exist([folder 'lrdinputs.mat'],'file')
    matobj = matfile([folder 'lrdinputs.mat']);
    if isprop(matobj,'stimstart')%~isempty(matobj.stimstart)
         eval(['load ' folder 'lrdinputs stimstart'])
    else
         stimstart = data.dt;
    end
    if isprop(matobj,'measurementindex')%~isempty(matobj.stimstart)
        eval(['load ' folder 'lrdinputs measurementindex'])
        feedbackflag = measurementindex;
    else
        feedbackflag = 0; 
    end
    if isprop(matobj,'fp_for_fbk')%~isempty(matobj.stimstart)
        eval(['load ' folder 'lrdinputs fp_for_fbk'])
         %this will only be used if feedbackflag ~= 0 
    end
    if isprop(matobj,'Llu')%~isempty(matobj.stimstart)
        eval(['load ' folder 'lrdinputs Llu'])
         %this will only be used if feedbackflag ~= 0 
    end
    if isprop(matobj,'Yfp_traj')%~isempty(matobj.stimstart)
        eval(['load ' folder 'lrdinputs Yfp_traj'])
        yfptrajflag = 1; 
         %this will only be used if feedbackflag ~= 0 
    end
end
data.stimstart = stimstart; % ms; wait this long before starting positive pulse (prior to this impose negative part of biphasic pulse)
%data.stimstart = 0; % ms; wait this long before starting positive pulse (prior to this impose negative part of biphasic pulse)
%Xsize=1;Xmesh=50;%% 1 cm=100 cells
%Xsize=0.02;Xmesh=1;%% 1 cm=100 cells

%x = linspace(0,Xsize,Xmesh);
%t = linspace(0,Time,Tmesh);
%deltat = .01; % ms % too long? Why?
%deltat = .005; % ms
deltat = data.dt; % simulation time step
finaltime = N*Time + deltat; %ms
numstep = round(finaltime/deltat);
%feedbackflag = 0; % zero to run model open-loop, nonzero to add observer feedback
% The current logic should be that the "measurement" time is the start of
% the interval, regardless of when the stimulus occurs relative to the
% start of the interval. 
if feedbackflag
    % Here feedback will be applied to the unscaled (dimensional) model, so
    % C matrices and gains should be chosen accordingly
    feedbackduration = data.dt;%0;%2*data.dt; %; %0.5; % ms; how long to apply observer feedback 
%    feedbacknumstep = round(feedbackduration/deltat); % number of timesteps over which to apply feedback
%    obsvgain = zeros(length(yin),1); % this should be the dimensional (with units) observer gain
    obsvgain = Llu';
%    obsvgain(1) = 0; 
%    obsvgain(10) = obsvgain(10)/2; 
%    obsvgain(10) = obsvgain(10)*1.5; 
    Cmatrix = zeros(1,length(yin)); % output matrix
    Cmatrix(measurementindex) = 1; 
    yfp = fp_for_fbk; % fixed point for this offset
end

fname = [folder 'lrddata_1cell_b' num2str(bcl)]; % save simulation data in this file
%save data data % no longer needed, since "data" structure is saved in file specified
%above

%  xs = linspace(0,Xsize,10*Xmesh);
%N=300;%% number of beats
%N=120;%% number of beats
%tp=linspace(0,Time*N,Tmesh*N);

% [Xp,Tp]=meshgrid(tp,x); % why is this backwards? Ans: Xp and Tp are only used for plotting and the error is corrected through thddddddde axis labels
%t0 = clock;


%TMP=zeros(length(t),length(x),N);
%Cai=zeros(length(t),length(x),N);
% Initial conditions
%   load llrd_n300
%
%  X0=ones(1,length(x))'*x0;
%  data.cs = spline(x,X0');% Initial conditions
% Initial conditions
%load  initial_alternans
%data.cs=cs_0;
%   New=ppval(xs,data.cs); % doesn't appear to be used for anything



%yall = []; tall = [];
if 0
    yall = zeros(round(N*Time/outstep),18);
    tall = zeros(round(N*Time/outstep),1);
    %yinit = ppval(cs_0,x(1));
    %load y15sb1000f4000
    %load y15sb1000f30000
    %load y15sb1000f300000
    %load y15sb1000f600000
    %yinit = ylast;
    yinit = yin;
    %    opt=odeset('MaxStep',0.25);
    %    opt=odeset('MaxStep',1);
    opt=[];
    tic
    for i=1:N
        %    [tout,yout]=ode15s(@cell_2009_mod2,[(i-1) i]*Time,yinit);
        [tout,yout]=ode15s(@cell_2009_mod2,[((i-1)*Time):outstep:(i*Time)],yinit);
        %    [tout,yout]=ode15s(@cell_2009_mod2,[((i-1)*Time):outstep:(i*Time)],yinit,opt); % Need to specify MaxStep to be sufficiently small to "catch" the stimulus pulse
        yinit=yout(end,:); % need to use data.cs instead if using PDE code
        %    yall = [yall;yout];
        %    tall = [tall;tout];
        %    yall = [yall;yout(1:end-1,:)];
        %    tall = [tall;tout(1:end-1,:)];
        trange = tout/outstep+1;
        yall(trange(1:end-1),:) = yout(1:end-1,:);
        tall(trange(1:end-1)) = tout(1:end-1);
    end
    toc
    yend = yall(end,:);
    %
end

if 1
    %finaltime = data.bcl+deltat;
    %finaltime = 800; %ms
    
%    Y=zeros(17,N*Time/outstep+1);
%    iext = zeros(1,N*Time/outstep+1);
    Y=zeros(17,N*subdiv_per_cyc+1);
    iext = zeros(1,N*subdiv_per_cyc+1);
    %iext = zeros(1,numstep);
    %Y(:,1) = [4.0, 0.8, 0.01, 0.2, 1.0, 2.0, 5.0e-5, 1.0, 140.0, 8.0, -87.0, 0.005, 1.0, 1.0, 0.01, 1.0];
    %Y(:,1) = ppval(cs_0,x(1));
    Y(:,1) = yin;
    yprev = Y(:,1);
    Yprevcyc = yprev; % use for feedback based on Y one cycle ago
    ynext = zeros(size(yprev));
    chargectr = 0;
    cyclectr = 1; 
    tic
    for k=1:numstep-1
                if ~mod(k*deltat,Time)
                   disp(['Cycle number = ' num2str(cyclectr)])
                   cyclectr = cyclectr + 1; 
                end
        %%%    Y(:,k+1) = Y(:,k) + deltat*dcl_dn_1993_extstim(bcl,finaltime,k*deltat,Y(:,k));
        %%%    temp = cell_2009_mod2(k*deltat,Y(:,k));
        %
        %    dydt = cell_2009_mod(k*deltat,Y(:,k),data);
        %    Y((1:end-1),k+1) = Y((1:end-1),k) + deltat*temp(1:end-1);
        %    Y(end,k+1) = temp(end);    % stimulus current
        
        % Biphasic stimulus
        %    trel =  k*deltat-data.bcl*floor(k*deltat/data.bcl);
        trel =  k*deltat-data.bcl*floor((k-1)*deltat/data.bcl);
        %%%%    data.In=data.Is*(trel>=data.stimstart & trel<=data.fnsh+data.stimstart)- data.Is*data.fnsh./(data.bcl-data.fnsh)*( trel>=data.fnsh+data.stimstart | trel < data.stimstart) ;
        %%%  data.In=data.Is*(trel>=data.stimstart & trel<=data.fnsh+data.stimstart)- data.Is*data.fnsh./(data.bcl-data.fnsh)*( trel>data.fnsh+data.stimstart | trel < data.stimstart) ;   % If using biphasic stimulus (only works for stimstart = 0)
        if data.stimflag % 0 if through V, otherwise through K+
            data.In=data.Is*(trel>=data.stimstart & trel<=data.fnsh+data.stimstart); % If using monophasic stimulus
        else
            data.In=data.Is*(trel>=data.stimstart & trel<data.fnsh+data.stimstart)- data.Is*data.fnsh./(data.bcl-data.fnsh)*( trel>=data.fnsh+data.stimstart | trel < data.stimstart) ;   % If using biphasic stimulus
        end
        
        dydt = cell_2009_mod(k*deltat,yprev,data);
        ynext = yprev + deltat*dydt;
        
        % Experiment with relaxation term to hasten convergence of Na and K
%         nakdirs = zeros(17,17);
%         nakdirs(9,9) = 1; 
%         nakdirs(10,10) = 1; 
%         ynext = yprev + deltat*dydt -deltat*nakdirs*(yprev-Yfp_traj(:,k)); 

        if feedbackflag 
            % The "beginning" of the simulation is the measurement time for
            % observer feedback tests, so the feedback pulse should occur
            % then (and periodically afterward), which is not necessarily
            % the same as the stimulus time interval.
%%%            if trel>=data.dt && trel<=feedbackduration+data.dt %k <= feedbacknumstep %< or <=? % apply observer feedback temporarily
  %          if trel>=data.dt && trel<feedbackduration+data.dt %k <= feedbacknumstep %< or <=? % apply observer feedback temporarily
%%%            if trel>=data.stimstart && trel<data.stimstart+deltat
                % The interpretation here is that ynext is the estimated
                % state, whereas yfp is the "measurement". Granted, it's a
                % constant, noiseless measurement. More challenging tests
                % would involve replacing yfp with points sampled once per bcl from a
                % different, perhaps non-period-1 trajectory.
                %ynext = yprev + deltat*dydt - obsvgain*Cmatrix*(yprev-yfp);
%%%             elseif k == numstep - 1 
%%%                 ynext = yprev + deltat*dydt -0.5*obsvgain*Cmatrix*(yprev-yfp);
%            end
                if ~mod(k*deltat,Time)
                    if yfptrajflag % entire reference traj is available
%                        ynext = yprev + deltat*dydt - obsvgain*Cmatrix*(Yprevcyc-Yfp_traj(:,k-round(Time/deltat)+1))    -deltat*nakdirs*(yprev-Yfp_traj(:,k));
%                        ynext = yprev + deltat*dydt - obsvgain*Cmatrix*(Yprevcyc-Yfp_traj(:,k-round(Time/deltat)+1));
                        tempfbterm = -obsvgain*Cmatrix*(Yprevcyc-Yfp_traj(:,k-round(Time/deltat)+1));
%                         if norm(tempfbterm) > 1 %saturate feedback term
%                             fbterm = 1*tempfbterm/norm(tempfbterm);
%                         end
                         fbterm = tempfbterm; 
%                        ynext = yprev + deltat*dydt +
%                        fbterm-deltat*nakdirs*(yprev-Yfp_traj(:,k)); %saturation and relaxation 
                        ynext = yprev + deltat*dydt + fbterm;
                    else % only fixed point vector is available
                        ynext = yprev + deltat*dydt - obsvgain*Cmatrix*(Yprevcyc-yfp);
                    end
                    Yprevcyc = ynext; % or yprev? use for feedback based on Y one cycle ago
                end

        end
        yprev = ynext;
        %    iext(k) = data.In;
        if ~mod(k*deltat,outstep)
%            iext(1,k*(deltat/outstep)) = data.In;
%            Y(:,k*(deltat/outstep) + 1) = ynext;
            iext(1,round(k*(deltat/outstep))) = data.In;
            Y(:,round(k*(deltat/outstep)) + 1) = ynext;
        end
        chargectr = chargectr + data.In;
    end
    toc
    % LMM, 4/14/17: shouldn't this be the last simulated step, not the last sampled step?
%    yend = Y(:,end); % last sampled step
    yend = ynext; % last simulated step
    %save lrddata_1cell_b1000_solem12_start1_Kp0/1 Y
end

eval(['save ' fname ' data yin bcl ncyc subdiv_per_cyc deltat finaltime numstep outstep Y iext chargectr']);

if 0
    
    ca=[]; tmp=[];
    h0 = waitbar(0,' Matlab is working hard, Please wait about 5 minutes');
    
    
    for i=1:N
        
        sol=pd_LRD_strand(t,x,data);
        y=sol(end,:,:);
        y1=reshape(y,Xmesh,18);
        data.cs = spline(x,y1');
        t1= etime(clock,t0)/60
        
        TMP(:,:,i) = sol(:,:,1);
        
        %Cai(:,:,i)=sol(:,:,17);
        
        
        
        waitbar(i/N,h0)
        Cai(:,:,i)=conc_buffer(sol(:,:,8),data.trpnbar,data.kmtrpn,data.cmdnbar,data.kmcmdn);
        
        ca=[ca ; Cai(:,:,i)];
        
        tmp=[tmp ;TMP(:,:,i)];
        
        
        
    end
    close(h0)
    cs_0=data.cs;
    %  save  init_alternans cs_0
    % save strand_web3  sol Tp Xp
    
    
    figure
    h=mesh(Tp',Xp',tmp);
    ylabel(['Time [msec]'],'FontSize',18)
    xlabel(['Fiber'],'FontSize',18)
    view(-37,70)
    zlabel('tmp','FontSize',18)
    set(gca,'FontSize',18)
    axis tight
    box off
    grid off
    
    ylabel(['Time [msec]'],'FontSize',18)
    xlabel(['Fiber'],'FontSize',18)
    view(-37,70)
    zlabel('V_m','FontSize',18)
    
    figure
    h1=mesh(Tp',Xp',1000*ca);
    ylabel(['Time [msec]'],'FontSize',18)
    xlabel(['Fiber'],'FontSize',18)
    
    zlabel('[Ca_i] \mu M','FontSize',18)
    set(gca,'FontSize',18)
    view(-37,70)
    axis tight
    %colorbar
    box off
    grid off
    
    % function [cai]=conc_buffer(ca_t,a1,b1,a2,b2)
    %
    %
    %  	alp2 = a1+a2+b1+b2-ca_t;
    % 	alp1 = b1*b2 -ca_t.*(b1+b2)+a1*b2+a2*b1;
    % 	alp0 = -b1*b2*ca_t;
    %       Q=(3*alp1-alp2.^2)/9;
    %       R=(9*alp2.*alp1-27*alp0-2*alp2.^3)/54;
    %       T=(R + (Q.^3 + R.^2).^0.5).^(1/3) -Q./((R + (Q.^3 + R.^2).^0.5).^(1/3));
    %
    %      cai =abs(T-alp2/3);
end