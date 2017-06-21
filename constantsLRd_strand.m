% LMM: added simulation time step
data.dt = 0.005; % simulation time step in ms
% LMM: added stimulus flag
data.stimflag = 1; % 0 to apply biphasic stimulus through V
% Nonzero to apply monophasic stimulus through K 

data.F=96485; %% Faraday constant
data.R=8314;  %% Rey
data.Temp=310;  %% Absolute temperature K
l = 0.01;       % Length of the cell (cm)
a = 0.0011;     % Radius of the cell (cm)

vcell = 1000*pi*a*a*l;     %   3.801e-5 uL   % Cell volume (uL)
ageo = 2*pi*a*a+2*pi*a*l;  %   7.671e-5 cm^2    % Geometric membrane area (cm^2)
Acap = ageo*2;             %   1.534e-4 cm^2    % Capacitive membrane area (cm^2)
data.vmyo = vcell*0.68;    % Myoplasm volume (uL)
vmito = vcell*0.24;  % Mitochondria volume (uL)
data.vsr = vcell*0.06;    % SR volume (uL)
data.vnsr = vcell*0.0552;   % NSR volume (uL)
data.vjsr=  vcell*0.0048;   % JSR volume (uL)
data.vss= vcell*0.02;

data.AF=Acap/data.F;

data.frt=data.F/data.Temp/data.R;

   data.k_o=4.5;
%   data.k_o=4.5*(1+1e-7); % pert
   data.na_o=140;
%   data.na_o=140*(1+1e-7); % pert;
   data.ca_o=1.8;
 





% data.st=0; % commented out since it didn't appear to be used anywhere

data.Is=80.0;
%data.Is=0.0;
data.fnsh=0.5;




data.IKsCa_max=0.6;
data.IKsCa_Kd=38e-6;
data.sqrt=(data.k_o/5.4).^(1/2);
data.Qb=1.0;

data.tau=4.75;
data.Qa=data.tau*0.125;%*1.2

data.qn=9;%9
data.tautr=120;
 



%%%% SODIUM COMPARTMENT
data.GNa=16;                          % mS/cm^2
data.GNab=0.004; 







%%%% CALCIUM COMPARTMENT


%% L-type channel
data.gacai=1;         % Activity coefficient of Ca
data.gacao=0.341;     % Activity coefficient of Ca
data.hLca=1;% Hill coefficient of Ca inactivation
data.kmca=6e-4;     % Half-saturation concentration of Ca channel (mM)
data.pca= 5.4e-4;     % Permiability of membrane to Ca (cm/s)
data.gacai=1;         % Activity coefficient of Ca
data.gacao=0.341;     % Activity coefficient of Ca
data.pna=6.75e-7;     % Permiability of membrane to Na (cm/s)
data.ganai=0.75;      % Activity coefficient of Na
data.ganao=0.75;      % Activity coefficient of Na
data.pk=1.93e-7;       % Permiability of membrane to K (cm/s)
data.gaki=0.75;       % Activity coefficient of K
data.gako=0.75;       % Activity coefficient of K


data.gcat = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TauF Default Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data.taufc1= .0197;
%data.taufc2= .02;
%data.taufc3= .0337;
%data.taufc_thresh= 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TauF Adjusted Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.taufc1= .0394; %initial = .0197
data.taufc2= .04;   %initial = .02
data.taufc3= .0674; %initial = .0337
data.taufc_thresh= 10; %initial = 10



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcium in sarcoplasmic
%%%%%%%%%%%%%%%%%%%%%% reticulum

data.hLeak=1;
data.kmup = 0.00092;    % Half-saturation concentration of iup (mM)
data.iupbar = 0.00875/1.0;  % Max. current through iup channel (mM/ms)
data.nsrbar = 15;       % Max. [Ca] in NSR (mM)
data.Jsrbar = 15;       % Max. [Ca] in NSR (mM)
data.ibarpca = 1.15; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
data.kmpca = 0.5e-3; % Half-saturation concentration of sarcolemmal Ca pump (mM)
data.cmdnbar = 0.050;   % Max. [Ca] buffered in CMDN (mM)
data.trpnbar = 0.070;   % Max. [Ca] buffered in TRPN (mM)
data.kmcmdn = 0.00238;  % Equilibrium constant of buffering for CMDN (mM)

 data.kmtrpn = 0.0005;   % Equilibrium constant of buffering for TRPN (mM)
 data.trpnf = 40;   % forward  buffered in TRPN (mM)
 data.trpnb = 0.02;   % backward  TRPN (mM)
  data.cmdnf = 100;   % forward  buffered in TRPN (mM)
 data.cmdnb = 0.238;   % backward  TRPN (mM)

data.csqnbar = 10;      % Max. [Ca] buffered in CSQN (mM)
data.kmcsqn = 0.8;      % Equilibrium constant of buffering for CSQN (mM)
% 
% data.csqnf = 100;  
% data.csqnb = 80; 
data.gcab=0.003016;




 data.taudiff=0.15;
% data.KmCaMK=0.15;
% 
% 
% data.dKmPLBmax=0.00017;
% data.dJupmax=0.75;
% 
% data.iupmax=0.004375; 
% data.Kmup=0.00092;
% data.nsrmax=15.0;
% 
% 
% 

data.BSRmax=0.047; data.KmBSR=0.00087; 
data.BSLmax=1.124; data.KmBSL=0.0087;
% data.slf=115;
% data.slb=1.0;
% 
% data.srf=115;
% data.srb=0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CaMKII%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data.CaMK0=0.05;
% data.KmCam=0.0015;
% data.betaCamK=6.8e-4;
% data.alphaCamK=0.05;
%%%%%%%%%% SodiumCalcium exchanger


data.c1 =0.00025;   % Scaling factor for inaca (uA/uF)
data.c2 = 0.0001;   % Half-saturation concentration of NaCa exhanger (mM)
data.gammas = 0.15;  % Position of energy barrier controlling voltage dependance of inaca






%%%% POTASSIUM COMPARTMENT



data.GKsmax = 0.433;

data.gkrmax =  0.02614;
% data.GK1 = 0.5;% IK1    K 2004  Time-independent potassium current

% 
% data.GKp =0.00276 ; % K 2004

%Calculates Slowly Activating K Current

data.prnak=0.01833;  


data.GK1max=0.75;
data.GKpmax= 0.00552;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sodium-Potassium Pump */



%        inak;    % NaK pump current (uA/uF)
%        fnak;    % Voltage-dependance parameter of inak
%        sigma;   % [Na]o dependance factor of fnak

data.kmnai = 10;    % Half-saturation concentration of NaK pump (mM)
data.kmko = 1.5;    % Half-saturation concentration of NaK pump (mM)
data.ibarnak = 2.25*1.00; % Max. current through Na-K pump (uA/uF)



%data.gitodv = 0.19;

%%%%%%%%%%%%%%%%%%%%%%%%%%  CHLOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% data.Clo=100;
% data.Kmto2=0.1502;
% data.GClb=2.25e-4;
% 
% 
% 
% data.CTKClmax=7.0756e-6;
% data.CTNaClmax=9.8443e-6;
% 
% 
% 
% data.PCl=4e-7;
% 
% 
% 
% 
