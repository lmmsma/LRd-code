function dy=cell_2009_mod(t,y,data)
% LMM: altered to output stimulus as 18th variable. Altered to stimulate more
% than once. Later removed 18th variable, since it makes more sense to have
% it as an input. 
% April 2017: Added data.stimflag switch to choose between stimulating through V or K+ 
 
V=y(1)  ;
H = y(2);m = y(3);  J=y(4); d=y(5);f=y(6); xr=y(7);
ca_T=y(8); na_i=y(9);k_i=y(10);jsr_T=y(11); nsr=y(12);xs=y(13);B=y(14);G=y(15);xs2=y(16);
 Rel=y(17);%Stimul=y(18);

% % Biphasic stimulus
% %In=data.Is*(t>0 & t<=data.fnsh)- data.Is*data.fnsh./(data.bcl-data.fnsh)*( t>=data.fnsh) ;
% trel =  t-data.bcl*floor(t/data.bcl);
% %In=data.Is*(trel>0 & trel<=data.fnsh)- data.Is*data.fnsh./(data.bcl-data.fnsh)*( trel>=data.fnsh) ;
% In=data.Is*(trel>=data.stimstart & trel<=data.fnsh+data.stimstart)- data.Is*data.fnsh./(data.bcl-data.fnsh)*( trel>=data.fnsh+data.stimstart | trel < data.stimstart) ;
In = data.In;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[jsr]=conc_buffer(jsr_T,data.csqnbar,data.kmcsqn, 0, 0);
[ca_i]=conc_buffer(ca_T,data.trpnbar,data.kmtrpn,data.cmdnbar,data.kmcmdn);
%[ca_i]=conc_cai_buffer(ca_T,data);

[ina,inab,am,bm,aH,bH,aj,bj] =comp_ina2005(V,m,H,J,na_i,data);

[ilca,ilcana,ilcak,taud,dss,tauf,fss]=comp_ical2005(V,d,f,ca_i,na_i,k_i,data);
[inaca]=comp_inaca2000(V,ca_i,na_i,data);
[inak]=comp_inak2000(V,na_i,data);

[iks,xss,tauxs]=comp_iks2000(V,xs,xs2,ca_i,na_i,k_i,data);
[icat,bss,gss,taub,taug]=comp_icat2000(V,B,G,ca_i,data);
[ikr,xrss,tauxr]= comp_ikr95(V,xr,k_i,data);
[IK1] = comp_IK1(V,k_i,data);%  
[ikp] = comp_ikp(V,k_i,data); %  

[iserca,ipca,icab,itr]=calcium_2009(V,nsr,jsr,ca_i,data);

caiont = ilca+icab+ipca-2*inaca+icat; %


naiont = ina+inab+3*inaca+ ilcana+3*inak;

if data.stimflag % 0 if through V, otherwise through K+
    kiont = ikr+iks+IK1+ikp+ ilcak-2*inak-In; % If stimulating with K+ 
else
    kiont = ikr+iks+IK1+ikp+ ilcak-2*inak; 
end
%% Derivatives of state variables 

if data.stimflag % 0 if through V, otherwise through K+
    dV=  -(naiont+kiont+caiont);
else
    dV=  -(naiont+kiont+caiont-In); % If stimulating through V
end
 dH=aH*(1-H)-bH*H;
 dm=am*(1-m)-bm*m; 
 dJ=aj*(1-J)-bj*J;



dD=(dss-d)/taud;
df=(fss-f)/tauf;
dxr=(xrss-xr)/tauxr;
dxs=(xss-xs)/tauxs;
dxs2=(xss-xs2)/tauxs/4;

dnai=-naiont*data.AF/(data.vmyo);
dki=-kiont*data.AF/(data.vmyo);
dB=(bss-B)/taub;
dG=(gss-G)/taug;

dOver=0;

Qss=data.Qa/(1+(data.Qb./jsr).^data.qn); 
tauro=data.tau./(1+0.0123./jsr);
dRel=-(ilca.*Qss + Rel)./tauro;

dcai =-caiont*data.AF/(data.vmyo*2)-iserca*data.vnsr/data.vmyo+Rel*data.vjsr/data.vmyo;

dnsr = iserca-itr*data.vjsr./data.vnsr;%;

djsr = itr-Rel;

% RETURN DERIVATIVES
 
% dy = [dV;dH; dm;dJ;dD;df;dxr; dcai;dnai;dki ;djsr;dnsr;dxs;dB;dG;dxs2;dRel;dOver];
% dy = [dV;dH; dm;dJ;dD;df;dxr; dcai;dnai;dki ;djsr;dnsr;dxs;dB;dG;dxs2;dRel;In];
% dy = [dV;dH; dm;dJ;dD;df;dxr; dcai;dnai;dki ;djsr;dnsr;dxs;dB;dG;dxs2;dRel;0];
 dy = [dV;dH; dm;dJ;dD;df;dxr; dcai;dnai;dki ;djsr;dnsr;dxs;dB;dG;dxs2;dRel];


%% L-type calcium channel
function [ilca,ilcana,ilcak,taud,dss,tauf,fss]=comp_ical2005(v,d,f,cai,nai,ki,data)
% Calculates Currents through L-Type Ca Channel

dss=1./(1+exp(-(v+10)/6.24));
taud=dss.*(1-exp(-(v+10)/6.24))./(0.035*(v+10));
  dss1=1./(1+exp(-(v+60)/0.024));
  dss=dss*dss1;
fss=1./(1+exp((v+32)/8))+(0.6)./(1+exp((50-v)/20));

tauf=1./(0.0197*exp(-(0.0337*(v+10))^2)+0.02);
ibarca= data.pca*4*v*data.F*data.frt*((data.gacai*cai*exp(2*v*data.frt)-data.gacao*data.ca_o)/(exp(2*v*data.frt)-1));

ibarna= data.pna*(v*data.F*data.frt).*((data.ganai*nai*exp(v*data.frt)-data.ganao*data.na_o)./(exp(v*data.frt)-1));
ibark= data.pk*(v*data.F*data.frt).*((data.gaki*ki*exp(v*data.frt)-data.gako*data.k_o)./(exp(v*data.frt)-1));

fca =1./(1+(cai./data.kmca)^1);
ilca   = d.*f.*fca.*ibarca;
ilcana = d.*f*fca*ibarna;
ilcak = d.*f*fca*ibark;
%%

%% Fast and background Sodium channels


function [In,inab,am,bm,ah,bh,aj,bj]=comp_ina2005(V,m,H,J,Na_i,data)


ENa =log(data.na_o./Na_i)/data.frt;       % Nernst potential of Na, mV
                    

gNa =data.GNa*m*m*m*H*J;
In = gNa.*(V-ENa);


inab = data.GNab*(V-ENa);

a=1-1./(1+exp(-(V+40)/0.024));
ah= a.*0.135.*exp((80+V)./(-6.8));
bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));

aj =  a.*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V)).*(V+37.78)./(1+exp(0.311*(V+79.23)));

bj= (1-a).*(0.3*exp(-2.535e-7*V)./(1+exp(-0.1*(V+32))))+(a).*(0.1212*exp(-0.01052*V)./(1+exp(-0.1378*(V+40.14))));


am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
bm = 0.08*exp(-V/11);


%% Transient calcium channel


function [icat,bss,gss,taub,taug]=comp_icat2000(V,b,g,cai,data)
%Calculates Currents through T-Type Ca Channel

bss = 1./(1+exp(-(V+14.0)/10.8));
taub = 3.7+6.1/(1+exp((V+25.0)/4.5));
gss = 1./(1+exp((V+60.0)/5.6));

a=1-1./(1+exp(-V/0.0024));
taug = a.*(-0.875*V+12.0)+12.0*(1-a);


ECa = log(data.ca_o/cai)/2/data.frt;

icat = data.gcat*b*b*g*(V-ECa);


%% Time-independent and plato potassium current


function [IK1] = comp_IK1(V,K_i,data)
% IK1    Time-independent potassium current

GK1 = data.GK1max*sqrt(data.k_o/5.4);
EK = log(data.k_o/K_i)/data.frt;
x=(1+exp(0.6987.*(V-EK+11.724)))./(1+exp(0.6168.*(V-EK+4.872))).*0.004;
     K1_inf=1./(1+x);
gK1 = GK1*K1_inf;
IK1 = gK1*(V-EK);

function [ikp] = comp_ikp(V,K_i,data)



EK = log(data.k_o/K_i)/data.frt;

ikp = data.GKpmax*(V-EK)./(1+exp((7.488-V)./5.98)); % plato K 95


%% Slow Activating potassium Current


function [iks,xss,tauxs]=comp_iks2000(v,xs1,xs2,cai,nai,ki,data);

gks = data.GKsmax*(1+0.6/(1+(3.8e-5/cai)^1.4));
eks = log((data.k_o+data.prnak*data.na_o)/(ki+data.prnak*nai))/data.frt;

xss = 1./(1+exp(-(v-1.5)/16.7));

tauxs = 10000./(0.719*(v+30)./(1-exp(-0.148*(v+30)))+1.31*(v+30)./(exp(0.0687*(v+30))-1));


iks = gks*xs1*xs2*(v-eks);
%%
%% Rapidly Activating Potassium Current


function [ikr,xrss,tauxr]= comp_ikr95(v,xr,ki,data)

%Calculates Rapidly Activating K Current


gkr =data.gkrmax*(data.k_o/5.4).^(1/2);
ekr = log(data.k_o/ki)/data.frt;
r = 1/(1+exp((v+9)/22.4));

ikr = gkr*xr*r*(v-ekr);
xrss = 1/(1+exp(-(v+21.5)/7.5));
tauxr = 1/(0.00138*(v+14.2)/(1-exp(-0.123*(v+14.2)))+0.00061*(v+38.9)/(exp(0.145*(v+38.9))-1));
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%  Intracellular Calcium subsystem currents and buffers
function [iserca,ipca,icab,itr]=calcium_2009(v,nsr,jsr,ca_i,data)



 

ipca = (data.ibarpca*ca_i)/(data.kmpca+ca_i);	 % sarcolema pump Ca SERCA
icab =data.gcab*(v- log(data.ca_o/ca_i)/2/data.frt); % background Ca

iserca =data.iupbar.*ca_i/(ca_i+data.kmup)- data.iupbar/ data.nsrbar*nsr;

itr = (nsr-jsr)./data.tautr;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculates Na-Ca Exchanger Current

function [Inaca]=comp_inaca2000(v,cai,nai,data)
%Calculates Na-Ca Exchanger Current


%    inaca;               % NaCa exchanger current (uA/uF)
Inaca = data.c1*exp(( data.gammas-1)*v* data.frt).*((exp(v* data.frt).*nai.^3*data.ca_o- data.na_o^3*cai)./...
    (1+ data.c2*exp(( data.gammas-1)*v* data.frt).*(exp(v* data.frt).*nai.^3*data.ca_o+ data.na_o^3*cai)));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sodium-Potassium Pump
function [inak]=comp_inak2000(v,nai,data)
% Sodium-Potassium Pump */



sigma = (exp(data.na_o/67.3)-1)/7;

fnak = 1/(1+0.1245*exp((-0.1*v*data.frt)) + 0.0365*sigma*exp(-v*data.frt));

inak = data.ibarnak*fnak./(1+(data.kmnai./nai).^2)./(1+data.kmko./data.k_o);

%%

function [cai]=conc_buffer(ca_t,a1,b1,a2,b2)
%% Online supplement Biophys 2009 Livshitz and Rudy

 	alp2 = a1+a2+b1+b2-ca_t;
	alp1 = b1*b2 -ca_t.*(b1+b2)+a1*b2+a2*b1;
	alp0 = -b1*b2*ca_t;
      Q=(3*alp1-alp2.^2)/9;
      R=(9*alp2.*alp1-27*alp0-2*alp2.^3)/54;
      T=(R + (Q.^3 + R.^2).^0.5).^(1/3) -Q./((R + (Q.^3 + R.^2).^0.5).^(1/3));
       
     cai =abs(T-alp2/3);
