function [Vstn, spkhststn, Vgpe, spkhstgpe,Vgpi,spkhstgpi,I_strD1] = stn_gpev6(Vstn, Vgpe,Istn, Igpe, niter,wstr,Vgpi,Igpi,wsgpi,DA,time1,time2,dt)

%Izhikevich Parameters
astn=0.005; bstn=0.265;cstn=-65;dstn=1.5;%'(C) tonic bursting'
agpe=0.1;bgpe=0.2;cgpe=-65;dgpe=2;
agpi=0.1;bgpi=0.2;cgpi=-65;dgpi=2;

% Izhikevich variable 'U' intialization 

Ustn = zeros(size(Vstn));
Ugpe = zeros(size(Vgpe));
Ugpi = zeros(size(Vgpe));

% dt in the 'V' ODE 
taustn = 0.1;
taugpe = 0.05;
taugpi=0.1;

% Membrane capacitances
Cstn = 1;
Cgpe = 1;

% Peak voltage
vpeak = 40;

[m,n] = size(Vstn);
% stn ,gpe initilization
stn_zeros = zeros(m,n);
gpe_zeros = zeros(m,n);
 
% spk variable intialization
stncol_zeros = zeros(m*n/2,9);
% gpecol_zeros = zeros(m*n/2,1);
stncurr_spk = stn_zeros;
gpicurr_spk=stn_zeros;
gpecurr_spk = gpe_zeros;
spkhststn1 = stncol_zeros;
spkhststn2 = stncol_zeros;
spkhstgpe1 = stncol_zeros;
spkhstgpe2 = stncol_zeros;
spkhstgpi1 = stncol_zeros;
spkhstgpi2 = stncol_zeros;



% psp variable initilization

h_strD1=zeros(50,50);
h_strD2=zeros(50,50);
h_strD11=zeros(25,50);
h_strD12=h_strD11;
h_strD21=zeros(25,50);
h_strD22=h_strD21;
h_gs = stn_zeros;
h_nmdastn=stn_zeros;
h_ampastn=stn_zeros;
h_nmdagpe=gpe_zeros;
h_ampagpe=gpe_zeros;
hlat_gaba_gpe=gpe_zeros;
hlat_gaba_gpi=gpe_zeros;
h_nmdagpi=gpe_zeros;
h_ampagpi=gpe_zeros;

% decay constants(ms)
taunmda=160; 
tauampa=6;
taugaba=4;
taunmdagpi=60;

% dt/T in PSP
lam_nmda = dt/taunmda; 
lam_ampa = dt/tauampa;
lam_gaba= dt/taugaba;
lam_nmdagpi=dt/taunmdagpi;

% RMP of receptors
Egaba = -60;
Enmda = 0;
Eampa = 0;

mg0=1; % magnesium conc.

% weight cahnging parameters
CD1=1;
CD2=1;

%  Cd1str= 1/(exp(-DA))
% Cd1str= 1/(exp(-DA/0.9))-1;
% Cd1str=2.5*tanh(2.5*DA);

Cd1str=10./(1+exp(-7.5*(DA-1)));
Cd2str=7./(1+exp(7.5*DA));
% Calculation of DBS current
amp=1000;
dur_mono=0.1;
dur_bi=0.2;
f=130;
% Itemp=zeros(1,niter);
Idbs=zeros(50,50);
%     Itemp=monophasic(amp,dur_mono,niter,f); % 
  Itemp=biphasic(amp,dur_bi,niter,f); % 
 
 nlat=50;
 radius=25;
for j = 1:niter
     Idbs=zeros(50,50);
%   Idbs=onepoint(Idbs,Itemp(j),nlat,radius,j);
%--------------------------------------striatum current-------------------------------%  
 if j>=time1 && j<=time2 
        fstr1=4;
        fstr2=8;
        D2=DA; %  current dopamine level
        D1=D2;
        [wsg,wgs,wlatstn,wlatgpe]= weightcal(D2);
        
   

    spkD11 = spkgen(j, 1, fstr1, 0); % random spike input for D1 striatum
    spkD12 = spkgen(j, 1, fstr2, 0);
       
    SpkD11=repmat(spkD11,n/2,n);
    SpkD12=repmat(spkD12,n/2,n);
    
    spkD21 = spkgen(j, 1, fstr1, 0); % random spike input for D2 striatum
    spkD22 = spkgen(j, 1, fstr2, 0);
    
    SpkD21=repmat(spkD21,n/2,n);
    SpkD22=repmat(spkD22,n/2,n);
       
    % psp variable D1 striatum 
    h_strD11= (1-lam_gaba).* h_strD11 + lam_gaba.*SpkD11;
    h_strD12= (1-lam_gaba).* h_strD12 + lam_gaba.*SpkD12;
    h_strD1=[h_strD11; h_strD12];
    % D1 striatum current
    I_strD1=wstr(2).*h_strD1.*(Egaba-Vgpi);
    
%      subplot(2,2,1)
%     imagesc(SpkD11)
%     subplot(2,2,2)
%     imagesc(SpkD12)
%     subplot(2,2,3)
%     imagesc(SpkD21)
%     subplot(2,2,4)
%     imagesc(SpkD22)
    
   % psp variable D2 striatum
    h_strD21= (1-lam_gaba).* h_strD21 + lam_gaba.*SpkD21;
    h_strD22= (1-lam_gaba).* h_strD22 + lam_gaba.*SpkD22;
    h_strD2=[h_strD21; h_strD22];
    % D2 striatum current
    I_strD2=wstr(1).*h_strD2.*(Egaba-Vgpe);
        
 else
        D2=0.1; % low dopamine level
        D1=D2;
        fstr1=1;
        fstr2=1;
        [wsg,wgs,wlatstn,wlatgpe]= weightcal(D2);
     
   spkD11 = spkgen(j, 1, fstr1, 0); % random spike input for D1 striatum
    spkD12 = spkgen(j, 1, fstr2, 0);
       
    SpkD11=repmat(spkD11,n/2,n);
    SpkD12=repmat(spkD12,n/2,n);
    
    spkD21 = spkgen(j, 1, fstr1, 0); % random spike input for D2 striatum
    spkD22 = spkgen(j, 1, fstr2, 0);
    
    SpkD21=repmat(spkD21,n/2,n);
    SpkD22=repmat(spkD22,n/2,n);
%     
%     subplot(2,2,1)
%     imagesc(SpkD11)
%     subplot(2,2,2)
%     imagesc(SpkD12)
%     subplot(2,2,3)
%     imagesc(SpkD21)
%     subplot(2,2,4)
%     imagesc(SpkD22)

    % psp variable D1 striatum 
    h_strD11= (1-lam_gaba).* h_strD11 + lam_gaba.*SpkD11;
    h_strD12= (1-lam_gaba).* h_strD12 + lam_gaba.*SpkD12;
    h_strD1=[h_strD11; h_strD12];
    % D1 striatum current
    I_strD1=wstr(2).*h_strD1.*(Egaba-Vgpi);
    
   % psp variable D2 striatum
    h_strD21= (1-lam_gaba).* h_strD21 + lam_gaba.*SpkD21;
    h_strD22= (1-lam_gaba).* h_strD22 + lam_gaba.*SpkD22;
    h_strD2=[h_strD21; h_strD22];
    % D2 striatum current
    I_strD2=wstr(1).*h_strD2.*(Egaba-Vgpe);

 end
    
%  -------------------------------------GPE-------------------------------------------%
    
    %--------------------- Input from stn to gpe----------------------------------%
    % psp variable
    h_nmdagpe = (1-lam_nmda).* h_nmdagpe + lam_nmda.*stncurr_spk;
    h_ampagpe = (1-lam_ampa).* h_ampagpe + lam_ampa.*stncurr_spk;
    
    % nmda and ampa currents
    Inmdagpe = wsg.*h_nmdagpe.*(Enmda - Vgpe);
    Iampagpe = wsg.*h_ampagpe.*(Eampa - Vgpe);
    %-------------------------- Input from gpe to gpe(laterals)-------------------%
    % psp variable
    hlat_gaba_gpe = (1-lam_gaba).* hlat_gaba_gpe + lam_gaba.*gpecurr_spk;
    tmplat_gaba_gpe = hlat_gaba_gpe.*(Egaba - Vgpe);
    
    % lateral currents
    Ilatgpe = conv2(tmplat_gaba_gpe, wlatgpe, 'same');
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vgpe));
    
    % total curent recieved
%     Itmpgpe = Ilatgpe + B.*Inmdagpe+Iampagpe + Igpe+(1-CD2*D2).*I_strD2; % linear D2 function
    Itmpgpe = Ilatgpe + B.*Inmdagpe+Iampagpe + Igpe+(Cd2str.*I_strD2);
    
    % V ,U updated
    dvgpe = taugpe.*((0.04.*Vgpe.*Vgpe)+5.*Vgpe+140 - Ugpe + Itmpgpe)./Cgpe;
    dugpe = taugpe.*agpe.*(bgpe.*(Vgpe) - Ugpe);
    Vgpe_nxt = Vgpe + dvgpe;
    Ugpe_nxt = Ugpe + dugpe;
    
    % replacing criteria
    indg = find(Vgpe_nxt > vpeak);
    Vgpe_nxt(indg) = cgpe.*ones(size(indg));
    Ugpe_nxt(indg) = Ugpe(indg) + dgpe.*ones(size(indg));
    gpecurr_spk = gpe_zeros;
    gpecurr_spk(indg) = ones(size(indg));
    
    tmpg1 = reshape(gpecurr_spk(1:m/2,:), m/2*n,1);
    spkhstgpe1= [spkhstgpe1 tmpg1];
    
    tmpg2 = reshape(gpecurr_spk(m/2+1:m,:), m/2*n,1);
    spkhstgpe2= [spkhstgpe2 tmpg2];
    
    spkhstgpe=[spkhstgpe1 ; spkhstgpe2];
    
    Vgpe = Vgpe_nxt;
    Ugpe = Ugpe_nxt;
%---------------------------------------------------STN----------------------------------------%
    %----------------------------------------Input from gpe to stn---------------------%
    % psp variable
    h_gs = (1-lam_gaba).* h_gs + lam_gaba.*gpecurr_spk;% input from gpe to nmda stn
    % gaba current
    I_gs = wgs.*h_gs.*(Egaba - Vstn);

    %---------------------------------------Input from stn to stn(laterals)--------------%
    % psp variable
    h_nmdastn = (1-lam_nmda).* h_nmdastn + lam_nmda.*stncurr_spk; %psp nmda stn lat
    h_ampastn = (1-lam_ampa).* h_ampastn + lam_ampa.*stncurr_spk; %psp ampa stn lat
    
    tmplat_nmda_stn = h_nmdastn.*(Enmda - Vstn);
    tmplat_ampa_stn = h_ampastn.*(Eampa - Vstn);
    
    
    Ilat_nmda_stn = conv2(tmplat_nmda_stn, wlatstn, 'same'); % lat nmda stn
    Ilat_ampa_stn = conv2(tmplat_ampa_stn, wlatstn, 'same');% lat ampa stn
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vstn));
    
    % total current stn recieves
    Itmpstn =  B.*Ilat_nmda_stn + Ilat_ampa_stn + I_gs + Istn+Idbs; % total currents
%     Itmpstn =  I_gs + Istn; % w/o laterals total currents
 
    % V,U updated
    dvstn = taustn.*((0.04.*Vstn.*Vstn)+5.*Vstn - Ustn + Itmpstn +140)./Cstn;
    dustn = taustn.*astn.*(bstn.*(Vstn) - Ustn);
    Vstn_nxt = Vstn + dvstn;
    Ustn_nxt = Ustn + dustn;
    
    inds = find(Vstn_nxt > vpeak);
    
    Vstn_nxt(inds) = cstn.*ones(size(inds));
    Ustn_nxt(inds) = Ustn(inds) + dstn.*ones(size(inds));
    stncurr_spk = stn_zeros;
    stncurr_spk(inds) = ones(size(inds));
    
    tmps1 = reshape(stncurr_spk(1:m/2,:), m/2*n,1);
    spkhststn1= [spkhststn1 tmps1];
 
    tmps2 = reshape(stncurr_spk(m/2+1:m,:), m/2*n,1);
    spkhststn2= [spkhststn2 tmps2];
    
    spkhststn=[spkhststn1 ; spkhststn2];
     
    Vstn = Vstn_nxt;
    Ustn = Ustn_nxt;

%------------------------------------GPI--------------------------------------%
    %--------------------- Input from stn to gpi----------------------------------%
    % psp variable 
    h_nmdagpi = (1-lam_nmdagpi).* h_nmdagpi + lam_nmdagpi.*stncurr_spk;
    h_ampagpi = (1-lam_ampa).* h_ampagpi + lam_ampa.*stncurr_spk;
   %nmda and ampa currents   
    Inmdagpi = wsgpi.*h_nmdagpi.*(Enmda - Vgpi);
    Iampagpi = wsgpi.*h_ampagpi.*(Eampa - Vgpi);
    
    % %-------------------------- inpur from gpi to gpi laterals---------------------%
    %psp variable
    hlat_gaba_gpi = (1-lam_gaba).* hlat_gaba_gpi + lam_gaba.*gpicurr_spk;
    tmplat_gaba_gpi = hlat_gaba_gpi.*(Egaba - Vgpi);
    
    Ilatgpi = conv2(tmplat_gaba_gpi, wlatgpe, 'same');
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vgpi));
    
  % total current gpi recieves  
%  Itmpgpi =  B.*Inmdagpi+Iampagpi+Igpi+Ilatgpi+I_strD1.*Cd1str;  %with laterals
 Itmpgpi =  B.*Inmdagpi+Iampagpi+Igpi+I_strD1.*Cd1str; % w/o laterals
%  Itmpgpi =  Igpi+I_strD1.*Cd1str; % w/o laterals and IP
%  Itmpgpi =  Igpi+Ilatgpi+I_strD1.*Cd1str;  %with laterals and no IP

   % V,U updated 
    dvgpi = taugpi.*((0.04.*Vgpi.*Vgpi)+5.*Vgpi+140 - Ugpi + Itmpgpi)./Cgpe;
    dugpi = taugpi.*agpi.*(bgpi.*(Vgpi) - Ugpi);
    Vgpi_nxt = Vgpi + dvgpi;
    Ugpi_nxt = Ugpi + dugpi;
    
    
    indgi = find(Vgpi_nxt > vpeak);
    
    Vgpi_nxt(indgi) = cgpi.*ones(size(indgi));
    Ugpi_nxt(indgi) = Ugpi(indgi) + dgpi.*ones(size(indgi));
    gpicurr_spk = gpe_zeros;
    gpicurr_spk(indgi) = ones(size(indgi));
    
    tmpgi1 = reshape(gpicurr_spk(1:m/2,:), m/2*n,1);
    spkhstgpi1= [spkhstgpi1 tmpgi1];
    tmpgi2 = reshape(gpicurr_spk(m/2+1:m,:), m/2*n,1);
    spkhstgpi2= [spkhstgpi2 tmpgi2];
    
    spkhstgpi=[spkhstgpi1 ; spkhstgpi2];
    
    Vgpi = Vgpi_nxt;
    Ugpi = Ugpi_nxt;
    
 %-------------------------------------------------------------end-----------------------%   
% pause()    
end
%  figure()
%  colormap('gray')
%  
%  subplot(2,2,1); 
% imagesc(I_strD1.*Cd1str)
% title('D1 striatum');
% 
%  subplot(2,2,2); 
% imagesc(spkhstgpe)
% title('gpespiking');
% 
% subplot(2,2,3), 
% imagesc(spkhststn)
% title('stnspiking');
% 
% subplot(2,2,4); 
% imagesc(spkhstgpi)
% title('gpispiking');

% fignum3=fignum3+1;

% figure()
% subplot(2,1,1)
% plot(sum(spkhststn1))
% hold on
% 
% subplot(2,1,2)
% plot(sum(spkhststn2))
% hold on
% 
% subplot(2,1,1)
% plot(sum(spkhstgpe1),'r')
% subplot(2,1,2)
% plot(sum(spkhstgpe2),'r')
% 
% subplot(2,1,1)
% plot(sum(spkhstgpi1),'g')
% subplot(2,1,2)
% plot(sum(spkhstgpi2),'g')
% hold off

% pause
end



% pause

