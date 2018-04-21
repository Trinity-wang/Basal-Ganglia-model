%sim stn gpe
clear all
clc
close all

n = 50; % number of neurons
dt=0.1;% step size
fignum=10;
fignum2=100;
fignum3=2000;
 
% pulse duration
%  time1=1000; % start of pulse
%  time2=2000; % end of pulse
niter =2500; % number of iterations Time=niter*dt
% actionselected=zeros(runs,9); 
D2=0.1:0.1:0.9;

%%  Parallel Map initialization
PARALLEL=0 ; % Parallelizing code
Amap=1;      % Create a serial Map
N=1; % Number of times to be averaged
Np=2; % no. of processors to be used
if PARALLEL
  Amap = map([Np 1],{},0:Np-1);   % Create parallel map.
end
A = zeros(N,numel(D2),Amap); 
Aloc=local(A);

myI = global_ind(A,1);  % Get local i indices.
trials=1:numel(Aloc);
%%
for trials=1:N
% initial();
    
 i=1; %frequency and actionselection variable   

  
 for D2=0.1:0.1:0.9
 %% Pulse bin
 time1=1000; % start of pulse
 time2=2000; % end of pulse
 
 

%%  Initialization Module
% random initialization of voltage
 Vstn = -60*(rand(n,n)-0.5*ones(n,n));
 Vgpe = -60*(rand(n,n)-0.5*ones(n,n));
 Vgpi = -60*(rand(n,n)-0.5*ones(n,n));
  
 % izhikevich currents
 Istn = 30*ones(n,n);  
 Igpe = 10*ones(n,n);
 Igpi = 10*ones(n,n);
 
% weights from striatum and stn to gpi
 wsgpi=1;  % from stn to gpi
 wstrgpe=1; % from striatum to gpe
 wstrgpi=.6;% from striatum to gpi
 wstr=[wstrgpe wstrgpi];
 %%
 % [Vstn, spkstn, Vgpe, spkgpe] = stn_gpev5(Vstn, Vgpe, Istn, Igpe,niter,wstr,D2,D2spikes); % STN-GPe module
 [Vstn, spkstn, Vgpe, spkgpe,Vgpi,spkgpi] = stn_gpev6(Vstn, Vgpe,Istn, Igpe, niter,wstr,Vgpi,Igpi,wsgpi,D2,time1,time2,dt); % STN-GPe with GPi
 base=n*n/2;
 bin_size=40;
 pulse = time1:time2;
 pul_dur = length(pulse);
 k = round(pul_dur/bin_size);
 bin = zeros(1,k);

 for p = 1:k
     bin(p)=time1+bin_size;
     if time1>=time2;
         break;
     end
      Ihop1(p) = sum(sum(spkgpi(1:base,time1:bin(p)-1)))/(base*dt*1e-3*bin_size);
      Ihop2(p) = sum(sum(spkgpi(base+1:n*n,time1:bin(p)-1)))/(base*dt*1e-3*bin_size);
     time1 = bin(p);
 end
 Ihop=[Ihop1;Ihop2];
 I=Ihop./max(max(Ihop));
 I = 1-I;
%  figure(fignum)
%  plot(I)
%  fignum=fignum+1;

 Aloc(trials,i)=selection(I);
 i=i+1;
 
% % plots 
%  plotting(spkstn,spkgpe,spkgpi,niter,fignum,dt)
%  fignum3=fignum3+1;
%  figure(fignum2)
%  plot(I(1,:));
%  hold on;
%  plot(I(2,:),'r');
%  hold off
%  
 figure(fignum3)
 colormap('gray')
subplot(3,1,1), imagesc(spkstn)
title('stnspiking');
title(num2str(time1:100:time2));

subplot(3,1,2); 
imagesc(spkgpe)
title('gpespiking');

subplot(3,1,3); 
imagesc(spkgpi)
title('gpispiking');


fignum3=fignum3+1;

%% pulse duration
 time1=1000; % start of pulse
 time2=2000; % end of pulse
%% Frequency calculation
FREQ= freqcal(i,base,time1,time2,spkstn,spkgpe,spkgpi,dt,niter);
 
%% 
%   [V act_sel]= action_select1(I(:,i-1),fignum3,2);

 
  end
%% pulse and background
 
[FREQT]=freq_trial(FREQ,trials);
 
%% 
Aloc
A = put_local(A,Aloc);
%%
%   actionselected
%  stn_pulse_block1=stnfrequency_pulse1
%  stn_pulse_block2=stnfrequency_pulse2
%  
%  stn_bck1= (stn_beforepulse_block1+stn_afterpulse_block1)/2
%  stn_bck2= (stn_beforepulse_block2+stn_afterpulse_block2)/2
%  
%  gpe_pulse_block1=gpefrequency_pulse1
%  gpe_pulse_block2=gpefrequency_pulse2
%   
%  gpe_bck1= (gpe_beforepulse_block1+gpe_afterpulse_block1)/2
%  gpe_bck2= (gpe_beforepulse_block2+gpe_afterpulse_block2)/2
%   
%  gpi_pulse_block1=gpifrequency_pulse1
%  gpi_pulse_block2=gpifrequency_pulse2
% 
%   gpi_bck1= (gpi_beforepulse_block1+gpi_afterpulse_block1)/2
%   gpi_bck2= (gpi_beforepulse_block2+gpi_afterpulse_block2)/2
 end
gpe_rad=0.5
gpe_amp=2

save 30thdec_as_(gpe_amp=2)_thresh(0,27)_10trials;

%% action selection percentage calculation
% D2=0.1:0.1:0.9;
% ZERO=zeros(runs,numel(D2));
% TWO=zeros(runs,numel(D2));
% ONE=zeros(runs,numel(D2));
% 
% zero=0;two=0;one=0;
% for i=1:numel(D2)
% for j=1:trials
% if actionselected(j,i)==0
% zero=zero+1;
% elseif actionselected(j,i)==2
% two=two+1;
% else
% one=one+1;
% end
% end
% ZERO(i)=zero;zero=0;TWO(i)=two;two=0;ONE(i)=one;one=0;
% end
% 
% zeropercent=(ZERO./trials).*100;
% onepercent=(ONE./trials).*100;
% twopercent=(TWO./trials).*100;
% 
% figure(1)
% plot(D2,zeropercent,'r--*');
% hold on
% plot(D2,twopercent,'g--*');
% hold on
% plot(D2,onepercent);
% hold off
% title('percentage of action selection :GO in green,NOGO in red and explore in blue');
% xlabel('Dopamine');
% ylabel('Percentage of actionselected');
% axis([0 1 0 100]);