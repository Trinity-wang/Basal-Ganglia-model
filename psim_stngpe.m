addpath('C:\Users\Alekhya\Documents\MATLAB\pMatlab-v2.0.15\parallel code for action selection_gen0ld');
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
tic()
%%  Parallel Map initialization
PARALLEL=0 ; % Parallelizing code
Amap=1;      % Create a serial Map
N=100; % Number of times to be averaged
Np=3; % no. of processors to be used

% if PARALLEL
%   Amap = map([Np 1],{},0:Np-1);   % Create parallel map.
% end
Aloc = zeros(N,numel(D2)); 
steps= zeros(N,numel(D2));
% Aloc=local(A);
FREQ=[];
FREQT=[];
% trials=1:numel(Aloc);
%%
 for trials=1:N
 initial(N,D2);
    
 i=1; %frequency and actionselection variable   

  
 for D2=0.1:0.1:0.9
     D2;
     
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
 wsgpi=1.15;  % from stn to gpi
 wstrgpe=1; % from striatum to gpe
 wstrgpi=.8;% from striatum to gpi
 wstr=[wstrgpe wstrgpi];
 %%
 % [Vstn, spkstn, Vgpe, spkgpe] = stn_gpev5(Vstn, Vgpe, Istn, Igpe,niter,wstr,D2,D2spikes); % STN-GPe module
 [Vstn, spkstn, Vgpe, spkgpe,Vgpi,spkgpi,IstrD1] = stn_gpev6(Vstn, Vgpe,Istn, Igpe, niter,wstr,Vgpi,Igpi,wsgpi,D2,time1,time2,dt); % STN-GPe with GPi
  
 base=n*n/2;
 bin_size=30;
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
 Ihop=[Ihop1./max(Ihop1);Ihop2./max(Ihop2)];
%  I=Ihop./max(max(Ihop));
I=Ihop;
I = 1-I;
%  figure(fignum) 
%  plot(I)
%  fignum=fignum+1; 
 
 [Aloc(trials,i),steps(trials,i)]=selection(I);
 
   Aloc(trials,i)
% % plots 
%    plotting(spkstn,spkgpe,spkgpi,niter,fignum,dt)  
%  fignum3=fignum3+1;
%  figure(fignum2)
%  plot(I(1,:));
%  hold on;
%  plot(I(2,:),'r');
%  hold off
%  
% fignum3=fignum3+1;

%% pulse duration
 time1=1000; % start of pulse
 time2=2000; % end of pulse
%% Frequency calculation
FREQ= freqcal(i,base,time1,time2,spkstn,spkgpe,spkgpi,dt,niter,FREQ);
 
%% 
%   [V act_sel]= action_select1(I(:,i-1),fignum3,2);
i=i+1;
  
  end
%% pulse and background
 
[FREQT]=freq_trial(FREQ,trials,FREQT); 
 
%% 
Aloc(trials);
% A = put_local(A,Aloc);
toc();
  end
% gpe_rad=0.5
% gpe_amp=2 

 save 13aug_es(0.15)_rs(6)_rg(25)thres(0.15)_5trials;

%% action selection percentage calculation
runs=N;
D2=0.1:0.1:0.9;
ZERO=zeros(1,numel(D2));
TWO=zeros(1,numel(D2));
ONE=zeros(1,numel(D2));
actionselected=Aloc;
zero=0;two=0;one=0;
for i=1:numel(D2)
for j=1:trials
if actionselected(j,i)==0
zero=zero+1;
elseif actionselected(j,i)==2
two=two+1;
else
one=one+1;
end
end
ZERO(i)=zero;zero=0;TWO(i)=two;two=0;ONE(i)=one;one=0;
end

zeropercent=(ZERO./trials).*1;
onepercent=(ONE./trials).*1;
twopercent=(TWO./trials).*1;

figure(1)
plot(D2,zeropercent,'r','linewidth',2);
hold on
plot(D2,twopercent,'g','linewidth',2);
hold on
plot(D2,onepercent,'m','linewidth',2);
hold off
% title('percentage of action selection :GO in green,NOGO in red and explore in blue');
xlabel('Dopamine');
ylabel('Probability of Selection');
axis([0.05 0.95 0 1.015]);