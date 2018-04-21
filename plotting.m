function plotting(spkstn,spkgpe,spkgpi,niter,fignum3,dt)
figure(fignum3)
 for i=1:niter
subplot(3,1,2)
plot((1:niter)*dt,i*spkstn(i,1:niter),'.');
title('stnspiking');
xlabel('time');
ylabel('spatial 2500 neurons');
axis([0 niter*dt 0 2500]);
hold on;

subplot(3,1,1)
plot((1:niter)*dt,i*spkgpe(i,1:niter),'.');
title('gpespiking');
xlabel('time');
ylabel('spatial 2500 neurons');
axis([0 2500*dt 0 2500]);
hold on;

subplot(3,1,3)
plot((1:niter)*dt,i*spkgpi(i,1:niter),'.');
title('gpispiking');
xlabel('time');
ylabel('spatial 2500 neurons');
axis([0 2500*dt 0 2500]);
hold on;
 end
end