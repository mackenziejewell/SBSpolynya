% tsview.m
% run profview first

format compact
echo off
time=onrts(:,1);
days=ceil(time(length(time))/(24*3600))
hours = days*24;
taus = onrts(:,2);
fprintf('Modrun: %s\n',modrun);
fprintf('Mean surface stress: %f\n',mean(taus));
Hs = onrts(:,3);
fprintf('Mean wave height: %f m\n',mean(Hs));
Td = onrts(:,4);
fprintf('Mean wave period: %f s\n',mean(Td));
ub = onrts(:,5);
fprintf('Mean wave-orbital vel.: %f m/s\n',mean(ub));
Cd = onrts(:,6);
fprintf('Drag coefficient stats: \n')
stats(Cd);
taub = onrts(:,7:8);
fprintf('Bottom shear-stress stats: \n')
stats(sqrt(taub(:,1).^2+taub(:,2).^2));
tausfmx = onrts(:,9);
fprintf('Tausf stats: \n')
stats(tausfmx);
Qw = onrts(:,10);
fprintf('Mean heat flux: %f; std. dev. %f\n',mean(Qw),std(Qw));
xtran = onrts(:,11);
fprintf('Cross-shelf trans. stats: \n');
stats(xtran);



if(0)
subplot(4,1,1)
h1=plot(time/3600,taus,'-w');
set(h1,'linewidth',2);
hold
h1=plot(time/3600,taub,'--r');
set(h1,'linewidth',1.2)


subplot(5,1,1)
h1= plot(tw(2,:),'-r');
set(h1,'linewidth',1.2)
ax=[0 hours -2 0];
axis(ax)
hold on
h1= plot(tf(2,:),'--r');
set(h1,'linewidth',1.2)
set(gca,'xticklabel',' ')
ylabel('\circC')
text(ax(1),ax(4)+.12*(ax(4)-ax(3)),...
    'a) Surface temperature (-) and freezing point (--)')
end

figure(1)
clf

subplot(4,1,1)
h1= plot(ice1(1,:),'-b');
set(h1,'linewidth',2)
axis([0 hours 0 .1])
ax= axis;
ax(2) = hours;
axis(ax);
hold on
set(gca,'xticklabel',' ')
ylabel(' ')
text(ax(1),ax(4)+.12*(ax(4)-ax(3)),...
    'a) Ice volume concentration at surface')

subplot(4,1,2)
dz = z(2)-z(3);
dztb = 2*z(nz);
totalice1 =sum(ice1(2:nz-1,:)*dz)+(ice1(1,:)+ice1(nz,:))*dztb; 
h1=plot(totalice1,'-b');
axis([0 hours 0 .1])
set(h1,'linewidth',2)
set(gca,'xticklabel',' ')
ylabel('meters')
ax= axis;
ax(2) = hours;
axis(ax);
text(ax(1),1.12*ax(4),'b) Total volume of ice per m^2')


subplot(4,1,3)
totalsed1 =2650*(sum(sed1(2:nz-1,:)*dz)+(sed1(1,:)+sed1(nz,:))*dztb); 
h1=plot(totalsed1./h,'-r');
set(h1,'linewidth',2);
axis([ 0 hours 0 .5 ])
hold on
totalsed2 =2650*(sum(sed2(2:nz-1,:)*dz)+(sed2(1,:)+sed2(nz,:))*dztb); 
%h1=plot(totalsed2./h,'-r');
%set(h1,'linewidth',2);
ylabel('kg m^{-3}')
ax= axis;;
ax(2) = hours;
axis(ax);
set(gca,'xticklabel',' ')
text(ax(1),1.12*ax(4),'c) Depth-mean sediment concentration')

subplot(4,1,4)
ice1sed1 = (sum(ice1(2:nz-1,:).*sed1(2:nz-1,:)*2650*dz)+...
            (ice1(1,:).*sed1(1,:)*2650*dztb)+...
            (ice1(nz,:).*sed1(nz,:)*2650*dztb))...
	    ./(sum(ice1(2:nz-1,:)*dz)+ice1(1,:)*dztb+ice1(nz,:)*dztb);
h1=plot(ice1sed1,'-k')
set(h1,'linewidth',2);
hold on
ax= axis;
ax(2) = hours;
axis(ax);
ylabel('kg m^{-3}')
text(ax(1),1.1*ax(4),'d) Sediment concentration in ice')
xlabel('Time (hours)')

fprintf('h: %f\n',h)
%assumes delta-t = 1 hour:
if(any(ice1(1,:)>.01)),
  ice_time=min(find(ice1(1,:)>.01))-1;
  fprintf('Ice formation time: %f\n',ice_time);
else
  ice_time = 99;
  fprintf('Ice conc. always < 0.01\n');
end
meansed1 = totalsed1(nt)/h;

subplot(4,1,1)
hold on
text(ice_time,0.017,'\downarrow')
text(ice_time-2.6,.047,'{\itt_{ice}} ({\it{C_I}} > 0.01)')
subplot(4,1,4)
hold on
mix_time = input('Enter mix time: ')
text(mix_time-.5,.09,'\downarrow')
text(mix_time-5,.08,'{\itt_{mix}}')
%plot([mix_time mix_time],[0 1],'--k');

figure(2)
clf
subplot(2,1,1)
h1=plot(mean(K(:,2:nt)),'-k');
hold on
set(h1,'linewidth',2);
axis([0 hours 0 .025])
set(gca,'xticklabel','')
ylabel('K (m^2/s)')
plot([mix_time mix_time],[0 .025],'--k');
text(0,1.08*.025,'a) Depth-Mean Eddy Viscosity')

subplot(2,1,2)
h2=plot( max(rho)-min(rho),'-r');
hold on
set(h2,'linewidth',2);
axis([0 hours 0 12])
ylabel('\Delta\rho (kg/m^3)')
xlabel('Time (hours)')
plot([mix_time mix_time],[0 12],'--k');
text(0,1.08*12,'b) Maximum Density Contrast')

fprintf(' tmix: %f \n tice: %f \n tmix/tice: %f \n', ...
   mix_time,ice_time,mix_time/ice_time);
fprintf(' Cmax: %f \n Max co-conc: %f \n',max(totalsed1/h),max(ice1sed1))

