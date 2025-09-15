% profview.m
% Process and display output from onr1d.f
%
if(exist('onrprof')~=1),
   onrload
end

nice = 1;
nsed = 2;

nz= min(find(diff(onrprof(:,1))>0)); 
nt = length(onrprof)/nz;
tindex = 1:nz:length(onrprof);
icol = 1;
z=onrprof(1:nz,icol);
h=ceil(z(1,1));
ha=5*ceil(h/5);

icol = icol+1;
tw=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
s = reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
K=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
u=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
v=reshape(onrprof(:,icol),nz,nt);
speed=sqrt(u.^2+v.^2);
icol = icol+1;
ice1=reshape(onrprof(:,icol),nz,nt);
if(nice==2),
  icol = icol+1;
  ice2=reshape(onrprof(:,icol),nz,nt);
end
icol = icol+1;
sed1=reshape(onrprof(:,icol),nz,nt);
if(nsed==2),
  icol = icol+1;
  sed2=reshape(onrprof(:,icol),nz,nt);
end
icol = icol+1;
rho=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
tf=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
stress=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
q=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
l=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
Ri=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
Sm=reshape(onrprof(:,icol),nz,nt);
icol = icol+1;
Sh=reshape(onrprof(:,icol),nz,nt);

i=1;
figure(1)
clf
k=tindex(i):(tindex(i)+(nz-1) );
subplot(2,2,1)
h1a=plot(u(:,i),z,'-r','erasemode','none');
set(h1a,'linewidth',1.2)
xlabel('m/s')
ylabel('Elevation (m)')
axis([-1 1 0 ha])
text(-.5,1.08*ha,'a) Cross-shelf velocity')
hold on
plot([0 0],[0 ha],'-w')

subplot(2,2,2)
h1b=plot(v(:,i),z,'-r','erasemode','none');
set(h1b,'linewidth',1.2)
set(gca,'yticklabel',' ')
axis([ 0 3 0 ha])
text(0,1.08*ha,'b) Alongshore velocity')
xlabel('m/s')
hold on

subplot(2,2,3)
h1c=plot(K(:,i),z,'-b','erasemode','none');
set(h1c,'linewidth',1.2)
axis([0 .04 0 ha])
text(0.,1.08*ha,'c) Eddy viscosity')
xlabel('m2 s-1')
ylabel('Elevation (m)')
hold on

subplot(2,2,4)
h1d=plot(stress(:,i),z,'-r','erasemode','none');
set(h1d,'linewidth',1.2)
axis([ 0 1 0 ha])
set(gca,'yticklabel',' ')
text(0,1.08*ha,'d) Stress')
xlabel('N m-2')
hold on

figure(2)
subplot(2,2,1)
h2a=plot(tw(:,i),z,'-r','erasemode','none');
set(h2a,'linewidth',1.2)
xlabel('deg C')
ylabel('Elevation (m)')
ax = [-2 2 0 ha];
axis(ax);
text(ax(1),1.08*ax(4),'a) Temperature')
hold on

subplot(2,2,2)
h2b=plot(s(:,i),z,'-k','erasemode','none');
set(h2b,'linewidth',1.2);
set(gca,'yticklabel',' ')
xlabel('psu')
ax = [15 35 0 ha];
axis(ax)
text(ax(1),1.08*ax(4),'b) Salinity')
hold on

subplot(2,2,3)
h2c=loglog(ice1(:,i),z,'-b','erasemode','none');
set(h2c,'linewidth',1.2)
axis([1e-6, 1e-1, .1, ha])
xlabel('Concentration')
ylabel('Elevation (m)')
text(1.e-6,2*ha,'c) Ice')
hold on
plot([0 .8],[h h],'-k')

subplot(2,2,4)
h2d1=loglog(sed1(:,i)+eps,z,'-g','erasemode','none');
hold on
h2d2=loglog(sed1(:,i)+eps,z,'--g','erasemode','none');
set(h2d1,'linewidth',1.2)
set(h2d2,'linewidth',1.2)
axis([1e-8, 1e-2, .1, ha])
set(gca,'yticklabel',' ')
xlabel('Concentration')
text(1.e-8,2*ha,'d) Sediment')
plot([0 .8],[h h],'-k')

for i=12:12:length(tindex),
   k=tindex(i):(tindex(i)+(nz-1) );
   set(h1a,'xdata',u(:,i),'ydata',z);
   set(h1b,'xdata',v(:,i),'ydata',z);
   set(h1c,'xdata',K(:,i),'ydata',z);
   set(h1d,'xdata',stress(:,i),'ydata',z);
   set(h2a,'xdata',tw(:,i),'ydata',z);
   set(h2b,'xdata',s(:,i),'ydata',z);
   set(h2c,'xdata',ice1(:,i),'ydata',z);
   set(h2d1,'xdata',sed1(:,i),'ydata',z);
   set(h2d2,'xdata',sed2(:,i),'ydata',z);
   drawnow
   pause(.3)
 end





