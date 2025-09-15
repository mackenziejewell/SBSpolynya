% shademaps.m

time = 0:nt-1;
figure(1)
clf
subplot(3,1,1)
surf(time,z,tw)
shading interp
axis([ 0. 72 0. ha  -1.7 0])
view(0,90)
set(gca,'xticklabel',' ')
ylabel('Elevation (m)')
text(0,1.08*ha,'a) Temperature (\circC)')
h1=colorbar('vert');

subplot(3,1,2)
surf(time,z,s)
axis([ 0. 72 0. ha  25 35])
view(0,90)
shading interp
%xlabel('Time (h)')
set(gca,'xticklabel',' ')
ylabel('Elevation (m)')
text(0,1.08*ha,'b) Salinity (psu)')
h1=colorbar('vert');

subplot(3,1,3)
surf(time,z,rho)
axis([ 0. 72 0. ha  0 50 ])
view(0,90)
shading interp
xlabel('Time (h)')
ylabel('Elevation (m)')
text(0,1.08*ha,'c) Density (-1000 kg/m^{3})')
h1=colorbar('vert');

figure(2)
clf
subplot(2,1,1)
surf(time,z,log10(ice1+1.e-9))
%pcolor(time,z,log10(ice1+1.e-9))
axis([ 0. 72 0. ha  -8 -1])
view(0,90)
shading interp
hold on
%contour(time,z,log10(ice1+1.e-9),'-k')
set(gca,'xticklabel',' ')
ylabel('Elevation (m)')
text(0,1.08*ha,'a) Ice volume concentration')
h1=colorbar('vert')

subplot(2,1,2)
surf(time,z,log10(sed1+1.e-9))
axis([ 0. 72 0. ha  -8 -1])
view(0,90)
shading interp
%set(gca,'xticklabel',' ')
xlabel('Time (h)')
ylabel('Elevation (m)')
text(0,1.08*ha,'b) Sediment volume concentration')
h1=colorbar('vert')

figure(3)
clf
subplot(3,1,1)
surf(time(2:73),z,K(:,2:73))
axis([ 0. 72 0. ha  0 .02])
view(0,90)
shading interp
set(gca,'xticklabel',' ')
ylabel('Elevation (m)')
%xlabel('Time (h)')
text(0,1.08*ha,'a) Eddy viscosity (m^2/s)')
h1=colorbar('vert')

subplot(3,1,2)
surf(time,z,l)
axis([ 0. 72 0. ha 0 1.5])
view(0,90)
shading interp
set(gca,'xticklabel',' ')
ylabel('Elevation (m)')
%xlabel('Time (h)')
text(0,1.08*ha,'b) Length scale l (m)')
h1=colorbar('vert')

subplot(3,1,3)
surf(time,z,q)
axis([ 0. 72 0. ha  0 .15])
view(0,90)
shading interp
%set(gca,'xticklabel',' ')
ylabel('Elevation (m)')
xlabel('Time (h)')
text(0,1.08*ha,'c) TKE velocity scale q (m/s)')
h1=colorbar('vert')



