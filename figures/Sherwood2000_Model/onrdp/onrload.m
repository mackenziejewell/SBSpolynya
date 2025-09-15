% onrload.m
% Load output from onr1d.f
%
clear
format compact;
direc = '/home1/sherwood/proj/arctic/newcases/';
%direc = '/home3/sherwood/cdroms/cd0/arctic/newcases/1dcases/';
modrun= input('Enter case for output with SINGLE QUOTES: ')
eval(['load ',direc,modrun,'/onrprof.dat'])
eval(['load ',direc,modrun,'/onrts.dat'])
if(exist('nice')~=1),
  nice = 1;
end
if(exist('nsed')~=1),
  nsed = 2;
end
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




