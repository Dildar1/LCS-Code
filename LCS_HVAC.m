clear all
clc 

display('Enter 1 for Stable Manifolds')
display('Enter 2 for UnStable Manifolds')

sou=input(' Enter Value = ');
runframe=1;
dframe=1;
frame=74;
Nframe=1;
for k=1:runframe

%% ************ FORWARD INTEGRATION *********************************

if sou==1
firstframe=frame;
lastframe=frame+Nframe;

frameread=firstframe:dframe:lastframe;
dt=3*dframe;
h=dt*dframe;
t0=firstframe*dt;        % At this time step this manifold is constructed
T=(lastframe-firstframe)*dt;
t=0:h:T;
end
%% ********************************************************************

%% ************ BACKWARD INTEGRATION *********************************

if sou==2  
firstframe=frame;
lastframe=frame+Nframe;
dt=3*dframe;

frameread=lastframe:-dframe:firstframe;
h=-dt*dframe;
t0=firstframe*dt;       % At this time step this manifold is constructed
T=(lastframe-firstframe)*dt;
t=T:h:0;
end

%% ********************************************************************

%% **********************************************************************
% *********** File Input From CFD software ******************************
p='velocity (';
l=strcat(p,num2str(frameread(1)),').txt');
%xx=dlmread(l,'',1,0);

xx=readtable(l);
xx=table2array(xx);

x0=[xx(:,1)];
y0=[xx(:,2)];
%% **********************************************************************
% x=linspace(min(x0),max(x0),1200);  % domain in x
% y=linspace(min(y0),max(y0),150);  % domain in y

x=linspace(0,9,301);       % domain in x
y=linspace(0,3,99);      % domain in y
[zeta,eta] = meshgrid(x, y);
[ny nx]=size(zeta);
zz=zeta;
ee=eta;

%% ********************** INTEGRATION ********************************   

for i=1:length(t)
  
l=strcat(p,num2str(frameread(i)),').txt');
%xx=dlmread(l,'',1,0);
xx=readtable(l);
xx=table2array(xx);

vx=[xx(:,3)];
vy=[xx(:,4)];

u2=TriScatteredInterp(x0,y0,vx);
v2=TriScatteredInterp(x0,y0,vy);

Zeta=zeta+h*u2(zeta,eta);
Eta=eta+h*v2(zeta,eta);
indx=find(isnan(Zeta)==1);
Zeta(indx)=zeta(indx);
Eta(indx)=eta(indx);
indy=find(isnan(Eta)==1);
Zeta(indy)=zeta(indy);
Eta(indy)=eta(indy);
zeta=Zeta;
eta=Eta;
waitbar(i/length(t))
close
end
x=zz';
y=ee';
X=Zeta';
Y=Eta';
%% **********************************************************************

clear sigma
%% ************* Calculation of flow map and FTLE field *****************
for j=2:ny-1
    for i=2:nx-1
        
        J=[(X(i+1,j)-X(i-1,j))/(x(i+1,j)-x(i-1,j))   (X(i,j+1)-X(i,j-1))/(y(i,j+1)-y(i,j-1)); ...
           (Y(i+1,j)-Y(i-1,j))/(x(i+1,j)-x(i-1,j))   (Y(i,j+1)-Y(i,j-1))/(y(i,j+1)-y(i,j-1))];
       
        Delta=J'*J;
        sigma(i,j)=log(sqrt(max(eig(Delta))))/abs(T);

    end
end
        sigma=abs(sigma)';
%% **********************************************************************

%% *************** Plotting *********************************************

sig=sigma/max(max(sigma));
for i=2:ny-1
    for j=2:nx-1
    if 2*sqrt(zz(i,j)^2+ee(i,j)^2)<=1
        sig(i,j)=0;
    end
    end
end

surf(sig)
axis equal
axis off
view(0,90)
lightangle(180,60)
set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.3,'DiffuseStrength',0.9,...
    'SpecularStrength',1,'SpecularExponent',25,...
    'BackFaceLighting','reverselit')


shading interp
material shiny
colormap hot

end