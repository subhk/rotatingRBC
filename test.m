clear all; clc;
clf;

cd ~/Documents/rb_;

% parameters
N=10; % number of gridpoints
rho=1; % fluid density
mu=0.001; % fluid viscosity
g=1; % gravity
alpha=2.5; % wavenumber in x
d=1; % thermal dilatation
k=1; % thermal diffusivity
Ty=-1; % vertical gradient of temperature
L=1; % domain height

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(N,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
y=y/scale; 
Z=zeros(N,N); I=eye(N); 

% renaming the differentiation matrices
dy=D; dyy=DD;
dx=i*alpha*I; 
dxx=-alpha^2*I;
Delta=dxx+dyy;


% system matrices
A=[mu*Delta, Z, -dx, Z; ...
   Z, mu*Delta, -dy, rho*g*d*I;  ...
   dx, dy, Z, Z;  ...
   Z, -Ty*I, Z, k*Delta];

E=blkdiag(rho*I,rho*I,Z,I);

% boundary conditions
II=eye(4*N); ddy=blkdiag(dy,dy,dy,dy);
u0=1; uL=N; v0=N+1; vL=2*N; T0=3*N+1; TL=4*N;
loc=[u0,uL,v0,vL,T0,TL]; 
C=[ddy([u0,uL],:); II([v0,vL,T0,TL],:)];
A(loc,:)=C;
E(loc,:)=0; 


% computing eigenmodes
[U, S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

plot(alpha,real(s(1:5)),'b*')
xlabel('wavenumber \alpha'); ylabel('exponential growth rate')
grid on; hold on;


% validation against the theory
Pr=mu/(rho*k);
Ra=rho*g*(-Ty*L)*L^3/(mu*k);

j=1:5;
alphavec=linspace(0,8,100);
Stheo=zeros(length(j),length(alphavec));
for ind=1:length(alphavec);
    alpha=alphavec(ind);
    Stheo(:,ind)=-0.5*(1+Pr)*(j.^2*pi^2+alpha^2)+sqrt(0.25*(Pr-1)^2*(j.^2*pi^2+alpha^2).^2+alpha^2*Ra*Pr./(j.^2*pi^2+alpha^2));
end
plot(alphavec,Stheo,'k-')
print('-dpng','-r100','rayleigh_benard.png');







