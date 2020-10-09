clear all; clc;

cd ~/Documents/rb_;

% parameters
N=80; % number of gridpoints
E = 1e-5;
invE = 1/E;

% differentiation matrices
% scale=2; %-2/L;
% [DM, z] = chebdif(N,4); 
% D=DM(:,:,1)*scale;    
% DD=DM(:,:,2)*scale^2;  
% D4=DM(:,:,4)*scale^2; 
% z=z/scale; 
% Z=zeros(N,N); I=eye(N); 

%%%%%
[DM, z] = cheb(N-1); 
D=DM;    
DD=DM*DM;   
D4=DD*DD;  
Z=zeros(N,N); I=eye(N); 

alpha = 20:0.01:80;

for i=1:1:length(alpha)
    
    k = alpha(i)
    
    kx=k/sqrt(2);
    ky=kx;
    
    % renaming the differentiation matrices
    dz=D; dzz=DD;
    dx=i*kx*I; 
    dy=i*ky*I; 
    dxx=-kx^2*I;
    dyy=-ky^2*I;
    
%     Delta=dxx+dyy+dzz;
    
    
    Delta=-k^2*I+dzz;

    % system matrices
    A=[-Delta, -invE*dz, Z, Z;
       invE*dz, 2*k^2*dzz-k^4*I, -dzz, Z; 
       Z, dzz, -I, Z;
       Z, -I, Z, -Delta];

    E=[Z, Z, Z, Z; 
        Z, Z, Z, -k^2*I; 
        Z, Z, Z, Z; 
        Z, Z, Z, Z];

    % boundary conditions
    II=eye(4*N); 
    
%     ddzz=blkdiag(dzz,dzz,dzz, dzz); 
    
    ddz=blkdiag(dz,dz,dz,dz); 
    
    w0=1; wL=N; 
    u0=N+1; uL=2*N; 
    phi0=2*N+1; phiL=3*N; 
    T0=3*N+1; TL=4*N;
    
    loc=[u0,uL, w0,wL, phi0,phiL, T0,TL]; 
    
    C=[ddz([w0,wL],:);  II([u0,uL,phi0,phiL,T0,TL],:)];  %II([phi0,phiL],:);  II([T0,TL],:)];
    
    A(loc,:)=C;
    E(loc,:)=0; 


    % computing eigenmodes
    [U, S]=eig(A,E);
    
    lm = sort(diag(S));
    lmda=lm(1);


    lmda=lm(real(lm)>0);
    rac(i)=real(min(lmda));
    
%     s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
%     rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

    
end


% plot(alpha,real(s(1:5)),'b*')
% xlabel('wavenumber \alpha'); ylabel('exponential growth rate')
% grid on; hold on;
% 
% 
% % validation against the theory
% Pr=mu/(rho*k);
% Ra=rho*g*(-Ty*L)*L^3/(mu*k);
% 
% j=1:5;
% alphavec=linspace(0,8,100);
% Stheo=zeros(length(j),length(alphavec));
% for ind=1:length(alphavec);
%     alpha=alphavec(ind);
%     Stheo(:,ind)=-0.5*(1+Pr)*(j.^2*pi^2+alpha^2)+sqrt(0.25*(Pr-1)^2*(j.^2*pi^2+alpha^2).^2+alpha^2*Ra*Pr./(j.^2*pi^2+alpha^2));
% end
% plot(alphavec,Stheo,'k-')
% print('-dpng','-r100','rayleigh_benard.png');
