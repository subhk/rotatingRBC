clear all; clc;

cd /home/subha/rrb;

% parameters
N=80; % number of gridpoints
Ek = 5e-6;
els = 0.1;
q = 1;

ii = sqrt(-1); 

%%%%%
[DM, z] = cheb(N-1); 
D=DM;    
DD=DM*DM;
D2=DD; 
D4=DD*DD;  
Z=zeros(N,N); I=eye(N); 
%%%%

mag = 1;
BC = 2;

%%%%%
% BC - strss free, conducting for current
%%%%

%%%%
% Horizonal magnetic field B_x
if mag==1
    d = 0.3;
    exp1 = exp(-z.^2/d^2);

    bx  = z.*exp1;
    aa = max(bx);
    bb = 1/aa;

    bx   = bb*bx;
    dbx  = bb*(exp1 - 2*z.^2.*exp1/d^2); 
    d2bx = bb*(4*z.^3/d^2.*exp1 - 6*z/d^2.*exp1);

    bx_   = diag(bx);
    dbx_  = diag(dbx);
    d2bx_ = diag(d2bx);
end

if mag==2
    bx_   = 1*I;
    dbx_  = 0*I;
    d2bx_ = 0*I;
end



alpha = 0.1:0.1:80;

for i=1:1:length(alpha)
    
    k = alpha(i);
    
    kx=k/sqrt(2);
    ky=kx;
    
    % renaming the differentiation matrices
    dz=D; dzz=DD;

    
%     Delta=dxx+dyy+dzz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BC==1 % stress-free && conducting
        
            %%% A-matrix %%%
        A11 = I;                           A12 = Z;
        A13 = q*(dzz - k^2*I);             A14 = Z;               
        A15 = Z;                           A16 = Z;

        A21 = dz;                          A22 = Z;
        A23 = Z;                           A24 = Ek*(dzz - k^2*I);    
        A25 = -els*dbx_*ii*ky;             A26 = els*bx_*ii*kx; 

        A31 = Ek*(k^4*I - 2*dzz*k^2);      A32 = Ek*dzz;
        A33 = Z;                           A34 = -dz;
        A35 = els*(-d2bx_*ii*kx + bx_*ii*kx*(dzz - k^2*I)); 
        A36 = Z;

        A41 = dzz;                         A42 = -I;
        A43 = Z;                           A44 = Z;
        A45 = Z;                           A46 = Z;

        A51 = bx_*ii*kx;                   A52 = Z;
        A53 = Z;                           A54 = Z;
        A55 = (dzz - k^2*I);               A56 = Z;

        A61 = dbx_*ii*ky;                  A62 = Z;
        A63 = Z;                           A64 = bx_*ii*kx;
        A65 = Z;                           A66 = (dzz-k^2*I);
    

        % system matrices
        A=[A11, A12, A13, A14, A15, A16;
           A21, A22, A23, A24, A25, A26;
           A31, A32, A33, A34, A35, A36;
           A41, A42, A43, A44, A45, A46;
           A51, A52, A53, A54, A55, A56;
           A61, A62, A63, A64, A65, A66];

        E=[0*I, 0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, Ek*k^2*q*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I, 0*I];

        % boundary conditions
        II=eye(6*N); 

    %     ddzz=blkdiag(dzz,dzz,dzz, dzz); 

        ddz=blkdiag(dz,dz,dz,dz,dz,dz); 

        u0=1; uL=N; 
        phi0=N+1; phiL=2*N; 
        T0=2*N+1; TL=3*N; 
        w0=3*N+1; wL=4*N;
        b0=4*N+1; bL=5*N;
        j0=5*N+1; jL=6*N;

        loc=[u0,uL, phi0,phiL, T0,TL, w0,wL, b0,bL, j0,jL]; 

        C=[II([u0,uL,phi0,phiL,T0,TL],:); ddz([w0,wL],:); II([b0,bL],:); ddz([j0,jL],:)];
        
        
        A(loc,:)=C;
        E(loc,:)=0; 
        
        
       % computing eigenmodes
       [U, S]=eig(A,E);

       lm = sort(diag(S));
       lmda=lm(1);


       lmda=lm(real(lm)>0);
       rac(i)=real(min(lmda));
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if BC==2 % noslip && insulating
        
        DN = D;
        DN(1,:) = 0;
        DN(N,:) = 0;
        
        D2N=DD;
        
        S = diag([0;1./(1-z(2:N-1).^2);0]);

        D4N = (diag(1-z.^2)*D^4 - 8*diag(z)*D^3 - 12*D^2)*S;
        
        for p=1:N-2
            D4N(1,p+2)=D4N(1,p+2)+(-1*D4N(1,2)*DN(1,p+2)/DN(1,2));
            D4N(N,p)=D4N(N,p)+(-1*D4N(N,N-1)*DN(N,p)/DN(N,N-1));
        end
        
%         D4N([1,2,N-1,N],:)=0;



        
%         IN = I;
%         c1 = IN([1,N],:);
%         c2 = D([1,N],:);
%         C=[c1;c2];             % the constraint matrix
%         
%         r=[1,2,N-1,N];         % removed degrees of freedom
%         l=3:N-2;               % kept degrees of freedom
%         
%         G=-C(:,r)\C(:,l);      % give-back matrix
%         H=inv(C(:,r));         % non-homogeneous contribution
%         
%         A=D4;                  % operator
%         D4N=A(l,l)+A(l,r)*G;   %implement boundary conditions
%         
%         B=D2;
%         D2N=B(l,l)+B(l,r)*G;   %implement boundary conditions
        
        
        
        
        
        for p=1:N-2
            D2N(1,p+2)=D2N(1,p+2)+(-1*D2N(1,2)*DN(1,p+2)/DN(1,2));
            D2N(N,p)=D2N(N,p)+(-1*D2N(N,N-1)*DN(N,p)/DN(N,N-1));
        end
        
%         D2N([1,2,N-1,N],:)=0;

        D2N(1,2)=0;
        D2N(N,N-1)=0;
        
        
        D4N(1,2)=0;
        D4N(N,N-1)=0;
       
        
        
                    %%% A-matrix %%%
        A11 = I;                          
        A12 = q*(dzz - k^2*I);             A13 = Z;               
        A14 = Z;                           A15 = Z;

        A21 = DN;                          
        A22 = Z;                           A23 = Ek*(dzz - k^2*I);    
        A24 = -els*dbx_*ii*ky;             A25 = els*bx_*ii*kx; 

        A31 = Ek*(k^4*I - 2*D2N*k^2 + D4N);      
        A32 = Z;                           A33 = -dz;
        A34 = els*(-d2bx_*ii*kx + bx_*ii*kx*(dzz - k^2*I)); 
        A35 = Z;

        A41 = bx_*ii*kx;                   
        A42 = Z;                           A43 = Z;
        A44 = (dzz - k^2*I);               A45 = Z;

        A51 = dbx_*ii*ky;                  
        A52 = Z;                           A53 = bx_*ii*kx;
        A54 = Z;                           A55 = (dzz-k^2*I);
    

        % system matrices
        A=[A11, A12, A13, A14, A15;
           A21, A22, A23, A24, A25;
           A31, A32, A33, A34, A35;
           A41, A42, A43, A44, A45;
           A51, A52, A53, A54, A55];

        E=[0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, Ek*k^2*q*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I; 
           0*I, 0*I, 0*I, 0*I, 0*I];
       
       
               % boundary conditions
        II=eye(5*N); 
       
       
        ddz=blkdiag(dz,dz,dz,dz,dz); 

        u0=1; uL=N;  
        T0=N+1; TL=2*N; 
        w0=2*N+1; wL=3*N;
        b0=3*N+1; bL=4*N;
        j0=4*N+1; jL=5*N;

        loc=[u0,uL, T0,TL, w0,wL, b0,bL, j0,jL];

        
        C=[II([u0,uL],:); II([T0,TL],:);  II([w0,wL],:); II([b0,bL],:); ddz([j0,jL],:)]; 
        
        A(loc,:)=C;
        E(loc,:)=0; 
        
        
       % computing eigenmodes
       [U, S]=eig(A,E);

       lm = sort(diag(S));
       lmda=lm(1);


       lmda=lm(real(lm)>0);
       rac(i)=real(min(lmda));
        
    end     
    
    
    
%         % computing eigenmodes
%     [U, S]=eig(A,E);
%     
%     lm = sort(diag(S));
%     lmda=lm(1);
% 
% 
%     lmda=lm(real(lm)>0);
%     rac(i)=real(min(lmda));


    
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
