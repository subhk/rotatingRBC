clear all;
clc;

cd /home/subha/ShearInstability/;

N=60;
Re=150;
Ro=1;

ii = sqrt(-1); 

%%%%%
[z, DM] = chebdif(N, 4); 
D=DM(:,:,1);    
DD=DM(:,:,2);
D4=DM(:,:,4);  
Z=zeros(N,N); I=eye(N); 
%%%%

%%%
d   = 0.4;

U   = tanh(z/d);
dU  = (1-tanh(z/d).*tanh(z/d))/d;
d2U = 2*tanh(z/d).*(tanh(z/d).*tanh(z/d)-1)/d^2;

U_   = diag(U);
dU_  = diag(dU);
d2U_ = diag(d2U);
%%%

kx = 0.01:0.01:1;
ky = 0.01:0.01:1;

for i=1:1:length(kx)
    
    kx_ = kx(i)
    
%     printf('one outer loop is over')
    
    for j=1:1:length(ky)
    
        ky_ = ky(j);
        
        k2 = kx_*kx_ + ky_*ky_;
 
        k4 = k2*k2;
    % renaming the differentiation matrices
        dz=D; dzz=DD;
        
        
    % construction of matrix
    
       A11 = U_*(dzz-k2*I) - 1/(ii*kx_*Re)*(k4*I - 2*k2*dzz) - d2U_;
       A12 = dzz;
       A13 = 1/(ii*kx_*Ro)*I;
       
       A21 = dzz;
       A22 = -I;
       A23 = Z;
       
       A31 = ky_*ky_/(ii*kx_*Ro)*(I + Ro*dU_);
       A32 = Z;
       A33 = U_ - 1/(ii*kx_*Re)*(dzz - k2*I);
       
       B11 = 1/(kx_)*(dzz - k2*I);
       B12 = Z;
       B13 = Z;
       
       B21 = Z;
       B22 = Z;
       B23 = Z;
       
       B31 = Z;
       B32 = Z;
       B33 = 1/kx_*I;
    
       
       A = [A11, A12, A13; A21, A22, A23; A31, A32, A33];
       
       B = [B11, B12, B13; B21, B22, B23; B31, B32, B33];
       
       
               % boundary conditions
        II=eye(3*N); 

    %     ddzz=blkdiag(dzz,dzz,dzz, dzz); 

        ddz=blkdiag(dz,dz,dz); 

        u0=1; uL=N; 
        phi0=N+1; phiL=2*N; 
        w0=2*N+1; wL=3*N; 
        

        loc=[u0,uL, phi0,phiL, w0,wL]; 

        C = [II([u0,uL,phi0,phiL],:); II([w0,wL],:)];
        
        
        A(loc,:)=C;
        B(loc,:)=0;
        
       [V, S]=eig(A,B);

       lm = sort(diag(S));
       lmda=lm(1);


       lmda=lm(real(lm)>0);
       growth_rate=real(min(lmda))
       
       omega(i,j) = growth_rate;

    end
end