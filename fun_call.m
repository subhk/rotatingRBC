function [lmda] = fun_call(U, Re, E, Pr, strt, N, kx, ky)

    Ra = strt;

    % differentiation matrices
    scale = 2; %-2/L;
    [z, DM] = chebdif(N,4); 
    D = DM(:,:,1)*scale;    
    DD = DM(:,:,2)*scale^2;  
    % D4=DM(:,:,4)*scale^2; 
    z = z./scale; 
    Z = zeros(N,N); 
    I = eye(N); 
    
    invE = 1/E;
    k2 = kx*kx + ky*ky;

    ii = sqrt(-1);

    dz = D; dzz = DD;
  
    A11 = -dz*invE - ii*ky*Re*I;
    A12 = k2*I - dzz + ii*kx*Re*diag(U); 
    A13 = Z;
    A14 = Z;
    
    A21 = -k2*k2*I + 2*k2*dzz - ii*k2*kx*Re*diag(U) + ii*kx*Re*diag(U)*dzz;
    A22 = dz*invE;
    A23 = -dzz;
    A24 = k2*Ra*I;
    
    A31 = -I - (ii*Pr*Re/(E*Ra*k2))*ky*dz;
    A32 = (ii*Pr*Re/(E*Ra*k2))*kx*I;
    A33 = Z;
    A34 = k2*I - dzz + ii*kx*Pr*Re*diag(U);

    A41 = dzz;
    A42 = Z;
    A43 = -I;
    A44 = Z;
    
    A = [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44];
    
    B11 = Z;
    B12 = -I;
    B13 = Z;
    B14 = Z;
    
    B21 = k2*I - dzz;
    B22 = Z;
    B23 = Z;
    B24 = Z;

    B31 = Z;
    B32 = Z;
    B33 = Z;
    B34 = -I*Pr;

    B41 = Z;
    B42 = Z;
    B43 = Z;
    B44 = Z;
    
    B = [B11 B12 B13 B14; B21 B22 B23 B24; B31 B32 B33 B34; B41 B42 B43 B44];

    % boundary conditions
    II = eye(4*N); 
    
%    ddzz=blkdiag(dzz,dzz,dzz, dzz); 
    
    ddz = blkdiag(dz,dz,dz,dz); 
    
    u0 = 1; uL = N; 
    w0 = N+1; wL = 2*N; 
    phi0 = 2*N+1; phiL = 3*N; 
    T0 = 3*N+1; TL = 4*N;
    
    loc = [u0,uL, w0,wL, phi0,phiL, T0,TL]; 
    
    C = [II([u0,uL],:); ddz([w0,wL],:); II([phi0,phiL,T0,TL],:)];  
    
    A(loc,:) = C;
    B(loc,:) = 0; 
    
%     Am1 = (full(Am1));
%     Bm1 = (full(Bm1));
    
    % Solving the eigenvalue problem
    [V, L] = eig(A,B);
    lm = sort(diag(L));
    lmda=lm(1);
    for i=1:size(lm)
        if lm(i)>0 && real(lm(i))~=Inf 
        lmda=lm(i);        
        end
    end
    
    lmda;
end

