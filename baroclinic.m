clear all;
clc;
% =========================================================
% This code solves rotating plane layer convection 
% with imposed zonal flow - Teed et al GAFD-2011.
% ========================================================= 

% parameters
N = 41; % number of gridpoints

% differentiation matrices
scale = 2; %-2/L;
[z, DM] = chebdif(N,2); 
D = DM(:,:,1)*scale;    
DD = DM(:,:,2)*scale^2;  
% D4=DM(:,:,4)*scale^2; 
z = z./scale; 
Z = zeros(N,N); 
I = eye(N); 


% imposed shear profile
U = z;

Re = 1;          % Reynolds number     
E  = 1e-4;
Pr = 1;
count = 1;

kx_ = 20:0.2:34;
ky_ = 0.0;

for i=1:length(kx_)
    kx = kx_(i);
    ky = ky_;
    
    strt = 1e4;  % Least guess
    endd = 1e12; % Max guess

    diff = endd - strt;
    
    while diff>0.5

        lmda1 = fun_call(U, Re, E, Pr, strt, N, kx, ky); 
        lmda2 = fun_call(U, Re, E, Pr, endd, N, kx, ky);

        if lmda1<0 && lmda2>0 
            midd = (strt+endd)/2;
        end
            
        midd = (strt+endd)/2;
        lmda3 = fun_call(U, Re, E, Pr, midd, N, kx, ky);

        if lmda3<0
           strt = midd;
           endd = endd;
           
        elseif lmda3>0 
           strt = strt;
           endd = midd;
        end

        diff = endd - strt;

        if diff<0.5
           Rac = endd
           rac(i,1) = Rac;
           frequencies(i,1) = lmda3
        end

        count = count+1;

    end
    
end

plot(kx_,rac)
[minn,j] = min(rac);
k_crit = kx_(j) % Critical wavenumber k
Ra_crit = min(rac) % Critical Raleigh number 




