% -------------------------------------------------------------------------
% Compressible and Non Dimensional Reynolds equation 2D -------------------
% -------------------------------------------------------------------------
%
% Compressible and non dimensional solver
%
% Reynolds equation + Energy equation
% 
% 18/03/2016 Andrea Codrignani
% 21/03/2016
% -------------------------------------------------------------------------
clc; close all; clear all;


% Data --------------------------------------------------------------------

% switches
geo_type  = 0;      % kind of geometry used: 0) channel/slider/dimple, 1) Pin 2D, 2) Pin 3D
an        = 0;      % 0 no analytical solution, 1 analytical solution for incompressible Slide Bearing 

eos       = 2;      % Equation of state: 0) incompressible, 1) rho(P), 2) rho(P,T)
mu_T      = 1;      % Constitutive law:  0) viscosity always constant, 1) viscosity depends on Temperature

BC_Pin    = 1;      % Boundary condition for the Pressure inlet:     0 Dirichlet, 1 Neumann
BC_Pout   = 0;      % Boundary condition for the Pressure outlet:    0 Dirichlet, 1 Neumann
BC_Tout   = 1;      % Boundary condition for the Temperature outlet: 0 Dirichlet, 1 Neumann, 2 Second Derivative null
BC_side   = 1;      % Boundary condition for T and P on the side:    0 Dirichlet, 1 Neumann

scheme_T  = 0;      % Numerical scheme for the Energy Ttransport term: 0) central difference, 1) upwind

it_max    = 50;     % max number of iteration
toll      = 1e-6;   % tollerance for the iterative loop

cfr_sparc = 0;      % enable the comparison with sparc:  0) none, 1) channel, 2) slide bearing 
file_name = 'F';    % name of the sparc simulation in the folder ~/sparc/data/Box_3D/ or ~/sparc/data/Patin/

% Flow properties
P1        = 104755;
P2        = 1*P1;
T1        = 323.15;
T2        = 1*T1;
U         = 0.5;
W         = 0;

fluid = 0;          % type of fluid:  0) Oil,  1) Air

% Fluid properties
if fluid == 0     % PAO 2
    k_0       = 0.154;                                  % heat conductivity J/(s*m*K)
    fun_mu    = @(T) 0.0180 .* exp(-0.0248.*(T-293));   % Pierce formula
    mu_0      = fun_mu(T1);                             % reference viscosity
    Psg       = 316450000;                              % Pressure coeff for stiffened gas
    gamma_sg  = 2.8;                                    % gamma for stiffened gas
    Cv_0      = 2000;                                   % specific heat coefficient at constant volume
    Rsg       = (gamma_sg-1)*Cv_0;                      % Gas constant for stiffened gas
    fun_rho   = @(P,T) (P+Psg*gamma_sg)/(Rsg*T);        % Stiffened-Ideal Gas Equation of State
    rho_0     = fun_rho(P1,T1);                         % reference density
    
elseif fluid == 1  % Air
    k_0       = 0.0257;                                 % heat conductivity J/(s*m*K)
    fun_mu    = @(T) 1.458e-6.*(T.^(3/2))./(T+110.4);   % Suttherland formula
    mu_0      = fun_mu(T1);                             % reference viscosity
    Psg       = 0;                                      % Pressure coeff for stiffened gas
    gamma_sg  = 1.4;                                    % gamma for stiffened gas
    Cv_0      = 717.9;                                  % specific heat coefficient at constant volume
    Rsg       = (gamma_sg-1)*Cv_0;                      % Gas constant for stiffened gas
    fun_rho   = @(P,T) (P+Psg*gamma_sg)/(Rsg*T);        % Stiffened-Ideal Gas Equation of State
    rho_0     = fun_rho(P1,T1);                         % reference density
end

% Geometry -------------------------------------------------------------------------------
if geo_type == 0;     % Slide Bearing with dimple
    
    D   = 100e-6;      % [m]
    H   = 0e-6;        % [m]
    x_m = 500e-6;      % x offset of center of the dimple
    y_m = 500e-6;      % y offset of center of the dimple
    R   = D/2;
    
    if cfr_sparc == 0
        Lx  = 1.0000e-03;  % [m]
        Lz  = 0.001;       % [m]
        h2  = 1.e-07;      % height at the point (Lx,0)
        h1  = h2;
        h3  = h1;          % height at the point (0,Lz)
    elseif cfr_sparc == 1
        Lx  = 1.0000e-03;  % [m]
        Lz  = 0.001;       % [m]
        h2  = 1.e-06;      % height at the point (Lx,0)
        h1  = h2;
        h3  = h1;          % height at the point (0,Lz)
    elseif cfr_sparc == 2
        Lx  = 1.0000e-03;  % [m]
        Lz  = 0.001;       % [m]
        h2  = 1e-07;      % height at the point (Lx,0)
        h1  = 2.2*h2;
        h3  = h1;          % height at the point (0,Lz)
    end
    
    % Discretization
    Nx = 2000;
    Nz = 3;
    [x,z,hh] = dimp_2D(Nx,Nz,D,H,x_m,y_m,Lx,Lz,h1,h2,h3);
    
elseif geo_type == 1;   % Pin 2D
    
    Nx      = 1001; % must be an odd number in order to have the optimal precision
    Nz      = 3;    % must be 3
    load(sprintf('Pin_surf_1D_Nx%s.mat',num2str(Nx)));
    
    h_gap  = 1e-3;
    
    hh_max = max(max(h_1D))+h_gap;
    hh_3D = hh_max - h_1D;
    hh = [hh_3D',hh_3D',hh_3D'];
    h2 = min(hh_3D);              % takes the minimum gap height as reference for adimensionalisation
    
    x = x_1D;
    Lx = max(x)-min(x);
    Lz = 1e-3;
    z = [-0.5 0 0.5]*Lz;
    
elseif geo_type == 2;   % Pin 3D
    
    Nx = 50; 
    Nz = 50;
    load(sprintf('Pin_surf_Nx%s_Nz%s.mat',num2str(Nx),num2str(Nz)));
    h_gap  = 100e-6;
    hh_max = max(max(h_pin))+h_gap;
    hh = hh_max - h_pin;
    x = -x;
    Lx = max(x)-min(x);
    Lz = max(z)-min(z);
    
end

% Non-Dimensionalisation --------------------------------------------------

% other dimensional quantities
rho     = rho_0*ones(Nx,Nz);
mu      = mu_0*ones(Nx,Nz);

% reference quantities
L_r   = Lx;                                 % ref lenght
U_r   = sqrt(U^2+W^2);                      % ref velocity
mu_r  = mu_0;                               % ref dimanic viscosity
rho_r = rho_0;                              % ref density
h_r   = h2;                                 % ref height
Cv_r  = Cv_0;                               % ref heat capacity at cosnstant volume
k_r   = L_r*rho_r*Cv_r*U_r/2;               % ref thermal conductivity coefficient
P_r   = 6*mu_r*U_r*L_r/(h_r^2);             % ref pressure
T_r   = 2*mu_r*U_r*L_r/(rho_r*Cv_r*h_r^2);  % ref temperature

% non dimensional quantities
x_nd    = x/L_r;
z_nd    = z/L_r;
hh_nd   = hh/h_r;
U_nd    = U/U_r;
W_nd    = W/U_r;
Cv_nd   = Cv_0/Cv_r;
mu_nd   = mu/mu_r;
rho_nd  = rho/rho_r;
k_nd    = k_0/k_r;
P1_nd   = P1/P_r;
P2_nd   = P2/P_r;
T1_nd   = T1/T_r;
T2_nd   = T2/T_r;


% Grid definition (constant discretization)
dx_nd = x_nd(2)-x_nd(1);
dz_nd = z_nd(2)-z_nd(1);


% Iterative loop ----------------------------------------------------------
res_P   =  1;
res_T   =  1;
res_max =  1;
iter    =  0;
Fn0     = -1;
Tr0     = -1;
it_res  = zeros(it_max,3);

% allocation
hpx     = zeros(Nx,Nz);
hpz     = zeros(Nx,Nz);
drho_dx = zeros(Nx,Nz);
drho_dz = zeros(Nx,Nz);
dmu_dx  = zeros(Nx,Nz);
dmu_dz  = zeros(Nx,Nz);
dP_x    = zeros(Nx,Nz); 
dP_z    = zeros(Nx,Nz);
P_nd    = zeros(Nx,Nz);
T_nd    = zeros(Nx,Nz);

du_eeP_dx    = zeros(Nx,Nz);
du_ee_rc_dx  = zeros(Nx,Nz);
dw_eeP_dz    = zeros(Nx,Nz);
dw_ee_rc_dz  = zeros(Nx,Nz);

AP = spalloc(Nx*Nz,Nx*Nz,5*Nx*Nz); % allocation of a sparse matrix (the last input is the maximum number of non zero elements)
AT = spalloc(Nx*Nz,Nx*Nz,5*Nx*Nz);
bP = ones(Nx*Nz,1);
bT = ones(Nx*Nz,1); 

ti = tic;

while res_max > toll && iter < it_max
    
    iter = iter + 1;
    
    % differentiation of the known variables ------------------------------
    
    % gap height hh
    for j = 1 : Nz
        hpx(:,j) = diff_2o_generic(x_nd,hh_nd(:,j));
    end
    for i = 1 : Nx
        hpz(i,:) = diff_2o_generic(z_nd,hh_nd(i,:));
    end
    % density rho
    for j = 1 : Nz
        drho_dx(:,j) = diff_2o_generic(x_nd,rho_nd(:,j));
    end
    for i = 1 : Nx
        drho_dz(i,:) = diff_2o_generic(z_nd,rho_nd(i,:));
    end
    % viscosity mu
    for j = 1 : Nz
        dmu_dx(:,j) = diff_2o_generic(x_nd,mu_nd(:,j));
    end
    for i = 1 : Nx
        dmu_dz(i,:) = diff_2o_generic(z_nd,mu_nd(i,:));
    end
    % K
    K = rho_nd.*(hh_nd.^3)./(mu_nd);
    dK_dx = ( drho_dx.*(hh_nd.^3)./mu_nd + 3*(hh_nd.^2).*hpx.*rho_nd./mu_nd - rho_nd.*(hh_nd.^3).*dmu_dx./(mu_nd.^2) );
    dK_dz = ( drho_dz.*(hh_nd.^3)./mu_nd + 3*(hh_nd.^2).*hpz.*rho_nd./mu_nd - rho_nd.*(hh_nd.^3).*dmu_dz./(mu_nd.^2) );
    
    
    % Assembling of Momentum Equation Matrix ---------------------------------------
    for i = 1 : Nx
        for j = 1 : Nz
            
            % BC inlet
            if i == 1
                
                if BC_Pin == 0
                    AP((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bP((i-1)*Nz+j)            = P1_nd;
                elseif BC_Pin == 1
                    AP((i-1)*Nz+j,(i-1)*Nz+j)    = -1/dx_nd;
                    AP((i-1)*Nz+j,(i-1)*Nz+j+Nz) =  1/dx_nd;
                    bP((i-1)*Nz+j)               =  0;
                end
            
            % BC outlet    
            elseif i == Nx  
                
                if BC_Pout == 0
                    AP((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bP((i-1)*Nz+j)            = P2_nd;
                elseif BC_Pout == 1
                    AP((i-1)*Nz+j,(i-1)*Nz+j)    =  1/dx_nd;
                    AP((i-1)*Nz+j,(i-1)*Nz+j-Nz) = -1/dx_nd;
                    bP((i-1)*Nz+j)               =  0;
                end
            
            % BC side 1    
            elseif j == 1  
                
                if BC_side == 0     % Dirichlet
                    AP((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bP((i-1)*Nz+j)            = P2_nd;
                elseif BC_side == 1 % Neumann
                    AP((i-1)*Nz+j,(i-1)*Nz+j)   = -1/dz_nd;
                    AP((i-1)*Nz+j,(i-1)*Nz+j+1) =  1/dz_nd;
                    bP((i-1)*Nz+j)              = 0;
                end
                
            % BC side 2    
            elseif j == Nz  
                
                if BC_side == 0     % Dirichlet
                    AP((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bP((i-1)*Nz+j)            = P2_nd;
                elseif BC_side == 1 % Neumann
                    AP((i-1)*Nz+j,(i-1)*Nz+j)   =  1/dz_nd;
                    AP((i-1)*Nz+j,(i-1)*Nz+j-1) = -1/dz_nd;
                    bP((i-1)*Nz+j)              = 0;
                end
               
            % Domain
            else 
                
                AP((i-1)*Nz+j,(i-1)*Nz+j)    = -2*K(i,j)/(dx_nd^2) -2*K(i,j)/(dz_nd^2);       % P i j
                
                AP((i-1)*Nz+j,(i-1)*Nz+j-1)  = -1./(2*dz_nd)*dK_dz(i,j) + K(i,j)/(dz_nd^2);   % P i j-1
                AP((i-1)*Nz+j,(i-1)*Nz+j+1)  = +1./(2*dz_nd)*dK_dz(i,j) + K(i,j)/(dz_nd^2);   % P i j+1
                
                AP((i-1)*Nz+j,(i-1)*Nz+j-Nz) = -1./(2*dx_nd)*dK_dx(i,j) + K(i,j)/(dx_nd^2);   % P i-1 j
                AP((i-1)*Nz+j,(i-1)*Nz+j+Nz) = +1./(2*dx_nd)*dK_dx(i,j) + K(i,j)/(dx_nd^2);   % P i+1 j
                
                bP((i-1)*Nz+j) = rho_nd(i,j)*U_nd*hpx(i,j) + hh_nd(i,j)*U_nd*drho_dx(i,j) + ...
                                 rho_nd(i,j)*W_nd*hpz(i,j) + hh_nd(i,j)*W_nd*drho_dz(i,j);
                
            end
            
        end
    end
    
    % Solving Momentum Equation -------------------------------------------
    Pv = AP\bP;
    
    % remapping the pressure field
    for i = 1 : Nx 
            P_nd(i,:) = Pv(((i-1)*Nz+1):(i*Nz));
    end
    
    % Pressure gradient
    for j = 1 : Nz
        dP_x(:,j) = diff_2o_generic(x_nd,P_nd(:,j));
    end
    for i = 1 : Nx
        dP_z(i,:) = diff_2o_generic(z_nd,P_nd(i,:));
    end
    
    % Terms of Energy Equation Matrix -------------------------------------
    u_ee = -(hh_nd.^3).*(dP_x./(mu_nd)) + U_nd*hh_nd;
    w_ee = -(hh_nd.^3).*(dP_z./(mu_nd)) + W_nd*hh_nd;
    u_ee_rc = Cv_nd*rho_nd.*u_ee;
    w_ee_rc = Cv_nd*rho_nd.*w_ee;
    for j = 1 : Nz
        du_eeP_dx(:,j)    = diff_2o_generic(x_nd,3*P_nd(:,j).*u_ee(:,j));
        du_ee_rc_dx(:,j)  = diff_2o_generic(x_nd,u_ee_rc(:,j));
    end
    for i = 1 : Nx
        dw_eeP_dz(i,:)    = diff_2o_generic(z_nd,3*P_nd(i,:).*w_ee(i,:));
        dw_ee_rc_dz(i,:)  = diff_2o_generic(z_nd,w_ee_rc(i,:));
    end
    
    lap_kh = k_nd*hh_nd;
    
    diss_x = 3*(hh_nd.^3).*(dP_x.^2)./(mu_nd) + (U_nd^2)*mu_nd./hh_nd;
    diss_z = 3*(hh_nd.^3).*(dP_z.^2)./(mu_nd) + (W_nd^2)*mu_nd./hh_nd;
    
    % Assembling Energy Equation Matrix -----------------------------------
    for i = 1 : Nx
        for j = 1 : Nz
            
            % BC inlet
            if i == 1
                
                AT((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                bT((i-1)*Nz+j)            = T1_nd;
                
            % BC outlet Dirichlet or Neumann
            elseif i == Nx
                
                AT((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                bT((i-1)*Nz+j)            = T2_nd;
                
                if BC_Tout == 0     % Dirichlet
                    AT((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bT((i-1)*Nz+j)            = T2_nd;
                elseif BC_Tout == 1 % Neumann
                    AT((i-1)*Nz+j,(i-1)*Nz+j)    =  1/dx_nd;
                    AT((i-1)*Nz+j,(i-1)*Nz+j-Nz) = -1/dx_nd;
                    bT((i-1)*Nz+j)               =  0;
                elseif BC_Tout == 2 % BC outlet second derivative null
                    AT((i-1)*Nz+j,(i-1)*Nz+j)      =  1/(dx_nd^2);   % T Nx    j
                    AT((i-1)*Nz+j,(i-1)*Nz+j-Nz)   = -2/(dx_nd^2);   % T Nx-1  j
                    AT((i-1)*Nz+j,(i-1)*Nz+j-2*Nz) =  1/(dx_nd^2);   % T Nx-2  j
                    bT((i-1)*Nz+j)                 =  0;
                end
                
            % BC side 1    
            elseif j == 1  
                
                if BC_side == 0     % Dirichlet
                    AT((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bT((i-1)*Nz+j)            = T2_nd;
                elseif BC_side == 1 % Neumann
                    AT((i-1)*Nz+j,(i-1)*Nz+j)   = -1/dz_nd;
                    AT((i-1)*Nz+j,(i-1)*Nz+j+1) =  1/dz_nd;
                    bT((i-1)*Nz+j)              = 0;
                end
                
            % BC side 2    
            elseif j == Nz  
                
                if BC_side == 0     % Dirichlet
                    AT((i-1)*Nz+j,(i-1)*Nz+j) = 1;
                    bT((i-1)*Nz+j)            = T2_nd;
                elseif BC_side == 1 % Neumann
                    AT((i-1)*Nz+j,(i-1)*Nz+j)   =  1/dz_nd;
                    AT((i-1)*Nz+j,(i-1)*Nz+j-1) = -1/dz_nd;
                    bT((i-1)*Nz+j)              = 0;
                end
            
            % Domain    
            else  
                
                if scheme_T == 0    % central difference (for transport terms)
                    
                    AT((i-1)*Nz+j,(i-1)*Nz+j)    = + 2*lap_kh(i,j)/(dx_nd^2) + 2*lap_kh(i,j)/(dz_nd^2)...
                                                   + du_ee_rc_dx(i,j)        + dw_ee_rc_dz(i,j);            % T i j
                    
                    AT((i-1)*Nz+j,(i-1)*Nz+j-1)  = - lap_kh(i,j)/(dz_nd^2)   - w_ee_rc(i,j)/(2*dz_nd);      % T i j-1
                    AT((i-1)*Nz+j,(i-1)*Nz+j+1)  = - lap_kh(i,j)/(dz_nd^2)   + w_ee_rc(i,j)/(2*dz_nd);      % T i j+1
                    
                    AT((i-1)*Nz+j,(i-1)*Nz+j-Nz) = - lap_kh(i,j)/(dx_nd^2)   - u_ee_rc(i,j)/(2*dx_nd);      % T i-1 j
                    AT((i-1)*Nz+j,(i-1)*Nz+j+Nz) = - lap_kh(i,j)/(dx_nd^2)   + u_ee_rc(i,j)/(2*dx_nd);      % T i+1 j
                    
                    bT((i-1)*Nz+j) = diss_x(i,j) + diss_z(i,j) - du_eeP_dx(i,j) - dw_eeP_dz(i,j);
                    
                elseif scheme_T == 1    % upwind (for transport terms)
                    
                    AT((i-1)*Nz+j,(i-1)*Nz+j)    = + 2*lap_kh(i,j)/(dx_nd^2) + 2*lap_kh(i,j)/(dz_nd^2)...
                                                   + du_ee_rc_dx(i,j)        + dw_ee_rc_dz(i,j)...
                                                   + u_ee_rc(i,j)/(dx_nd)    + w_ee_rc(i,j)/(dz_nd);        % T i j
                    
                    AT((i-1)*Nz+j,(i-1)*Nz+j-1)  = - lap_kh(i,j)/(dz_nd^2)   - w_ee_rc(i,j)/(dz_nd);        % T i j-1
                    AT((i-1)*Nz+j,(i-1)*Nz+j+1)  = - lap_kh(i,j)/(dz_nd^2);                                 % T i j+1
                    
                    AT((i-1)*Nz+j,(i-1)*Nz+j-Nz) = - lap_kh(i,j)/(dx_nd^2)   - u_ee_rc(i,j)/(dx_nd);        % T i-1 j
                    AT((i-1)*Nz+j,(i-1)*Nz+j+Nz) = - lap_kh(i,j)/(dx_nd^2);                                 % T i+1 j
                    
                    bT((i-1)*Nz+j) = diss_x(i,j) + diss_z(i,j) - du_eeP_dx(i,j) - dw_eeP_dz(i,j);
                    
                end
                
            end
            
        end
    end 
    
    % Solving Energy Equation ---------------------------------------------
    Tv = AT\bT;
    
    % remapping the temperature field
    for i = 1 : Nx
             T_nd(i,:) = Tv(((i-1)*Nz+1):(i*Nz));
    end
    
    % Equation of state (update density) ----------------------------------
    if     eos == 1 
        rho_nd = (P_nd*P_r + Psg*gamma_sg)./(Rsg.*T1_nd*T_r) / rho_r;    % rho(P)
    elseif eos == 2
        rho_nd = (P_nd*P_r + Psg*gamma_sg)./(Rsg.*T_nd*T_r)  / rho_r;    % rho(P,T)
    end
    
    % Update viscosity ----------------------------------------------------
    if mu_T == 1;
        for i = 1 : Nx
            for j = 1 : Nz
                mu_nd(i,j) = fun_mu(T_nd(i,j)*T_r) / mu_r;
            end
        end
    end
    
    % Residual update (Norm L2 of Pressure and Temperature) ---------------
    Fn = sum(sum(P_nd));         % pseudo normal force
    res_P = abs(Fn-Fn0);
    Fn0 = Fn;
    
    Tr = sum(sum(T_nd));         % pseudo mean temperature
    res_T = abs(Tr-Tr0);
    Tr0 = Tr;
    
    res_max = max([res_P res_T]);
    
    it_res(iter,1) = iter;
    it_res(iter,2) = res_P;
    it_res(iter,3) = res_T;
    
    fprintf('\nit:%d\n',iter)
    fprintf('Log residual P: %f\n',log10(res_P))
    fprintf('Log residual T: %f\n',log10(res_T))

    condest(AP)
    condest(AT)
    
end

tm = toc(ti);
fprintf('\ntime for the while cycle: tm = %f [s]\n\n',tm)

% residual plot
figure(80)
plot(it_res(1:iter,1),log10(it_res(1:iter,2)),'b'); hold on
plot(it_res(1:iter,1),log10(it_res(1:iter,3)),'r'); grid on
xlabel('iteration'); ylabel('residual'); title('Residual')
legend('Pressure','Temperature')

if sum(sum(isnan(P_nd))) >= 1 % if there is a NaN in the solution stop the code
    error('NaN in the solution! function stopped'); return
end

% Dimensionalization ------------------------------------------------------
mu   = mu_nd*mu_r;
P    = P_nd*P_r;
T    = T_nd*T_r;
rho  = rho_nd*rho_r;

dP_x = dP_x*(P_r/L_r);
dP_z = dP_z*(P_r/L_r);


% Post processing ---------------------------------------------------------

fprintf('Pressure: P_max = %f [Pa],  P_max/P_min = %f [Pa]\n',max(max(P)),max(max(P))/min(min(P)));

% geometry figure
figure(1)
surf(z,x,hh); hold on
surf(z,x,zeros(Nx,Nz))
ylabel('x'); xlabel('z'); zlabel('h')
title('Geometry')

figure(2)
surf(z,x,P)
ylabel('x'); xlabel('z'); zlabel('P')
title('Pressure distribution')

figure(7)
surf(z,x,rho)
ylabel('x'); xlabel('z'); zlabel('\rho')
title('Density distribution')

figure(9)
surf(z,x,T)
ylabel('x'); xlabel('z'); zlabel('T')
title('Temperature distribution')

% Velocity Profile --------------------------------------------------------
% iz_slice = 1; % index of z for the xy plane cut
iz_slice = floor(Nz/2); % index of z for the xy plane cut
P_xy = P(:,iz_slice);
h_xy = hh(:,iz_slice);
Ny   = 20;
y_xy = zeros(Nx,Ny);
for i = 1 : Nx 
    y_xy(i,:) = linspace(0,h_xy(i),Ny);
end
dP_xy = diff_2o_generic(x,P(:,iz_slice));
u_xy = zeros(Nx,Ny);
for i = 1 : Nx
    u_xy(i,:) = 1/(2*mu_0)*dP_xy(i)*(y_xy(i,:).^2-y_xy(i,:)*h_xy(i)) + U/h_xy(i).*y_xy(i,:);
end

if geo_type == 0 || geo_type == 2; % Dimple o Patin
    figure(3)
    ix_profile = 1;
    plot(u_xy(ix_profile,:),y_xy(ix_profile,:),'b')
    hold on
    ix_profile = floor(Nx/2);
    plot(u_xy(ix_profile,:),y_xy(ix_profile,:),'r')
    ix_profile = Nx;
    plot(u_xy(ix_profile,:),y_xy(ix_profile,:),'k')
    grid on
    title('Velocity Profile'); xlabel('u(y)'); ylabel('y')
    legend('inlet','middle','outlet')
    
elseif geo_type == 1; % Pin
    R = 0.004; % pin radius [m]
    ix_prof_1 = floor((max(x)-R)/(max(x))*(Nx/2)) + 2; % index of the first node of the pin surface
    ix_prof_2 = floor(Nx/2);
    ix_prof_3 = floor((max(x)+R)/(max(x))*(Nx/2));     % index of the last node of the pin surface
    
    figure(3)
    plot(u_xy(ix_prof_1,:),y_xy(ix_prof_1,:),'b'); hold on    
    plot(u_xy(ix_prof_2,:),y_xy(ix_prof_2,:),'r')
    plot(u_xy(ix_prof_3,:),y_xy(ix_prof_3,:),'k')
    grid on
    title('Velocity Profile over the Pin surface'); xlabel('u(y)'); ylabel('y')
    legend('Pin inlet','Pin middle','Pin outlet')
end


% Shear Stress ------------------------------------------------------------
% Slice
tau_sl = zeros(Nx,Ny);   % shear stress of the iz_slice
for i = 1 : Nx
    tau_sl(i,:) = mu_0*( 1/(2*mu_0)*dP_xy(i)*(2.*y_xy(i,:)-h_xy(i)) + U/h_xy(i) );
end 

figure(4)
plot(tau_sl(1,:),y_xy(1,:))
title('tau slice')
grid on

% Surface
tau_up = zeros(Nx,Nz);   % shear stress of the upper wall 
tau_lo = zeros(Nx,Nz);   % shear stress of the lower wall

for i = 1 : Nx
    for j = 1 : Nz
        tau_up(i,j) = mu_0*( 1/(2*mu_0)*dP_x(i,j)*( hh(i,j)) + U/hh(i,j) );
        tau_lo(i,j) = mu_0*( 1/(2*mu_0)*dP_x(i,j)*(-hh(i,j)) + U/hh(i,j) );
    end
end

figure(5)
surf(z,x,tau_up)
title('tau_u_p')
ylabel('x'); xlabel('z'); zlabel('tau_u_p')

figure(6)
surf(z,x,tau_lo)
title('tau_l_o')
ylabel('x'); xlabel('z'); zlabel('tau_l_o')

% Normal Force ------------------------------------------------------------
Int = zeros(Nz,1);
for i = 1 : Nx
    Int(i) = trapez_int_disc(z,P(i,:)-P1);
end
Fn_int = trapez_int_disc(x,Int);

fprintf('\nNormal force: % f [N]',Fn_int)

% Tangential Force --------------------------------------------------------
int1 = zeros(Nz,1); int2 = zeros(Nz,1); 
for j = 1 : Nz
    int1(j) = trapez_int_disc(x,tau_up(:,j));
    int2(j) = trapez_int_disc(x,tau_lo(:,j));
end
Ft_up = trapez_int_disc(z,int1(:));
Ft_lo = trapez_int_disc(z,int2(:));
 
Surf  = Lx*Lz;
Cf_up = Ft_up/(0.5*Surf*rho_0*U^2);
Cf_lo = Ft_lo/(0.5*Surf*rho_0*U^2);

fprintf('\nintegrated upper wall, Ft = %f [N]; skin Cf = %f\n',Ft_up, Cf_up)
fprintf(  'integrated lower wall, Ft = %f [N]; skin Cf = %f\n',Ft_lo, Cf_lo)

Cf_sys_up = Ft_up/Fn_int;
Cf_sys_lo = Ft_lo/Fn_int;

fprintf('\nFriction coefficient of the whole system (Ft_up/Fn) Cf_up: %f\n',Cf_sys_up)
fprintf(  'Friction coefficient of the whole system (Ft_lo/Fn) Cf_lo: %f\n',Cf_sys_lo)

% Analytical  Solution ----------------------------------------------------
% only for patin
if an == 1
h_an = hh(:,1)./h2;
rh = h1/h2;             % ratio h1/h2
P_adim = rh/(1-rh^2)*(1./(h_an.^2)-1/(rh^2)) - 1/(1-rh).*(1./(h_an)-1/rh);
P_ref = P1;
P_an = P_adim *6*mu_0*U*Lx/(h2^2) + P_ref;

err_an = (P_an-P(:,1))./P_an;

figure(2)
hold on
plot3(zeros(Nx,1),x,P_an,'r*')

figure(8)
plot(x,err_an)
xlabel('x'); ylabel('error'); grid on
title('error (analytic - numeric)/analytical')

end

figure(10)
surf(z,x,mu)
title('dinamic viscosity \mu')
ylabel('x'); xlabel('z'); zlabel('\mu')


figure(2)
figure(7)
figure(9)
figure(80)

% figure(10); spy(AT)

% Pin 2D visualizations ---------------------------------------------------
if geo_type  == 1

    figure(20)
    plot(x,hh(:,1),'b'); grid on; holn on
    plot(x,zeros(Nx,1),'b');
    xlabel('x'); ylabel('h')
    title('Pin Geometry')
    
    figure(21)
    plot(x,T(:,1),'b'); grid on
    xlabel('x'); ylabel('T')
    title('Temperature Distribution')
    
    figure(22)
    plot(x,rho(:,1),'b'); grid on
    xlabel('x'); ylabel('\rho')
    title('Density Distribution')
    
    figure(23)
    plot(x,P(:,1),'b'); grid on
    xlabel('x'); ylabel('P')
    title('Pressure Distribution')
    
end


% comparison with SPARC ---------------------------------------------------

% Channel
if cfr_sparc == 1
    
    % the sparc post processing data are created and saved in ~/Tribology/channel_box/Post_proc/
    % with the script Box_2D_post
    load(sprintf('/net/istmauriga/localhome/hi217/Tribology/channel_box/Post_proc/sim_channel_SPARC_%s',file_name))
    
    figure(21)
    plot(x,T(:,1),'b')
    hold on; grid on
    plot(x_mez+max(x_mez),T_mez,'r')
    xlabel('x'); ylabel('T')
    title('Comparison of Temperature Distribution')
    legend('Rey comp','SPARC')
    
    figure(20)
    plot(x,rho(:,1),'b')
    hold on; grid on
    plot(x_mez+max(x_mez),rho_mez,'r')
    xlabel('x'); ylabel('\rho')
    title('Comparison of Density Distribution')
    legend('Rey comp','SPARC')
    
    figure(22)
    plot(x,P(:,1),'b')
    hold on; grid on
    plot(x_mez+max(x_mez),P_mez,'r')
    xlabel('x'); ylabel('P')
    title('Comparison of Pressure Distribution')
    legend('Rey comp','SPARC')
    
%     % for presentation
%     figure(30)
%     plot(x,T(:,1),'b','linewidth',3)
%     hold on; grid on
%     plot(x_mez+max(x_mez),T_mez,'r','linewidth',3)
% %     xlabel('x'); ylabel('T')
% %     title('Comparison of Temperature Distribution')
%     legend('Rey comp','SPARC')
%     
%     figure(31)
%     plot(x,rho(:,1),'b','linewidth',3)
%     hold on; grid on
%     plot(x_mez+max(x_mez),rho_mez,'r','linewidth',3)
% %     xlabel('x'); ylabel('\rho')
% %     title('Comparison of Density Distribution')
%     legend('Rey comp','SPARC')
%     
%     figure(32)
%     plot(x,P(:,1),'b','linewidth',3)
%     hold on; grid on
%     plot(x_mez+max(x_mez),P_mez,'r','linewidth',3)
% %     xlabel('x'); ylabel('P')
% %     title('Comparison of Pressure Distribution')
%     legend('Rey comp','SPARC')
    
end

% Slide Bearing
if cfr_sparc == 2
    
    % the sparc post processing data are created and saved in ~/Tribology/channel_box/Post_proc/
    % with the script Box_2D_post
    load(sprintf('/net/istmauriga/localhome/hi217/Tribology/Inclined_plate/Comp_alpha_%s',file_name))
    
    Pref = P2;
    alpha_Rey = mu_0*Lx*U/(Pref*h2^2);
    fprintf('\nalpha Reynolds: %f, alpha Sparc: %f\n\n',alpha_Rey,alpha_num)
    
    figure(21)
    plot(x,T(:,1),'b')
    hold on; grid on
    plot(x_mez+max(x_mez),T_mez,'r')
    xlabel('x'); ylabel('T')
    title('Comparison of Temperature Distribution')
    legend('Rey comp','SPARC')
    
    figure(20)
    plot(x,rho(:,1),'b')
    hold on; grid on
    plot(x_mez+max(x_mez),rho_mez,'r')
    xlabel('x'); ylabel('\rho')
    title('Comparison of Density Distribution')
    legend('Rey comp','SPARC')
    
    figure(22)
    plot(x,P(:,1),'b')
    hold on; grid on
    plot(x_mez+max(x_mez),P_mez,'r')
    plot(x,P_an,'k')
    xlabel('x'); ylabel('P')
    title('Comparison of Pressure Distribution')
    legend('Rey comp','SPARC','incompressible')
    
    
    
end
