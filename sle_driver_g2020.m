% SLE calculation driver
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

clear

% Choose test configuration
config = 1; % 1 = path dependence | 2 = perturbation | 3 = pump | 
            % 4 = transition | 5 = no iso | 6,7,8 = deglaciation

% define constants
params.rho_ice = 917; % kg/m^3 
params.rho_ocean = 1028; % kg/m^3 
params.rho_water = 1000; % kg/m^3 
%params.Aoc = 3.618e14; % m^2
% Schematic case
params.Aoc = 1; 

if config==1
    % 1) define configuration with path dependence
    THICK = [1, 1, 0, 0];
    BED = [0, -0.892, -0.892, 0];
elseif config==2
    % 2) define configuration  with minor perturbation
    THICK = [1, 1, 1, 1];
    BED = [0, -0.893, -0.892, 0];
elseif config==3
    % 3) define configuration with sea-level pump
    THICK = [1, 1, 0, 0, 1, 1, 0, 0, 1];
    BED = [0, -0.892, -0.892, 0, 0, -0.892, -0.892, 0, 0];
elseif config==4
    % 4) define configuration with gr-fl transtion
    THICK = [1, 1, 0];
    %BED = [0, -0.893, -0.893];
    BED = [0, -2, -2];
elseif config==5
    % 5) define configuration with constant bed
    THICK = [2, 0.5, 1.5];
    BED = [-1, -1, -1];
elseif config==6
    % 6) define configuration with deglaciation
    THICK = [3, 2.5, 1.8, 0.7,  0];
    BED = [-2, -2, -1.5, -1, -0.5 ];
elseif config==7
    % 7) define configuration with bed transitioning to above
    THICK = [3, 2.5, 1.8, 0.7,  0.5, 0.2];
    BED = [-2, -2, -1.5, -1, 0.5, 1 ];
elseif config==8
    % 8) define configuration with bed below
    THICK = [3, 2.5, 1.8, 0.7,  0.5, 0.2];
    BED = [-2, -2, -1.5, -1, -0.7, -0.5 ];
else
    disp('unknown config')
    return
end

% Constructiong consistent configuration 
SURFACEg = BED+THICK;
SURFACEf = THICK*(1-params.rho_ice/params.rho_ocean);
% Grounded ice masks. Equivalent to Equation 5.  
F = THICK + params.rho_ocean/params.rho_ice*BED;
GROUND_MASK = F; 
GROUND_MASK(GROUND_MASK<0) = 0; 
GROUND_MASK(GROUND_MASK>0) = 1;
% Surface and base 
SURFACE = SURFACEg.*GROUND_MASK + SURFACEf.*(1-GROUND_MASK); 
BASE = SURFACE-THICK;

nt = length(BED)-1;
sle_g2020 = zeros(1,nt);
% step through problem
for n = 1:nt
    sle_g2020(1,n) = g2020_func(BED(1,n:(n+1)),BASE(1,n:(n+1)),SURFACE(1,n:(n+1)),params);
end
sle_step_g2020 = sum(sle_g2020);

% leap through problem from t0 to tend
sle_leap_g2020 = g2020_func(BED(1,[1,end]),BASE(1,[1,end]),SURFACE(1,[1,end]),params);

SURFACE
THICK
BASE
BED
GROUND_MASK

% Stepwise change
sle_g2020

% compare step and leap
[sle_step_g2020, sle_leap_g2020]

% plot configuration
plot_columns
