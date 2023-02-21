% SLE calculation driver
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

clear

% Choose test configuration
config = 1; % 1 = path dependence | 2 = perturbation | 3 = pump | 
            % 4 = no iso

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
    F = THICK + params.rho_ocean/params.rho_ice*BED;
elseif config==2
    % 2) define configuration  with minor perturbation
    THICK = [1, 1, 0, 0];
    BED = [0, -0.893, -0.893, 0];
    F = THICK + params.rho_ocean/params.rho_ice*BED;
elseif config==3
    % 3) define configuration with sea-level pump
    THICK = [1, 1, 0, 0, 1, 1, 0, 0, 1];
    BED = [0, -0.892, -0.892, 0, 0, -0.892, -0.892, 0, 0];
    F = THICK + params.rho_ocean/params.rho_ice*BED;
elseif config==4
    % 4) define configuration with constant bed
    THICK = [2, 1, 1.5];
    BED = [-1, -1, -1];
    F = THICK + params.rho_ocean/params.rho_ice*BED;
else
    disp('unknown config')
    return
end

% Constructiong consistent configuration 
SURFACEg = BED+THICK;
SURFACEf = THICK*(1-params.rho_ice/params.rho_ocean);
% Grounded ice masks. Equivalent to Equation 5.  
GROUND_MASK = F; 
GROUND_MASK(GROUND_MASK<0) = 0; 
GROUND_MASK(GROUND_MASK>0) = 1;
% Surface and base 
SURFACE = SURFACEg.*GROUND_MASK + SURFACEf.*(1-GROUND_MASK); 
BASE = SURFACE-THICK;

nt = length(BED)-1;
sle_a2020 = zeros(1,nt);
% step through problem
for n = 1:nt
    sle_a2020(1,n) = a2020_func(BED(1,n:(n+1)),BASE(1,n:(n+1)),SURFACE(1,n:(n+1)),params);
end
sle_step = sum(sle_a2020);

% leap through problem from t0 to tend
sle_leap = a2020_func(BED(1,[1,end]),BASE(1,[1,end]),SURFACE(1,[1,end]),params);

SURFACE
THICK
BASE
BED
GROUND_MASK

% Stepwise change
sle_a2020

% compare step and leap
[sle_step, sle_leap]

