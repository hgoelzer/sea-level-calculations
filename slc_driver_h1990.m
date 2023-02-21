% SL calculation driver for H1990
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

clear

% Choose test configuration
config = 9;
% 1 = path dependence | 2 = perturbation | 3 = pump |
% 4 = transition | 5 = no iso | 6,7,8 = deglaciation |
% 9 = 2d case (n,t)

% Verbose mode
verbflg = 0;
% Plotting mode
pltflg = 1; 
select = 1; % Select grid cell to plot 
% Output mode. Sum over grid cells, e.g. for complex examples
sumflg = 1;
% Scaling for 0=schematic vs 1=real world setups
sclflg = 0;

% define constants
params.rho_ice = 917; % kg/m^3 
params.rho_ocean = 1028; % kg/m^3 
params.rho_water = 1000; % kg/m^3 
if sclflg
    % Real world case
    params.Aoc = 3.625e14; % m^2
else
    % Schematic case
    params.Aoc = 1; 
end

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
elseif config==9
    % 9) define 2D configuration (grid index,time) 
    THICK = [3, 2.5, 1.8, 0.7,  0.5, 0.2; 3, 2.5, 1.8, 0.7,  0.5, 0.2];
    BED = [-2, -2, -1.5, -1, -0.7, -0.5; -2, -2, -1.5, -1, 0.5, 1  ];
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

nt = size(BED,2)-1;
nc = size(BED,1);
slc_h1990 = zeros(nc,nt);
% step through problem
for n = 1:nt
    slc_h1990(:,n) = h1990_func(BED(:,n:(n+1)),BASE(:,n:(n+1)),SURFACE(:,n:(n+1)),params);
end
slc_step_h1990 = sum(slc_h1990,2);

% leap through problem from t0 to tend
slc_leap_h1990 = h1990_func(BED(:,[1,end]),BASE(:,[1,end]),SURFACE(:,[1,end]),params);

% output configuration
if verbflg
    SURFACE
    THICK
    BASE
    BED
    GROUND_MASK
end

% Output results
if sumflg
    % sum over all grid cells
    % compare step and leap
    h1990_step_leap = [sum(slc_step_h1990,1), sum(slc_leap_h1990,1)]
    % compare all transitions
    sum(slc_h1990,1)
else
    % compare step and leap
    h1990_step_leap = [slc_step_h1990, slc_leap_h1990]
    % compare all transitions
    slc_h1990
end

% plot configuration
if pltflg
    plot_columns(BED,BASE,SURFACE,select)
end
