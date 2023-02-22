% SL calculation driver to compare Vaf, H1990, A2020 and G2020
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

clear

% Choose test configuration
config = 1;
% 1 = path dependence | 2 = perturbation | 3 = pump |
% 4 = transition | 5 = no iso | 6,7,8 = deglaciation |
% 9 = 2d case (n,t) | 10 = vaf test | 11 = ISMIP6

%% Settings
% Verbose mode
verbflg = 0;
% Plotting mode
pltflg = 0; 
select = 1; % Select grid cell to plot 

%% More settings, may be overwritten with config specific settings
% Output mode. Sum over grid cells, e.g. for real world examples
sumflg = 0;
res = 1;
% Scaling for 0=schematic vs 1=real world setups
sclflg = 0;

% define constants
params.rho_ice = 917; % kg/m^3 
params.rho_ocean = 1028; % kg/m^3 
params.rho_water = 1000; % kg/m^3 
AF2 = 1;

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
elseif config==10
    % 10) define configuration for vaf test
    THICK = [3, 3, 0.7,  0.5, 0.2];
    BED = [-2, -1.5, -1, 0.5, 1 ];
elseif config==11
    % 11) define configuration from ISMIP6 run
    prep_ISMIP6
    sumflg = 1;
    sclflg = 1;
    res = 16000;
else
    disp('unknown config')
    return
end

if sclflg
    % Real world case
    params.Aoc = 3.625e14; % m^2
    params.res = res;
else
    % Schematic case
    params.Aoc = 1;
    params.res = 1;
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
slc_vaf = zeros(nc,nt);
slc_h1990 = zeros(nc,nt);
slc_a2020 = zeros(nc,nt);
slc_g2020 = zeros(nc,nt);
% step through problem
for n = 1:nt
    slc_vaf(:,n) = vaf_func(BED(:,n:(n+1)),BASE(:,n:(n+1)),SURFACE(:,n:(n+1)),params);
    slc_h1990(:,n) = h1990_func(BED(:,n:(n+1)),BASE(:,n:(n+1)),SURFACE(:,n:(n+1)),params);
    slc_a2020(:,n) = a2020_func(BED(:,n:(n+1)),BASE(:,n:(n+1)),SURFACE(:,n:(n+1)),params);
    slc_g2020(:,n) = g2020_func(BED(:,n:(n+1)),BASE(:,n:(n+1)),SURFACE(:,n:(n+1)),params);
end
slc_step_vaf = sum(slc_vaf,2);
slc_step_h1990 = sum(slc_h1990,2);
slc_step_a2020 = sum(slc_a2020,2);
slc_step_g2020 = sum(slc_g2020,2);

% leap through problem from t0 to tend
slc_leap_vaf = vaf_func(BED(:,[1,end]),BASE(:,[1,end]),SURFACE(:,[1,end]),params);
slc_leap_h1990 = h1990_func(BED(:,[1,end]),BASE(:,[1,end]),SURFACE(:,[1,end]),params);
slc_leap_a2020 = a2020_func(BED(:,[1,end]),BASE(:,[1,end]),SURFACE(:,[1,end]),params);
slc_leap_g2020 = g2020_func(BED(:,[1,end]),BASE(:,[1,end]),SURFACE(:,[1,end]),params);

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
    vaf_step_leap = [params.res^2*sum(slc_step_vaf.*AF2,1), params.res^2*sum(slc_leap_vaf.*AF2,1)]
    h1990_step_leap = [params.res^2*sum(slc_step_h1990.*AF2,1), params.res^2*sum(slc_leap_h1990.*AF2,1)]
    a2020_step_leap = [params.res^2*sum(slc_step_a2020.*AF2,1), params.res^2*sum(slc_leap_a2020.*AF2,1)]
    g2020_step_leap = [params.res^2*sum(slc_step_g2020.*AF2,1), params.res^2*sum(slc_leap_g2020.*AF2,1)]
    % compare all transitions
    sum(slc_vaf.*AF2,1)*params.res^2
    sum(slc_h1990.*AF2,1)*params.res^2
    sum(slc_a2020.*AF2,1)*params.res^2
    sum(slc_g2020.*AF2,1)*params.res^2
else
    % compare step and leap
    vaf_step_leap = [slc_step_vaf, slc_leap_vaf]
    h1990_step_leap = [slc_step_h1990, slc_leap_h1990]
    a2020_step_leap = [slc_step_a2020, slc_leap_a2020]
    g2020_step_leap = [slc_step_g2020, slc_leap_g2020]
    % compare all transitions
    slc_vaf
    slc_h1990
    slc_a2020
    slc_g2020
end

% plot configuration
if pltflg
    plot_columns(BED,BASE,SURFACE,select)
end
