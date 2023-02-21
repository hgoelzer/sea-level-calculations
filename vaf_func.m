% SL calculation function
% Implementinig Vaf method 
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

function [slc] = vaf (BED,BASE,SURFACE,params)
% Expect pairwise definitions of geometry for t0 and t1
% dim(vars) = [:,2]
% BASE => Base of the ice.
% BED => Bedrock elevation. 
% SURFACE => Surface of the ice. 

% Check input
if (size(BED,2) ~= 2 )
    disp('Error in variable size')
    return
end

% use provided constants.
rho_ice = params.rho_ice; % kg/m^3 
rho_ocean = params.rho_ocean; % kg/m^3 
rho_water = params.rho_water; % kg/m^3 
Aoc = params.Aoc; % m^2

% compute ice thickness 
THICK = SURFACE - BASE; 

% Level set function.
F = THICK + rho_ocean/rho_ice*BED; % Equation 1. 

% Grounded ice masks. Equivalent to Equation 5.  
GROUND_MASK = F; 
GROUND_MASK(GROUND_MASK<0) = 0; 
GROUND_MASK(GROUND_MASK>0) = 1;	

% Floatation height and the height above flotation (HAF) 
H0 = -rho_ocean/rho_ice * min(BED,0); % Equation 7. 
HF = GROUND_MASK.*(THICK - H0); % Equation 8.  

% The HAF method. 
% Change in HAF 
DEL_HF = HF(:,2) - HF(:,1); % Equation 9. 

% Calculate sea-level contribution
slc = -(rho_ice/rho_water) * DEL_HF / Aoc; 
