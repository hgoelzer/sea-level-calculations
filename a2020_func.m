% SLE calculation function
% Implementinig Adhikari et al., 2020
% https://tc.copernicus.org/articles/14/2819/2020/
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

function [sle] = a2020 (BED,BASE,SURFACE,params)
% Expect pairwise definitions of geometry for t0 and t1
% dim(vars) = [1,2]
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

% Change in ice thickness relative to t0
DEL_H = THICK - THICK(:,1); 

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
DEL_HF = HF - HF(:,1); % Equation 9. 

% A2020 method of computing sea-level contribution from ice sheets. 
GG = GROUND_MASK(:,1).*GROUND_MASK;
DEL_HM = DEL_H.*GG + DEL_HF.*(1-GG); % Equation 11. 
DEL_HV = (1-rho_water/rho_ocean)*(DEL_H - DEL_HF).*(1-GG); % Equation 12. 
DEL_HS = DEL_HM(:,2) + DEL_HV(:,2); % Equation 10. 

% Calculate sea-level contribution
sle = -(rho_ice/rho_water) * DEL_HS / Aoc; % text end of 2.3 
