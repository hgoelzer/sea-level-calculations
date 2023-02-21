% SL calculation function
% Implementinig Goelzer et al., 2020,
% https://tc.copernicus.org/articles/14/833/2020/
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

function [slc_corr] = g2020 (BED,BASE,SURFACE,params)
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

% no RSL in this example
SLR = BED*0;

THICK = SURFACE-BASE;

% Volume above floatation Equ 1 in corrigendum !!  
Vaf = max(0,THICK+min(BED-SLR,0)*(rho_ocean/rho_ice));

% All following 3 quantities in ocean water volume
Vaf_owv = Vaf*(rho_ice/rho_ocean); % Like Equ 2, but not devided by Aoc
% Potential ocean volume
Vpov = max(0,-BED); % Equ 8.
% Density correction for transformation from ice to freshwater 
Vden = THICK*((rho_ice/rho_water)-(rho_ice/rho_ocean)); % Equ 10.

% SL components
slc_af =  -(Vaf_owv(:,2) - Vaf_owv(:,1))/Aoc; % Equ 2 and Equ 3
slc_pov = -(Vpov(:,2) - Vpov(:,1))/Aoc; % Equ 9.
slc_den = -(Vden(:,2) - Vden(:,1))/Aoc; % Equ 11.

% Corrected ocean water volume contribution
slc_corr = slc_af + slc_pov + slc_den;
