% SL calculation function
% Implementinig Huybrechts 1990
% Heiko Goelzer (heig@norceresearch.no), Feb 2023

function [slc] = h1990(BED,BASE,SURFACE,params)
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


% Level set function.
F = THICK + rho_ocean/rho_ice*BED; % Equation 1. 

% Grounded ice masks. Equivalent to Equation 5.  
%GROUND_MASK = F; 
GROUND_MASK = (-rho_ice/rho_ocean*THICK) < BED;
GROUND_MASK(GROUND_MASK<0) = 0; 
GROUND_MASK(GROUND_MASK>0) = 1;	

nn = size(BED,1);
VOLTOT = zeros(nn,1);
for n = 1:nn
    if (GROUND_MASK(n,1) ~= 0)
        VOLTOT(n) = THICK(n,1);
    end
end
%VOLTOT = GROUND_MASK(:,2).*THICK;
VOL2 = VOLTOT;
VOLTOT = zeros(nn,1);
VOLCORR = zeros(nn,1);
for n = 1:nn
    if (GROUND_MASK(n,2) ~= 0)
        VOLTOT(n) = THICK(n,2);
    end
    VOLCOR2 = -BED(n,1)*rho_ocean/rho_ice;
    if (BED(n,1) > 0)
        VOLCOR2 = 0;
    end
    if (GROUND_MASK(n,1)==0 && GROUND_MASK(n,2)~=0) 
        VOLCORR(n)=-VOLCOR2;
    elseif (GROUND_MASK(n,1)~=0 && GROUND_MASK(n,2)==0)
        VOLCORR(n)=VOLCOR2;
    end
end

%VOLTOT = GROUND_MASK(:,2).*THICK;
%
%VOLCOR2 = -BED(:,1)*rho_ocean/rho_ice;
%if(HBEDIN(I,J).GT.0)
%    VOLCOR2=0
%end
%if(IGRORIG(I,J).EQ.0.AND.IGR(I,J).NE.0)THEN
%    VOLCORR=VOLCORR-VOLCOR2
%elseif(IGRORIG(I,J).NE.0.AND.IGR(I,J).EQ.0)THEN
%    VOLCORR=VOLCORR+VOLCOR2
%end     

% Calculate sea-level contribution
%sle = (VOLTOT+VOLCORR)*(-2.512e-15)
slc = -(VOLTOT+VOLCORR-VOL2)*rho_ice/rho_water*Aoc;
