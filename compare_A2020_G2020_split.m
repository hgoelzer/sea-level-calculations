% Compare SL calculation for A2020 and G2020
% Violaine Coulon March 2023 (violaine.coulon@ulb.be)
clear 

% Read ISMIP6-AIS-2300 output. In this case with bedrock changes
filename=['AIS_ULB_fETISh-KoriBU2_expAE14']; % ISMIP6 result to load
H_snapshots0 = double(ncread(['../Data/ISMIP6/ULB/lithk_',filename,'.nc'],'lithk')); % adapt path if needed
B_snapshots0 = double(ncread(['../Data/ISMIP6/ULB/topg_',filename,'.nc'],'topg')); % adapt path if needed
H_snapshots0(isnan(H_snapshots0))=0; % remove potential nan values
B_snapshots0(isnan(B_snapshots0))=0; % remove potential nan values
time=2015:1:2300;
delta=16.e3; % spatial resolution from experiment -- adapt if needed

% Define split time
tsplit=2150;

% Reference time
tref = 1;

% define some constants.
rho_ice = 917;		% kg/m^3 
rho_ocean = 1027; % kg/m^3 
rho_water = 1000; % kg/m^3 
Aoc=3.618e14;

%% SPLIT 1

time1=time(1):tsplit;
H_snapshots=H_snapshots0(:,:,1:find(time==tsplit));
B_snapshots=B_snapshots0(:,:,1:find(time==tsplit));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADHIKARI 2020 METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute ice thickness + define other variables
THICK = H_snapshots;
BED = B_snapshots;

% Change in ice thickness relative to initial ice thickness.
DEL_H = THICK - THICK(:,:,tref);

% Level set function. In the present context, it is sufficent to consider equations 1 and 2...
% ...(rather than equations 3 and 4) to define the ocean and land domains, and their interfaces.
F = THICK + rho_ocean/rho_ice*BED; % Eq. 1 of the paper.

% Grounded ice masks. Equivalent to Equation 5.
GROUND_MASK = F;
GROUND_MASK(GROUND_MASK<0) = 0;
GROUND_MASK(GROUND_MASK>0) = 1;

% Floatation height and the height above flotation (HAF)
H0 = -rho_ocean/rho_ice * min(BED,0);	% Equation 7.
HF = GROUND_MASK.*(THICK - H0);			% Equation 8.

% The HAF method.
% Change in HAF relative to AD 2000
DEL_HF = HF - HF(:,:,tref);	% Equation 9.

% New method of computing sea-level contribution from ice sheets.
GG = GROUND_MASK(:,:,tref).*GROUND_MASK;
DEL_HM = DEL_H.*GG + DEL_HF.*(1-GG);	% Equation 11.
DEL_HV = (1-rho_water/rho_ocean)*(DEL_H - DEL_HF).*(1-GG);	% Equation 12.
DEL_HS = DEL_HM + DEL_HV;	% Equation 10.

%%% Conversion in SLE m:

DEL_HM_sle = -squeeze(sum(DEL_HM,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HV_sle = -squeeze(sum(DEL_HV,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HS_sle = -squeeze(sum(DEL_HS,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HF_sle = -squeeze(sum(DEL_HF-DEL_HF(:,:,1),[1 2])*delta^2*rho_ice/(rho_ocean*Aoc));

%%% Calculate sea-level contribution

del_hm_sle1 = (DEL_HM_sle-DEL_HM_sle(1));
del_hv_sle1 = (DEL_HV_sle-DEL_HV_sle(1));
del_hs_sle1 = (DEL_HS_sle-DEL_HS_sle(1));
del_hf_sle1 = (DEL_HF_sle-DEL_HF_sle(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GOELZER 2020 METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SLref=0;

% VAF variation (SL equivalent in ocean water)
VAF=max(0,THICK+min(BED-SLref,0)*(rho_ocean/rho_ice))*rho_ice/(rho_ocean);
% Potential Ocean volume variation (SL equivalent in ocean water)
POV=max(0,SLref-BED);
% Density correction for transformation from ice to freshwater and not ocean water (SL equivalent in ocean water)
DEN=THICK*((rho_ice/rho_water)-(rho_ice/rho_ocean));
SLC=VAF+POV+DEN; 

%%% Conversion in SLE m:
SLC_AF=-squeeze(sum(VAF-VAF(:,:,tref),[1 2])*delta^2)./Aoc;
SLC_POV=-squeeze(sum(POV-POV(:,:,tref),[1 2])*delta^2.)./Aoc;
SLC_DEN=-squeeze(sum(DEN-DEN(:,:,tref),[1 2])*delta^2.)./Aoc;
SLC_CORR=SLC_AF+SLC_POV+SLC_DEN; 

%%% Sea-level contribution from modelled ice-sheet
slc_af1=SLC_AF-SLC_AF(1);
slc_pov1=SLC_POV-SLC_POV(1);
slc_den1=SLC_DEN-SLC_DEN(1);
slc_corr1=SLC_CORR-SLC_CORR(1);

%% SPLIT 2

time2=tsplit:time(end);
H_snapshots=H_snapshots0(:,:,find(time==tsplit):end);
B_snapshots=B_snapshots0(:,:,find(time==tsplit):end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADHIKARI 2020 METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute ice thickness + define other variables
THICK = H_snapshots;
BED = B_snapshots;

% Change in ice thickness relative to initial ice thickness.
DEL_H = THICK - THICK(:,:,tref);

% Level set function. In the present context, it is sufficent to consider equations 1 and 2...
% ...(rather than equations 3 and 4) to define the ocean and land domains, and their interfaces.
F = THICK + rho_ocean/rho_ice*BED; % Eq. 1 of the paper.

% Grounded ice masks. Equivalent to Equation 5.
GROUND_MASK = F;
GROUND_MASK(GROUND_MASK<0) = 0;
GROUND_MASK(GROUND_MASK>0) = 1;

% Floatation height and the height above flotation (HAF)
H0 = -rho_ocean/rho_ice * min(BED,0);	% Equation 7.
HF = GROUND_MASK.*(THICK - H0);			% Equation 8.

% The HAF method.
% Change in HAF relative to AD 2000
DEL_HF = HF - HF(:,:,tref);	% Equation 9.

% New method of computing sea-level contribution from ice sheets.
GG = GROUND_MASK(:,:,tref).*GROUND_MASK;
DEL_HM = DEL_H.*GG + DEL_HF.*(1-GG);	% Equation 11.
DEL_HV = (1-rho_water/rho_ocean)*(DEL_H - DEL_HF).*(1-GG);	% Equation 12.
DEL_HS = DEL_HM + DEL_HV;	% Equation 10.

%%% Conversion in SLE m:

DEL_HM_sle = -squeeze(sum(DEL_HM,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HV_sle = -squeeze(sum(DEL_HV,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HS_sle = -squeeze(sum(DEL_HS,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HF_sle = -squeeze(sum(DEL_HF-DEL_HF(:,:,1),[1 2])*delta^2*rho_ice/(rho_ocean*Aoc));

%%% Calculate sea-level contribution

del_hm_sle2 = (DEL_HM_sle-DEL_HM_sle(1));
del_hv_sle2 = (DEL_HV_sle-DEL_HV_sle(1));
del_hs_sle2 = (DEL_HS_sle-DEL_HS_sle(1));
del_hf_sle2 = (DEL_HF_sle-DEL_HF_sle(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GOELZER 2020 METHOD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SLref=0;

% VAF variation (SL equivalent in ocean water)
VAF=max(0,THICK+min(BED-SLref,0)*(rho_ocean/rho_ice))*rho_ice/(rho_ocean);
% Potential Ocean volume variation (SL equivalent in ocean water)
POV=max(0,SLref-BED);
% Density correction for transformation from ice to freshwater and not ocean water (SL equivalent in ocean water)
DEN=THICK*((rho_ice/rho_water)-(rho_ice/rho_ocean));
SLC=VAF+POV+DEN; 

%%% Conversion in SLE m:
SLC_AF=-squeeze(sum(VAF-VAF(:,:,tref),[1 2])*delta^2)./Aoc;
SLC_POV=-squeeze(sum(POV-POV(:,:,tref),[1 2])*delta^2.)./Aoc;
SLC_DEN=-squeeze(sum(DEN-DEN(:,:,tref),[1 2])*delta^2.)./Aoc;
SLC_CORR=SLC_AF+SLC_POV+SLC_DEN; 

%%% Sea-level contribution from modelled ice-sheet
slc_af2=SLC_AF-SLC_AF(1);
slc_pov2=SLC_POV-SLC_POV(1);
slc_den2=SLC_DEN-SLC_DEN(1);
slc_corr2=SLC_CORR-SLC_CORR(1);

%% ADD SPLIT 1 and SPLIT 2

disp(['tsplit: ',int2str(tsplit)])
disp(['VAF A2020: ',num2str(del_hf_sle1(end)+del_hf_sle2(end),6),' m'])
disp(['DEL_HS A2020: ',num2str(del_hs_sle1(end)+del_hs_sle2(end),6),' m'])
disp(['VAF G2020: ',num2str(slc_af1(end)+slc_af2(end),6),' m'])
disp(['SLC G2020: ',num2str(slc_corr1(end)+slc_corr2(end),6),' m'])
