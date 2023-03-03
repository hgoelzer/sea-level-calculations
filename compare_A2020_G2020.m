% Compare SL calculation for A2020 and G2020
% Violaine Coulon March 2023 (violaine.coulon@ulb.be)
clear 

% Read ISMIP6-AIS-2300 output. In this case with bedrock changes
filename=['AIS_ULB_fETISh-KoriBU2_expAE14']; % ISMIP6 result to load
H_snapshots = double(ncread(['../Data/ISMIP6/ULB/lithk_',filename,'.nc'],'lithk')); % adapt path if needed
B_snapshots = double(ncread(['../Data/ISMIP6/ULB/topg_',filename,'.nc'],'topg')); % adapt path if needed
H_snapshots(isnan(H_snapshots))=0; % remove potential nan values
B_snapshots(isnan(B_snapshots))=0; % remove potential nan values
time=2015:1:2300;
delta=16.e3; % spatial resolution from experiment -- adapt if needed

% Reference time (change reference time to test path dependency)
tref = 1;
% tref = 100;
% tref = 200;

% define some constants.
rho_ice = 917; % kg/m^3 
rho_ocean = 1027; % kg/m^3 
rho_water = 1000; % kg/m^3 
Aoc=3.618e14;

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
DEL_HM_sle = squeeze(sum(DEL_HM,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HV_sle = squeeze(sum(DEL_HV,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HS_sle = squeeze(sum(DEL_HS,[1 2])*delta^2*rho_ice/(rho_water*Aoc));
DEL_HF_sle = squeeze(sum(DEL_HF-DEL_HF(:,:,1),[1 2])*delta^2*rho_ice/(rho_ocean*Aoc));

%%% Sea-level contribution from modelled ice-sheet
del_hm_sle = -(DEL_HM_sle-DEL_HM_sle(1));
del_hv_sle = -(DEL_HV_sle-DEL_HV_sle(1));
del_hs_sle = -(DEL_HS_sle-DEL_HS_sle(1));
del_hf_sle = -(DEL_HF_sle-DEL_HF_sle(1));

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
% G2020 corrected sea-level contribution
SLC=VAF+POV+DEN; 

%%% Conversion in SLE m:
SLC_AF=-squeeze(sum(VAF-VAF(:,:,tref),[1 2])*delta^2)./Aoc;
SLC_POV=-squeeze(sum(POV-POV(:,:,tref),[1 2])*delta^2.)./Aoc;
SLC_DEN=-squeeze(sum(DEN-DEN(:,:,tref),[1 2])*delta^2.)./Aoc;
SLC_CORR=SLC_AF+SLC_POV+SLC_DEN; 

%%% Sea-level contribution from modelled ice-sheet
slc_af=SLC_AF-SLC_AF(1);
slc_pov=SLC_POV-SLC_POV(1);
slc_den=SLC_DEN-SLC_DEN(1);
slc_corr=SLC_CORR-SLC_CORR(1);

disp(['tref: ',int2str(tref)])
disp(['VAF A2020: ',num2str(del_hf_sle(end),6),' m'])
disp(['DEL_HS A2020: ',num2str(del_hs_sle(end),6),' m'])
disp(['VAF G2020: ',num2str(slc_af(end),6),' m'])
disp(['SLC G2020: ',num2str(slc_corr(end),6),' m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT DIFFERENCE MAP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute maps of sea-level contributions for A2020 and G2020 in ocean water equivalent
DEL_HS_map = -DEL_HS*rho_ice/(rho_water);
DEL_HS_map = DEL_HS_map-DEL_HS_map(:,:,1);
SLC_map=-(SLC-SLC(:,:,tref));
SLC_map=SLC_map-SLC_map(:,:,1);

t=286; % change to test another year
figure;
imagesc(SLC_map(:,:,t)'-DEL_HS_map(:,:,t)',[0 20])
%cmap = crameri('bilbao');
%colormap(cmap)
cm = colormap('pink');
colormap(flipud(cm))
hold on
contour(GROUND_MASK(:,:,end)',1,'k')
title(['SLC (G2020) - DEL_{HS} (A2020) at year ',int2str(time(t))])
axis off
axis xy
axis equal
c=colorbar;
c.Label.String = 'meters ocean equivalent';
c.FontSize = 12;
c.Label.FontSize = 12;

%function cmap = crameri(ColormapName,varargin) 
%% crameri returns perceptually-uniform scientific colormaps created
%% by Fabio Crameri. 
%% 
%%% Syntax 
%% 
%%  crameri 
%%  cmap = crameri('ColormapName') 
%%  cmap = crameri('-ColormapName') 
%%  cmap = crameri(...,NLevels)
%%  cmap = crameri(...,'pivot',PivotValue) 
%%  crameri(...)
%% 
%%% Description 
%% 
%% crameri without any inputs displays the options for colormaps. 
%% 
%% cmap = crameri('ColormapName') returns a 256x3 colormap.  For a visual
%% depiction of valid colormap names, type |crameri|. 
%%
%% cmap = crameri('-ColormapName') a minus sign preceeding any ColormapName flips the
%% order of the colormap. 
%%
%% cmap = crameri(...,NLevels) specifies a number of levels in the colormap.  Default
%% value is 256. 
%%
%% cmap = crameri(...,'pivot',PivotValue) centers a diverging colormap such that white 
%% corresponds to a given value and maximum extents are set using current caxis limits. 
%% If no PivotValue is set, 0 is assumed. 
%%
%% crameri(...) without any outputs sets the current colormap to the current axes.  
%% 
%%% Examples 
%% For examples, type: 
%% 
%%  showdemo crameri_documentation
%%
%%% Author Info 
%% This function was written by Chad A. Greene of the University of Texas
%% Institute for Geophysics (UTIG), August 2018, using Fabio Crameri's 
%% scientific colormaps, version 4.0. http://www.fabiocrameri.ch/colourmaps.php
%% 
%%% Citing this colormap: 
%% Please acknowledge the free use of these colormaps by citing
%% 
%% Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
%% 
%% Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and 
%% StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018.
%% 
%% For more on choosing effective and accurate colormaps for science, be sure
%% to enjoy this fine beach reading: 
%% 
%% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True 
%% colors of oceanography: Guidelines for effective and accurate colormap selection. 
%% Oceanography 29(3):9-13, http://dx.doi.org/10.5670/oceanog.2016.66.
%% 
%% See also colormap and caxis.  
%
%%% Display colormap options: 
%
%if nargin==0
%   figure('menubar','none','numbertitle','off','Name','crameri options:')
%   
%   if license('test','image_toolbox')
%      imshow(imread('crameri6.0.png')); 
%   else
%      axes('pos',[0 0 1 1])
%      image(imread('crameri6.0.png')); 
%      axis image off
%   end
%   
%   return
%end
%
%%% Error checks: 
%
%assert(isnumeric(ColormapName)==0,'Input error: ColormapName must be a string.') 
%
%%% Set defaults: 
%
%NLevels = 256; 
%autopivot = false; 
%PivotValue = 0; 
%InvertedColormap = false; 
%
%%% Parse inputs: 
%
%% Does user want to flip the colormap direction? 
%dash = strncmp(ColormapName,'-',1); 
%if any(dash) 
%   InvertedColormap = true; 
%   ColormapName(dash) = []; 
%end
%
%% Standardize all colormap names to lowercase: 
%ColormapName = lower(ColormapName); 
%
%% Oleron's too hard for me to remember, so I'm gonna use dem or topo. 
%if ismember(ColormapName,{'dem','topo'})
%   ColormapName = 'oleron'; 
%end
%
%% Does the user want to center a diverging colormap on a specific value? 
%% This parsing support original 'zero' syntax and current 'pivot' syntax. 
%tmp = strncmpi(varargin,'pivot',3); 
%if any(tmp) 
%   autopivot = true; 
%   try
%      if isscalar(varargin{find(tmp)+1})
%         PivotValue = varargin{find(tmp)+1}; 
%         tmp(find(tmp)+1) = 1; 
%      end
%   end
%   varargin = varargin(~tmp); 
%end
%
%% Has user requested a specific number of levels? 
%tmp = isscalar(varargin); 
%if any(tmp) 
%   NLevels = varargin{tmp}; 
%end
%
%%% Load RGB values and interpolate to NLevels: 
%
%try
%   S = load('CrameriColourMaps6.0.mat',ColormapName); 
%   cmap = S.(ColormapName); 
%catch
%   error(['Unknown colormap name ''',ColormapName,'''. Try typing crameri with no inputs to check the options and try again.'])
%end
%
%% Interpolate if necessary: 
%if NLevels~=size(cmap,1) 
%   cmap = interp1(1:size(cmap,1), cmap, linspace(1,size(cmap,1),NLevels),'linear');
%end
%
%%% Invert the colormap if requested by user: 
%
%if InvertedColormap
%   cmap = flipud(cmap); 
%end
%
%%% Adjust values to current caxis limits? 
%
%if autopivot
%   clim = caxis; 
%   maxval = max(abs(clim-PivotValue)); 
%   cmap = interp1(linspace(-maxval,maxval,size(cmap,1))+PivotValue, cmap, linspace(clim(1),clim(2),size(cmap,1)),'linear');
%end
%
%%% Clean up 
%
%if nargout==0
%   colormap(gca,cmap) 
%   clear cmap  
%end
%
%end
