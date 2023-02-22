% prepare ISMIP6 input
sel = [1, 86,186,286];
nx = 381;

ncl = ncload('../Data/ISMIP6/NORCE/lithk_AIS_NORCE_CISM5-MAR364-ERA-t1_expAE03.nc');
nct = ncload('../Data/ISMIP6/NORCE/topg_AIS_NORCE_CISM5-MAR364-ERA-t1_expAE03.nc');
nca = ncload('../Data/ISMIP6/Grids/af2_ISMIP6_AIS_16000m.nc');
THICK = reshape(ncl.lithk(:,:,sel),[nx*nx,4]);
BED = reshape(nct.topg(:,:,sel),[nx*nx,4]);
AF2 = reshape(nca.af2(:,:),[nx*nx,1]);

