
addpath C:/cygwin64/home/Takashi/OTPS/TMD
addpath C:/cygwin64/home/Takashi/OTPS/TMD/FUNCTIONS
addpath C:/cygwin64/home/Takashi/OTPS/t_tide_v1.3beta

gfile='D:/ROMS/Data/Yaeyama/Yaeyama1_grd_v10.nc'
base_date=datenum(2000,1,1);
pred_date=datenum(2000,1,1);
ofile='tidetest.nc';
% model_file='Model_tpxo9_atlas_30';
model_file='Model_tpxo9';

otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,'ESPRESSO')