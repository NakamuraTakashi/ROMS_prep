
addpath 'D:/cygwin64/home/Takashi/ROMS_prep/OTPS/TMD'
addpath 'D:/cygwin64/home/Takashi/ROMS_prep/OTPS/TMD/FUNCTIONS'
addpath 'D:/cygwin64/home/Takashi/ROMS_prep/OTPS/t_tide_v1.3beta'

% gfile='D:/ROMS/Data/Yaeyama/Yaeyama1_grd_v10.nc'
% ofile='Yaeyama1_tide_tpxo9_atlas30_20000101.nc';
% gfile='C:\cygwin64\home\Takashi\ROMS_prep\output\CT_0.08_grd_v2.nc'
% ofile='CT_tide_tpxo9_atlas30_20000101.nc';
gfile='D:/cygwin64/home/Takashi/ROMS_prep/Projects/TokyoBay/TokyoBay1_grd_v1.0.nc';
ofile='TokyoBay1_tide_tpxo9_atlas30_20000101.nc';
% gfile='D:/cygwin64/home/Takashi/ROMS_prep/Projects/RedSea/RedSea1_grd_v0.nc';
% ofile='RedSea1_tide_tpxo9_atlas30_20000101.nc';

base_date=datenum(2000,1,1);
pred_date=datenum(2000,1,1);

model_dir='D:/cygwin64/home/Takashi/ROMS_prep/OTPS/DATA';

% otps2frc_v5_tpxo9_atlas30(gfile,base_date,pred_date,ofile,model_dir,'Yaeyama1')
% otps2frc_v5_tpxo9_atlas30(gfile,base_date,pred_date,ofile,model_dir,'Coral Triangle')
otps2frc_v5_tpxo9_atlas30(gfile,base_date,pred_date,ofile,model_dir,'TokyoBay1')
% otps2frc_v5_tpxo9_atlas30(gfile,base_date,pred_date,ofile,model_dir,'RedSea1')