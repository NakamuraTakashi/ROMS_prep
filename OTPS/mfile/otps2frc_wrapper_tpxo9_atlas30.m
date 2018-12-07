
addpath C:/cygwin64/home/Takashi/OTPS/TMD
addpath C:/cygwin64/home/Takashi/OTPS/TMD/FUNCTIONS
addpath C:/cygwin64/home/Takashi/OTPS/t_tide_v1.3beta

% gfile='D:/ROMS/Data/Yaeyama/Yaeyama1_grd_v10.nc'
% ofile='Yaeyama1_tide_tpxo9_atlas30_20000101.nc';
gfile='C:\cygwin64\home\Takashi\ROMS_prep\output\CT_0.08_grd_v2.nc'
ofile='CT_tide_tpxo9_atlas30_20000101.nc';

base_date=datenum(2000,1,1);
pred_date=datenum(2000,1,1);

model_dir='C:/cygwin64/home/Takashi/OTPS/DATA';

% otps2frc_v5_tpxo9_atlas30(gfile,base_date,pred_date,ofile,model_dir,'Yaeyama1')
otps2frc_v5_tpxo9_atlas30(gfile,base_date,pred_date,ofile,model_dir,'Coral Triangle')