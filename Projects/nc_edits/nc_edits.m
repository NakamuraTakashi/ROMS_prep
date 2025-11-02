% === Copyright (c) 2024-2025 Takashi NAKAMURA  =====

%% Check ncfile contents
% NC_FILE = 'shiraho_roms_grd_JCOPET_v18.0.nc';
% NC_FILE = '\\wsl.localhost\ubuntu-24.04\home\nakamulab2\COAWST\Projects\Shiraho_reef2_eco\SR_veg_eco2_sg_his_20231007.nc';
% NC_FILE = './tokyobay_grid/TokyoBay3_grd_v3.1.nc';
% NC_FILE = './tokyobay_grid/TokyoBay3_river_NZ30_1990_2020v4.nc';
% NC_FILE = './tokyobay_grid/test2.nc';
% NC_FILE = './palau_grid/Palau2_grd_v0.0.nc';
% NC_FILE = '../Palau/Palau1_grd_v1.0.nc';
NC_FILE = './Kushimoto_grid/Kushimoto_grd_v0.0.nc';
% NC_FILE = 'D:\COAWST_DATA\Yaeyama\Yaeyama1\Grid\Yaeyama1_grd_v10.nc';
% NC_FILE = 'D:\COAWST_DATA\Yaeyama\Yaeyama2\Grid\Yaeyama2_grd_v11.3.nc';
ncdisp(NC_FILE)

%% Read 2D (x,y) ncdata
NC_VAR = 'h';
% NC_VAR = 'mask_rho';
% NC_VAR = 'aquaculture_01';
% NC_VAR = 'aquaculture_02';
% NC_VAR = 'aquaculture_03';
% NC_VAR = 'aquaculture_04';
% NC_VAR = 'p_sand';
% NC_VAR = 'sgd_src';
% NC_VAR = 'p_sgrass_01';
% NC_VAR = 'p_coral_01';
% NC_VAR = 'p_coral_02';
% NC_VAR = 'p_algae_01';
ncdata = ncread(NC_FILE,NC_VAR);

%% Read 3D (x,y,time) ncdata 
NC_VAR = 'mudmass_02';
itime = 1;
ncdata = ncread(NC_FILE,NC_VAR,[1,1,1,itime],[Inf,Inf,1,1]);

%% Read 4D (x,y,z,time) ncdata 
NC_VAR = 'mudmass_02';
itime = 1;
iz = 1;
ncdata = ncread(NC_FILE,NC_VAR,[1,1,iz,itime],[Inf,Inf,1,1]);

%% write nc data as a CSV file
ncdata2=flipud(ncdata');
% CSV_FILENAME = "./palau_grid/Palau2_grd_v1.0_"+NC_VAR+"_2.csv";
% CSV_FILENAME = "./Yaeyama_grid/Yaeyama2_grd_v11.3_"+NC_VAR+".csv";
CSV_FILENAME = "./Kushimoto_grid/Kushimoto_grd_v0.0_"+NC_VAR+".csv";

writematrix(ncdata2, CSV_FILENAME);
% ncdata2= zeros(size(ncdata));

%% read CSV file and overwrite NC_VAR to ncfile
% NC_OUTFILE = 'shiraho_roms_grd_JCOPET_v18.1_nosg.nc';
% NC_OUTFILE = './tokyobay_grid/TokyoBay3_grd_v3.1.nc';
% NC_OUTFILE = './palau_grid/Palau2_grd_v0.1.nc';
% NC_OUTFILE = '../Palau/Palau1_grd_v1.0.nc';
NC_OUTFILE = './Kushimoto_grid/Kushimoto_grd_v0.1.nc';
% NC_OUTFILE = 'D:\COAWST_DATA\Yaeyama\Yaeyama1\Grid\Yaeyama1_grd_v10.1.nc';

% CSV_INFILE = "mask_rho_2.csv";

% CSV_INFILE = "p_sand_fin.csv";
% CSV_INFILE = CSV_FILENAME;
% NC_VAR = 'p_sand';

% CSV_INFILE = "./palau_grid/Palau2_grd_h_2.csv";
% CSV_INFILE = "./palau_grid/Palau1_grd_v1.0_h_2.csv";
CSV_INFILE = "./Kushimoto_grid/Kushimoto_grd_v0.1_h.csv";
NC_VAR = 'h';

% CSV_INFILE = "sgd_src_fin.csv";
% NC_VAR = 'sgd_src';

% CSV_INFILE = "p_sgrass_nosg.csv";
% NC_VAR = 'p_sgrass_01';

csvdata = readmatrix(CSV_INFILE);
ncdata3=flipud(csvdata)';
ncwrite(NC_OUTFILE,NC_VAR,ncdata3);
%% 
NC_OUTFILE = 'D:\COAWST_DATA\Yaeyama\Yaeyama2\Grid\Yaeyama2_grd_v11.3.nc';
csvdata = ones(size(ncdata));
% csvdata = zeros(size(ncdata));
ncdata3=flipud(csvdata)';
ncwrite(NC_OUTFILE,NC_VAR,ncdata3);