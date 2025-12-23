% === Copyright (c) 2025 Takashi NAKAMURA  =====
%
% STEP 1
%  Trim time slice from the his file
%   ncks -d ocean_time,2 ocean_his.nc ocean_ini.nc
%  example:
%   ncks -d ocean_time,24 ../../../COAWST_OUTPUT/Shizugawa/SZ3_eco_oys_test1/SZ3_veg_eco_his_20230701.nc Shizugawa3_ini_ROMS_Nz15_20230702.00.nc

% STEP 2
%  Run this matlab scriot 

%% Check ncfile contents
ini = 'Shizugawa3_ini_ROMS_Nz15_20230702.00.nc';
grd = 'D:/COAWST_DATA/Shizugawa/Shizugawa3/Grid/Shizugawa3_grd_v0.4c.nc';
plgn = 'Shizugawa/input/Shizugawa_Bay_area/Shizugawa_Bay_area.shp';
ncdisp(ini)

%% Read polygon data

pg = shaperead(plgn);
Npg = size(pg.X,2)-1;
lon_pg = pg.X(1:Npg);
lat_pg = pg.Y(1:Npg);

%% Read grid data
h     = ncread(grd,'h');
rmask = ncread(grd,'mask_rho');
latr = ncread(grd,'lat_rho');
lonr = ncread(grd,'lon_rho');
pm = ncread(grd,'pm');
pn = ncread(grd,'pn');


%% Read ini data
temp = ncread(ini,'temp');
L=size(temp,2);
M=size(temp,1);
N=size(temp,3);
time = ncread(ini,'ocean_time');
dye=zeros(M,L,N);

%% Dye area
Vol = 0;
for i=1:L
    for j=1:M
        if (inpolygon(lonr(j,i),latr(j,i),lon_pg,lat_pg) && rmask(j,i)==1)
            dye(j,i,:) = 1;     % Set dye concentration to one in water cells
            Vol = Vol + h(j,i)/pm(j,i)/pn(j,i);
        else
            dye(j,i,:) = 0;     % Set dye concentration to zero in land cells
        end
    end
end
disp ('Water volume (m^3): ')
disp (Vol)

%% Add dye to the netcdf file

NC_VARNAME = 'dye_01';
NC_LONGNAME = 'dye concentration, type 01';
NC_UNIT = 'kilogram meter-3';
NC_TIME = 'ocean_time';
% dye_01
nccreate(ini,NC_VARNAME,...
          'Dimensions',{'xi_rho',M, 'eta_rho',L, 's_rho',N, 'ocean_time',time },...
          'Datatype','double')
ncwriteatt(ini,NC_VARNAME,'long_name',NC_LONGNAME);
ncwriteatt(ini,NC_VARNAME,'units',NC_UNIT);
ncwriteatt(ini,NC_VARNAME,'time', NC_TIME);

ncwrite(ini,NC_VARNAME,dye);
