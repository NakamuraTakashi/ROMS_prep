
#gfortran mod_calendar.f90 mod_interpolation.f90 mod_roms_netcdf.f90 set_scoord.f90 pzcon.f potmp.f bryHYCOM2ROMSv4.F90 -O2 -I/usr/include -L/usr/lib -lnetcdff -o bryHYCOM2ROMS.exe
gfortran mod_calendar.f90 mod_interpolation.f90 mod_roms_netcdf.f90 set_scoord.f90 pzcon.f potmp.f bryHYCOM2ROMSv4.F90 -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o bryHYCOM2ROMS.exe

export OMP_NUM_THREADS=12

./bryHYCOM2ROMS.exe
