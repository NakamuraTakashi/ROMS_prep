gfortran -fbounds-check -fno-align-commons mod_calendar.f90 mod_interpolation.f90 set_scoord.f90 naotidej.f pzcon.f potmp.f iniHYCOMpNAO2ROMSv3.F90 -O2 -I/usr/include -L/usr/lib -lnetcdff -o iniHYCOMpNAO2ROMS.exe
./iniHYCOMpNAO2ROMS.exe
