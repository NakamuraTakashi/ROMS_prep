gfortran -fbounds-check -fno-align-commons mod_calendar.f90 mod_interpolation.f90 set_scoord.f90 naotidej.f pzcon.f potmp.f bryHYCOMpNAO2ROMSv4.F90 -O2 -I/usr/include -L/usr/lib -lnetcdff -o bryHYCOMpNAO2ROMS.exe
./bryHYCOMpNAO2ROMS.exe
