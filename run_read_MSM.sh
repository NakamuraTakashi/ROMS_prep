gfortran mod_calendar.f90 mod_interpolation.f90 read_MSM.f90 -O2 -I/usr/include -L/usr/lib -lnetcdff -lkriging -llapack -lblas -o read_MSM.exe
./read_MSM.exe
