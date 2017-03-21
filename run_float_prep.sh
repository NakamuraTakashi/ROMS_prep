#gfortran mod_calendar.f90 float_prep.f90 -I/usr/include -L/usr/lib -lnetcdff -lkriging -llapack -lblas -o float_prep.exe
gfortran mod_calendar.f90 float_prep.f90 -I/usr/include -L/usr/lib -lnetcdff -o float_prep.exe
./float_prep.exe
