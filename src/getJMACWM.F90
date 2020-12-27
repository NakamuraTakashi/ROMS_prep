
!!!=== Copyright (c) 2018-2021 Takashi NAKAMURA  =====

PROGRAM getJMACWM
  use netcdf
  use eccodes
  use mod_calendar
 
  implicit none
      
! SETTINGS of input and output files  -----------------------------------

  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  integer :: Ryear, Rmonth, Rday
  
  real(8), parameter :: Slat = 24.45d0   ! Sampling Latitude
  real(8), parameter :: Slon = 124.28d0  ! Sampling Longitude

  integer, parameter :: N_Param = 3
      
!----------------------------------------------------------------------
!--- JMA_CWM parameter setting -----------------
  character(256) :: CWM_dir
  character(256) :: SWAN_prefix
 
  character(len=*), parameter :: GRIB_prefix  = "Z__C_RJTD_"
  character(len=*), parameter :: GRIB_suffix  =   &
        "0000_CWM_GPV_Rjp_Gll0p05deg_FD0000-0300_grib2.bin"
  character(11) :: OUT_suffix    = "_201701.dat"
  character(10) :: GRIB_yyyymmddhh = "2017010100"
  character(256) :: GRIB_FILE
  character(256) :: OUT_FILE
  
  character(256) :: GRIB_NAME(N_Param) = (/  &
     "swh  "                                 &
    ,"perpw"                                 &
    ,"dirpw"                                 &
    /)
     
  integer :: GRIB_STEP(2) = (/ 0, 3 /)

  real(8), allocatable :: out_data(:,:) ! output forcing data
       
  real(8), allocatable :: lat(:), lon(:)
  real(8), allocatable :: grib_data(:, :)

  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(2) :: hh
  
  integer :: i,  j
  integer :: Im, Jm
  integer :: Is, Js
  
  integer :: istart, iend
  integer :: itime, iparam
  integer :: ifile,idx,iret,igrib
  integer :: iyear, imonth, iday
  integer :: ihour, imin
  integer :: idays,ihours,ijdate,Sjdate
  integer :: YYYYMMDD, hhmm, hhmmss
  real(8), allocatable :: values(:)
  real(8) :: Hs, Tp, Dp
  real(8) :: dL1, dL2
  
  integer :: itmp
  character(6), allocatable :: ctmp(:)
  integer, allocatable :: tmp(:)

  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/wave_cwm/CWM_dir
  namelist/wave_cwm/SWAN_prefix

!  read (*, nml=grd)
  read (5, nml=sdate)
  rewind(5)
  read (5, nml=edate)
  rewind(5)
  read (5, nml=refdate)
  rewind(5)
  read (5, nml=wave_cwm)
  rewind(5)

!---- Read JMA CWM GRIB2 file --------------------------------

  write (YYYY, "(I4.4)") Syear
  write (MM, "(I2.2)") Smonth
  write (DD, "(I2.2)") Sday
  GRIB_yyyymmddhh = YYYY//MM//DD//"00"

  GRIB_FILE = trim(CWM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
                GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix
  !Open GRIB file
  write(*,*) "OPEN: ", trim( GRIB_FILE )
  call codes_grib_multi_support_on	(	iret	)	
  call codes_open_file(ifile, GRIB_FILE,'r')
  call codes_grib_new_from_file(ifile,igrib, iret)

  ! Get dimension data
  call codes_get(igrib,'Nj', Jm)
  call codes_get(igrib,'Ni', Im)
      
  write(*,*) Im, Jm
  
!  ! Allocate variable
!  allocate(lat(Jm))
!  allocate(lon(Im))
  
  call codes_get(igrib,'distinctLatitudes',lat)
  call codes_get(igrib,'distinctLongitudes',lon)
  
  call codes_release(igrib)
  call codes_close_file(ifile)
  
  allocate(values(Im*Jm))
  allocate(grib_data(Im, Jm))
  
  dL1= abs(lat(1)-slat)
  Js = 1
  do j=2,Jm
    dL2= abs(lat(j)-slat)
    if(dL2<dL1) then
      Js = j
      dL1 = dL2
    end if
  end do
  
  dL1= abs(lon(1)-slon)
  Is = 1
  do i=2,Im
    dL2= abs(lon(i)-slon)
    if(dL2<dL1) then
      Is = i
      dL1 = dL2
    end if
  end do
  write(*,*) Is, Js, lon(Is), lat(Js)

  
!-Create the Output file --------------------------------

  OUT_suffix(2:5)=YYYY
  OUT_suffix(6:7)=MM
  OUT_FILE = trim( SWAN_prefix )//OUT_suffix
  open(unit=20, file= trim( OUT_FILE ), status='replace')
  write(20, "(a4)") 'TPAR'
  
!---- LOOP set up --------------------------------
  itime = 1
  ihours = 0
  iyear = Syear
  imonth = Smonth
  iday = Sday
  call jd(iyear, imonth, iday, Sjdate)
!---- LOOP1 START --------------------------------  
  DO
    ihour = mod(ihours,24)
    ijdate = Sjdate + int(ihours/24)
    call cdate( ijdate, iyear, imonth, iday )
    ! Check end date
    if(iyear==Eyear .and. imonth==Emonth .and. iday==Eday) then
      write(*,*) "Completed!!!"
      STOP
    endif

    write (YYYY, "(I4.4)") iyear
    write (MM, "(I2.2)") imonth
    write (DD, "(I2.2)") iday ! 1+int((itime-1)*1/24)
    write (hh, "(I2.2)") ihour ! mod((itime-1)*1,24)

    GRIB_yyyymmddhh = YYYY//MM//DD//hh

    ihours = ihours + 6  !!! Files exist 6 hourly interval
  
    GRIB_FILE = trim(CWM_dir)//YYYY//"/"//MM//"/"//DD//"/"// &
                  GRIB_prefix//GRIB_yyyymmddhh//GRIB_suffix    
    !Open GRIB file
    write(*,*) "OPEN: ", trim( GRIB_FILE )
    call codes_index_create(idx,GRIB_FILE,'shortName:s,stepRange:i',iret)
    if (iret /= CODES_SUCCESS) then
      write(*,*) "CANNOT OPEN: ", trim( GRIB_FILE )
      exit
    end if
    
!    call codes_index_get_size(idx,'shortName',itmp) !!!!!!!
!    allocate(ctmp(itmp)) !!!
!    call codes_index_get(idx,'shortName',ctmp) !!!!!!
!    write(*,*) 'shortName: ', itmp, ctmp !!!!!!!!!
!    deallocate(ctmp) !!!!!!
!    call codes_index_get_size(idx,'stepRange',itmp) !!!!!!!!
!    allocate(tmp(itmp)) !!!!!!
!    call codes_index_get(idx,'stepRange',tmp) !!!!!!!!!!
!    write(*,*) 'stepRange: ', itmp, tmp !!!!!!!!
!    deallocate(tmp) !!!!!!!!
    
!-  LOOP2 START --------------------------------
    DO iparam=1,2
      write(*,*) "READ: ", trim( GRIB_FILE )
      call codes_index_select(idx,'stepRange',GRIB_STEP(iparam))
      call codes_index_select(idx,'shortName',GRIB_NAME(1))
      call codes_new_from_index(idx,igrib, iret)
      call codes_get(igrib,'validityDate',YYYYMMDD)
      write(*,*) 'validityDate=', YYYYMMDD
      call codes_get(igrib,'validityTime',hhmm)
      write(*,*) 'validityTime=', hhmm
      
      call codes_get(igrib,'values', values)
      call codes_release(igrib)
      do i=1, Jm
        istart = 1 + Im*(i-1)
        iend   = Im*i
        grib_data(:,i) = values(istart:iend)
      end do
      Hs = grib_data(Is,Js)
      
      call codes_index_select(idx,'stepRange',GRIB_STEP(iparam))
      call codes_index_select(idx,'shortName',GRIB_NAME(2))
      call codes_new_from_index(idx,igrib, iret)
      call codes_get(igrib,'values', values)
      call codes_release(igrib)
      do i=1, Jm
        istart = 1 + Im*(i-1)
        iend   = Im*i
        grib_data(:,i) = values(istart:iend)
      end do
      Tp = grib_data(Is,Js)
      
      call codes_index_select(idx,'stepRange',GRIB_STEP(iparam))
      call codes_index_select(idx,'shortName',GRIB_NAME(3))
      call codes_new_from_index(idx,igrib, iret)
      call codes_get(igrib,'values', values)
      call codes_release(igrib)
      do i=1, Jm
        istart = 1 + Im*(i-1)
        iend   = Im*i
        grib_data(:,i) = values(istart:iend)
      end do
      Dp = grib_data(Is,Js)

      hhmmss = hhmm*100
      
      write(* , "(I8.8 '.' I6.6 f9.3 f9.3 f9.3 f9.3)") YYYYMMDD, hhmmss, Hs, Tp, Dp,30.0
      write(20, "(I8.8 '.' I6.6 f9.3 f9.3 f9.3 f9.3)") YYYYMMDD, hhmmss, Hs, Tp, Dp,30.0
      
    END DO
!-LOOP2 END --------------------------------
    call codes_index_release(idx)

    itime = itime + 1

  END DO
!-LOOP1 END --------------------------------

  deallocate(values)
  deallocate(lat)
  deallocate(lon)
  deallocate(grib_data)
  close(20)
  write(*,*) 'FINISH!!'

!-End of Main program --------------------------------------------
  
END PROGRAM getJMACWM
      
