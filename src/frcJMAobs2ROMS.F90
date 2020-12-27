
!!!=== Copyright (c) 2020-2021 Takashi NAKAMURA  =====

PROGRAM frcJMAobs2ROMS
  use netcdf
  use mod_roms_netcdf
  use mod_calendar
 
  implicit none
  
!-------------------------------------------------------------------------------
  integer :: Ryear, Rmonth, Rday
  integer, parameter :: N_InPar  = 8 
  integer, parameter :: N_OutPar = 9
  character(256) :: FRC_prefix, FRC_suffix
  character(9) :: FRC_yyyymmdd = "_20020701" 
  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"
  character(256) :: FRC_FILE(2)

  character(256) :: NC_NAME(N_OutPar) = (/    &
     "Pair      "                             &
    ,"Uwind     "                             &
    ,"Vwind     "                             &
    ,"Tair      "                             &
    ,"Qair      "                             &
    ,"cloud     "                             &
    ,"rain      "                             &
    ,"swrad     "                             &
    ,"lwrad_down"                             &
    /)
  character(256) :: NC_LNAME(N_OutPar) = (/   &
     "surface air pressure               "    &
    ,"surface u-wind component           "    &
    ,"surface v-wind component           "    &
    ,"surface air temperature            "    &
    ,"surface air relative humidity      "    &
    ,"cloud fraction                     "    &
    ,"rain fall rate                     "    &
    ,"solar shortwave radiation flux     "    &
    ,"downwelling longwave radiation flux"    &
    /)
  character(256) :: NC_UNIT(N_OutPar) = (/    &
     "millibar                 "              &
    ,"meter second-1           "              &
    ,"meter second-1           "              &
    ,"Celsius                  "              &
    ,"percentage               "              &
    ,"0 to 1                   "              &
    ,"kilogram meter-2 second-1"              &
    ,"watt meter-2             "              &
    ,"watt meter-2             "              &
    /)

  real(8), parameter :: PI = 3.141592653589793d0
  integer :: airvar(8)
  character(256) :: Pair_NAME, wind_NAME, Tair_NAME, Qair_NAME, cloud_NAME
  character(256) :: rain_NAME, swrad_NAME
  character(256) :: LW_prefix, LW_suffix
  character(256) :: FILENAME(8)
  character(1000) :: buff
  character(10), allocatable :: cdata(:)
  integer :: is,iyr
  integer :: iqt1,iqt2
  real(8) :: t
  integer :: Nt
  integer, allocatable :: iyear(:), imonth(:), iday(:), ihour(:)
  real(8), allocatable :: rdata(:,:)
  character(10) :: cwdir,ccloud
  character(10) :: cangle(16)
  real(8) :: rangle(16)
  real(8) :: rwspd, angle
  real(8), allocatable :: time(:), out_data(:,:) ! output forcing data
  integer :: ios
  integer :: data_num, header_num
       
  integer :: time_dimid
  integer :: dimids(1)
  integer :: start1D(1), count1D(1)

  integer :: i,j,k
  integer :: idays,ihours,ijdate,Sjdate
  integer :: itime

  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(2) :: hh

  integer :: ncid,var_id
  
  integer :: iparam,iOutPar,nnc
  real(8) :: u, v
  logical :: file_exists
!
!-------------------------------------------------------------------------------
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/frc_jmaobs/airvar
  namelist/frc_jmaobs/Pair_NAME, wind_NAME, Tair_NAME, Qair_NAME, cloud_NAME
  namelist/frc_jmaobs/rain_NAME, swrad_NAME
  namelist/frc_jmaobs/LW_prefix, LW_suffix
  namelist/frc_jmaobs/FRC_prefix
  ! Read parameters in namelist file
  read (5, nml=refdate)
  rewind(5)
  read (5, nml=frc_jmaobs)
  rewind(5)

!---- Modify time-unit description for NetCDF output ---------------------------------
      
  write (YYYY, "(I4.4)") Ryear
  write (MM, "(I2.2)") Rmonth
  write (DD, "(I2.2)") Rday
  
  TIME_ATT(12:15)=YYYY
  TIME_ATT(17:18)=MM
  TIME_ATT(20:21)=DD

!---- Set input file name ---------------------------------

  FILENAME(1) = trim( Pair_NAME )
  FILENAME(2) = trim( wind_NAME )
  FILENAME(3) = trim( Tair_NAME )
  FILENAME(4) = trim( Qair_NAME )
  FILENAME(5) = trim( cloud_NAME )
  FILENAME(6) = trim( rain_NAME )
  FILENAME(7) = trim( swrad_NAME )

!---- Read JMA observation data ------------------------------------
  DO iparam=1, N_InPar
    if( airvar(iparam)==0 ) cycle

    iOutPar = iparam
    if( iparam >= 3) iOutPar = iOutPar + 1

!---- Read data file(s) ------------------------------------------    

    if( iparam == 8 ) then
      Nt = 0
      allocate( iyear(366*24) )
      allocate( imonth(366*24) )
      allocate( iday(366*24) )
      allocate( ihour(366*24) )
      allocate( rdata(1,366*24) )
      allocate( cdata(24) )

    ! Read Long-wave radiation data
      DO k=1,12
        write (MM, "(I2.2)") k
        FILENAME(8) = trim( LW_prefix )//MM//trim( LW_suffix )
        ! Check file
        write(*,*) "CHECK: ", trim( FILENAME(8) )
        inquire(FILE=trim( FILENAME(8) ), EXIST=file_exists)
    
        if( file_exists ) then
          write(*,*) "PASSED"
        else
          write(*,*) "Not found..."
          cycle
        endif

        write(*,*) 'Read: ', trim( FILENAME(iparam) )
        open(50, file = trim( FILENAME(iparam) ) )
  
        data_num = 59
        header_num = 18
        data_num = data_num - header_num ! Remove header lines
  
        do i=1, header_num 
          read(50,'(a)') buff
!          write(*,*) trim(buff)
          if(i==15) read(buff(38:41),*) iyr
        end do
  
        do i=1, data_num
          read(50,'(a)') buff
          if( buff(15:16) == '--' ) exit
!          write(*,*) trim(buff)

          do j=1,24
            is = 21 + (j-1)*5
            cdata(j) = buff(is:is+4)
            if(cdata(j)==' XXX') cycle
  
            Nt = Nt + 1
            iyear(Nt) = iyr
            read(buff(15:16),*) imonth(Nt)
            read(buff(17:18),*) iday(Nt)
            ihour(Nt) = j
            read(cdata(j),*) rdata(1,Nt)
            rdata(1,Nt) = rdata(1,Nt)*0.01d0*1.0d6/3600.0d0  ! 0.01 MJ m-2 h-1 -> W m-2

!            write(*,*) iyear(Nt),imonth(Nt),iday(Nt),ihour(Nt),rdata(1,Nt)
          enddo        
        enddo
        close(50)
      END DO

    else
    ! Read JMA observation data
      write(*,*) 'Read: ', trim( FILENAME(iparam) )
      open(50, file = trim( FILENAME(iparam) ) )
    
      data_num = 0
      do
        read(50,*,iostat=ios)
        if(ios==-1) exit
        data_num = data_num+1
      end do
  
      if( iparam == 2 ) then !! for Wind data
        header_num = 6
        nnc = 1
        allocate( cdata(9) )
        allocate( rdata(1:1+nnc,data_num) )

        open(51, file = ANGLE_DIR//'/wind_angle.csv' )
        do i=1,16
          read(51,*) cangle(i), rangle(i)
!          write(*,*) cangle(i), rangle(i)
        end do
        close(51)

      elseif( iparam == 6 ) then !! for rain data
        header_num = 5
        nnc = 0
        allocate( cdata(8) )
        allocate( rdata(1,data_num) )
      else
        header_num = 5
        nnc = 0
        allocate( cdata(7) )
        allocate( rdata(1,data_num) )
      endif
      data_num = data_num - header_num ! Remove header lines
  
      allocate( iyear(data_num) )
      allocate( imonth(data_num) )
      allocate( iday(data_num) )
      allocate( ihour(data_num) )
      
      rewind(50)
    
      do i=1, header_num 
        read(50,'(a)') buff
!        write(*,*) trim(buff)
      end do
      Nt = 0
      do i=1, data_num
        read(50,'(a)') buff
!        write(*,*) trim(buff)
        read(buff,*) cdata

        if(iparam == 2) then  ! for wind
          read( cdata(6),*) iqt1
          read( cdata(8),*) iqt2
          if( iqt1 >= 5 .and. iqt2 >= 5 ) then
            Nt = Nt + 1
            read( cdata(1),*) iyear(Nt)
            read( cdata(2),*) imonth(Nt)
            read( cdata(3),*) iday(Nt)
            read( cdata(4),*) ihour(Nt)
            read( cdata(5),*) rwspd
            read( cdata(7),*) cwdir

            if(trim(cwdir)==trim(cangle(1))) then
              angle = rangle(1)
            elseif(trim(cwdir)==trim(cangle(2))) then
              angle = rangle(2)
            elseif(trim(cwdir)==trim(cangle(3))) then
              angle = rangle(3)
            elseif(trim(cwdir)==trim(cangle(4))) then
              angle = rangle(4)
            elseif(trim(cwdir)==trim(cangle(5))) then
              angle = rangle(5)
            elseif(trim(cwdir)==trim(cangle(6))) then
              angle = rangle(6)
            elseif(trim(cwdir)==trim(cangle(7))) then
              angle = rangle(7)
            elseif(trim(cwdir)==trim(cangle(8))) then
              angle = rangle(8)
            elseif(trim(cwdir)==trim(cangle(9))) then
              angle = rangle(9)
            elseif(trim(cwdir)==trim(cangle(10))) then
              angle = rangle(10)
            elseif(trim(cwdir)==trim(cangle(11))) then
              angle = rangle(11)
            elseif(trim(cwdir)==trim(cangle(12))) then
              angle = rangle(12)
            elseif(trim(cwdir)==trim(cangle(13))) then
              angle = rangle(13)
            elseif(trim(cwdir)==trim(cangle(14))) then
              angle = rangle(14)
            elseif(trim(cwdir)==trim(cangle(15))) then
              angle = rangle(15)
            elseif(trim(cwdir)==trim(cangle(16))) then
              angle = rangle(16)
            else
              angle = 0.0d0
            endif

            rdata(1,Nt) = rwspd*cos(angle/180.0d0*PI)
            rdata(2,Nt) = rwspd*sin(angle/180.0d0*PI)

!            write(*,*) iyear(Nt),imonth(Nt),iday(Nt),ihour(Nt),rwspd,cwdir,rdata(1,Nt),rdata(2,Nt)
          endif

        elseif(iparam == 5) then ! for cloud
          read( cdata(6),*) iqt1
          if( iqt1 >= 5 ) then
            Nt = Nt + 1
            read( cdata(1),*) iyear(Nt)
            read( cdata(2),*) imonth(Nt)
            read( cdata(3),*) iday(Nt)
            read( cdata(4),*) ihour(Nt)
            read( cdata(5),*) ccloud
            if(trim(ccloud)=='10-') then
              rdata(1,Nt) = 0.95d0
            elseif(trim(ccloud)=='0+') then
              rdata(1,Nt) = 0.05d0
            else
              read( ccloud,*) rdata(1,Nt)
              rdata(1,Nt) = rdata(1,Nt)/10.0d0
            endif
!            write(*,*) iyear(Nt),imonth(Nt),iday(Nt),ihour(Nt),rdata(1,Nt)
          endif
        
        elseif(iparam == 6) then ! for rain
          read( cdata(7),*) iqt1
          if( iqt1 >= 5 ) then
            Nt = Nt + 1
            read( cdata(1),*) iyear(Nt)
            read( cdata(2),*) imonth(Nt)
            read( cdata(3),*) iday(Nt)
            read( cdata(4),*) ihour(Nt)
            read( cdata(5),*) rdata(1,Nt)
            rdata(1,Nt) = rdata(1,Nt)/3600.0d0    ! mm h-1 -> kg m-2 s-1
!            write(*,*) iyear(Nt),imonth(Nt),iday(Nt),ihour(Nt),rdata(1,Nt)
          endif

        else
          read( cdata(6),*) iqt1
          if( iqt1 >= 5 ) then
            Nt = Nt + 1
            read( cdata(1),*) iyear(Nt)
            read( cdata(2),*) imonth(Nt)
            read( cdata(3),*) iday(Nt)
            read( cdata(4),*) ihour(Nt)
            read( cdata(5),*) rdata(1,Nt)
            if(iparam == 7) then  ! for swrad
              rdata(1,Nt) = rdata(1,Nt)*1.0d6/3600.0d0  ! MJ m-2 h-1 -> W m-2
            endif 
!            write(*,*) iyear(Nt),imonth(Nt),iday(Nt),ihour(Nt),rdata(1,Nt)
          endif         
        endif
  
      end do
      close(50)
  
    endif

    allocate( time(Nt) )
    allocate( out_data(iOutPar:iOutPar+nnc,Nt) )

    do i=1, Nt
      call ndays(imonth(i), iday(i), iyear(i), Rmonth, Rday, Ryear, idays)     
      t = dble(idays)+dble(ihour(i))/24.0d0
      if(iparam <= 6) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        time(i) = t
      else
        time(i) = t - 0.5d0/24.0d0
      endif
#if !defined LOCAL_TIME
      time(i) = time(i) - 9.0d0/24.0d0
#endif
      out_data(iOutPar:iOutPar+nnc,i) = rdata(1:1+nnc,i)
    enddo

!---- Create the forcing netCDF file --------------------------------

    if( iparam == 2 ) then !! for Wind data
      FRC_FILE(1) = trim(FRC_prefix)//'_wind.nc'
    else     
      FRC_FILE(1) = trim(FRC_prefix)//'_'//trim( NC_NAME(iOutPar) )//'.nc'
      nnc = 0
    endif
    
    write(*,*) "CREATE: ", trim( FRC_FILE(1) )
    call check( nf90_create(trim( FRC_FILE(1) ), nf90_clobber, ncid) )
    call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid) )
    dimids = (/ time_dimid /)
    call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, var_id) )
    call check( nf90_put_att(ncid, var_id, 'long_name', 'atmospheric forcing time') )
    call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )
    DO i=iOutPar,iOutPar+nnc
      call check( nf90_def_var(ncid, trim( NC_NAME(i) ), NF90_REAL, dimids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', trim( NC_LNAME(i) )) )
      call check( nf90_put_att(ncid, var_id, 'units',     trim( NC_UNIT(i) ) ) )
      call check( nf90_put_att(ncid, var_id, 'time',      'time') )
    END DO

! End define mode.
    call check( nf90_enddef(ncid) )
    call check( nf90_close(ncid) ) 
  
!---- Write forcing netCDF file --------------------------------

    write(*,*) time(1),TIME_ATT
    start1D = (/ 1 /)
    count1D = (/ Nt /)
    call writeNetCDF_1d( 'time', trim( FRC_FILE(1) )                  &
        , 1, time, start1D, count1D )
    DO i=iOutPar,iOutPar+nnc
      call writeNetCDF_1d(trim( NC_NAME(i) ), trim( FRC_FILE(1))      &
        , 1, out_data(i,:), start1D, count1D )
    END DO

    deallocate( cdata )
    deallocate( iyear, imonth, iday, ihour )
    deallocate( rdata, time, out_data )
    
  END DO
  
  write(*,*) 'FINISH!!'

!---- End of Main program --------------------------------------------
      
END PROGRAM frcJMAobs2ROMS
      
