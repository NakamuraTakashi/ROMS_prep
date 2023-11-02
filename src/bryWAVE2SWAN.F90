
!!!=== Copyright (c) 2018-2023 Takashi NAKAMURA  =====

PROGRAM bryWAVE2SWAN
  use netcdf
  use eccodes
  use mod_interpolation 
  use mod_calendar
  use mod_roms_netcdf
!  use mod_swan

  implicit none
      
! SETTINGS of input and output files  -----------------------------------

  integer :: Syear, Smonth, Sday
  integer :: Eyear, Emonth, Eday
  integer :: Ryear, Rmonth, Rday
  character(256) :: GRID_FILE
  character(256) :: BRY_prefix
  integer :: SNWE(4)

  integer, parameter :: N_Param = 3
      
!----------------------------------------------------------------------
!--- JMA_CWM parameter setting -----------------
  character(256) :: CWM_dir
  character(256) :: SWAN_prefix
 
  character(len=*), parameter :: GRIB_prefix  = "Z__C_RJTD_"
  character(len=*), parameter :: GRIB_suffix  =   &
        "0000_CWM_GPV_Rjp_Gll0p05deg_FD0000-0300_grib2.bin"
  character(10) :: GRIB_yyyymmddhh = "2017010100"
  character(256) :: GRIB_FILE
  character(256), allocatable :: OUT_FILE(:)
  
  character(256) :: GRIB_NAME(N_Param) = (/  &
     "swh  "                                 &
    ,"perpw"                                 &
    ,"dirpw"                                 &
    /)
     
  integer :: GRIB_STEP(2) = (/ 0, 3 /)

  integer :: p1,p2,p3
  character(256) :: p4

  real(8), allocatable :: out_data(:,:) ! output forcing data
       
  real(8), allocatable :: lat(:), lon(:)
  real(8), allocatable :: in_data(:,:,:)
  real(8), allocatable :: rmask_dg(:,:)  ! land mask of donor grid
  real(8), allocatable :: latr_dg(:,:)
  real(8), allocatable :: lonr_dg(:,:)
  integer :: Irdg_min, Irdg_max, Jrdg_min, Jrdg_max

  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(2) :: hh
  
  integer :: i,  j
  integer :: Nxr_dg, Nyr_dg
  integer :: Ldg, Mdg

  integer :: Nout
  integer :: Idg(1000), Jdg(1000)
  integer :: Irgs(1000), Irge(1000), Jrgs(1000), Jrge(1000)
  integer :: Idgt, Jdgt
  integer :: iout, ibry
  integer :: corner(1000)=0
  character(4) :: IXXXX, JXXXX
  
  integer :: istart, iend
  integer :: itime, iparam, istep
  integer :: ifile,idx,iret,igrib
  integer :: iyear, imonth, iday
  integer :: ihour, imin
  integer :: idays,ihours,ijdate,Sjdate
  integer :: YYYYMMDD, hhmm, hhmmss
  real(8), allocatable :: values(:)
  real(8) :: Hs, Tp, Dp
  real(8) :: dL1, dL2
  
  integer :: Nxr, Nyr
  integer :: L, M  
  real(8), allocatable :: yr(:, :)
  real(8), allocatable :: xr(:, :)
  real(8), allocatable :: angler(:, :)
  real(8), allocatable :: rmask(:, :)

  integer :: ncid,var_id

!-------------------------------------------------------------------------------
  namelist/grd/GRID_FILE
  namelist/sdate/Syear, Smonth, Sday
  namelist/edate/Eyear, Emonth, Eday
  namelist/refdate/Ryear, Rmonth, Rday
  namelist/bry/BRY_prefix, SNWE
  namelist/wave_cwm/CWM_dir
  namelist/wave_cwm/SWAN_prefix

  read (5, nml=grd)
  rewind(5)
  read (5, nml=sdate)
  rewind(5)
  read (5, nml=edate)
  rewind(5)
  read (5, nml=refdate)
  rewind(5)
  read (*, nml=bry)
  rewind(5)
  read (5, nml=wave_cwm)
  rewind(5)

!---- Read ROMS grid netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( GRID_FILE )
  
  ! Open NetCDF grid file
  call check( nf90_open(trim( GRID_FILE ), nf90_nowrite, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)
  L = Nxr-1
  M = Nyr-1
  allocate( yr(0:L, 0:M) )
  allocate( xr(0:L, 0:M) )
  allocate( rmask(0:L, 0:M) )
  allocate( angler(0:L, 0:M) )
#if defined UTM_COORD
  ! Get variable id
  call check( nf90_inq_varid(ncid, 'x_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, xr) )
  call check( nf90_inq_varid(ncid, 'y_rho', var_id) )
  call check( nf90_get_var(ncid, var_id, yr) )
#else
 call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, yr) )
  call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
  call check( nf90_get_var(ncid, var_id, xr) )
#endif  
  call check( nf90_inq_varid(ncid, 'mask_rho', var_id) ) 
  call check( nf90_get_var(ncid, var_id, rmask) )
  call check( nf90_inq_varid(ncid, 'angle', var_id) )
  call check( nf90_get_var(ncid, var_id, angler) )
  ! Close NetCDF file
  call check( nf90_close(ncid) )

!---- Read JMA CWM GRIB2 file --------------------------------
#if defined JMA_CWM
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
  call codes_get(igrib,'Nj', Nyr_dg)
  call codes_get(igrib,'Ni', Nxr_dg)

  Ldg = Nxr_dg-1
  Mdg = Nyr_dg-1      
  write(*,*) Ldg, Mdg
  
  call codes_get(igrib,'distinctLatitudes',lat)
  call codes_get(igrib,'distinctLongitudes',lon)

  allocate(values(Nxr_dg*Nyr_dg))
  allocate(in_data(0:Ldg, 0:Mdg, 3))
  allocate(latr_dg(0:Ldg, 0:Mdg))
  allocate(lonr_dg(0:Ldg, 0:Mdg))
  allocate(rmask_dg(0:Ldg, 0:Mdg))

  DO WHILE (iret /= CODES_END_OF_FILE)
    call codes_get(igrib,'stepRange',p1)
    call codes_get(igrib,'shortName', p4)
    if (p1==GRIB_STEP(1) .and.             &
        trim(p4)==trim(GRIB_NAME(3))  ) exit
    call codes_release(igrib)
    call codes_grib_new_from_file(ifile,igrib, iret)
  END DO
  call codes_get(igrib,'values', values)
  call codes_release(igrib)
  call codes_close_file(ifile)

  do j=0, Mdg
    istart = 1 + Nxr_dg*j
    iend   = Nxr_dg*(j+1)
    in_data(:,j,3) = values(istart:iend)
  end do

  do j=0,Mdg
    do i=0,Ldg
      latr_dg(i,j) = lat(j+1)
      lonr_dg(i,j) = lon(i+1)
      if(in_data(i,j,3)<0) then
        rmask_dg(i,j) = 0.0d0
      else
        rmask_dg(i,j) = 1.0d0
      endif
    enddo
!    write(9,*) in_data(:,j,3)
  enddo
!  write(*,*) in_data(0,0,3), rmask_dg(0,0)

#elif defined ERA5

#elif defined SWAN

#elif defined ROMS

#endif

!---- Check IJ indices for extracting wave data --------------------------------

  write(*,*) '***************** Start extracting IJ indices *****************'
  call seek_IJrange(                                 &
        0, Ldg, 0, Mdg, lonr_dg, latr_dg             & 
      , 0, L,   0, M,   xr,      yr                  &
      , Irdg_min, Irdg_max, Jrdg_min, Jrdg_max)

  Nout = 0
    
  if( SNWE(1)==1 ) then ! south
    do i=0,L
      if( rmask(i,0)==0 ) cycle

      call nearest_id(                                  &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
        , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
        , xr(i,0), yr(i,0)                              &
        , Idgt, Jdgt )

!      write(*,*) i,xr(i,0), yr(i,0),Idgt, Jdgt
      
      if(Nout == 0) then
        Nout =  Nout+1
        Idg(Nout) = Idgt
        Jdg(Nout) = Jdgt
        Irgs(Nout) = i
        Jrgs(Nout) = 0
      else
        Irge(Nout) = i
        Jrge(Nout) = 0
        if(Idg(Nout)/=Idgt .or. Jdg(Nout)/=Jdgt) then
          Nout =  Nout+1
          Idg(Nout) = Idgt
          Jdg(Nout) = Jdgt 
          Irgs(Nout) = i
          Jrgs(Nout) = 0
        endif
      endif
      if( i==0 ) corner(Nout)=-1
      if( i==L ) corner(Nout)=-2
    enddo
  endif
  if (SNWE(4) == 1) then ! East
    do j=0,M
      if( rmask(L,j)==0 ) cycle

      call nearest_id(                                  &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
        , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
        , xr(L,j), yr(L,j)                              &
        , Idgt, Jdgt )

!      write(*,*) j,xr(L,j), yr(L,j),Idgt, Jdgt
      
      if(Nout == 0) then
        Nout =  Nout+1
        Idg(Nout) = Idgt
        Jdg(Nout) = Jdgt
        Irgs(Nout) = L
        Jrgs(Nout) = j
      else
        Irge(Nout) = L
        Jrge(Nout) = j
        if(Idg(Nout)/=Idgt .or. Jdg(Nout)/=Jdgt) then
          Nout =  Nout+1
          Idg(Nout) = Idgt
          Jdg(Nout) = Jdgt 
          Irgs(Nout) = L
          Jrgs(Nout) = j
        else
          if( corner(Nout)==-2 ) corner(Nout)=2
        endif
      endif
      if( j==M ) corner(Nout)=-3
    enddo
  endif
  if( SNWE(2)==1 ) then ! North
    do i=L,0,-1
      if( rmask(i,M)==0 ) cycle

      call nearest_id(                                  &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
        , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
        , xr(i,M), yr(i,M)                              &
        , Idgt, Jdgt )

!      write(*,*) i,xr(i,M), yr(i,M),Idgt, Jdgt
      
      if(Nout == 0) then
        Nout =  Nout+1
        Idg(Nout) = Idgt
        Jdg(Nout) = Jdgt
        Irgs(Nout) = i
        Jrgs(Nout) = M
      else
        Irge(Nout) = i
        Jrge(Nout) = M
        if(Idg(Nout)/=Idgt .or. Jdg(Nout)/=Jdgt) then
          Nout =  Nout+1
          Idg(Nout) = Idgt
          Jdg(Nout) = Jdgt 
          Irgs(Nout) = i
          Jrgs(Nout) = M
        else
          if( corner(Nout)==-3 ) corner(Nout)=3
        endif
      endif
      if( i==0 ) corner(Nout)=-4
    enddo
  endif
  if( SNWE(3)==1 ) then ! West
    do j=M,0,-1
      if( rmask(0,j)==0 ) cycle

      call nearest_id(                                  &
          Irdg_min, Irdg_max, Jrdg_min, Jrdg_max        &
        , lonr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , latr_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max)  &
        , rmask_dg(Irdg_min:Irdg_max,Jrdg_min:Jrdg_max) & 
        , xr(0,j), yr(0,j)                              &
        , Idgt, Jdgt )
      
!      write(*,*) j,xr(0,j), yr(0,j),Idgt, Jdgt

      if(Nout == 0) then
        Nout =  Nout+1
        Idg(Nout) = Idgt
        Jdg(Nout) = Jdgt
        Irgs(Nout) = 0
        Jrgs(Nout) = j
      else
        Irge(Nout) = 0
        Jrge(Nout) = j
        if(Idg(Nout)/=Idgt .or. Jdg(Nout)/=Jdgt) then
          Nout =  Nout+1
          Idg(Nout) = Idgt
          Jdg(Nout) = Jdgt 
          Irgs(Nout) = 0
          Jrgs(Nout) = j
        else
          if( corner(Nout)==-4 ) corner(Nout)=4
          if( corner(1)==-1 .and. j==0 ) then
            corner(Nout)=1
            Nout =  Nout-1
          endif
        endif
      endif
    enddo
  endif
  write(*,*) Nout
!-Create the Output file --------------------------------
  open(unit=10, file= 'swan_in_BOUNDSPEC_SEGMENT_IJ.txt', status='replace')
!  && BOUNDARY FORCING &&
!  BOUND SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGREES
!  BOUNDSPEC SEGMENT IJ 30 0 63 0 63 191 50 191 CONSTANT FILE '../../Data/Shiraho_reef2/Wave/Shiraho_wave2_201901.dat'
  write(10, '( "&& BOUNDARY FORCING &&" )' )
  write(10, '( "BOUND SHAPESPEC JONSWAP MEAN DSPR DEGREES" )' )
  allocate(OUT_FILE(Nout))
  do iout=1,Nout
    write (IXXXX, "(I4.4)") Idg(iout)
    write (JXXXX, "(I4.4)") Jdg(iout)
    OUT_FILE(iout) = trim( SWAN_prefix )//'_'//YYYY//MM//DD//'_I'//IXXXX//'_J'//JXXXX//'.dat'
    open(unit=20+iout, file= trim( OUT_FILE(iout) ), status='replace')
    write(20+iout, "(a4)") 'TPAR'
    if(corner(iout)==1) then
      write(10, '( "BOUNDSPEC SEGMENT IJ", 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, " VARIABLE FILE 0 ", A )' )  &
        Irge(Nout), Jrge(Nout), Irgs(iout), Jrgs(iout), Irge(iout), Jrge(iout), "'"//trim(OUT_FILE(iout))//"'"
    elseif(corner(iout)==2) then
      write(10, '( "BOUNDSPEC SEGMENT IJ", 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, " VARIABLE FILE 0 ", A )' )  &
        Irgs(iout), Jrgs(iout), Irge(iout), Jrgs(iout), Irge(iout), Jrge(iout), "'"//trim(OUT_FILE(iout))//"'"
    elseif(corner(iout)==3) then
      write(10, '( "BOUNDSPEC SEGMENT IJ", 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, " VARIABLE FILE 0 ", A )' )  &
        Irgs(iout), Jrgs(iout), Irgs(iout), Jrge(iout), Irge(iout), Jrge(iout), "'"//trim(OUT_FILE(iout))//"'"
    elseif(corner(iout)==4) then
      write(10, '( "BOUNDSPEC SEGMENT IJ", 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, 1x,I4, " VARIABLE FILE 0 ", A )' )  &
        Irgs(iout), Jrgs(iout), Irge(iout), Jrgs(iout), Irge(iout), Jrge(iout), "'"//trim(OUT_FILE(iout))//"'"
    else
      write(10, '( "BOUNDSPEC SEGMENT IJ", 1x,I4, 1x,I4, 1x,I4, 1x,I4, " VARIABLE FILE 0 ", A )' )  &
        Irgs(iout), Jrgs(iout), Irge(iout), Jrge(iout), "'"//trim(OUT_FILE(iout))//"'"
    endif
  enddo

  close(10)
  
!---- Write SWAN BOUNDSPEC SEGMENT IJinput script --------------------------------

  
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
    call codes_open_file(ifile, GRIB_FILE,'r')
    call codes_grib_new_from_file(ifile,igrib, iret)
    if (iret /= CODES_SUCCESS) then
      write(*,*) "CANNOT OPEN: ", trim( GRIB_FILE )
      exit
    end if  
    
!-  LOOP2 START --------------------------------
    DO istep=1,2
      DO iparam=1,3
        write(*,*) "READ: ", trim(GRIB_NAME(iparam))
        DO WHILE (iret /= CODES_END_OF_FILE)
          call codes_get(igrib,'stepRange',p1)
          call codes_get(igrib,'shortName', p4)
          if (p1==GRIB_STEP(istep) .and.             &
              trim(p4)==trim(GRIB_NAME(iparam))  ) exit
          call codes_release(igrib)
          call codes_grib_new_from_file(ifile,igrib, iret)
        END DO
        call codes_get(igrib,'validityDate',YYYYMMDD)
        write(*,*) 'validityDate=', YYYYMMDD
        call codes_get(igrib,'validityTime',hhmm)
        write(*,*) 'validityTime=', hhmm
        call codes_get(igrib,'values', values)

        call codes_release(igrib)
        call codes_grib_new_from_file(ifile,igrib, iret)
        do j=0, Mdg
          istart = 1 + Nxr_dg*j
          iend   = Nxr_dg*(j+1)
          in_data(:,j,iparam) = values(istart:iend)
        end do
      END DO


      hhmmss = hhmm*100

! -----  Write SWAN TPAR input files -------------------------------------------- 
      do iout=1,Nout

        Hs = in_data(Idg(iout),Jdg(iout),1)
        Tp = in_data(Idg(iout),Jdg(iout),2)
        Dp = in_data(Idg(iout),Jdg(iout),3)
        
!        write(* , "(I8.8 '.' I6.6 f9.3 f9.3 f9.3 f9.3)") YYYYMMDD, hhmmss, Hs, Tp, Dp, 20.0
        write(20+iout, "(I8.8 '.' I6.6 f9.3 f9.3 f9.3 f9.3)") YYYYMMDD, hhmmss, Hs, Tp, Dp, 20.0
!        write(* , "(I8.8 '.' I4.4 f9.3 f9.3 f9.3 f9.3)") YYYYMMDD, hhmm, Hs, Tp, Dp, 20.0
!        write(20+iout, "(I8.8 '.' I4.4 f9.3 f9.3 f9.3 f9.3)") YYYYMMDD, hhmm, Hs, Tp, Dp, 20.0

      enddo

    END DO
!-LOOP2 END --------------------------------
!    call codes_index_release(idx)
    call codes_release(igrib)
    call codes_close_file(ifile) 

    itime = itime + 1

  END DO
!-LOOP1 END --------------------------------

  deallocate(values)
  deallocate(lat)
  deallocate(lon)
  deallocate(in_data)

  do iout=1,Nout
    close(20+iout)
  enddo

  write(*,*) 'FINISH!!'

!-End of Main program --------------------------------------------
  
END PROGRAM bryWAVE2SWAN
      
