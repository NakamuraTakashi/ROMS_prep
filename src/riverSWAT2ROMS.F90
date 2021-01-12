
!!!=== Copyright (c) 2021 Takashi NAKAMURA  ===== 

PROGRAM riverSWAT2ROMS
  use netcdf
  use mod_utility
  use mod_roms_netcdf
  use mod_calendar
 
  implicit none
      
! ---------------------------------------------------------------------
 
  character(30) :: TIME_ATT  = "days since 2000-01-01 00:00:00"

  real(8), allocatable :: time_all(:)
  integer :: start1D(1), count1D(1)
  integer :: start2D(2), count2D(2)
  integer :: start3D(3), count3D(3)
  integer :: start4D(4), count4D(4)
  
  integer :: Nzr, Nsrc

  character(4) :: YYYY
  character(2) :: MM
  character(2) :: DD
  character(11) :: YYYYMMDDpHH
  integer :: jdate, jdate_Ref
  real(8) :: d_jdate, d_jdate_Ref

  integer :: ncid,var_id
  integer :: ncid2,var_id2
  character(256) :: varname
  character(2) :: varnum
  
  character(2000) :: headerline
  integer :: iyear, imonth, iday
  integer :: ihour, imin
  integer :: jday, iunit, gis_id
  integer :: itime
  character(10) :: name
  real(8), allocatable :: river_data(:)
  real(8), allocatable :: data(:,:,:), river_Vshape(:,:)
  real(8), allocatable :: river_flow(:,:), river_flow2(:,:)
  real(8) :: river_time(1)
  character(20), allocatable :: label(:), label2(:)
  integer :: ios
  integer :: Ncol,Ntype
  integer :: i,j,k
  !
  integer :: Ryear, Rmonth, Rday
  character(256) :: ROMS_HISFILE
  integer :: romsvar(N_var)
  character(256) :: OUT_FILE, GLOBAL_ATT
  integer :: Isrc, Jsrc, river_dir, river_dir2, NCS, NNS, NCS2, NNS2

  character(256) :: RIVER_FILE
  character(20) :: cha_name
  logical :: Yield
  integer :: LOCAL_TIME
  real(8) :: Flow0, T0, S0, TAlk0, TIC_0, Oxyg0
  real(8) :: DOC_0(Ndom), POC_0(Npom), Phyt_0(Nphy), Zoop_0(Nzoo), PIC_0(Npim)
  real(8) :: d13C_TIC0, NO3_0, NH4_0, PO4_0, d15N_NO3_0, d15N_NH4_0
  real(8) :: d13C_DOC_0(Ndom), d13C_POC_0(Npom), d13C_Phyt_0(Nphy), d13C_Zoop_0(Nzoo)
  real(8) :: d13C_PIC_0(Npim)
  real(8) :: DON_0(Ndom), PON_0(Npom), DOP_0(Ndom), POP_0(Npom)
  real(8) :: d15N_DON_0(Ndom), d15N_PON_0(Npom), d15N_Phyt_0(Nphy), d15N_Zoop_0(Nzoo)
  real(8), allocatable :: MUD0(:), SAND0(:) 
  character(20) :: Flow_Label, T_Label, S_Label, TAlk_Label, TIC_Label, Oxyg_Label
  character(20) :: DOC_Label(Ndom), POC_Label(Npom), Phyt_Label(Nphy), Zoop_Label(Nzoo), PIC_Label(Npim)
  character(20) :: d13C_TIC_Label, NO3_Label, NH4_Label, PO4_Label, d15N_NO3_Label, d15N_NH4_Label
  character(20) :: d13C_DOC_Label(Ndom), d13C_POC_Label(Npom), d13C_Phyt_Label(Nphy), d13C_Zoop_Label(Nzoo)
  character(20) :: d13C_PIC_Label(Npim)
  character(20) :: DON_Label(Ndom), PON_Label(Npom), DOP_Label(Ndom), POP_Label(Npom)
  character(20) :: d15N_DON_Label(Ndom), d15N_PON_Label(Npom), d15N_Phyt_Label(Nphy), d15N_Zoop_Label(Nzoo)
  character(20), allocatable :: MUD_Label(:), SAND_Label(:) 
  real(8) :: Flow_sf, T_sf, S_sf, TAlk_sf, TIC_sf, Oxyg_sf
  real(8) :: DOC_sf(Ndom), POC_sf(Npom), Phyt_sf(Nphy), Zoop_sf(Nzoo), PIC_sf(Npim)
  real(8) :: d13C_TIC_sf, NO3_sf, NH4_sf, PO4_sf, d15N_NO3_sf, d15N_NH4_sf
  real(8) :: d13C_DOC_sf(Ndom), d13C_POC_sf(Npom), d13C_Phyt_sf(Nphy), d13C_Zoop_sf(Nzoo)
  real(8) :: d13C_PIC_sf(Npim)
  real(8) :: DON_sf(Ndom), PON_sf(Npom), DOP_sf(Ndom), POP_sf(Npom)
  real(8) :: d15N_DON_sf(Ndom), d15N_PON_sf(Npom), d15N_Phyt_sf(Nphy), d15N_Zoop_sf(Nzoo)
  real(8), allocatable :: MUD_sf(:), SAND_sf(:) 
  real(8) :: Flow_off, T_off, S_off, TAlk_off, TIC_off, Oxyg_off
  real(8) :: DOC_off(Ndom), POC_off(Npom), Phyt_off(Nphy), Zoop_off(Nzoo), PIC_off(Npim)
  real(8) :: d13C_TIC_off, NO3_off, NH4_off, PO4_off, d15N_NO3_off, d15N_NH4_off
  real(8) :: d13C_DOC_off(Ndom), d13C_POC_off(Npom), d13C_Phyt_off(Nphy), d13C_Zoop_off(Nzoo)
  real(8) :: d13C_PIC_off(Npim)
  real(8) :: DON_off(Ndom), PON_off(Npom), DOP_off(Ndom), POP_off(Npom)
  real(8) :: d15N_DON_off(Ndom), d15N_PON_off(Npom), d15N_Phyt_off(Nphy), d15N_Zoop_off(Nzoo)
  real(8), allocatable :: MUD_off(:), SAND_off(:) 

  integer :: iFlow, iT, iS, iTAlk, iTIC, iOxyg
  integer :: iDOC(Ndom), iPOC(Npom), iPhyt(Nphy), iZoop(Nzoo), iPIC(Npim)
  integer :: id13C_TIC, iNO3, iNH4, iPO4, id15N_NO3, id15N_NH4
  integer :: id13C_DOC(Ndom), id13C_POC(Npom), id13C_Phyt(Nphy), id13C_Zoop(Nzoo)
  integer :: id13C_PIC(Npim)
  integer :: iDON(Ndom), iPON(Npom), iDOP(Ndom), iPOP(Npom)
  integer :: id15N_DON(Ndom), id15N_PON(Npom), id15N_Phyt(Nphy), id15N_Zoop(Nzoo)
  integer, allocatable :: iMUD(:), iSAND(:) 

  integer :: spherical
  integer :: N_s_rho
  integer :: Vtransform, Vstretching
  real(8) :: THETA_S , THETA_B, TCLINE, DCRIT

  namelist/refdate/Ryear, Rmonth, Rday
  namelist/roms2roms/ROMS_HISFILE, romsvar
  namelist/river1/OUT_FILE, Isrc, Jsrc, river_dir, river_dir2, GLOBAL_ATT, NCS, NNS

  namelist/river2/RIVER_FILE, cha_name
  namelist/river2/Yield
  namelist/river2/LOCAL_TIME
  namelist/river2/Flow0, T0, S0, TAlk0, TIC_0, Oxyg0
  namelist/river2/DOC_0, POC_0, Phyt_0, Zoop_0, PIC_0
  namelist/river2/d13C_TIC0, NO3_0, NH4_0, PO4_0, d15N_NO3_0, d15N_NH4_0
  namelist/river2/d13C_DOC_0, d13C_POC_0, d13C_Phyt_0, d13C_Zoop_0
  namelist/river2/d13C_PIC_0
  namelist/river2/DON_0, PON_0, DOP_0, POP_0
  namelist/river2/d15N_DON_0, d15N_PON_0, d15N_Phyt_0, d15N_Zoop_0
  namelist/river2/MUD0, SAND0 
  namelist/river2/Flow_Label, T_Label, S_Label, TAlk_Label, TIC_Label, Oxyg_Label
  namelist/river2/DOC_Label, POC_Label, Phyt_Label, Zoop_Label, PIC_Label
  namelist/river2/d13C_TIC_Label, NO3_Label, NH4_Label, PO4_Label, d15N_NO3_Label, d15N_NH4_Label
  namelist/river2/d13C_DOC_Label, d13C_POC_Label, d13C_Phyt_Label, d13C_Zoop_Label
  namelist/river2/d13C_PIC_Label
  namelist/river2/DON_Label, PON_Label, DOP_Label, POP_Label
  namelist/river2/d15N_DON_Label, d15N_PON_Label, d15N_Phyt_Label, d15N_Zoop_Label
  namelist/river2/MUD_Label, SAND_Label 
  namelist/river2/Flow_sf, T_sf, S_sf, TAlk_sf, TIC_sf, Oxyg_sf
  namelist/river2/DOC_sf, POC_sf, Phyt_sf, Zoop_sf, PIC_sf
  namelist/river2/d13C_TIC_sf, NO3_sf, NH4_sf, PO4_sf, d15N_NO3_sf, d15N_NH4_sf
  namelist/river2/d13C_DOC_sf, d13C_POC_sf, d13C_Phyt_sf, d13C_Zoop_sf
  namelist/river2/d13C_PIC_sf
  namelist/river2/DON_sf, PON_sf, DOP_sf, POP_sf
  namelist/river2/d15N_DON_sf, d15N_PON_sf, d15N_Phyt_sf, d15N_Zoop_sf
  namelist/river2/MUD_sf, SAND_sf 
  namelist/river2/Flow_off, T_off, S_off, TAlk_off, TIC_off, Oxyg_off
  namelist/river2/DOC_off, POC_off, Phyt_off, Zoop_off, PIC_off
  namelist/river2/d13C_TIC_off, NO3_off, NH4_off, PO4_off, d15N_NO3_off, d15N_NH4_off
  namelist/river2/d13C_DOC_off, d13C_POC_off, d13C_Phyt_off, d13C_Zoop_off
  namelist/river2/d13C_PIC_off
  namelist/river2/DON_off, PON_off, DOP_off, POP_off
  namelist/river2/d15N_DON_off, d15N_PON_off, d15N_Phyt_off, d15N_Zoop_off
  namelist/river2/MUD_off, SAND_off 

  namelist/zcoord/N_s_rho
  namelist/zcoord/Vtransform, Vstretching
  namelist/zcoord/THETA_S, THETA_B, TCLINE, DCRIT

  ! Read parameters in namelist file
  
  read (*, nml=refdate)
  rewind(5)
  read (*, nml=roms2roms)
  rewind(5)
  read (*, nml=river1)

  NCS2 = NCS
  NNS2 = NNS
  if(NCS == 0 ) NCS2 = 1
  if(NNS == 0 ) NNS2 = 1
  allocate(  MUD0(NCS2), SAND0(NNS2) )
  allocate(  MUD_Label(NCS2), SAND_Label(NNS2) )
  allocate(  MUD_sf(NCS2), SAND_sf(NNS2) )
  allocate(  MUD_off(NCS2), SAND_off(NNS2) )
  allocate(  iMUD(NCS2), iSAND(NNS2) )

  rewind(5)
  read (*, nml=river2)
  rewind(5)
  read (*, nml=zcoord)


  Nzr = N_s_rho
  Nsrc = 1   !!!! 1 <- river number
  allocate( river_Vshape(Nsrc,Nzr) )
  allocate( data(Nsrc,Nzr,1) )
  allocate( river_flow(Nsrc,1) )
  allocate( river_flow2(Nsrc,1) )

  river_Vshape(1,:) = 1.0d0/dble(Nzr)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  data(:,:,:) = 0.0d0

!-Modify time-unit description ---------------------------------
  
  write (YYYY, "(I4.4)") Ryear
  write (MM, "(I2.2)") Rmonth
  write (DD, "(I2.2)") Rday
  TIME_ATT(12:15)=YYYY
  TIME_ATT(17:18)=MM
  TIME_ATT(20:21)=DD

  call jd(Ryear, Rmonth, Rday, jdate_Ref)
  d_jdate_Ref = dble(jdate_Ref)
!  write(*,*) d_jdate_Ref

  
!  OUT_FILE = trim( OUT_prefix ) //trim( OUT_suffix )

  call createNetCDFriver( trim( OUT_FILE ), trim( GLOBAL_ATT )     &
        , TIME_ATT, 1, N_s_rho, romsvar, NCS, NNS  )

  write(*,*) "WRITE: ", trim( OUT_FILE )
  call check( nf90_open(trim( OUT_FILE ), NF90_WRITE, ncid) )

  write(*,*)  'Write: river'
  call check( nf90_inq_varid(ncid, 'river', var_id) )
  call check( nf90_put_var(ncid, var_id, 1 ) )  !!!! river number
  write(*,*)  'Write: river_direction'
  call check( nf90_inq_varid(ncid, 'river_direction', var_id) )
  call check( nf90_put_var(ncid, var_id, river_dir) )
  write(*,*)  'Write: river_Xposition'
  call check( nf90_inq_varid(ncid, 'river_Xposition', var_id) )
  call check( nf90_put_var(ncid, var_id, Isrc ) )
  write(*,*)  'Write: river_Eposition'
  call check( nf90_inq_varid(ncid, 'river_Eposition', var_id) )
  call check( nf90_put_var(ncid, var_id, Jsrc ) )
  write(*,*)  'Write: river_Vshape'
  call check( nf90_inq_varid(ncid, 'river_Vshape', var_id) )
  call check( nf90_put_var(ncid, var_id, river_Vshape) )


!-Read SWAT+ file ---------------------------------

  write(*,*) 'Read SWAT+ output file'
  write(*,*) 'Open: ', trim( RIVER_FILE )
  open(unit=20, file=trim( RIVER_FILE ), action='read')

 !----- SWAT+ channel_day.txt output file ------------
#if defined CHANNEL_DAY
 !----- SWAT+ Rev 2019.59.2/3; channel_day.txt/channl_sd_day.txt output ------------
  ! read SWAT+ file header
  write(*,*) 'Read Header: ', trim( RIVER_FILE )
  read(20,'(a)') headerline
  write(*,*) trim( headerline )
  ! Check data label
  read(20,'(a)') headerline
  write(*,*) trim( headerline )

  ! Read label data
  read(headerline,*) label 
!  write(*,*) label

# if defined SWAT_PLUS_REV_2019_59_2 || defined SWAT_PLUS_OKMT

  Ncol = 59  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate( label(Ncol) )
  allocate( label2(7) )

  ! Read label data
  read(headerline,*) label2, label 
!  write(*,*) label2
!  write(*,*) label

# elif defined SWAT_PLUS_REV_2019_59_3

  Ncol = 58  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate( label(Ncol) )

  ! Read label data
  read(headerline,*) label 
!  write(*,*) label

# endif

  ! Read unit line
  read(20,'(a)') headerline
  write(*,*) trim( headerline )

  do i=1, Ncol
    write(*,*) i,trim(label(i))    
  enddo


!----- SWAT+ subdaity output ------------

#elif defined CHANNEL_SUBDAY
# if defined SWAT_PLUS_REV_2019_59_2 || defined SWAT_PLUS_REV_2019_59_3

  write(*,*) "STOP: SWAT_PLUS_REV_2019_59_2/3 with CHANNEL_SUBDAY is not supported..."
  close(20)
  stop

# elif defined SWAT_PLUS_OKMT
 !----- SWAT+ Mr. Okamoto version; subdaity output ------------
  ! read SWAT+ file header
  write(*,*) 'Read Header: ', trim( RIVER_FILE )
  read(20,'(a)') headerline
  write(*,*) trim( headerline )

  Ncol = 20  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate( label(Ncol) )

  do i=1, Ncol
    label(i) = trim( adjustl(headerline( (31+(i-1)*16):(31+i*16)) ) )
    write(*,*) i,trim(label(i))    
  enddo


# endif
#endif

  allocate( river_data(Ncol) )

! ==== Check columun number of the parameter ======================
  ! Flow
  if (trim(Flow_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(Flow_Label) == trim( label(i) ) ) then
        iFlow = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(Flow_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iFlow: ", iFlow
  endif

  ! temp
  if (trim(T_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(T_Label) == trim( label(i) ) ) then
        iT = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(T_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iT: ", iT
  endif

  ! Salt
  if (trim(S_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(S_Label) == trim( label(i) ) ) then
        iS = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(S_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iS: ", iS
  endif
  
  ! mud_
  do j=1,NCS
    if ( trim(MUD_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(MUD_Label(j)) == trim( label(i) ) ) then
          iMUD(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(MUD_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iMUD: ",j, iMUD(j)
    endif
  enddo

  ! sand_
  do j=1,NNS
    if ( trim(SAND_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(SAND_Label(j)) == trim( label(i) ) ) then
          iSAND(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(SAND_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iSAND: ",j, iSAND(j)
    endif
  enddo
!  write(*,*) "iSAND(:): ", iSAND

  ! TIC
  if (trim(TIC_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(TIC_Label) == trim( label(i) ) ) then
        iTIC = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(TIC_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iTIC: ", iTIC
  endif

  ! alkalinity
  if (trim(TAlk_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(TAlk_Label) == trim( label(i) ) ) then
        iTAlk = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(TAlk_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iTAlk: ", iTAlk
  endif

  ! oxygen
  if (trim(Oxyg_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(Oxyg_Label) == trim( label(i) ) ) then
        iOxyg = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(Oxyg_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iOxyg: ", iOxyg
  endif

  ! DOC_
  do j=1,Ndom
    if ( trim(DOC_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(DOC_Label(j)) == trim( label(i) ) ) then
          iDOC(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(DOC_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iDOC: ",j, iDOC(j)
    endif
  enddo

  ! POC_
  do j=1,Npom
    if ( trim(POC_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(POC_Label(j)) == trim( label(i) ) ) then
          iPOC(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(POC_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iPOC: ",j, iPOC(j)
    endif
  enddo

  ! phytoplankton_
  do j=1,Nphy
    if ( trim(Phyt_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(Phyt_Label(j)) == trim( label(i) ) ) then
          iPhyt(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(Phyt_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iPhyt: ",j, iPhyt(j)
    endif
  enddo

  ! zooplankton_
  do j=1,Nzoo
    if ( trim(Zoop_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(Zoop_Label(j)) == trim( label(i) ) ) then
          iZoop(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(Zoop_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iZoop: ",j, iZoop(j)
    endif
  enddo

  ! PIC_
  do j=1,Npim
    if ( trim(PIC_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(PIC_Label(j)) == trim( label(i) ) ) then
          iPIC(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(PIC_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iPIC: ",j, iPIC(j)
    endif
  enddo

  ! NO3
  if (trim(NO3_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(NO3_Label) == trim( label(i) ) ) then
        iNO3 = i
        exit
      endif
     if(i==Ncol) then
        write(*,*) "Could not find "//trim(NO3_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iNO3: ", iNO3
  endif

  ! NH4
  if (trim(NH4_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(NH4_Label) == trim( label(i) ) ) then
        iNH4 = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(NH4_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iNH4: ", iNH4
  endif

  ! PO4
  if (trim(PO4_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(PO4_Label) == trim( label(i) ) ) then
        iPO4 = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(PO4_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "iPO4: ", iPO4
  endif

  ! DON_
  do j=1,Ndom
    if ( trim(DON_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(DON_Label(j)) == trim( label(i) ) ) then
          iDON(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(DON_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iDON: ",j, iDON(j)
    endif
  enddo

  ! PON_
  do j=1,Npom
    if ( trim(PON_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(PON_Label(j)) == trim( label(i) ) ) then
          iPON(j) = i
          exit
        endif
      enddo
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(PON_Label(j))//" in "//trim( RIVER_FILE )
      endif
      write(*,*) "iPON: ",j, iPON(j)
    endif
  enddo

  ! DOP_
  do j=1,Ndom
    if ( trim(DOP_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(DOP_Label(j)) == trim( label(i) ) ) then
          iDOP(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(DOP_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iDOP: ",j, iDOP(j)
    endif
  enddo

  ! POP_
  do j=1,Npom
    if ( trim(POP_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(POP_Label(j)) == trim( label(i) ) ) then
          iPOP(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(POP_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "iPOP: ",j, iPOP(j)
    endif
  enddo

  ! TI13C
  if (trim(d13C_TIC_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(d13C_TIC_Label) == trim( label(i) ) ) then
        id13C_TIC = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(d13C_TIC_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "id13C_TIC: ", id13C_TIC
  endif

  ! DO13C_
  do j=1,Ndom
    if ( trim(d13C_DOC_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d13C_DOC_Label(j)) == trim( label(i) ) ) then
          id13C_DOC(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d13C_DOC_Label(j))//" in "//trim( RIVER_FILE )
        endif
      enddo
      write(*,*) "id13C_DOC: ",j, id13C_DOC(j)
    endif
  enddo

  ! PO13C_
  do j=1,Npom
    if ( trim(d13C_POC_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d13C_POC_Label(j)) == trim( label(i) ) ) then
          id13C_POC(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d13C_POC_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id13C_POC: ",j, id13C_POC(j)
    endif
  enddo

  ! phyt13C_
  do j=1,Nphy
    if ( trim(d13C_Phyt_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d13C_Phyt_Label(j)) == trim( label(i) ) ) then
          id13C_Phyt(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d13C_Phyt_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id13C_Phyt: ",j, id13C_Phyt(j)
    endif
  enddo

  ! zoop13C_
  do j=1,Nzoo
    if ( trim(d13C_Zoop_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d13C_Zoop_Label(j)) == trim( label(i) ) ) then
          id13C_Zoop(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d13C_Zoop_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id13C_Zoop: ",j, id13C_Zoop(j)
    endif
  enddo

  ! PI13C_
  do j=1,Npim
    if ( trim(d13C_PIC_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d13C_PIC_Label(j)) == trim( label(i) ) ) then
          id13C_PIC(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d13C_PIC_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id13C_PIC: ",j, id13C_PIC(j)
    endif
  enddo

  ! 15NO3
  if (trim(d15N_NO3_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(d15N_NO3_Label) == trim( label(i) ) ) then
        id15N_NO3 = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(d15N_NO3_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "id15N_NO3: ", id15N_NO3
  endif

  ! 15NH4
  if (trim(d15N_NH4_Label) /= "Null" ) then
    do i=1, Ncol
      if( trim(d15N_NH4_Label) == trim( label(i) ) ) then
        id15N_NH4 = i
        exit
      endif
      if(i==Ncol) then
        write(*,*) "Could not find "//trim(d15N_NH4_Label)//" in "//trim( RIVER_FILE )
        stop
      endif
    enddo
    write(*,*) "id15N_NH4: ", id15N_NH4
 endif

  ! DO15N_
  do j=1,Ndom
    if ( trim(d15N_DON_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d15N_DON_Label(j)) == trim( label(i) ) ) then
          id15N_DON(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d15N_DON_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id15N_DON: ",j, id15N_DON(j)
    endif
  enddo

  ! PO15N_
  do j=1,Npom
    if ( trim(d15N_PON_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d15N_PON_Label(j)) == trim( label(i) ) ) then
          id15N_PON(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d15N_PON_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id15N_PON: ",j, id15N_PON(j)
    endif
  enddo

  ! phyt15N_
  do j=1,Nphy
    if ( trim(d15N_Phyt_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d15N_Phyt_Label(j)) == trim( label(i) ) ) then
          id15N_Phyt(j) = i
          exit
        endif
        if(i==Ncol) then
          write(*,*) "Could not find "//trim(d15N_Phyt_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id15N_Phyt: ",j, id15N_Phyt(j)
    endif
  enddo

  ! zoop15N_
  do j=1,Nzoo
    if ( trim(d15N_Zoop_Label(j)) /= "Null" ) then
      do i=1, Ncol
        if( trim(d15N_Zoop_Label(j)) == trim( label(i) ) ) then
          id15N_Zoop(j) = i
          exit
        endif
         if(i==Ncol) then
          write(*,*) "Could not find "//trim(d15N_Zoop_Label(j))//" in "//trim( RIVER_FILE )
          stop
        endif
      enddo
      write(*,*) "id15N_Zoop: ",j, id15N_Zoop(j)
    endif
  enddo

!--- Preparation for Loop ---------------------------------------
  start1D = (/ 1 /)
  count1D = (/ 1 /)
  start3D = (/ 1, 1,   itime /)
  count3D = (/ 1, Nzr, 1     /)
  itime = 0
! === LOOP start ==================================================================
  DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined SWAT_PLUS_REV_2019_59_3 || defined SWAT_PLUS_REV_2019_59_2
# if defined CHANNEL_DAY
  !----- SWAT+ Rev 2019.59.2/3; channel_day.txt/channl_sd_day.txt output ------------
    read(20,*,iostat=ios) jday, imonth, iday, iyear, iunit, gis_id, name, river_data
    if(ios==-1) exit  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!    write(*,*) label
!    write(*,*) river_data
    if( trim( name ) /= trim( cha_name )) cycle !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) trim( name )

# elif defined CHANNEL_SUBDAY



# endif

#elif defined SWAT_PLUS_OKMT
# if defined CHANNEL_DAY


# elif defined CHANNEL_SUBDAY
  !----- SWAT+ Mr. Okamoto version; subdaily output ------------
    read(20,*,iostat=ios) name, iyear, imonth, iday, ihour, river_data
    if(ios==-1) exit
    if( trim( name ) /= trim( cha_name )) cycle
    write(*,*) trim( name )
# endif
#endif

    itime = itime + 1

    call jd(iyear, imonth, iday, jdate)
    d_jdate = dble(jdate)

  ! river_time
    river_time(1) = d_jdate - d_jdate_Ref - dble(LOCAL_TIME)/24.0d0

#if defined CHANNEL_DAY
    river_time(1) = river_time(1) - 0.5d0
#elif defined CHANNEL_SUBDAY
    river_time(1) = river_time(1) + ( dble(ihour) - 0.5d0 )/24.0d0
#endif

    start1D = (/ itime /)
    count1D = (/ 1     /) 
    write(*,*)  'Write: river_time', river_time, TIME_ATT
    call check( nf90_inq_varid(ncid, 'river_time', var_id) )
    call check( nf90_put_var(ncid, var_id, river_time, start = start1D, count = count1D) )

  ! river_transport
    start2D = (/ 1, itime /)
    count2D = (/ 1, 1     /)
    if( trim(Flow_Label) == "Null" ) then
      river_flow(1,1) = Flow0
    else
      river_flow(1,1) = river_data(iflow) * Flow_sf + Flow_off
    endif
    river_flow2(1,1) = river_flow(1,1) * dble( river_dir2 )

    write(*,*)  'Write: river_transport', river_flow2(1,1)
    call check( nf90_inq_varid(ncid, 'river_transport', var_id) )
    call check( nf90_put_var(ncid, var_id, river_flow2, start = start2D, count = count2D) )

    if(river_flow(1,1)<1.0d-20) river_flow(1,1)=1.0d-20 !!! Error handling


    start3D = (/ 1, 1,   itime /)
    count3D = (/ 1, Nzr, 1     /)
    
  ! temp
    i=6
    if (romsvar( i ) == 1 ) then

      if( trim(T_Label) == "Null" ) then
        data(1,:,:) = T0
      else
        data(1,:,:) = river_data(iT) * T_sf + T_off
      endif
  
      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1)        
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! Salt
    i=7
    if (romsvar( i ) == 1 ) then

      if( trim(S_Label) == "Null" ) then
        data(1,:,:) = S0
      else
        data(1,:,:) = river_data(iS) * S_sf + S_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif
      endif
  
      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1)       
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! mud_
    i=8
    if (romsvar( i ) == 1 ) then
      do j=1,NCS

        if( trim(MUD_Label(j)) == "Null" ) then
          data(1,:,:) = MUD0(j)
        else
          data(1,:,:) = river_data(iMUD(j)) * MUD_sf(j) + MUD_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif 
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! sand_
    i=9
    if (romsvar( i ) == 1 ) then
      do j=1,NNS

        if( trim(SAND_Label(j)) == "Null" ) then
          data(1,:,:) = SAND0(j)
        else
          data(1,:,:) = river_data(iSAND(j)) * SAND_sf(j) + SAND_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif 
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! TIC
    i=10
    if (romsvar( i ) == 1 ) then

      if( trim(TIC_Label) == "Null" ) then
        data(1,:,:) = TIC_0
      else
        data(1,:,:) = river_data(iTIC) * TIC_sf + TIC_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif
      endif
  
      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! alkalinity
    i=11
    if (romsvar( i ) == 1 ) then

      if( trim(TAlk_Label) == "Null" ) then
        data(1,:,:) = TAlk0
      else
        data(1,:,:) = river_data(iTAlk) * TAlk_sf + TAlk_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! oxygen
    i=12
    if (romsvar( i ) == 1 ) then

      if( trim(Oxyg_Label) == "Null" ) then
        data(1,:,:) = Oxyg0
      else
        data(1,:,:) = river_data(iOxyg) * Oxyg_sf + Oxyg_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! DOC_
    i=13
    if (romsvar( i ) == 1 ) then
      do j=1,Ndom

        if( trim(DOC_Label(j)) == "Null" ) then
          data(1,:,:) = DOC_0(j)
        else
          data(1,:,:) = river_data(iDOC(j)) * DOC_sf(j) + DOC_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! POC_
    i=14
    if (romsvar( i ) == 1 ) then
      do j=1,Npom

        if( trim(POC_Label(j)) == "Null" ) then
          data(1,:,:) = POC_0(j)
        else
          data(1,:,:) = river_data(iPOC(j)) * POC_sf(j) + POC_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! phytoplankton_
    i=15
    if (romsvar( i ) == 1 ) then
      do j=1,Nphy

        if( trim(Phyt_Label(j)) == "Null" ) then
          data(1,:,:) = Phyt_0(j)
        else
          data(1,:,:) = river_data(iPhyt(j)) * Phyt_sf(j) + Phyt_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! zooplankton_
    i=16
    if (romsvar( i ) == 1 ) then
      do j=1,Nzoo

        if( trim(Zoop_Label(j)) == "Null" ) then
          data(1,:,:) = Zoop_0(j)
        else
          data(1,:,:) = river_data(iZoop(j)) * Zoop_sf(j) + Zoop_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! PIC_
    i=17
    if (romsvar( i ) == 1 ) then
      do j=1,Npim

        if( trim(PIC_Label(j)) == "Null" ) then
          data(1,:,:) = PIC_0(j)
        else
          data(1,:,:) = river_data(iPIC(j)) * PIC_sf(j) + PIC_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! NO3
    i=18
    if (romsvar( i ) == 1 ) then

      if( trim(NO3_Label) == "Null" ) then
        data(1,:,:) = NO3_0
      else
        data(1,:,:) = river_data(iNO3) * NO3_sf + NO3_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif  
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! NH4
    i=19
    if (romsvar( i ) == 1 ) then

      if( trim(NH4_Label) == "Null" ) then
        data(1,:,:) = NH4_0
      else
        data(1,:,:) = river_data(iNH4) * NH4_sf + NH4_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif  
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! PO4
    i=20
    if (romsvar( i ) == 1 ) then

      if( trim(PO4_Label) == "Null" ) then
        data(1,:,:) = PO4_0
      else
        data(1,:,:) = river_data(iPO4) * PO4_sf + PO4_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif  
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! DON_
    i=21
    if (romsvar( i ) == 1 ) then
      do j=1,Ndom

        if( trim(DON_Label(j)) == "Null" ) then
          data(1,:,:) = DON_0(j)
        else
          data(1,:,:) = river_data(iDON(j)) * DON_sf(j) + DON_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! PON_
    i=22
    if (romsvar( i ) == 1 ) then
      do j=1,Npom

        if( trim(PON_Label(j)) == "Null" ) then
          data(1,:,:) = PON_0(j)
        else
          data(1,:,:) = river_data(iPON(j)) * PON_sf(j) + PON_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! DOP_
    i=23
    if (romsvar( i ) == 1 ) then
      do j=1,Ndom

        if( trim(DOP_Label(j)) == "Null" ) then
          data(1,:,:) = DOP_0(j)
        else
          data(1,:,:) = river_data(iDOP(j)) * DOP_sf(j) + DOP_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! POP_
    i=24
    if (romsvar( i ) == 1 ) then
      do j=1,Npom

        if( trim(POP_Label(j)) == "Null" ) then
          data(1,:,:) = POP_0(j)
        else
          data(1,:,:) = river_data(iPOP(j)) * POP_sf(j) + POP_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif


  ! TI13C
    i=25
    if (romsvar( i ) == 1 ) then

      if( trim(d13C_TIC_Label) == "Null" ) then
        data(1,:,:) = d13C_TIC0
      else
        data(1,:,:) = river_data(id13C_TIC) * d13C_TIC_sf + d13C_TIC_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif  
      endif
  
      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! DO13C_
    i=26
    if (romsvar( i ) == 1 ) then
      do j=1,Ndom

        if( trim(d13C_DOC_Label(j)) == "Null" ) then
          data(1,:,:) = d13C_DOC_0(j)
        else
          data(1,:,:) = river_data(id13C_DOC(j)) * d13C_DOC_sf(j) + d13C_DOC_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! PO13C_
    i=27
    if (romsvar( i ) == 1 ) then
      do j=1,Npom

        if( trim(d13C_POC_Label(j)) == "Null" ) then
          data(1,:,:) = d13C_POC_0(j)
        else
          data(1,:,:) = river_data(id13C_POC(j)) * d13C_POC_sf(j) + d13C_POC_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! phyt13C_
    i=28
    if (romsvar( i ) == 1 ) then
      do j=1,Nphy

        if( trim(d13C_Phyt_Label(j)) == "Null" ) then
          data(1,:,:) = d13C_Phyt_0(j)
        else
          data(1,:,:) = river_data(id13C_Phyt(j)) * d13C_Phyt_sf(j) + d13C_Phyt_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! zoop13C_
    i=29
    if (romsvar( i ) == 1 ) then
      do j=1,Nzoo

        if( trim(d13C_Zoop_Label(j)) == "Null" ) then
          data(1,:,:) = d13C_Zoop_0(j)
        else
          data(1,:,:) = river_data(id13C_Zoop(j)) * d13C_Zoop_sf(j) + d13C_Zoop_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! PI13C_
    i=30
    if (romsvar( i ) == 1 ) then
      do j=1,Npim

        if( trim(d13C_PIC_Label(j)) == "Null" ) then
          data(1,:,:) = d13C_PIC_0(j)
        else
          data(1,:,:) = river_data(id13C_PIC(j)) * d13C_PIC_sf(j) + d13C_PIC_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! 15NO3
    i=31
    if (romsvar( i ) == 1 ) then

      if( trim(d15N_NO3_Label) == "Null" ) then
        data(1,:,:) = d15N_NO3_0
      else
        data(1,:,:) = river_data(id15N_NO3) * d15N_NO3_sf + d15N_NO3_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif  
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ) , data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! 15NH4
    i=32
    if (romsvar( i ) == 1 ) then

      if( trim(d15N_NH4_Label) == "Null" ) then
        data(1,:,:) = d15N_NH4_0
      else
        data(1,:,:) = river_data(id15N_NH4) * d15N_NH4_sf + d15N_NH4_off
        if(Yield) then
          data(1,:,:) = data(1,:,:)/river_flow(1,1)
        endif  
      endif

      varname = 'river_'//trim( VAR_NAME(i) )
      write(*,*) 'Write: ', trim( varname  ), data(1,1,1) 
      call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
    endif

  ! DO15N_
    i=33
    if (romsvar( i ) == 1 ) then
      do j=1,Ndom

        if( trim(d15N_DON_Label(j)) == "Null" ) then
          data(1,:,:) = d15N_DON_0(j)
        else
          data(1,:,:) = river_data(id15N_DON(j)) * d15N_DON_sf(j) + d15N_DON_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! PO15N_
    i=34
    if (romsvar( i ) == 1 ) then
      do j=1,Npom

        if( trim(d15N_PON_Label(j)) == "Null" ) then
          data(1,:,:) = d15N_PON_0(j)
        else
          data(1,:,:) = river_data(id15N_PON(j)) * d15N_PON_sf(j) + d15N_PON_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! phyt15N_
    i=35
    if (romsvar( i ) == 1 ) then
      do j=1,Nphy

        if( trim(d15N_Phyt_Label(j)) == "Null" ) then
          data(1,:,:) = d15N_Phyt_0(j)
        else
          data(1,:,:) = river_data(id15N_Phyt(j)) * d15N_Phyt_sf(j) + d15N_Phyt_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif
    
        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

  ! zoop15N_
    i=36
    if (romsvar( i ) == 1 ) then
      do j=1,Nzoo

        if( trim(d15N_Zoop_Label(j)) == "Null" ) then
          data(1,:,:) = d15N_Zoop_0(j)
        else
          data(1,:,:) = river_data(id15N_Zoop(j)) * d15N_Zoop_sf(j) + d15N_Zoop_off(j)
          if(Yield) then
            data(1,:,:) = data(1,:,:)/river_flow(1,1)
          endif  
        endif

        write(varnum,'(I2.2)') j
        varname = 'river_'//trim( VAR_NAME(i) )//varnum
        write(*,*) 'Write: ', trim( varname ), data(1,1,1) 
        call check( nf90_inq_varid(ncid, trim( varname ), var_id) )
        call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      enddo
    endif

!    stop !!!!!!!!!!!!!!!!! For debug

  END DO  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  close(20)
  call check( nf90_close(ncid) )

  write(*,*) 'FINISH!!'    
      
END PROGRAM riverSWAT2ROMS
      
