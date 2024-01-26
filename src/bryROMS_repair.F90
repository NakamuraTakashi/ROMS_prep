
!!!=== Copyright (c) 2024 Takashi NAKAMURA  =====

PROGRAM bryROMS_repair
  use netcdf
  use mod_utility
  use mod_roms_netcdf
 
  implicit none
  
! -------------------------------------------------------------------------
  character(256) :: NC_FILE
  integer :: romsvar(7)
  integer :: SNWE(4)
  integer :: itime
  real(8) :: StrangeValue
!  integer :: icoord  ! 1: rho point, 2: u point, 3: v point, 4: psi point
! -------------------------------------------------------------------------
  character(256) :: BRY_VARNAME
  integer :: start2D(2), count2D(2)
  integer :: start3D(3), count3D(3)
  
  real(8), allocatable :: data1(:,:)
  real(8), allocatable :: data2(:,:,:)
  
  integer :: i,j,k
  integer :: ibry, ivar
  integer :: is

  integer :: Nxr, Nyr
  integer :: Nxp, Nyp
  integer :: Nxu, Nyu
  integer :: Nxv, Nyv
  integer :: Nzr
  integer :: L, M, N
  
  integer :: ncid,var_id

  namelist/bry_repair/NC_FILE,romsvar,SNWE,itime,StrangeValue

  read (5, nml=bry_repair)
  rewind(5)
!-Read ROMS netCDF file --------------------------------
  write(*,*) "OPEN: ", trim( NC_FILE )

  call check( nf90_open(trim( NC_FILE ), NF90_WRITE, ncid) )
  ! Get dimension data
  call get_dimension(ncid, 'xi_rho',  Nxr)
  call get_dimension(ncid, 'eta_rho', Nyr)
  call get_dimension(ncid, 's_rho',   Nzr)

  ! Open NetCDF grid file
  Nxu = Nxr-1
  Nyu = Nyr
  Nxv = Nxr
  Nyv = Nyr-1
  Nxp = Nxr-1
  Nyp = Nyr-1
  N = Nzr

! ==== Loop start =========================================
  DO ibry=1, 4

    if(SNWE(ibry)==0) cycle

    if(ibry==1 .or. ibry==2  ) then
      L = Nxr-1
    else
      L = Nyr-1
    endif

    DO ivar=1,7
      if(romsvar(ivar)==0) cycle

      is = 0
      if( (ibry==1 .or. ibry==2) .and. (ivar==2 .or. ivar==4) ) is = 1
      if( (ibry==3 .or. ibry==4) .and. (ivar==3 .or. ivar==5) ) is = 1

      BRY_VARNAME = trim(VAR_NAME(ivar))//'_'//trim(BRY_NAME(ibry))
    
      IF(ivar==1.or.ivar==4.or.ivar==5) THEN

        write(*,*) "CHECK1: ", trim( BRY_VARNAME ),is,L,N
  
        allocate( data1(is:L, 3) )
      
  !     Get variables
        start2D = (/ 1,     itime-1 /)
        count2D = (/ L+1-is, 3       /)
        call check( nf90_inq_varid(ncid, trim( BRY_VARNAME ), var_id) ) 
        call check( nf90_get_var(ncid, var_id, data1, start = start2D, count = count2D ) )
  
  !     Check and repair values
        do i=is, L
          if( StrangeValue==0.0d0 ) then
            data1(i,2) = 0.5d0*( data1(i,1) + data1(i,3) )
          else if(data1(i,2) <= -StrangeValue .or. data1(i,2) >= StrangeValue) then
            data1(i,2) = 0.5d0*( data1(i,1) + data1(i,3) )
          endif
        enddo
    
  !-  --- Write ROMS grid netCDF file --------------------------------
        write(*,*)  'Write: ', trim( BRY_VARNAME )
        call check( nf90_inq_varid(ncid, trim( BRY_VARNAME ), var_id ) )
        call check( nf90_put_var(ncid, var_id, data1, start = start2D, count = count2D) )
  
        deallocate( data1 )
        
      ELSE
  
        write(*,*) "CHECK2: ", trim( BRY_VARNAME ),is,L,N

        allocate( data2(is:L, 1:N, 3) )
      
  !     Get variables
        start3D = (/ 1,      1,   itime-1 /)
        count3D = (/ L+1-is, Nzr, 3       /)
        call check( nf90_inq_varid(ncid, trim( BRY_VARNAME ), var_id) ) 
        call check( nf90_get_var(ncid, var_id, data2, start = start3D, count = count3D ) )
  !     Check and repair values
        do k=1, N
          do i=is, L
            if( StrangeValue==0.0d0 ) then
              data2(i,k,2) = 0.5d0*( data2(i,k,1) + data2(i,k,3) )
            else if(data2(i,k,2) <= -StrangeValue .or. data2(i,k,2) >= StrangeValue) then
              data2(i,k,2) = 0.5d0*( data2(i,k,1) + data2(i,k,3) )
            endif
          enddo
        enddo
    
  !-  --- Write ROMS grid netCDF file --------------------------------
        write(*,*)  'Write: ', trim( BRY_VARNAME )
        call check( nf90_inq_varid(ncid, trim( BRY_VARNAME ), var_id ) )
        call check( nf90_put_var(ncid, var_id, data2, start = start3D, count = count3D) )
  
        deallocate( data2 )
  
      END IF

    END DO

  END DO
!=======================================================================
  call check( nf90_close(ncid) ) 
  write(*,*) 'FINISH!!'
      
END PROGRAM bryROMS_repair
      
