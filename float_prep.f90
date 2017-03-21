
!!!=== ver 2014/10/15   Copyright (c) 2014 Takashi NAKAMURA  =====

    PROGRAM float_prep
      use netcdf
      use mod_calendar
     
      implicit none
      
! Initial floats locations for all grids:
!
!   G      Nested grid number
!   C      Initial horizontal coordinate type (0: grid units, 1: spherical)
!   T      Float trajectory type (1: Lagrangian, 2: isobaric, 3: Geopotential)
!   N      Number floats to be released at (Fx0,Fy0,Fz0)
!   Ft0    Float release time (days) after model initialization
!   Fx0    Initial float X-location (grid units or longitude)
!   Fy0    Initial float Y-location (grid units or latitude)
!   Fz0    Initial float Z-location (grid units or depth)
!   Fdt    Float cluster release time interval (days)
!   Fdx    Float cluster X-distribution parameter
!   Fdy    Float cluster Y-distribution parameter
!   Fdz    Float cluster Z-distribution parameter
!
! POS = G, C, T, N,   Ft0,    Fx0,    Fy0,    Fz0,    Fdt,    Fdx,    Fdy,   Fdz
!
      integer, parameter :: G   = 1      ! Nested grid number
      integer, parameter :: C   = 0      ! Initial horizontal coordinate type (0: grid units, 1: spherical)
      integer, parameter :: T   = 1      ! Float trajectory type (1: Lagrangian, 2: isobaric, 3: Geopotential)
      integer, parameter :: N   = 1      ! Number floats to be released at (Fx0,Fy0,Fz0)
      real(8), parameter :: Ft0 = 0.0d0  ! Float release time (days) after model initialization
      real(8)            :: Fx0          ! Initial float X-location (grid units or longitude)
      real(8)            :: Fy0          ! Initial float Y-location (grid units or latitude)
      real(8), parameter :: Fz0 = 8.0d0 !!!0.01d0 ! Initial float Z-location (grid units or depth)
      real(8), parameter :: Fdt = 0.0d0  ! Float cluster release time interval (days)
      real(8), parameter :: Fdx = 0.0d0  ! Float cluster X-distribution parameter
      real(8), parameter :: Fdy = 0.0d0  ! Float cluster Y-distribution parameter
      real(8), parameter :: Fdz = 0.0d0  ! Float cluster Z-distribution parameter
      
      real(8), parameter :: cf  = 4.0d0   ! Float number per coverage (0-1): [float num.] = int ( cf * [coverage (0-1)] )
      
      integer :: NFLOATS
      integer :: iFLOATS
      
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "shiraho_reef_grid10.nc"
      character(len=*), parameter :: OUT_FILE  = "float_points.dat"

      real(8), allocatable :: lat_rho(:, :)
      real(8), allocatable :: lon_rho(:, :)
      real(8), allocatable :: x_rho(:, :)
      real(8), allocatable :: y_rho(:, :)
      real(8), allocatable :: p_coral(:, :)
      real(8) :: dx, dy
      real(8) :: rnd
      
      integer :: ncid,var_id
      integer :: N_xi_rho, N_eta_rho
      integer :: Im, Jm, Nt
      
      integer :: i,j,k

      
!---- Read ROMS grid netCDF file --------------------------------
      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      
      allocate(lat_rho(N_xi_rho, N_eta_rho))
      allocate(lon_rho(N_xi_rho, N_eta_rho))
      allocate(x_rho(N_xi_rho, N_eta_rho))
      allocate(y_rho(N_xi_rho, N_eta_rho))
      allocate(p_coral(N_xi_rho, N_eta_rho))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'lat_rho', var_id) ) ! latitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lat_rho) )
      call check( nf90_inq_varid(ncid, 'lon_rho', var_id) ) ! longitude at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, lon_rho) )
      call check( nf90_inq_varid(ncid, 'x_rho', var_id) ) ! x at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, x_rho) )
      call check( nf90_inq_varid(ncid, 'y_rho', var_id) ) ! y at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, y_rho) )
      call check( nf90_inq_varid(ncid, 'p_coral', var_id) ) ! coral coverage at RHO-points (degree_east)
      call check( nf90_get_var(ncid, var_id, p_coral) )
      
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "*** SUCCESS read GRID file!"
      
!---- LOOP START --------------------------------

      write(*,*) "CREATE: ", OUT_FILE

      open(50, file=OUT_FILE)
      open(51, file='check.dat')
      
      NFLOATS = 0

      do j=2, N_eta_rho-1
        do i=2, N_xi_rho-1
        
          iFLOATS = int( cf * p_coral(i,j) )
          
          NFLOATS =NFLOATS + iFLOATS
          
          dx = (x_rho(i+1,j)-x_rho(i-1,j))*0.5d0
          dy = (y_rho(i,j+1)-y_rho(i,j-1))*0.5d0
          
          do k=1, iFLOATS
          
            call random_number(rnd)
            
!            Fx0 = x_rho(i,j)+(rnd-0.5d0)*dx
!            Fy0 = y_rho(i,j)+(rnd-0.5d0)*dy
            Fx0 = dble(i)+(rnd-0.5d0)
            Fy0 = dble(j)+(rnd-0.5d0)
            
            write(50,100) G, C, T, N, Ft0, Fx0, Fy0, Fz0, Fdt, Fdx, Fdy, Fdz
            write(51,*) Fx0, Fy0
          end do
          
        end do
      end do
      
      write(50,*) 'NFLOATS = ', NFLOATS
      write(*,*)  'NFLOATS = ', NFLOATS

      close(50)
      close(51)
      
      write(*,*) "*** SUCCESS create float points file!"

!---- Formats --------------------------------------------

! POS =       G,     C,     T,     N,    Ft0,      Fx0,     Fy0,    Fz0,     Fdt,      Fdx,     Fdy,    Fdz
100   FORMAT(I2,1x, I2,1x, I2,1x, I2,1x, f8.2,1x, f8.1, 1x,f8.1,1x,f6.3,1x, f8.2,1x, f8.1,1x, f8.1,1x, f6.3 )

!---- End of Main program --------------------------------------------

    CONTAINS
    
!---- NetCDF utility -------------------------------------------------
     
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension


      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
    END PROGRAM float_prep

!------------------------------------------------------------------------------
! Initial float location KEYWORDS.  The model variable name is not used as
! keyword in some parameters.  However, it namce is provided in brackets.
!------------------------------------------------------------------------------
!
!  G         Nested grid number [ng].
!
!  C         Initial horizontal location coordinate type [Fcoor]:
!
!              Fcoor = 0,  rho grid units
!                             0.5 =< Fx0 =< Lm(ng)+0.5,
!                             0.5 =< Fy0 =< Mm(ng)+0.5
!
!              Fcoor = 1,  Fx0 is longitude (west values are negative).
!                          Fy0 is latitude (south values are negative).
!
!  T         Float trajectory type [Ftype]:
!
!              Ftype = 1,  3D Lagrangian floats.
!                          (flt_Lagran)
!
!              Ftype = 2,  Isobaric floats, p=g*(z+zeta)=constant.
!                          (flt_Isobar)
!
!              Ftype = 3,  Geopotential floats, constant depth.
!                          (flt_Geopot)
!
!  N         Number of floats [Fcount] to be released at the (Fx0,Fy0,Fz0)
!              location.  It must be equal or greater than one. If Fcount
!              is greater than one, a cluster distribution of floats
!              centered at (Fx0,Fy0,Fz0) is activated.
!
!              NOTE:  The total number of floats trajectories to compute
!              ====   must add to NFLOATS.
!
!  Ft0       Time (days) of float release after model initialization (real).
!
!  Fx0       Initial float x-location (real; grid units or longitude).
!
!  Fy0       Initial float y-location (real; grid units or latitude).
!
!  Fz0       Initial float z-location (real; vertical level or depth).
!              If Fz0 is less than or equal to zero, Fz0 is the initial
!              depth in meters.  Otherwise, if Fz0 is positive and
!
!                             0 < Fz0 =< N(ng),
!
!              Fz0 is the initial position relative to the W grid
!              (0 = bottom; N = surface).
!
!              If the float trajectory type is Isobaric (Ftype=2) or
!              Geopotential (Ftype=3), Fz0 must be a negative number
!              representing depth in meters. If the float type is
!              Lagrangian (Ftype=1), Fz0 can be a level (positive) or
!              a depth in meters (negative).
!
!              WARNING: If the depth in meters at particular horizontal
!                       is not bounded, the floats are released at the
!                       surface.
!
!  Fdt       Float cluster release time interval (real; days), only used
!              if Fcount > 1:
!
!              Fdt = 0,  Fcount floats will be deployed simultaneously
!                        with a distribution centered at (Fx0,Fy0,Fz0)
!                        and defined by (Fdx,Fdy,Fdz).
!
!              Fdt > 0,  a cluster of floats will be deployed from
!                        (Fx0,Fy0,Fz0) at Fdt intervals until Fcount
!                        floats are released.
!
!  Fdx       Cluster x-distribution parameter (real), only used if
!              Fcount > 1 and Fdt = 0.
!
!  Fdy       Cluster y-distribution parameter (real), only used if
!              Fcount > 1 and Fdt = 0.
!
!  Fdz       Cluster z-distribution parameter (real), only used if
!              Fcount > 1 and Fdt = 0.
!
!            NOTE:  the parameters (Fdx,Fdy,Fdz) can be used to specify
!            ====   any type of cluster distribution, for instance:
!
!                   * Lines of floats:
!
!                         Fdx > 0,  Fdy = 0,  Fdz = 0
!                         Fdx = 0,  Fdy > 0,  Fdz = 0
!                         Fdx = 0,  Fdy = 0,  Fdz > 0
!
!             The USER can use any of these parameters to design any
!             cluster distribution in routine "init_floats.F".

