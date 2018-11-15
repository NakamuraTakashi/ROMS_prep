
!!!=== Copyright (c) 2018 Takashi NAKAMURA  =====

!!!**** UTILITY MODULE ****************

      MODULE mod_utility

      implicit none
      real(8), parameter :: PI = 3.14159265359d0
      
      CONTAINS

! --- Compute coriolis parameter ---------------------
 
      real(8) function Coriolis(lat)
      real(8), intent(in) :: lat      ! Latitude (degree)
      real(8) :: omega
      
      omega = 2.0d0*PI*366.25d0/(24.0d0*3600.0d0*365.25d0)
      Coriolis = 2.0d0*omega * sin( PI/180.0d0*lat )

      return
      end function Coriolis
!
! --- Compute distance ---------------------
 
      real(8) function distance(lat1, lon1, lat2, lon2)
      real(8), intent(in) :: lat1      ! Latitude (degree) at point1
      real(8), intent(in) :: lon1      ! Longitude (degree) at point1
      real(8), intent(in) :: lat2      ! Latitude (degree) at point2
      real(8), intent(in) :: lon2      ! Longitude (degree) at point2

      distance = 0.0d0

      return
      end function distance
!
! --- Compute mask_u, mask_v, mask_psi ---------------------
 
      SUBROUTINE uvp_masks(N_xi_rho, N_eta_rho, mask_rho, mask_u, mask_v, mask_psi)
      integer, intent(in ) :: N_xi_rho
      integer, intent(in ) :: N_eta_rho
      real(8), intent(in ) :: mask_rho(N_xi_rho,N_eta_rho) ! Land/Sea mask on PHO-points.
      real(8), intent(out) :: mask_u(N_xi_rho-1,N_eta_rho) ! Land/Sea mask on U-points.
      real(8), intent(out) :: mask_v(N_xi_rho,N_eta_rho-1) ! Land/Sea mask on V-points.
      real(8), intent(out) :: mask_psi(N_xi_rho-1,N_eta_rho-1) ! Land/Sea mask on PSI-points.
      integer :: i,j
      
      do i=1, N_xi_rho-1
        do j=1, N_eta_rho
          mask_u(i,j) = mask_rho(i+1,j)*mask_rho(i,j)
        enddo
      enddo
      
      do i=1, N_xi_rho
        do j=1, N_eta_rho-1
          mask_v(i,j) = mask_rho(i,j+1)*mask_rho(i,j)
        enddo
      enddo

      return
      END SUBROUTINE uvp_masks

      END MODULE mod_utility

