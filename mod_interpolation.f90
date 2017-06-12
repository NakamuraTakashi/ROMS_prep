
!!!=== Copyright (c) 2014-2017 Takashi NAKAMURA  =====

!!!**** Interpolation MODULE ************************************

  MODULE mod_interpolation

  CONTAINS

!!! **********************************************************************
!!!  Linear interpolation on the 2D grid
!!! **********************************************************************

    SUBROUTINE LinearInterpolation2D_grid(Nx, Ny, X, Y, V, Im, Jm, Xgrd, Ygrd, Vgrd)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny
      real(8), intent( in) :: X(Nx), Y(Ny)
      real(8), intent( in) :: V(Nx,Ny)
      integer, intent( in) :: Im, Jm
      real(8), intent( in) :: Xgrd(Im,Jm), Ygrd(Im,Jm)
! output parameters
      real(8), intent(out) :: Vgrd(Im,Jm)
     
      integer :: i,j
      
      do i=1,Im
        do j=1,Jm
          call LinearInterpolation2D_point(Nx, Ny, X, Y, V, Xgrd(i,j), Ygrd(i,j), Vgrd(i,j))
        enddo
      enddo
      
      RETURN
    END SUBROUTINE LinearInterpolation2D_grid
    
! **********************************************************************

    SUBROUTINE LinearInterpolation2D_grid2(Nx, Ny, X, Y, V, vmin, vmax, Im, Jm, Xgrd, Ygrd, Vgrd)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny
      real(8), intent( in) :: X(Nx), Y(Ny)
      real(8), intent( in) :: V(Nx,Ny)
      real(8), intent( in) :: vmin, vmax
      integer, intent( in) :: Im, Jm
      real(8), intent( in) :: Xgrd(Im,Jm), Ygrd(Im,Jm)
! output parameters
      real(8), intent(out) :: Vgrd(Im,Jm)
      
      real(8) :: Vnew(Nx,Ny)
      real(8) :: mask(Nx,Ny)
      real(8) :: masksum
      integer :: i,j,k,loop
      integer :: loopflag
      
      Vnew(:,:)=V(:,:)
      
!$omp parallel
!$omp do private(i,j)
      do i=1,Nx
        do j=1,Ny
          if( (vmin > Vnew(i,j)) .or. (Vnew(i,j) > vmax) ) then
            mask(i,j)=0.0d0
          else
            mask(i,j)=1.0d0
          endif
        enddo
      enddo
!$omp end do
!$omp end parallel

      do
        loopflag = 0

        do i=2,Nx-1
          do j=2,Ny-1
            if( mask(i,j)==0.0d0 ) then
              masksum = mask(i-1,j)+mask(i+1,j)+mask(i,j-1)+mask(i,j+1)
              if(masksum == 0.0d0) then
                loopflag = 1
              else
                Vnew(i,j) =( mask(i-1,j)*Vnew(i-1,j)                   &
     &                      +mask(i+1,j)*Vnew(i+1,j)                   &
     &                      +mask(i,j-1)*Vnew(i,j-1)                   &
     &                      +mask(i,j+1)*Vnew(i,j+1) )                 &
     &                     /masksum
                mask(i,j) = 1.0d0
              endif
            endif
          enddo
        enddo
        if(loopflag == 0) exit
      enddo

!$omp parallel
!$omp do private(i)
      do i=1,Nx
        if( mask(i,1)==0.0d0 ) then
          Vnew(i,1) = Vnew(i,2)
        endif
      enddo
!$omp end do
!$omp do private(i)
      do i=1,Nx
        if( mask(i,Ny)==0.0d0 ) then
          Vnew(i,Ny) = Vnew(i,Ny-1)
        endif
      enddo
!$omp end do
!$omp do private(j)
      do j=1,Ny
        if( mask(1,j)==0.0d0 ) then
          Vnew(1,j) = Vnew(2,j)
        endif
      enddo
!$omp end do
!$omp do private(j)
      do j=1,Ny
        if( mask(Nx,j)==0.0d0 ) then
          Vnew(Nx,j) = Vnew(Nx-1,j)
        endif
      enddo
!$omp end do
!$omp do private(i,j)
      do i=1,Im
        do j=1,Jm
          call LinearInterpolation2D_point(Nx, Ny, X, Y, Vnew, Xgrd(i,j), Ygrd(i,j), Vgrd(i,j))
        enddo
      enddo
!$omp end do
!$omp end parallel
      
      RETURN
    END SUBROUTINE LinearInterpolation2D_grid2
    
!!! **********************************************************************
!!!  Linear interpolation on the 3D grid
!!! **********************************************************************
    SUBROUTINE LinearInterpolation3D_grid(Nx, Ny, Nz, X, Y, Z, V, Im, Jm, Km, Xgrd, Ygrd, Zgrd, Vgrd)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny, Nz
      real(8), intent( in) :: X(Nx), Y(Ny), Z(Nz)
      real(8), intent( in) :: V(Nx,Ny,Nz)
      integer, intent( in) :: Im, Jm, Km
      real(8), intent( in) :: Xgrd(Im,Jm), Ygrd(Im,Jm)
      real(8), intent( in) :: Zgrd(Im,Jm,Km)
! output parameters
      real(8), intent(out) :: Vgrd(Im,Jm,Km)
     
      integer :: i,j,k
      
!$omp parallel
!$omp do private(i,j,k)
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            call LinearInterpolation3D_point(Nx, Ny, Nz, X, Y, Z, V   &
   &                , Xgrd(i,j), Ygrd(i,j), Zgrd(i,j,k), Vgrd(i,j,k))
          enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
      
      RETURN
    END SUBROUTINE LinearInterpolation3D_grid
    
! **********************************************************************

    SUBROUTINE LinearInterpolation3D_grid2(Nx, Ny, Nz, X, Y, Z, V, vmin, vmax, Im, Jm, Km, Xgrd, Ygrd, Zgrd, Vgrd)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny, Nz
      real(8), intent( in) :: X(Nx), Y(Ny), Z(Nz)
      real(8), intent( in) :: V(Nx,Ny,Nz)
      real(8), intent( in) :: vmin, vmax
      integer, intent( in) :: Im, Jm, Km
      real(8), intent( in) :: Xgrd(Im,Jm), Ygrd(Im,Jm)
      real(8), intent( in) :: Zgrd(Im,Jm,Km)
! output parameters
      real(8), intent(out) :: Vgrd(Im,Jm,Km)
     
      real(8) :: Vnew(Nx,Ny,Nz)
      real(8) :: mask(Nx,Ny,Nz)
      real(8) :: masksum
      integer :: i,j,k
      integer :: loopflag
      
      Vnew(:,:,:)=V(:,:,:)
      
      do i=1,Nx
        do j=1,Ny
          do k=1,Nz
            if( (vmin > Vnew(i,j,k)) .or. (Vnew(i,j,k) > vmax) ) then
              mask(i,j,k)=0.0d0
            else
              mask(i,j,k)=1.0d0
            endif
          enddo
        enddo
      enddo
      
      
      do
        loopflag = 0

        do i=2,Nx-1
          do j=2,Ny-1
            do k=2,Nz-1
              if( mask(i,j,k)==0.0d0 ) then
                masksum = mask(i-1,j,k)+mask(i+1,j,k)  &
     &                   +mask(i,j-1,k)+mask(i,j+1,k)  &
     &                   +mask(i,j,k-1)+mask(i,j,k+1)
     
                if(masksum <= 1.0d0) then
                  loopflag = 1
                else
                  Vnew(i,j,k) =( mask(i-1,j,k)*Vnew(i-1,j,k)       &
     &                          +mask(i+1,j,k)*Vnew(i+1,j,k)       &
     &                          +mask(i,j-1,k)*Vnew(i,j-1,k)       &
     &                          +mask(i,j+1,k)*Vnew(i,j+1,k)       &
     &                          +mask(i,j,k-1)*Vnew(i,j,k-1)       &
     &                          +mask(i,j,k+1)*Vnew(i,j,k+1) )     &
     &                        /masksum
                  mask(i,j,k) = 1.0d0
                endif
              endif
            enddo
          enddo
        enddo
        if(loopflag == 0) exit
      enddo
      
      do i=1,Nx
        do j=1,Ny
          if( mask(i,j,1)==0.0d0 ) then
            Vnew(i,j,1) = Vnew(i,j,2)
          endif
        enddo
      enddo
      do i=1,Nx
        do j=1,Ny
          if( mask(i,j,Nz)==0.0d0 ) then
            Vnew(i,j,Nz) = Vnew(i,j,Nz-1)
          endif
        enddo
      enddo
      do i=1,Nx
        do k=1,Nz
          if( mask(i,1,k)==0.0d0 ) then
            Vnew(i,1,k) = Vnew(i,2,k)
          endif
        enddo
      enddo
      do i=1,Nx
        do k=1,Nz
          if( mask(i,Ny,k)==0.0d0 ) then
            Vnew(i,Ny,k) = Vnew(i,Ny-1,k)
          endif
        enddo
      enddo
      do j=1,Ny
        do k=1,Nz
          if( mask(1,j,k)==0.0d0 ) then
            Vnew(1,j,k) = Vnew(2,j,k)
          endif
        enddo
      enddo
      do j=1,Ny
        do k=1,Nz
          if( mask(Nx,j,k)==0.0d0 ) then
            Vnew(Nx,j,k) = Vnew(Nx-1,j,k)
          endif
        enddo
      enddo

      
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            call LinearInterpolation3D_point(Nx, Ny, Nz, X, Y, Z, Vnew &
   &                , Xgrd(i,j), Ygrd(i,j), Zgrd(i,j,k), Vgrd(i,j,k))
          enddo
        enddo
      enddo
      
      RETURN
    END SUBROUTINE LinearInterpolation3D_grid2
    
!!! **********************************************************************

    SUBROUTINE LinearInterpolation3D_grid3(Nx, Ny, Nz, X, Y, Z, V, vmin, vmax, Im, Jm, Km, Xgrd, Ygrd, Zgrd, Vgrd)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny, Nz
      real(8), intent( in) :: X(Nx), Y(Ny), Z(Nz)
      real(8), intent( in) :: V(Nx,Ny,Nz)
      real(8), intent( in) :: vmin, vmax
      integer, intent( in) :: Im, Jm, Km
      real(8), intent( in) :: Xgrd(Im,Jm), Ygrd(Im,Jm)
      real(8), intent( in) :: Zgrd(Im,Jm,Km)
! output parameters
      real(8), intent(out) :: Vgrd(Im,Jm,Km)
     
      real(8) :: Vnew(Nx,Ny,Nz)
      real(8) :: mask(Nx,Ny,Nz)
      real(8) :: masksum
      integer :: i,j,k,loop
      integer :: loopflag
      
      Vnew(:,:,:)=V(:,:,:)

!$omp parallel
!$omp do private(i,j,k)
      
      do i=1,Nx
        do j=1,Ny
          do k=1,Nz
            if( (vmin > Vnew(i,j,k)) .or. (Vnew(i,j,k) > vmax) ) then
              mask(i,j,k)=0.0d0
            else
              mask(i,j,k)=1.0d0
            endif
          enddo
        enddo
      enddo
!$omp end do
!$omp do private(i,j,k,loopflag,masksum)
      do k=1,Nz
!        do loop=1,500
        do
          loopflag = 0
          do i=2,Nx-1
            do j=2,Ny-1
              if( mask(i,j,k)==0.0d0 ) then
                masksum = mask(i-1,j,k)+mask(i+1,j,k)  &
     &                   +mask(i,j-1,k)+mask(i,j+1,k)  !&
!     &                   +mask(i,j,k-1)+mask(i,j,k+1)
     
                if(masksum == 0.0d0) then
                  loopflag = 1
                else
                  Vnew(i,j,k) =( mask(i-1,j,k)*Vnew(i-1,j,k)       &
     &                          +mask(i+1,j,k)*Vnew(i+1,j,k)       &
     &                          +mask(i,j-1,k)*Vnew(i,j-1,k)       &
     &                          +mask(i,j+1,k)*Vnew(i,j+1,k)  )    &
!     &                          +mask(i,j,k-1)*Vnew(i,j,k-1)       &
!     &                          +mask(i,j,k+1)*Vnew(i,j,k+1) )     &
     &                        /masksum
                  mask(i,j,k) = 1.0d0
                endif
              endif
            enddo
          enddo
          if(loopflag == 0) exit
        enddo
        if(loopflag == 1) then
!          write(*,*) 'OF', k
          do i=2,Nx-1
            do j=2,Ny-1
              if( mask(i,j,k)==0.0d0 ) then
                Vnew(i,j,k)=Vnew(i,j,k-1)
              endif
            enddo
          enddo
        endif
      enddo
!$omp end do
!$omp do private(i,k)
      do i=1,Nx
        do k=1,Nz
          if( mask(i,1,k)==0.0d0 ) then
            Vnew(i,1,k) = Vnew(i,2,k)
          endif
        enddo
      enddo
!$omp end do
!$omp do private(i,k)
      do i=1,Nx
        do k=1,Nz
          if( mask(i,Ny,k)==0.0d0 ) then
            Vnew(i,Ny,k) = Vnew(i,Ny-1,k)
          endif
        enddo
      enddo
!$omp end do
!$omp do private(j,k)
      do j=1,Ny
        do k=1,Nz
          if( mask(1,j,k)==0.0d0 ) then
            Vnew(1,j,k) = Vnew(2,j,k)
          endif
        enddo
      enddo
!$omp end do
!$omp do private(j,k)
      do j=1,Ny
        do k=1,Nz
          if( mask(Nx,j,k)==0.0d0 ) then
            Vnew(Nx,j,k) = Vnew(Nx-1,j,k)
          endif
        enddo
      enddo
!$omp end do
!$omp do private(i,j,k)
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            call LinearInterpolation3D_point(Nx,Ny,Nz,X, Y, Z, Vnew  &
   &                , Xgrd(i,j), Ygrd(i,j), Zgrd(i,j,k), Vgrd(i,j,k))
          enddo
        enddo
      enddo
!$omp end do
!$omp end parallel
      
      RETURN
    END SUBROUTINE LinearInterpolation3D_grid3
    
!!! **********************************************************************
!!! **********************************************************************


    SUBROUTINE LinearInterpolation2D_point(Nx, Ny, X, Y, V, xi, yi, vi)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny
      real(8), intent( in) :: X(Nx), Y(Ny)
      real(8), intent( in) :: V(Nx,Ny)
      real(8), intent( in) :: xi, yi
! output parameters
      real(8), intent(out) :: vi
     
      integer :: i,j,k
      integer :: iL, iR, jB, jT
      real(8) :: xL, xR, yB, yT
      real(8) :: f,g
      
!       vLT              vRT
!      YT |-------------|
!         |             |
!         |             |
!         | f           |
!         |---* (xi,yi) |
!         |   |g        |
!         |   |         |
!      YB |-------------|
!     vLB XL           XR vRB

! Check the cornar points.
      do i=2,Nx
        if( (X(i-1)-xi)*(X(i)-xi) <= 0.0d0) then
          iL = i-1
          iR = i
          xL = X(iL)
          xR = X(iR)
          exit
        end if
      enddo
      do j=2,Ny
        if( (Y(j-1)-yi)*(Y(j)-yi) <= 0.0d0) then
          jB = j-1
          jT = j
          yB = Y(jB)
          yT = Y(jT)
          exit
        end if
      enddo
      f = abs(xi-xL)/abs(xR-xL)
      g = abs(yi-yB)/abs(yT-yB)
      
      vi = V(iL,jB)*(1.0d0-f)*(1.0d0-g)                         &
     &    +V(iR,jB)* f *(1.0d0-g)                               &
     &    +V(iL,jT)*(1.0d0-f)* g                                &
     &    +V(iR,jT)* f * g
      RETURN
    END SUBROUTINE LinearInterpolation2D_point
    
! **********************************************************************


    SUBROUTINE LinearInterpolation3D_point(Nx, Ny, Nz, X, Y, Z, V, xi, yi, zi, vi)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny, Nz
      real(8), intent( in) :: X(Nx), Y(Ny), Z(Nz)
      real(8), intent( in) :: V(Nx,Ny,Nz)
      real(8), intent( in) :: xi, yi, zi
! output parameters
      real(8), intent(out) :: vi
     
      integer :: i,j,k
      integer :: iL, iR, jB, jT, kF, kB
      real(8) :: xL, xR, yB, yT, zF, zB
      real(8) :: f,g,h
      real(8) :: vLBF, vRBF, vLTF, vRTF, vLBB, vRBB, vLTB, vRTB
      
!          vLTB______________ vRTB
!             /             /| 
!            / |           / |
!           /  |          /  |
!      vLTF/____________ vRTF|
!      YT |    |        |    |
!         |  __f_ *(xi,yi)   |
!    vLBB | /|   /| ___ | __ | vRBB
!         |/____/h|g    |   /zB
!         |  |  | |     |  /
!         | /   | /     | /
!      YB |_____|_______|/ zF
!    vLBF XL           XR vRBF

! Check the cornar points.
      do i=2,Nx
        if( (X(i-1)-xi)*(X(i)-xi) <= 0.0d0) then
          iL = i-1
          iR = i
          xL = X(iL)
          xR = X(iR)
          exit
        end if
      enddo
      
      do j=2,Ny
        if( (Y(j-1)-yi)*(Y(j)-yi) <= 0.0d0) then
          jB = j-1
          jT = j
          yB = Y(jB)
          yT = Y(jT)
          exit
        end if
      enddo
      
      do k=2,Nz
        if( (Z(k-1)-zi)*(Z(k)-zi) <= 0.0d0) then
          kF = k-1
          kB = k
          zF = Z(kF)
          zB = Z(kB)
          exit
        end if
        if( k==Nz ) then
          if(ABS(Z(1)-zi) <= ABS(Z(Nz)-zi) ) then
            kF = 1
            kB = 1
            zF = zi
            zB = Z(kB)
          else
            kF = Nz
            kB = Nz
            zF = Z(kF)
            zB = zi
          endif
        end if
      enddo
      
      f = abs(xi-xL)/abs(xR-xL)
      g = abs(yi-yB)/abs(yT-yB)
      h = abs(zi-zF)/abs(zB-zF)
      
      vLBF=V(iL,jB,kF)
      vRBF=V(iR,jB,kF)
      vLTF=V(iL,jT,kF)
      vRTF=V(iR,jT,kF)
      vLBB=V(iL,jB,kB)
      vRBB=V(iR,jB,kB)
      vLTB=V(iL,jT,kB)
      vRTB=V(iR,jT,kB)
      
      vi =  vLBF *(1.0d0-f)*(1.0d0-g)*(1.0d0-h)   &
     &    + vRBF * f *(1.0d0-g)*(1.0d0-h)         &
     &    + vLTF *(1.0d0-f)* g *(1.0d0-h)         &
     &    + vRTF * f * g * (1.0d0-h)              &
     &    + vLBB * (1.0d0-f)*(1.0d0-g)*h          &
     &    + vRBB * f *(1.0d0-g)*h                 &
     &    + vLTB *(1.0d0-f)* g *h                 &
     &    + vRTB * f * g * h

      RETURN
    END SUBROUTINE LinearInterpolation3D_point
    
! **********************************************************************

    SUBROUTINE FillExtraPoints2D(Nx, Ny, V, vmin, vmax)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny
      real(8), intent(inout) :: V(Nx,Ny)
      real(8), intent( in) :: vmin, vmax
      
      real(8) :: mask(Nx,Ny)
      real(8) :: masksum
      integer :: i,j,k
      integer :: loopflag
      
      do i=1,Nx
        do j=1,Ny
          if( (vmin > V(i,j)) .or. (V(i,j) > vmax) ) then
            mask(i,j)=0.0d0
          else
            mask(i,j)=1.0d0
          endif
        enddo
      enddo
      
      
      do
        loopflag = 0

        do i=2,Nx-1
          do j=2,Ny-1
            if( mask(i,j)==0.0d0 ) then
              masksum = mask(i-1,j)+mask(i+1,j)+mask(i,j-1)+mask(i,j+1)
              if(masksum <= 1.0d0) then
                loopflag = 1
              else
                V(i,j) =( mask(i-1,j)*V(i-1,j)                   &
     &                      +mask(i+1,j)*V(i+1,j)                   &
     &                      +mask(i,j-1)*V(i,j-1)                   &
     &                      +mask(i,j+1)*V(i,j+1) )                 &
     &                     /masksum
                mask(i,j) = 1.0d0
              endif
            endif
          enddo
        enddo
        if(loopflag == 0) exit
      enddo
      
      do i=1,Nx
        if( mask(i,1)==0.0d0 ) then
          V(i,1) = V(i,2)
        endif
      enddo
      do i=1,Nx
        if( mask(i,Ny)==0.0d0 ) then
          V(i,Ny) = V(i,Ny-1)
        endif
      enddo
      do j=1,Ny
        if( mask(1,j)==0.0d0 ) then
          V(1,j) = V(2,j)
        endif
      enddo
      do j=1,Ny
        if( mask(Nx,j)==0.0d0 ) then
          V(Nx,j) = V(Nx-1,j)
        endif
      enddo
      
      RETURN
    END SUBROUTINE FillExtraPoints2D

    
  END MODULE mod_interpolation


