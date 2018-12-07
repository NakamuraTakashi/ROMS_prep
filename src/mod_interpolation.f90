
!!!=== Copyright (c) 2014-2018 Takashi NAKAMURA  =====

!!!**** Interpolation MODULE ************************************

  MODULE mod_interpolation

  CONTAINS

!!! **********************************************************************
!!!  Linear interpolation on the 2D grid
!!! **********************************************************************
    !
    SUBROUTINE weight2D_grid(NXd, NYd, Xd, Yd, nxr, nyr, xr, yr, Id_cont, w_cont)

      implicit none
! input parameters
      integer, intent( in) :: NXd, NYd
      real(8), intent( in) :: Xd(NXd), Yd(NYd)
      integer, intent( in) :: nxr, nyr
      real(8), intent( in) :: xr(nxr,nyr), yr(nxr,nyr)
! output parameters
      integer, intent(out) :: Id_cont(4, nxr*nyr)  ! ir,jr, Id1,Jd1
      real(8), intent(out) :: w_cont(4, nxr*nyr)    ! w1, w2, w3, w4
     
      real(8) :: f,g
      real(8) :: w(4)
      integer :: ir,jr,idc
      integer :: Idn,Jdn
      integer :: Id(4),Jd(4)

!
!          Xd(Id1)          Xd(Id1+1)
! Yd(Jd1+1) |----------------| Yd(Jd1+1) 
!           |                |
!           |                |
!           |   f            |
! yr(ir,jr) |------*         |
!           |      |g        |
!           |      |         |
!   Yd(Jd1) |----------------| Yd(Jd1) 
!         Xd(Id1) xr(ir,jr) Xd(Id1+1)
!
!

!$omp parallel
!$omp do private(ir,jr,Idn,Jdn,idc,Id,Jd,w,f,g)

      do ir=1,nxr
        do jr=1,nyr
          ! Check the cornar points.
          do Idn=2,Nxd
            if( (Xd(Idn-1)-xr(ir,jr))*(Xd(Idn)-xr(ir,jr)) <= 0.0d0 ) then
              Id(1) = Idn-1
              Id(2) = Idn
              Id(3) = Idn
              Id(4) = Idn-1
              exit
            end if
          enddo
          do Jdn=2,Nyd
            if( (Yd(Jdn-1)-yr(ir,jr))*(Yd(Jdn)-yr(ir,jr)) <= 0.0d0 ) then
              Jd(1) = Jdn-1
              Jd(2) = Jdn-1
              Jd(3) = Jdn
              Jd(4) = Jdn
              exit
            end if
          enddo
          f = abs(xr(ir,jr)-Xd(Id(1)))/abs(Xd(Id(2))-Xd(Id(1)))
          g = abs(yr(ir,jr)-Yd(Jd(1)))/abs(Yd(Jd(4))-Yd(Jd(1)))
          
          w(1) = (1.0d0-f)*(1.0d0-g)
          w(2) = f*(1.0d0-g)
          w(3) = f*g
          w(4) = (1.0d0-f)*g
          
          idc = nyr*(ir-1)+jr
          
          Id_cont(1, idc)=ir
          Id_cont(2, idc)=jr
          Id_cont(3, idc)=Id(1)
          Id_cont(4, idc)=Jd(1)
          
          w_cont(1, idc)=w(1)
          w_cont(2, idc)=w(2)
          w_cont(3, idc)=w(3)
          w_cont(4, idc)=w(4)
          
        enddo
      enddo
!$omp end do
!$omp end parallel

    END SUBROUTINE weight2D_grid
    
! **********************************************************************
    !
    SUBROUTINE weight2D_grid2(NXd, NYd, Xd, Yd, nxr, nyr, xr, yr, Id_cont, w_cont)

      implicit none
! input parameters
      integer, intent( in) :: NXd, NYd
      real(8), intent( in) :: Xd(NXd,NYd), Yd(NXd,NYd)
      integer, intent( in) :: nxr, nyr
      real(8), intent( in) :: xr(nxr,nyr), yr(nxr,nyr)
! output parameters
      integer, intent(out) :: Id_cont(8, nxr*nyr)  ! ir,jr, Id1,Jd1, Id2,Jd2, Id3,Jd3
      real(8), intent(out) :: w_cont(3, nxr*nyr)   ! w1, w2, w3
     
      real(8) :: w(3)
      integer :: ir,jr,idc
      integer :: Idn,Jdn
      integer :: Id(3),Jd(3)
      real(8) :: dmin(3)
      real(8) :: d
      real(8) :: a(3)
      integer :: is,js,ie,je,istep

!
!             Xd(Id4,Jd4)      Xd(Id3,Jd3)
!   Yd(Jd4,Jd4)|----------------|Yd(Id3,Jd3) 
!              |                |
!              |                |
!              |                |
!     yr(ir,jr)|------*         |
!              |      |         |
!              |      |         |
!   Yd(Jd1,Jd1)|----------------|Yd(Id2,Jd2) 
!                    xr(ir,jr) 
!             Xd(Id1,Jd1)       Xd(Id2,Jd2)
!


      do ir=1,nxr
        do jr=1,nyr
        
          istep = max(NXd/2,NYd/2)
          js = 1
          is = 1
          je = NYd
          ie = NXd

          do
            dmin = sqrt((Xd(NXd,NYd)-Xd(1,1))**2.0d0 + (Yd(NXd,NYd)-Yd(1,1))**2.0d0 )

            do Idn=is,ie, istep
              do Jdn=js,je, istep
                ! Check nearest points.
                d = sqrt((Xd(Idn,Jdn)-xr(ir,jr))**2.0d0 + (Yd(Idn,Jdn)-yr(ir,jr))**2.0d0 )
                if( d < dmin(1) ) then
                  dmin(3) = dmin(2)
                  dmin(2) = dmin(1)
                  dmin(1) = d
                  Id(3) = Id(2)
                  Id(2) = Id(1)
                  Id(1) = Idn
                  Jd(3) = Jd(2)
                  Jd(2) = Jd(1)
                  Jd(1) = Jdn
                else if( d < dmin(2) ) then
                  dmin(3) = dmin(2)
                  dmin(2) = d
                  Id(3) = Id(2)
                  Id(2) = Idn
                  Jd(3) = Jd(2)
                  Jd(2) = Jdn
                else if( d < dmin(3) ) then
                  dmin(3) = d
                  Id(3) = Idn
                  Jd(3) = Jdn
                end if
              enddo
            enddo
!            write(*,*) istep,Id,Jd,dmin
            if(istep == 1) then
              exit 
            else if (istep < 5 ) then
              istep = 1
              js = max(1,Jd(1)-3)
              is = max(1,Id(1)-3)
              je = min(NYd,Jd(1)+3)
              ie = min(NXd,Id(1)+3)
            else
              istep = istep/2
              js = max(1,Jd(1)-istep-2)
              is = max(1,Id(1)-istep-2)
              je = min(NYd,Jd(1)+istep+2)
              ie = min(NXd,Id(1)+istep+2)
            end if
            
          enddo
          
          a(1) = 0.5d0*abs((Xd(Id(2),Jd(2))-xr(ir,jr))*(Yd(Id(3),Jd(3))-yr(ir,jr)) &
     &                    -(Xd(Id(3),Jd(3))-xr(ir,jr))*(Yd(Id(2),Jd(2))-yr(ir,jr)))
          a(2) = 0.5d0*abs((Xd(Id(1),Jd(1))-xr(ir,jr))*(Yd(Id(3),Jd(3))-yr(ir,jr)) &
     &                    -(Xd(Id(3),Jd(3))-xr(ir,jr))*(Yd(Id(1),Jd(1))-yr(ir,jr)))
          a(3) = 0.5d0*abs((Xd(Id(2),Jd(2))-xr(ir,jr))*(Yd(Id(1),Jd(1))-yr(ir,jr)) &
     &                    -(Xd(Id(1),Jd(1))-xr(ir,jr))*(Yd(Id(2),Jd(2))-yr(ir,jr)))
          
          w(1) = a(1)/(a(1)+a(2)+a(3))
          w(2) = a(2)/(a(1)+a(2)+a(3))
          w(3) = a(3)/(a(1)+a(2)+a(3))
          
          idc = nyr*(ir-1)+jr
          
          Id_cont(1, idc)=ir
          Id_cont(2, idc)=jr
          Id_cont(3, idc)=Id(1)
          Id_cont(4, idc)=Jd(1)
          Id_cont(5, idc)=Id(2)
          Id_cont(6, idc)=Jd(2)
          Id_cont(7, idc)=Id(3)
          Id_cont(8, idc)=Jd(3)
          
          w_cont(1, idc)=w(1)
          w_cont(2, idc)=w(2)
          w_cont(3, idc)=w(3)
          
        enddo
      enddo

    END SUBROUTINE weight2D_grid2
    
! **********************************************************************

    SUBROUTINE weight3D_grid(NXd, NYd, NZd, Xd, Yd, Zd,                         &
    &            nxr, nyr, nzr, xr, yr, zr, Id_cont, w_cont)

      implicit none
! input parameters
      integer, intent( in) :: NXd, NYd, NZd
      real(8), intent( in) :: Xd(NXd), Yd(NYd), Zd(NZd)
      integer, intent( in) :: nxr, nyr, nzr
      real(8), intent( in) :: xr(nxr,nyr), yr(nxr,nyr)
      real(8), intent( in) :: zr(nxr,nyr,nzr)
! output parameters
      integer, intent(out) :: Id_cont(6, nxr*nyr*nzr)  ! ir,jr,kr, Id1,Jd1,Kd1
      real(8), intent(out) :: w_cont(8, nxr*nyr*nzr)    ! w1, w2, w3, w4, w5, w6, w7, w8

      real(8) :: f,g,h
      real(8) :: w(8)
      integer :: i,j,k
      integer :: Id,Jd,Kd
      integer :: idc

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

!$omp parallel
!$omp do private(i,j,k,idc,Id,Jd,Kd,w,f,g,h)

      do i=1,nxr
        do j=1,nyr
          ! Check the cornar points.
          do Id=2,NXd
            if( (Xd(Id-1)-xr(i,j))*(Xd(Id)-xr(i,j)) <= 0.0d0 ) then
              f = abs(xr(i,j)-Xd(Id-1))/abs(Xd(Id)-Xd(Id-1))
              exit
            end if
          enddo
          do Jd=2,NYd
            if( (Yd(Jd-1)-yr(i,j))*(Yd(Jd)-yr(i,j)) <= 0.0d0 ) then
              g = abs(yr(i,j)-Yd(Jd-1))/abs(Yd(Jd)-Yd(Jd-1))
              exit
            end if
          enddo
          
          do k=1,nzr
            do Kd=2,NZd
              h = 0.0d0
              if( (Zd(Kd-1)-zr(i,j,k))*(Zd(Kd)-zr(i,j,k)) <= 0.0d0 ) then
                h = abs(zr(i,j,k)-Zd(Kd-1))/abs(Zd(Kd)-Zd(Kd-1))
                exit
              end if
            enddo
          
            w(1) = (1.0d0-f)*(1.0d0-g)*(1.0d0-h)
            w(2) = f*(1.0d0-g)*(1.0d0-h)
            w(3) = f*g*(1.0d0-h)
            w(4) = (1.0d0-f)*g*(1.0d0-h)
            w(5) = (1.0d0-f)*(1.0d0-g)*h
            w(6) = f*(1.0d0-g)*h
            w(7) = f*g*h
            w(8) = (1.0d0-f)*g*h
          
            idc = nyr*nzr*(i-1)+nzr*(j-1)+k
          
            Id_cont(1, idc)=i
            Id_cont(2, idc)=j
            Id_cont(3, idc)=k
            Id_cont(3, idc)=Id-1
            Id_cont(4, idc)=Jd-1
            Id_cont(4, idc)=Kd-1
          
            w_cont(1, idc)=w(1)
            w_cont(2, idc)=w(2)
            w_cont(3, idc)=w(3)
            w_cont(4, idc)=w(4)
            w_cont(5, idc)=w(5)
            w_cont(6, idc)=w(6)
            w_cont(7, idc)=w(7)
            w_cont(8, idc)=w(8)
          enddo
        enddo
      enddo
!$omp end do
!$omp end parallel

    END SUBROUTINE weight3D_grid

! **********************************************************************
!
    SUBROUTINE interp2D_grid(NXd, NYd, Vd, nxr, nyr, vr, Id_cont, w_cont)

      implicit none
!     input parameters
      integer, intent( in) :: NXd, NYd
      real(8), intent( in) :: Vd(NXd,NYd)
      integer, intent( in) :: nxr, nyr
!     output parameters
      real(8), intent(out) :: vr(nxr,nyr)
      integer, intent(out) :: Id_cont(4, nxr*nyr)  ! ir,jr, Id1,Jd1
      real(8), intent(out) :: w_cont(4, nxr*nyr)    ! w1, w2, w3, w4
 
      real(8) :: w(4)
      integer :: ir,jr,idc
      integer :: Id(4),Jd(4)
      integer :: i,j
!
!$omp parallel
!$omp do private(i,j,ir,jr,idc,Id,Jd,w)
      
      do i=1,nxr
        do j=1,nyr
          ! Check the cornar points.
            !
            idc = nyr*(i-1)+j
            
            ir    = Id_cont(1, idc)
            jr    = Id_cont(2, idc)
            Id(1) = Id_cont(3, idc)
            Jd(1) = Id_cont(4, idc)
            Id(2) = Id(1)+1
            Jd(2) = Jd(1)
            Id(3) = Id(1)+1
            Jd(3) = Jd(1)+1
            Id(4) = Id(1)
            Jd(4) = Jd(1)+1
            
            w(1)= w_cont(1, idc)
            w(2)= w_cont(2, idc)
            w(3)= w_cont(3, idc)
            w(4)= w_cont(4, idc)
            
            vr(ir,jr) =  Vd(Id(1),Jd(1))*w(1) + Vd(Id(2),Jd(2))*w(2)     &
     &                 + Vd(Id(3),Jd(3))*w(3) + Vd(Id(4),Jd(4))*w(4)
          
        enddo
      enddo
!$omp end do
!$omp end parallel

    END SUBROUTINE interp2D_grid

! **********************************************************************
!
    SUBROUTINE interp2D_grid2(NXd, NYd, Vd, nxr, nyr, vr, Id_cont, w_cont)

      implicit none
!     input parameters
      integer, intent( in) :: NXd, NYd
      real(8), intent( in) :: Vd(NXd,NYd)
      integer, intent( in) :: nxr, nyr
!     output parameters
      real(8), intent(out) :: vr(nxr,nyr)
      integer, intent(out) :: Id_cont(8, nxr*nyr)  ! ir,jr, Id1,Jd1, Id2,Jd2, Id3,Jd3
      real(8), intent(out) :: w_cont(3, nxr*nyr)    ! w1, w2, w3
 
      real(8) :: w(3)
      integer :: ir,jr,idc
      integer :: Id(3),Jd(3)
      integer :: i,j
!
!$omp parallel
!$omp do private(i,j,ir,jr,idc,Id,Jd,w)
      
      do i=1,nxr
        do j=1,nyr
          ! Check the cornar points.
            !
            idc = nyr*(i-1)+j
            
            ir    = Id_cont(1, idc)
            jr    = Id_cont(2, idc)
            Id(1) = Id_cont(3, idc)
            Jd(1) = Id_cont(4, idc)
            Id(2) = Id_cont(5, idc)
            Jd(2) = Id_cont(6, idc)
            Id(3) = Id_cont(7, idc)
            Jd(3) = Id_cont(8, idc)
            
            w(1)= w_cont(1, idc)
            w(2)= w_cont(2, idc)
            w(3)= w_cont(3, idc)
            
            vr(ir,jr) =  Vd(Id(1),Jd(1))*w(1) + Vd(Id(2),Jd(2))*w(2)     &
     &                 + Vd(Id(3),Jd(3))*w(3)
          
        enddo
      enddo
!$omp end do
!$omp end parallel

    END SUBROUTINE interp2D_grid2

! **********************************************************************
! **********************************************************************
!
    SUBROUTINE interp3D_grid(NXd, NYd, NZd, Vd, nxr, nyr, nzr, vr, Id_cont, w_cont)

      implicit none
!     input parameters
      integer, intent( in) :: NXd, NYd, NZd
      real(8), intent( in) :: Vd(NXd,NYd,NZd)
      integer, intent( in) :: nxr, nyr, nzr
!     output parameters
      real(8), intent(out) :: vr(nxr,nyr,nzr)
      integer, intent(out) :: Id_cont(6, nxr*nyr*nzr)  ! ir,jr, Id1,Jd1
      real(8), intent(out) :: w_cont(8, nxr*nyr*nzr)    ! w1, w2, w3, w4
 
      real(8) :: w(8)
      integer :: idc
      integer :: Id(8),Jd(8),Kd(8)
      integer :: i,j,k
      integer :: ir,jr,kr
!
!$omp parallel
!$omp do private(i,j,k,ir,jr,kr,idc,Id,Jd,Kd,w)
      
      do i=1,nxr
        do j=1,nyr
          do k=1,nzr
          ! Check the cornar points.
            !
            idc = nyr*nzr*(i-1)+nzr*(j-1)+k
            
            ir    = Id_cont(1, idc)
            jr    = Id_cont(2, idc)
            kr    = Id_cont(3, idc)
            Id(1) = Id_cont(4, idc)
            Jd(1) = Id_cont(5, idc)
            Kd(1) = Id_cont(6, idc)
            Id(2) = Id(1)+1
            Jd(2) = Jd(1)
            Kd(2) = Kd(1)
            Id(3) = Id(1)+1
            Jd(3) = Jd(1)+1
            Kd(3) = Kd(1)
            Id(4) = Id(1)
            Jd(4) = Jd(1)+1
            Kd(4) = Kd(1)
            Id(5) = Id(1)
            Jd(5) = Jd(1)
            Kd(5) = Kd(1)+1
            Id(6) = Id(1)+1
            Jd(6) = Jd(1)
            Kd(6) = Kd(1)+1
            Id(7) = Id(1)+1
            Jd(7) = Jd(1)+1
            Kd(7) = Kd(1)+1
            Id(8) = Id(1)
            Jd(8) = Jd(1)+1
            Kd(8) = Kd(1)+1
            
            w(1)= w_cont(1, idc)
            w(2)= w_cont(2, idc)
            w(3)= w_cont(3, idc)
            w(4)= w_cont(4, idc)
            w(5)= w_cont(5, idc)
            w(6)= w_cont(6, idc)
            w(7)= w_cont(7, idc)
            w(8)= w_cont(8, idc)
            
            vr(ir,jr,kr) =  Vd(Id(1),Jd(1),Kd(1))*w(1) + Vd(Id(2),Jd(2),Kd(2))*w(2) &
     &                    + Vd(Id(3),Jd(3),Kd(3))*w(3) + Vd(Id(4),Jd(4),Kd(4))*w(4) &
     &                    + Vd(Id(5),Jd(5),Kd(5))*w(5) + Vd(Id(6),Jd(6),Kd(6))*w(6) &
     &                    + Vd(Id(7),Jd(7),Kd(7))*w(7) + Vd(Id(8),Jd(8),Kd(8))*w(8)
     
          enddo
        enddo
      enddo
!$omp end do
!$omp end parallel

    END SUBROUTINE interp3D_grid

! **********************************************************************
! **********************************************************************



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
!
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
    
! **********************************************************************

    SUBROUTINE LinearInterpolation2D_grid3(Nx, Ny, X, Y, V, vmin, vmax, Im, Jm, Xgrd, Ygrd, Vgrd)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny
      real(8), intent( in) :: X(Nx,Ny), Y(Nx,Ny)
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
          call LinearInterpolation2D_point2(Nx, Ny, X, Y, Vnew, Xgrd(i,j), Ygrd(i,j), Vgrd(i,j))
        enddo
      enddo
!$omp end do
!$omp end parallel
      
      RETURN
    END SUBROUTINE LinearInterpolation2D_grid3

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

    SUBROUTINE LinearInterpolation2D_point2(Nx, Ny, X, Y, V, xi, yi, vi)

      implicit none
! input parameters
      integer, intent( in) :: Nx, Ny
      real(8), intent( in) :: X(Nx,Ny), Y(Nx,Ny)
      real(8), intent( in) :: V(Nx,Ny)
      real(8), intent( in) :: xi, yi
! output parameters
      real(8), intent(out) :: vi
     
      real(8), parameter :: dmin  = 1.0d-6
      integer :: istep1
      integer :: istep2
      integer :: istep
      integer :: i,j
      integer :: is, ie, js, je

      integer :: in(4), jn(4)
      integer :: ic, jc
      real(8) :: d(4)
      real(8) :: dtmp
      real(8) :: sum_w
      
      
      js = 1
      is = 1
      je = Ny
      ie = Nx
      istep = max(Nx/2,Ny/2)
      do
        d = sqrt( (X(1,1)-X(Nx,Ny))**2.0d0+(Y(1,1)-Y(Nx,Ny))**2.0d0 )
        do j=js, je, istep
          do i=is, ie, istep
            dtmp = sqrt( (X(i,j)-xi)**2.0d0+(Y(i,j)-yi)**2.0d0 )
            if( dtmp < d(1)) then
              d(4) = d(3)
              d(3) = d(2)
              d(2) = d(1)
              d(1) = dtmp
              in(4) = in(3)
              in(3) = in(2)
              in(2) = in(1)
              in(1) = i
              jn(4) = jn(3)
              jn(3) = jn(2)
              jn(2) = jn(1)
              jn(1) = j
            else if( dtmp < d(2)) then
              d(4) = d(3)
              d(3) = d(2)
              d(2) = dtmp
              in(4) = in(3)
              in(3) = in(2)
              in(2) = i
              jn(4) = jn(3)
              jn(3) = jn(2)
              jn(2) = j
            else if( dtmp < d(3)) then
              d(4) = d(3)
              d(3) = dtmp
              in(4) = in(3)
              in(3) = i
              jn(4) = jn(3)
              jn(3) = j
            else if( dtmp < d(4)) then
              d(4) = dtmp
              in(4) = i
              jn(4) = j
            end if
          enddo
          
        enddo
        
        if(istep == 1) then
          exit 
        else if (istep < 5 ) then
          istep = 1
          js = max(1,jn(1)-3)
          is = max(1,in(1)-3)
          je = min(Ny,jn(1)+3)
          ie = min(Nx,in(1)+3)
        else
          istep = istep/2
          js = max(1,jn(1)-istep-2)
          is = max(1,in(1)-istep-2)
          je = min(Ny,jn(1)+istep+2)
          ie = min(Nx,in(1)+istep+2)
        end if
        

      enddo
      
!      sum_w = 1.0d0/d(1)+1.0d0/d(2)+1.0d0/d(3)+1.0d0/d(4)
!      
!      vi =( V(in(1),jn(1))*1.0d0/d(1)    &
!     &     +V(in(2),jn(2))*1.0d0/d(2)    &
!     &     +V(in(3),jn(3))*1.0d0/d(3)    &
!     &     +V(in(4),jn(4))*1.0d0/d(4) )/sum_w

      sum_w = 1.0d0/d(1)/d(1)+1.0d0/d(2)/d(2)+1.0d0/d(3)/d(3)+1.0d0/d(4)/d(4)
      
      vi =( V(in(1),jn(1))*1.0d0/d(1)/d(1)    &
     &     +V(in(2),jn(2))*1.0d0/d(2)/d(2)    &
     &     +V(in(3),jn(3))*1.0d0/d(3)/d(3)    &
     &     +V(in(4),jn(4))*1.0d0/d(4)/d(4) )/sum_w
     
!      sum_w = d(1)*d(1)+d(2)*d(2)+d(3)*d(3)+d(4)*d(4)
!      
!      vi =( V(in(1),jn(1))*(sum_w-d(1)*d(1))    &
!     &     +V(in(2),jn(2))*(sum_w-d(2)*d(2))    &
!     &     +V(in(3),jn(3))*(sum_w-d(3)*d(3))    &
!     &     +V(in(4),jn(4))*(sum_w-d(4)*d(4)))/sum_w

!      vi = V(in(1),jn(1))
      
    END SUBROUTINE LinearInterpolation2D_point2
    
! **********************************************************************
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


