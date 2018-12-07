!***********************************************************************
      SUBROUTINE set_depth (        &
     &              Im, Jm, N       &
     &            , h               &
# ifdef ICESHELF
     &            , zice            &
# endif
     &            , Vtransform      &
     &            , hc              &
     &            , sc_w, sc_r      &
     &            , Cs_w, Cs_r      &
!        output parameters
     &            , z_r, z_w)
!***********************************************************************
!
!  Imported variable declarations.
!
!    input parameters
      integer, intent( in) :: Im,Jm,N
      real(8), intent( in) :: h(Im,Jm)
#  ifdef ICESHELF
      real(8), intent( in) :: zice(Im,Jm)
#  endif
      integer, intent( in) :: Vtransform 
      real(8), intent( in) :: hc       
      real(8), intent( in) :: sc_w(0:N)       
      real(8), intent( in) :: sc_r(1:N)  
      real(8), intent( in) :: Cs_w(0:N)       
      real(8), intent( in) :: Cs_r(1:N)  
      real(8), intent(out) :: z_r(Im,Jm,N)
      real(8), intent(out) :: z_w(Im,Jm,0:N)

!  Local variable declarations.
!
      integer :: i, j, k

      real(8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(8) :: hinv, hwater, z_r0, z_w0

      real(8), parameter :: zeta0 = 0.0d0

!
!-----------------------------------------------------------------------
!  Original formulation: Compute time independent vertical depths
!                        (meters, negative) at RHO- and W-points.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = Zo_w + zeta(x,y,0) * [1.0 + Zo_w / h(x,y)]
!
!                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!  where zeta(x,y,0) = 0 for time independent depths.
!
!-----------------------------------------------------------------------
!
      IF (Vtransform.eq.1) THEN
        DO j=1,Jm
          DO i=1,Im
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N
            cff_r=hc*(sc_r(k)-Cs_r(k))
            cff_w=hc*(sc_w(k)-Cs_w(k))
            cff1_r=Cs_r(k)
            cff1_w=Cs_w(k)
            DO i=1,Im
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0d0/hwater
              z_w0=cff_w+cff1_w*hwater
              z_w(i,j,k)=z_w0+zeta0*(1.0d0+z_w0*hinv)
              z_r0=cff_r+cff1_r*hwater
              z_r(i,j,k)=z_r0+zeta0*(1.0d0+z_r0*hinv)
# ifdef ICESHELF
              z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
              z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  New formulation: Compute time independent vertical depths
!                   (meters, negative) at RHO- and W-points.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = zeta(x,y,0) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!  where zeta(x,y,0) = 0 for time independent depths.
!
!-----------------------------------------------------------------------
!
      ELSE IF (Vtransform.eq.2) THEN
        DO j=1,Jm
          DO i=1,Im
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N
            cff_r=hc*sc_r(k)
            cff_w=hc*sc_w(k)
            cff1_r=Cs_r(k)
            cff1_w=Cs_w(k)
            DO i=1,Im
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0d0/(hc+hwater)
              cff2_r=(cff_r+cff1_r*hwater)*hinv
              cff2_w=(cff_w+cff1_w*hwater)*hinv

              z_w(i,j,k)=zeta0+(zeta0+hwater)*cff2_w
              z_r(i,j,k)=zeta0+(zeta0+hwater)*cff2_r
# ifdef ICESHELF
              z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
              z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
            END DO
          END DO
        END DO
      END IF
      END SUBROUTINE set_depth
