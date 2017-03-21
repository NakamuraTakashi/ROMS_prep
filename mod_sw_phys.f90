
!!!=== ver 2015/11/11   Copyright (c) 2015 Takashi NAKAMURA  =====

!!!**** MODULE OF SEAWATER PHYSICS ****************

      MODULE mod_sw_phys

      implicit none

! --- convert  temperature to potential temperature (oC) ---------------------
 
      real(8) function sw_temp2potemp(temp,S,press,ref_press)
      real(8), intent(in) :: temp      ! Temperature (oC)
      real(8), intent(in) :: S         ! Salinity (psu)
      real(8), intent(in) :: press     ! Pressuer (Pa)
      real(8), intent(in) :: ref_press ! Reference pressuer (Pa)
      
      real(8) :: DP,P,Q
      real(8) :: R1,R2,R3,R4,R5
      real(8) :: S1,T,X
      integer :: i,j,N

      S1 = S-35.0d0
      P  = press
      T  = temp

      DP = ref_press - P
      N  = IFIX(ABS(DP)/1000.) + 1
      DP = DP/dble(N)

      DO 10 I=1,N
         DO 20 J=1,4

            R1 = ((-2.1687d-16*T+1.8676d-14)*T-4.6206d-13)*P
            R2 = (2.7759d-12*T-1.1351d-10)*S1
            R3 = ((-5.4481d-14*T+8.733d-12)*T-6.7795d-10)*T
            R4 = (R1+(R2+R3+1.8741d-8))*P+(-4.2393d-8*T+1.8932d-6)*S1
            R5 = R4+((6.6228d-10*T-6.836d-8)*T+8.5258Ed-6)*T+3.5803d-5

            X  = DP*R5

            GO TO (100,200,300,400),J

  100       CONTINUE
            T = T+.5*X
            Q = X
            P = P + .5*DP
            GO TO 20

  200       CONTINUE
            T = T + 0.29298322d0*(X-Q)
            Q = 0.58578644d0*X + 0.121320344d0*Q
            GO TO 20

  300       CONTINUE
            T = T + 1.707106781d0*(X-Q)
            Q = 3.414213562d0*X - 4.121320344d0*Q
            P = P + 0.5d0*DP
            GO TO 20

  400       CONTINUE
            T = T + (X-2.0d0*Q)/6.0d0
  20      CONTINUE
  10    CONTINUE

      sw_temp2potemp = T
      return
      end function sw_temp2potemp
    
! --- convert depth to pressure (Pa)---------------------------------   

      real(8) function sw_depth2press(depth)
      real(8), intent(in) :: depth   ! Depth (m)



      sw_temp2potemp = depth
      return
      end function sw_depth2press



      END MODULE mod_swphys

