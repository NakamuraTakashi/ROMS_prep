!****************************************************************************
!* System: Connectivity Modeling System (CMS)                               *
!* File : calendar.f90                                                      *
!* Last Modified: 2011-07-22                                                *
!* Code contributors: Judith Helgers, Ashwanth Srinivasan, Claire B. Paris  * 
!*                                                                          *
!* Copyright (C) 2011, University of Miami                                  *
!*                                                                          *
!* This program is free software: you can redistribute it and/or modify     *
!* it under the terms of the GNU Lesser General Public License as published *
!* by the Free Software Foundation, either version 3 of the License, or     *
!*(at your option) any later version.                                       *
!*                                                                          *
!* This program is distributed in the hope that it will be useful,          *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of           *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                     *
!* See the Lesser GNU General Public License for more details.              *
!*                                                                          *
!* You should have received a copy of the GNU Lesser General                *
!* Public License along with this program.                                  *
!* If not, see <http://www.gnu.org/licenses/>.                              *
!****************************************************************************

MODULE mod_calendar

IMPLICIT NONE

CONTAINS

!**************************************************************

!calculates how many days are in the given month.
SUBROUTINE DaysMonth(yyyy, mm, days)
 integer, intent(IN)  :: yyyy
 integer, intent(IN)  :: mm
 integer, intent(OUT) :: days

 IF (mm==0) days = 31
 IF (mm==1) days = 31
 IF (mm==2) then
    IF(mod(yyyy,100).ne.0.and.mod(yyyy,4).eq.0)then
      days = 29
    ELSEIF(mod(yyyy,400).eq.0) THEN
      days = 29
    ELSE
      days= 28
    ENDIF
 ENDIF
 IF (mm==3) days = 31
 IF (mm==4) days = 30
 IF (mm==5) days = 31
 IF (mm==6) days = 30
 IF (mm==7) days = 31
 IF (mm==8) days = 31
 IF (mm==9) days = 30
 IF (mm==10) days = 31
 IF (mm==11) days = 30
 IF (mm==12) days = 31

END SUBROUTINE DaysMonth

!**************************************************************

!converts calendar date to julian date number.
SUBROUTINE jd(yyyy, mm, dd, julian)

 integer, intent(IN)  :: yyyy
 integer, intent(IN)  :: mm
 integer, intent(IN)  :: dd
 integer, intent(OUT) :: julian
 integer:: A,B,C

!The algorithm is valid for all Gregorian calendar dates after 4800 BC
 A = (14-mm)/12
 B = yyyy + 4800 - A
 C = mm + (12*A) - 3
 julian = dd + ((153*C+2)/5) + 365*B + (B/4) - (B/100) + (B/400) - 32045
      
END SUBROUTINE jd

!**************************************************************

!converts calendar date to julian date number.
SUBROUTINE dble_jd(iyear, imonth, iday, ihour, imin, dble_julian)

 integer, intent(IN)  :: iyear
 integer, intent(IN)  :: imonth
 integer, intent(IN)  :: iday
 integer, intent(IN)  :: ihour
 integer, intent(IN)  :: imin
 real(8), intent(OUT) :: dble_julian
 integer:: julian

 Call jd(iyear, imonth, iday, julian)
 dble_julian = dble(julian) + (dble(ihour)/24.0d0) + (dble(imin)/1440.0d0)
      
END SUBROUTINE dble_jd

!**************************************************************

!converts ex. from 20240709 to 2024, 7, 9.
SUBROUTINE yyyymmdd2year_month_day(yyyymmdd, iyear, imonth, iday)

 integer, intent(IN )  :: yyyymmdd
 integer, intent(OUT)  :: iyear
 integer, intent(OUT)  :: imonth
 integer, intent(OUT)  :: iday

 iyear  = yyyymmdd/10000
 imonth = (yyyymmdd-iyear*10000)/100
 iday   = yyyymmdd-iyear*10000-imonth*100
      
END SUBROUTINE yyyymmdd2year_month_day

!**************************************************************

!converts calendar date to julian date number.
SUBROUTINE dble_jd2(yyyyMMdd, hhmm, dble_julian)

 integer, intent(IN)  :: yyyyMMdd
 integer, intent(IN)  :: hhmm
 real(8), intent(OUT) :: dble_julian
 integer :: iyear, imonth, iday, ihour, imin

 Call yyyymmdd2year_month_day(yyyyMMdd, iyear, imonth, iday)
 ihour  = hhmm/100
 imin   = hhmm-100*ihour
 Call dble_jd(iyear, imonth, iday, ihour, imin, dble_julian)
      
END SUBROUTINE dble_jd2

!**************************************************************

!converts julian date number to calendar date
SUBROUTINE cdate(jd, yyyy, mm, dd)
 integer, intent(IN)   :: jd
 integer, intent(OUT)  :: yyyy
 integer, intent(OUT)  :: mm
 integer, intent(OUT)  :: dd
 integer               :: j,g,dg,c,dc,b,db,a,da

 j = jd + 32044
 g = int(j/146097); dg = mod(j, 146097)
 c = int((int(dg/36524) + 1) * 3 / 4); dc = dg - c * 36524
 b = int(dc/1461); db = mod(dc, 1461)
 a = int((int(db/365) + 1) * 3 / 4); da = db - a * 365
 yyyy = g * 400 + c * 100 + b * 4 + a
 mm = int((da * 5 + 308)/ 153) - 2
 dd = int(da - (mm + 4) * 153 / 5) + 122
 yyyy = int(yyyy - 4800 + (mm + 2) / 12)
 mm = mod((mm + 2),12) + 1
 dd = dd + 1

END SUBROUTINE cdate

!**************************************************************

!converts julian date number to calendar date
SUBROUTINE dble_jd2cdate(dble_jd, iyear, imonth, iday, ihour, imin)
 real(8), intent(IN)   :: dble_jd
 integer, intent(OUT)  :: iyear
 integer, intent(OUT)  :: imonth
 integer, intent(OUT)  :: iday
 integer, intent(OUT)  :: ihour
 integer, intent(OUT)  :: imin

 integer :: jd

 jd = nint(dble_jd)  ! integer part of julian date number
 Call cdate(jd, iyear, imonth, iday)
 ihour = nint( (dble_jd - dble(jd))*24.0d0 )
 imin  = nint( (dble_jd - dble(jd) - dble(ihour)/24.0d0)*1440.0d0 )

END SUBROUTINE dble_jd2cdate

!**************************************************************

!converts julian date number to calendar date
SUBROUTINE dble_jd2yyyyMMMdd_hhmm(dble_jd, yyyyMMdd, hhmm)
 real(8), intent(IN)   :: dble_jd
 integer, intent(OUT)  :: yyyyMMdd
 integer, intent(OUT)  :: hhmm

 integer :: iyear, imonth, iday, ihour, imin

 Call dble_jd2cdate(dble_jd, iyear, imonth, iday, ihour, imin)
 yyyymmdd = iyear*10000 + imonth*100 + iday
 hhmm = ihour*100 + imin

END SUBROUTINE dble_jd2yyyyMMMdd_hhmm

!**************************************************************

!calculates number of days between two dates.
!days is positive if date1 is more recent than date2
SUBROUTINE ndays(mm1, dd1, yyyy1, mm2, dd2, yyyy2, days)
 
 integer, intent(IN)  :: mm1
 integer, intent(IN)  :: dd1
 integer, intent(IN)  :: yyyy1
 integer, intent(IN)  :: mm2
 integer, intent(IN)  :: dd2
 integer, intent(IN)  :: yyyy2
 integer, intent(OUT) :: days
 integer :: jd1, jd2

 CALL jd(yyyy1, mm1, dd1, jd1)
 CALL jd(yyyy2, mm2, dd2, jd2)
 days = jd1 - jd2

END SUBROUTINE ndays

!**************************************************************

!calculates number of months between two dates.
!nmonths is positive if date1 is more recent than date2
SUBROUTINE nmonths(mm1,yyyy1,mm2,yyyy2, months)
 integer, intent(IN)  :: mm1
 integer, intent(IN)  :: yyyy1
 integer, intent(IN)  :: mm2
 integer, intent(IN)  :: yyyy2
 integer, intent(OUT) :: months

 months = (yyyy1 - yyyy2)*12 + (mm1 - mm2)  

END SUBROUTINE nmonths

!**************************************************************

!converts ROMS ocean_time (sec) since (Ryear, Rmonth, Rday) to calendar date
SUBROUTINE oceantime2cdate(ocean_time, Ryear, Rmonth, Rday, YYYYMMDDpHH)
  real(8), intent(IN)  :: ocean_time  ! sec
  integer, intent(IN)  :: Ryear
  integer, intent(IN)  :: Rmonth
  integer, intent(IN)  :: Rday
  character(len=*), intent(OUT)  :: YYYYMMDDpHH
  integer :: yyyy, mm, dd, hh
  integer :: rjd, out_jd
  real(8) :: ocean_date

  ocean_date = ocean_time/24.0d0/60.0d0/60.0d0
  CALL jd(Ryear, Rmonth, Rday, rjd)
  out_jd = floor(ocean_date)+rjd
  CALL cdate(out_jd, yyyy, mm, dd)
  hh = int( ( ocean_date - dble( floor(ocean_date) ) )*24.0d0 )

  write (YYYYMMDDpHH, "(I4.4, I2.2, I2.2, '.', I2.2)") yyyy, mm, dd, hh

  write(*,*) "DEBUG :", yyyy, mm, dd, hh, YYYYMMDDpHH
   
 END SUBROUTINE oceantime2cdate
 
END MODULE mod_calendar

