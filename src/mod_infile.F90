
!!!=== Copyright (c) 2025-2026 Takashi NAKAMURA  =====

!!!**** INFILE MODULE ************************************

MODULE mod_infile

  implicit none

  TYPE T_INFILE
    character(256), pointer :: NAME(:)
    real(8), pointer :: time_all(:)
    integer :: Nt
    integer :: ItS, ItE
  END TYPE T_INFILE
  TYPE (T_INFILE), allocatable :: INFILE(:)

  CONTAINS

  SUBROUTINE infile_check_time( Nfile, t_Start, t_End, Nt, time, iNCt, idt )

  integer, intent( in) :: Nfile
  real(8), intent( in) :: t_Start, t_End
  integer, intent(out) :: Nt
  real(8), allocatable, intent(out) :: time(:)
  integer, allocatable, optional, intent(out) :: iNCt(:)  ! merged step -> file index
  integer, allocatable, optional, intent(out) :: idt(:)   ! merged step -> in-file index

  real(8) :: t
  integer :: i,j,k,l
  integer :: js,je

  write(*,*) INFILE(:)%Nt

  do j=1, Nfile-1
    do i=2,INFILE(j)%Nt
      if( INFILE(j+1)%time_all(1) <= INFILE(j)%time_all(i) ) then
        INFILE(j)%Nt = i-1
        exit
      end if
    end do
  end do
  
  write(*,*) INFILE(:)%Nt

  write(*,*) "******************************************************************"

  INFILE(:)%ItE = -1
  INFILE(:)%ItS = -1
  
  do j=Nfile,1,-1
    do i=INFILE(j)%Nt,1,-1
      t = INFILE(j)%time_all(i)
      if(t < t_End) then
        write(*,*) '*** FOUND: Ending point @ INFILE',j
        INFILE(j)%ItE=i
        exit
      endif
    end do
  end do
  write(*,*) INFILE(:)%ItE 
  
  do j=1,Nfile
    do i=INFILE(j)%ItE,1,-1
      t = INFILE(j)%time_all(i)
      if(t < t_Start) then
!        write(*,*) '*** FOUND: Starting point @ ATM_FILE',j
        exit
      endif
      INFILE(j)%ItS=i
    end do
  end do
  write(*,*) INFILE(:)%ItS 

  Nt = 0
  js = Nfile
  je = 1

  do j=1,Nfile
    if(INFILE(j)%ItS==-1) then
      cycle
    end if
    Nt = Nt + INFILE(j)%ItE - INFILE(j)%ItS + 1
    js = min(js,j)
    je = max(je,j)
  enddo

  allocate( time(Nt) )
  if( present(iNCt) .and. present(idt) ) allocate( iNCt(Nt), idt(Nt) )

  ! Build the merged time list (and optional step->file / step->in-file maps).
  ! NOTE: do NOT overwrite the loop variable j here (the earlier draft did
  ! `j = i + ...`, which is undefined behaviour for a DO index).
  i=1
  do j=js,je
    if( INFILE(j)%ItS==-1 ) cycle           ! skip files outside [t_Start,t_End]
    k = INFILE(j)%ItE - INFILE(j)%ItS        ! (number of steps from this file) - 1
    time(i:i+k) = INFILE(j)%time_all( INFILE(j)%ItS : INFILE(j)%ItE )
    if( present(iNCt) .and. present(idt) ) then
      do l=0,k
        iNCt(i+l) = j
        idt(i+l)  = INFILE(j)%ItS + l
      end do
    end if
    i = i+k+1
  enddo

  END SUBROUTINE infile_check_time
    
END MODULE mod_infile
      
! -------------------------------------------------------------------------
