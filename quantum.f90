module potential
  implicit none

  integer, parameter :: dp = selected_real_kind(16,300)

contains
  function potmod(r)
    real(kind=dp), intent(in) :: r
    real(kind=dp) :: potmod

    potmod = 0.15_dp*r**4 - 4.5_dp*r**2 - 0.77529_dp*r    
  end function potmod
end module potential

program quantum
  use potential
  use solver
  implicit none

  integer :: ierr, i, j, wfn_unit, en_unit
  integer :: numsteps, numlevels
  real(kind=dp) :: bound, dx, left, right

  real(kind=dp), dimension(:), allocatable :: eigenenergy, dens1, dens12
  real(kind=dp), dimension(:,:), allocatable :: wavefunc, wavefunc2, wavefunc12

  numlevels = 2
  bound = 6.0_dp
  numsteps = 10000
  dx = 2.0_dp * bound / real(numsteps,dp)

  allocate(eigenenergy(numlevels), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating eigenenergy in main.'
  allocate(wavefunc(numsteps+1,5), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating wavefunc in main.'
  allocate(wavefunc2(numsteps+1,5), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating wavefunc2 in main.'
  allocate(wavefunc12(numsteps+1,numsteps+1), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating wavefunc12 in main.'
  allocate(dens1(numsteps+1), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating dens1 in main.'
  allocate(dens12(numsteps+1), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating dens12 in main.'

  wfn_unit = 11
  open(wfn_unit, file='wfn.dat', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening file wfn.dat in main'

  en_unit = 12
  open(en_unit, file='energy.dat', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening file energy.dat in main'

  eigenenergy = 0.0_dp
  wavefunc = 0.0_dp

  do i = 1, numsteps+1
     wavefunc(i,1) = -bound + dx * real(i-1,dp)
     wavefunc(i,5) = potmod(wavefunc(i,1))
  end do
  wavefunc(1,2) = 0.0_dp
  wavefunc(1,3) = 1.0e-5_dp
  wavefunc(1,4) = tiny(1.0_dp)
  wavefunc(2,2) = wavefunc(1,3) * dx
  
  wavefunc2(:,:) = wavefunc(:,:)
  
  call bisection(wavefunc,eigenenergy,0.00025_dp,1.0e-14_dp)
  call solve_numerov(wavefunc,eigenenergy(1))
  call normalise(wavefunc,dx)
  call solve_numerov(wavefunc2,eigenenergy(2))
  call normalise(wavefunc2,dx)

  do j = 1, size(wavefunc2,1)
     do i = 1, size(wavefunc,1)
        wavefunc12(i,j) = wavefunc(i,2)*wavefunc2(j,2) - wavefunc2(i,2)*wavefunc(j,2)
     end do
  end do
  wavefunc12(:,:) = wavefunc12(:,:) / sqrt(2.0_dp)

  dens1(:) = wavefunc(:,2)**2
  do i = 1, size(wavefunc12,1)
     dens12(i) = dot_product(wavefunc12(i,:),wavefunc12(i,:))
  end do
  dens12(:) = 2.0_dp * dens12(:) * dx

  ! calculate left to right ratio
  left = 0.0_dp
  do i = 1, (numsteps+1)/2
     left = left + dens12(i)
  end do
  left = left / sum(dens12)
  right = 1.0_dp - left
  
  do i = 1, size(eigenenergy,1)
     if (eigenenergy(i).gt.maxval(wavefunc(:,5))) then
        write(en_unit,*) eigenenergy(i), '<---- Scatt. State!'
        print*, eigenenergy(i), '<---- Scatt. State!'
     else
        write(en_unit,*) eigenenergy(i)
        print*, eigenenergy(i)
     end if
  end do
  write(*, 16) 'Left-to-right ratio =', left, '-', right
  16 format (A,1X,F7.5,1X,A,1X,F7.5)
       
  do i = 1, numsteps+1
     write(wfn_unit,*) wavefunc(i,1), wavefunc(i,2), wavefunc2(i,2), dens12(i), wavefunc(i,5)
  end do
  
  deallocate(eigenenergy, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating eigenenergy in main.'
  deallocate(wavefunc, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating wavefunc in main.'
  deallocate(wavefunc2, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating wavefunc2 in main.'
  deallocate(wavefunc12, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating wavefunc12 in main.'
  deallocate(dens1, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating dens1 in main.'
  deallocate(dens12, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating dens12 in main.'
  close(wfn_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing file wfn.dat in main.'
  close(en_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing file energy.dat in main.'

end program quantum

module solver
  use potential
  implicit none

contains
  function schrodinger(wfn,v,e)
    implicit none

    real(kind=dp), intent(in) :: wfn, v, e
    real(kind=dp) :: schrodinger
    
    schrodinger = 2.0_dp * (v-e) * wfn
  end function schrodinger
    
  subroutine solve_verlet(array,e)
    implicit none
    
    real(kind=dp), dimension(:,:), allocatable :: array
    real(kind=dp) :: dx, e
    integer :: i
    
    if(.not.allocated(array)) stop 'Array not allocated in solve_verlet.'
    
    dx = array(2,1) - array(1,1)
    do i = 2, size(array,1)
       array(i,2) = array(i-1,2) + array(i-1,3) * dx + 0.5_dp * array(i-1,4) * dx**2
       array(i,4) = schrodinger(array(i,2), potmod(array(i,1)), e)
       array(i,3) = array(i-1,3) + 0.5_dp * (array(i-1,4)+array(i,4)) * dx
    end do
  end subroutine solve_verlet

  subroutine solve_numerov(array,e)
    implicit none
    real(kind=dp), parameter :: twelveth = 1.0_dp / 12.0_dp

    real(kind=dp), dimension(:,:), allocatable :: array
    real(kind=dp) :: dx, e, phi_j, phi_k, phi
    integer :: i
    
    if(.not.allocated(array)) stop 'Array not allocated in solve_verlet.'
    
    dx = array(2,1) - array(1,1)
    do i = 3, size(array,1)
       phi_j = array(i-2,2) * (1.0_dp-dx*twelveth*2.0_dp*(array(i-2,5)-e))
       phi_k = array(i-1,2) * (1.0_dp-dx*twelveth*2.0_dp*(array(i-1,5)-e))
       phi = 2.0_dp * phi_k - phi_j + dx**2 * (array(i-1,5)-e) * array(i-1,2)
       array(i,2) = phi / (1.0_dp-dx*twelveth*2.0_dp*(array(i,5)-e))
    end do
  end subroutine solve_numerov

  ! function dls_numerov(array,e) result(dls)
  !   implicit none
  !   real(kind=dp), parameter :: twelveth = 1.0_dp / 12.0_dp

  !   real(kind=dp), dimension(:,:), allocatable :: array
  !   real(kind=dp) :: dx, e, phi_j, phi_k, phi, dls
  !   integer :: i
    
  !   if(.not.allocated(array)) stop 'Array not allocated in solve_verlet.'
    
  !   dx = array(2,1) - array(1,1)
  !   array(size(array,1),2) = 0.0_dp
  !   array(size(array,1),3) = 1.0e-5_dp
  !   array(size(array,1),4) = tiny(1.0_dp)
  !   array(size(array,1)-1,2) = array(size(array,1),3) * dx
  !   do i = 3, size(array,1)/2
  !      phi_j = array(i-2,2) * (1.0_dp-dx*twelveth*2.0_dp*(array(i-2,5)-e))
  !      phi_k = array(i-1,2) * (1.0_dp-dx*twelveth*2.0_dp*(array(i-1,5)-e))
  !      phi = 2.0_dp * phi_k - phi_j + dx**2 * (array(i-1,5)-e) * array(i-1,2)
  !      array(i,2) = phi / (1.0_dp-dx*twelveth*2.0_dp*(array(i,5)-e))
  !      array(i,3) = (array(i,2)-array(i-1,2))/dx
  !      array(i,4) = 2.0_dp*(potmod(array(i,5)-e))
  !   end do
  !   do i = size(array,1)-2, size(array,1)/2+1, -1
  !      phi_j = array(i+2,2) * (1.0_dp-dx*twelveth*2.0_dp*(array(i+2,5)-e))
  !      phi_k = array(i+1,2) * (1.0_dp-dx*twelveth*2.0_dp*(array(i+1,5)-e))
  !      phi = 2.0_dp * phi_k - phi_j + dx**2 * (array(i+1,5)-e) * array(i+1,2)
  !      array(i,2) = phi / (1.0_dp-dx*twelveth*2.0_dp*(array(i,5)-e))
  !      array(i,3) = (array(i+1,2)-array(i,2))/dx
  !      array(i,4) = 2.0_dp*(potmod(array(i,5)-e))
  !   end do
  !   dls = min(abs(array(size(array,1)/2,2)/array(size(array,1)/2,3)-array(size(array,1)/2+1,2)/array(size(array,1)/2+1,3)), &
  !        & abs(array(size(array,1)/2,3)/array(size(array,1)/2,3)-array(size(array,1)/2+1,3)/array(size(array,1)/2+1,2)))
  ! end function dls_numerov
  
  subroutine bisection(wfn_array,en_array,de,tol)
    implicit none

    real(kind=dp), dimension(:,:), allocatable, intent(in) :: wfn_array
    real(kind=dp), dimension(:), allocatable :: en_array
    real(kind=dp), intent(in) :: de, tol
    real(kind=dp) :: low, mid, high, wl, wm, wh
    integer :: level

    low = minval(wfn_array(:,5))
    call solve_numerov(wfn_array,low)
    wl = wfn_array(size(wfn_array,1),2)

    high = low + 2.0_dp * de
    call solve_numerov(wfn_array,high)
    wh = wfn_array(size(wfn_array,1),2)

    mid = 0.5_dp * (low + high)
    call solve_numerov(wfn_array,mid)
    wm = wfn_array(size(wfn_array,1),2)

    do level = 1, size(en_array,1)
       do while (abs(wm).gt.tol)
          if (wl*wm.lt.0) then
             wh = wm
             high = mid
             mid = 0.5_dp * (low + high)
             call solve_numerov(wfn_array,mid)
             wm = wfn_array(size(wfn_array,1),2)
          else if (wm*wh.lt.0) then
             wl = wm
             low = mid
             mid = 0.5_dp * (low + high)
             call solve_numerov(wfn_array,mid)
             wm = wfn_array(size(wfn_array,1),2)
          else
             wm = wh
             wl = wm
             low = mid
             mid = high
             high = 2.0_dp*high - low
             call solve_numerov(wfn_array,high)
             wh = wfn_array(size(wfn_array,1),2)
          end if
          if (high-mid .lt. 1.0e-16_dp) exit
       end do
       en_array(level) = mid

       low = high + de
       call solve_numerov(wfn_array,low)
       wl = wfn_array(size(wfn_array,1),2)

       high = high + de * 2.0_dp
       call solve_numerov(wfn_array,high)
       wh = wfn_array(size(wfn_array,1),2)

       mid = 0.5_dp * (low + high)
       call solve_numerov(wfn_array,mid)
       wm = wfn_array(size(wfn_array,1),2)
    end do
  end subroutine bisection

  subroutine normalise(wfn_array, dx)
    implicit none

    real(kind=dp), intent(in) :: dx
    real(kind=dp), dimension(:,:), allocatable :: wfn_array

    real(kind=dp) :: normconst
    integer :: i

    if (.not.allocated(wfn_array)) stop 'wfn_array not allocated in normalise.'

    normconst = 1.0_dp / sqrt(dx * norm2(real(wfn_array(:,2),8))**2)
    do i = 2, 4
       wfn_array(:,i) = wfn_array(:,i) * normconst
    end do
  end subroutine normalise
    
end module solver

