!----------------------------------------------------------------------------
! EdgeDetectPFI: Fortran90 implementation of the algorithm proposed in the 
!                manuscript "EdgeDetectPFI: an algorithm for automatic edge
!                detection in potential field anomaly images - application to
!                dike-like magnetic structures" (SP Oliveira, FJF Ferreira, 
!                and J de Souza), published in Computers & Geosciences
!                http://dx.doi.org/

program EdgeDetectPFI

  integer :: n,ST_Mz,ST_minus,ierr,i
  double precision auxv(4),h,a,x,y,xv,xminus
  allocatable x(:),y(:),xv(:),xminus(:),ST_Mz(:),ST_minus(:)

  ! n       : total number of grid points
  ! ST_Mz   : signum transform of vertical derivative 
  ! STminus : signum transform of vertical derivative minus
  !           total horiz derivative
  ! ierr,i  : auxiliary integers
  ! auxv    : auxiliary I/O data vector
  ! h       : estimated depth
  ! a       : estimated half-width
  ! (x,y)   : grid coordinates
  ! xv      : radii within ST_Mz
  ! xminus  : radii within ST_minus
  
  !--------------------------------------------------------------------------
  ! Read input data

  ! detect dimensions and allocate variables

  open(unit=11,file='input.csv');
  ierr = 1
  n = 0
  do while(ierr.ne.-1)
     read(11,*,iostat=ierr) auxv
     n = n + 1
  end do
  n = n-1

  allocate( x(n),y(n),ST_Mz(n),ST_minus(n),xv(n),xminus(n) )

  ! read coordinates and signum-transformed data

  rewind(unit=11)
  do i = 1,n
     read(11,*) x(i),y(i),ST_Mz(i),ST_minus(i)
  end do
  close(unit=11)

  !--------------------------------------------------------------------------
  ! Compute radii and estimate depth & width (saving on CSV file)

  call radius(n,x,y,st_mz,xv)
  call radius(n,x,y,st_minus,xminus)

  open(unit=12,file='output.csv');
  do i = 1,n
     if ((xminus(i).gt.0.0D0)) then
        h = (xv(i)*xv(i) - xminus(i)*xminus(i))/(2.0D0*xminus(i))
        a = xv(i)*xv(i) - h*h
        if (a.ge.0.0D0) then
           a = dsqrt(a)
           write(12,99) x(i),y(i),h,2.0D0*a
        end if
     end if
  end do
  close(unit=12)

99 format(F12.2,',',F12.2,',',E16.7,',',E16.7)

  stop
end program EdgeDetectPFI

!----------------------------------------------------------------------------
subroutine radius(n,x,y,S,rad)

  ! radius: for each point (xi,yi) with S(xi,yi)=1, find the radius r of the 
  ! largest circle C with center in (xi,yi) such that S(x,y)=1 if (x,y) in C.
  !
  ! To find the circle, we find the point (xf,xf) such that S(x,y)=-1 that is
  ! closer to (xi,yi). Let dC be the distance from (xi,yi) to (xf,yf)
  !
  ! The radius must be less than dC. We select r = dC - delta, delta is half
  ! the grid length
  !
  ! Once r is found, all the points (xj,yj) within the circle C are assigned
  ! the radius rad_j = r, unless the current value rad_j > r.

  ! Input arguments: 
  ! n       : data dimensions
  ! x,y     : coordinates
  ! S       : signum-transformed data
  !
  ! Output arguments:
  ! rad     : radius of the largest circle C

  integer :: n,S(n)
  double precision :: x(n),y(n),rad(n)

  ! Internal variables:
  ! i,ii    : counters
  ! nS      : number of grid points (x,y) such that S(x,y) = 1
  ! ind     : vector with indices i such that S(xi,yi) = 1
  ! dist    : vector of distances from a point (xi,yi) to all grid points
  ! delta   : one half of the grid length
  ! penalty : auxliary vector that penalizes points with S(x,y)=1, so that
  !           they are not taken into account when searching (xf,yf)

  integer :: i,ii,nS,ind
  double precision :: dist,r,delta,penalty  
  allocatable ind(:),dist(:),penalty(:)
  allocate  ( ind(n),dist(n),penalty(n) )

  !--------------------------------------------------------------------------

  rad(1:n) = 0.0D0*x(1:n)
  penalty(1:n) = rad(1:n)
  delta = 0.5D0*max(x(2)-x(1),y(2)-y(1))

  nS = 0
  do i = 1,n
     if(S(i).eq.1) then
        nS = nS + 1
        ind(nS) = i
        penalty(i) = 1.0D10
     end if
  end do

  do i = 1,nS
     dist(1:n) = dsqrt((x(1:n)-x(ind(i)))**2.0D0 + (y(1:n)-y(ind(i)))**2.0D0)
     r = minval(dist + penalty) - delta
     if(r.lt.0.0D0) r = 0.0D0
     do j = 1,nS
        ii = ind(j)
        if((rad(ii).lt.r).and.(dist(ii).le.r)) rad(ii) = r
     end do
  end do

  return
end subroutine radius
