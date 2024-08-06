!#
!                                        *****************************
!                                        *         pbswi.f90         *
!                                        *       R. J. Purser        *
!                                        *       NOAA/NCEP/EMC       *
!                                        *       February 2022       *
!                                        *****************************
! NOAA/NCEP Environmental Modeling Center.
! jim.purser@noaa.gov
! Suite of routines to perform B-Spline-Weighted Interpllations (BSWI)
! using ! stencils of 4 pts, 6 pts, or 8 pts.
! The interpolations are performed at points within the central interval
! of ! each stencil laid on the uniform unit grid in 1D or 2D.
! The 4-point stencil scheme uses linearly-weighted quadratic interpolation,
! which is already widely used, but not generally recognized as being one
! member of the wider family of the interpolations by overlapping 
! Lagrange polynomials that are themselves weighted according to the values
! across the interpolation interval of segments of a B-spline of degree 
! one-less than that of the component Lagrnage polynomial. Thus, the 6-pt
! scheme weights three overlapping Lagrange cubics with weights from the
! three segments of the 2nd degree B-spline; the 8-pt scheme weights four
! overlapping quartic Lagrangen polynomials with weights from the four
! segments of the 3rd degree B-spline. Although this technique generalizes
! to non-uniformly spaced grids, the codes for uniform unit grids (where the
! assumed coordinates are also the indices in each dimension) are particularly
! simple since the coefficient matrices involved in constructing the 
! interpolation weights have rational components that have been precomputed.
! The interpolating polynomial of these schemes are of degree one-less than 
! the number of stencil points (for example, the 4-pt scheme is a cubic in
! the central interval, passing through the values at both sides of the 
! interval but, unlike the Lagrange cubic, does not generally fit the other
! values touched by the 4-pt stencil.) The interpolating polynomial is
! constructed explicitly by way of its power series, which seems to be the
! most efficient way to obtain the result.
! The most trivial members of the BSWI family, linear or bilinear schemes,
! with 2-pt stencils, are included for completeness (the relevant B-spline
! in this case is the single segment unit "hat" function, which weights by
! one the equally trivial linear Lagrange polynomial.)
!
! DEPENDENCIES
! Modules: pkind, pietc
!=============================================================================
module pbswi
!=============================================================================
private
public :: pown,getrp4,getrp6,getrp8,abswi2,abswi4,abswi6,abswi8

interface pown;   module procedure pown;              end interface
interface getrp4; module procedure getrp4;            end interface
interface getrp6; module procedure getrp6;            end interface
interface getrp8; module procedure getrp8;            end interface
interface abswi2; module procedure abswi2_1,abswi2_2; end interface
interface abswi4; module procedure abswi4_1,abswi4_2; end interface
interface abswi6; module procedure abswi6_1,abswi6_2; end interface
interface abswi8; module procedure abswi8_1,abswi8_2; end interface

contains

!============================================================================
subroutine pown(n,x,xps)!                                              [pown]
!============================================================================
! Return a vector, xps,  containing the powers of x up to x**n.
!============================================================================
use pkind, only: spi,dp
use pietc, only: u1
implicit none
integer(spi),           intent(in ):: n
real(dp),               intent(in ):: x
real(dp),dimension(0:n),intent(out):: xps
!----------------------------------------------------------------------------
real(dp)    :: xp
integer(spi):: i
!============================================================================
xp=u1
do i=0,n
   xps(i)=xp
   xp=xp*x
enddo
end subroutine pown

!============================================================================
subroutine getrp4(rp4)!                                              [getrp4]
!============================================================================
! Return the matrix rp4 that converts the 4 consecutive stencil values
! of a unit grid to the power-series coefficients in the central
! interval, relative to the origin at the lower side of it, of the
! B-spline-weighted interpolation polynomial.
!============================================================================
use pkind, only: spi,dp
use pietc, only: u2
implicit none
real(dp),dimension(0:3,0:3),intent(out):: rp4
!----------------------------------------------------------------------------
integer(spi),dimension(0:3,0:3):: irp4
data irp4/0,-1,2,-1, 2,0,-5,3, 0,1,4,-3, 0,0,-1,1/
!============================================================================
rp4=irp4/u2
end subroutine getrp4

!============================================================================
subroutine getrp6(rp6)!                                              [getrp6]
!============================================================================
! Return the matrix rp6 that converts the 6 consecutive stencil values
! of a unit grid to the power-series coefficients in the central
! interval, relative to the origin at the lower side of it, of the
! B-spline-weighted interpolation polynomial.
!============================================================================
use pkind, only: spi,dp
implicit none
real(dp),parameter:: u12=12
real(dp),dimension(0:5,0:5),intent(out):: rp6
!----------------------------------------------------------------------------
integer(spi),dimension(0:5,0:5):: irp6
data irp6/&
     0, 1, -2, 0,  2, -1, &
     0,-8, 14, 0,-11,  5, &
    12, 0,-24,-2, 24,-10, &
     0, 8, 14, 6,-26, 10, &
     0,-1, -2,-6, 14, -5, &
     0, 0,  0, 2, -3,  1/
!============================================================================
rp6=irp6/u12
end subroutine getrp6

!============================================================================
subroutine getrp8(rp8)!                                              [getrp8]
!============================================================================
! Return the matrix rp8 that converts the 8 consecutive stencil values
! of a unit grid to the power-series coefficients in the central
! interval, relative to the origin at the lower side of it, of the
! B-spline-weighted interpolation polynomial.
!============================================================================
use pkind, only: spi,dp
implicit none
real(dp),parameter:: u144=144
real(dp),dimension(0:7,0:7),intent(out):: rp8
!----------------------------------------------------------------------------
integer(spi),dimension(0:7,0:7):: irp8
data irp8/&
     0,  -2,   5, -1,  -6,   4,   1, -1, &
     0,  20, -36, -8,  48, -19, -12,  7, &
     0,-106, 171, 19,-138,  24,  51,-21, &
   144,   0,-280,  0, 186,  25,-110, 35, &
     0, 106, 171,-19,-114,-100, 135,-35, &
     0, -20, -36,  8,  12, 111, -96, 21, &
     0,   2,   5,  1,  18, -56,  37, -7, &
     0,   0,   0,  0,  -6,  11,  -6,  1/
!============================================================================
rp8=irp8/u144
end subroutine getrp8

!============================================================================
subroutine abswi2_1(lx,mx,ax,x, a,ff)!                               [abswi2]
!============================================================================
! Linear interpolation (trivial member of the BSWI family).
!============================================================================
use pkind, only: spi,dp
use pietc, only: u1
implicit none
integer(spi),             intent(in ):: lx,mx
real(dp),dimension(lx:mx),intent(in ):: ax
real(dp),                 intent(in ):: x
real(dp),                 intent(out):: a
logical,                  intent(out):: ff
!----------------------------------------------------------------------------
real(dp)    :: wx0,wx1
integer(spi):: ix
!============================================================================
ix=floor(x)
ff=ix<lx .or.ix>=mx
if(ff)return
wx1=x-ix
wx0=u1-wx1
a=wx0*ax(ix)+wx1*ax(ix+1)
end subroutine abswi2_1
!============================================================================
subroutine abswi2_2(lx,mx,ly,my,axy,x,y, a,ff)!                      [abswi2]
!============================================================================
! Bilinear interpolation (trivial member of the BSWI family).
!============================================================================
use pkind, only: spi,dp
use pietc, only: u1
implicit none
integer(spi),                   intent(in ):: lx,mx,ly,my
real(dp),dimension(lx:mx,ly:my),intent(in ):: axy
real(dp),                       intent(in ):: x,y
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp)    :: wx0,wx1,wy0,wy1
integer(spi):: ix,iy
!============================================================================
ix=floor(x); iy=floor(y)
ff=ix<lx .or.ix>=mx .or. iy<ly .or. iy>=my
if(ff)return
wx1=x-ix;   wy1=y-iy
wx0=u1-wx1; wy0=u1-wy1
a=wy0*(wx0*axy(ix,iy  )+wx1*axy(ix+1,iy  )) &
 +wy1*(wx0*axy(ix,iy+1)+wx1*axy(ix+1,iy+1))
end subroutine abswi2_2

!============================================================================
subroutine abswi4_1(lx,mx,rp4,ax,x, a,ff)!                           [abswi4]
!============================================================================
! Apply the B-Spline-Weighted Interpolation with 4-point stnecil to a 1D
! array of real values, ax(lx:mx), with a target point at (x) in
! the index unit coordinates. If the target is outside the region of the
! grid that allows for centered 4-point interpolation,
! the logical Failure Flag, FF, is returned .true. and no interpolation
! is attempted; otherwise, for a normal condition, the interpolated values
! is returned as "a".
! Real 4*4 array, rp4, is the constant coefficients matrix precomputed in
! subroutine getrp4
!===========================================================================
use pkind, only: spi,dp
implicit none
integer(spi),                   intent(in ):: lx,mx
real(dp),dimension(0:3,0:3),    intent(in ):: rp4
real(dp),dimension(lx:mx),      intent(in ):: ax
real(dp),                       intent(in ):: x
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp),dimension(-1:2):: xp,wx
real(dp)                :: rx
integer(spi)            :: i,ix
!============================================================================
ix=floor(x)
ff=ix<lx+1 .or. ix>=mx-1
if(ff)return
rx=x-ix
call pown(3,rx,xp)
wx=matmul(xp,rp4)
a=dot_product(wx,ax(ix-1:ix+2))
end subroutine abswi4_1
!============================================================================
subroutine abswi4_2(lx,mx,ly,my,rp4,axy,x,y, a,ff)!                  [abswi4]
!============================================================================
! Apply the B-Spline-Weighted Interpolation with 4-point stnecil to a 2D
! array of real values, axy(lx:mx,ly:my), with a target point at (x,y) in
! the index unit coordinates. If the target is outside the region of the
! grid that allows for centered 4-point interpolation in both directions,
! the logical Failure Flag, FF, is returned .true. and no interpolation
! is attempted; otherwise, for a normal condition, the interpolated values
! is returned as "a".
! Real 4*4 array, rp4, is the constant coefficients matrix precomputed in
! subroutine getrp4
!===========================================================================
use pkind, only: spi,dp
implicit none
integer(spi),                   intent(in ):: lx,mx,ly,my
real(dp),dimension(0:3,0:3),    intent(in ):: rp4
real(dp),dimension(lx:mx,ly:my),intent(in ):: axy
real(dp),                       intent(in ):: x,y
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp),dimension(-1:2):: ax,xp,yp,wx,wy
real(dp)                :: rx,ry
integer(spi)            :: i,ix,iy
!============================================================================
ix=floor(x); iy=floor(y)
ff=ix<lx+1 .or. ix>=mx-1 .or. iy<ly+1 .or. iy>=my-1
if(ff)return
rx=x-ix; ry=y-iy
call pown(3,rx,xp); call pown(3,ry,yp)
wx=matmul(xp,rp4); wy=matmul(yp,rp4)
do i=-1,2; ax(i)=dot_product(wy,axy(ix+i,iy-1:iy+2)); enddo
a=dot_product(wx,ax)
end subroutine abswi4_2

!============================================================================
subroutine abswi6_1(lx,mx,rp6,ax,x, a,ff)!                           [abswi6]
!============================================================================
! Apply the B-Spline-Weighted Interpolation with 6-point stnecil to a 1D
! array of real values, ax(lx:mx), with a target point at (x) in
! the index unit coordinates. If the target is outside the region of the
! grid that allows for centered 6-point interpolation,
! the logical Failure Flag, FF, is returned .true. and no interpolation
! is attempted; otherwise, for a normal condition, the interpolated values
! is returned as "a".
! Real 6*6 array, rp6, is the constant coefficients matrix precomputed in
! subroutine getrp6
!===========================================================================
use pkind, only: spi,dp
implicit none
integer(spi),                   intent(in ):: lx,mx
real(dp),dimension(0:5,0:5),    intent(in ):: rp6
real(dp),dimension(lx:mx),      intent(in ):: ax
real(dp),                       intent(in ):: x
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp),dimension(-2:3):: xp,wx
real(dp)                :: rx
integer(spi)            :: i,ix
!============================================================================
ix=floor(x)
ff=ix<lx+2 .or. ix>=mx-2
if(ff)return
rx=x-ix
call pown(5,rx,xp)
wx=matmul(xp,rp6)
a=dot_product(wx,ax(ix-2:ix+3))
end subroutine abswi6_1
!============================================================================
subroutine abswi6_2(lx,mx,ly,my,rp6,axy,x,y, a,ff)!                  [abswi6]
!============================================================================
! Apply the B-Spline-Weighted Interpolation with 6-point stnecil to a 2D
! array of real values, axy(lx:mx,ly:my), with a target point at (x,y) in
! the index unit coordinates. If the target is outside the region of the
! grid that allows for centered 6-point interpolation in both directions,
! the logical Failure Flag, FF, is returned .true. and no interpolation
! is attempted; otherwise, for a normal condition, the interpolated
! values is returned as "a".
! Real 6*6 array, rp6, is the constant coefficients matrix precomputed in
! subroutine getrp6
!===========================================================================
use pkind, only: spi,dp
implicit none
integer(spi),                   intent(in ):: lx,mx,ly,my
real(dp),dimension(0:5,0:5),    intent(in ):: rp6
real(dp),dimension(lx:mx,ly:my),intent(in ):: axy
real(dp),                       intent(in ):: x,y
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp),dimension(-2:3):: ax,xp,yp,wx,wy
real(dp)                :: rx,ry
integer(spi)            :: i,ix,iy
!============================================================================
ix=floor(x); iy=floor(y)
ff=ix<lx+2 .or. ix>=mx-2 .or. iy<ly+2 .or. iy>=my-2
if(ff)return
rx=x-ix; ry=y-iy
call pown(5,rx,xp); call pown(5,ry,yp)
wx=matmul(xp,rp6); wy=matmul(yp,rp6)
do i=-2,3; ax(i)=dot_product(wy,axy(ix+i,iy-2:iy+3)); enddo
a=dot_product(wx,ax)
end subroutine abswi6_2

!============================================================================
subroutine abswi8_1(lx,mx,rp8,ax,x, a,ff)!                           [abswi8]
!============================================================================
! Apply the B-Spline-Weighted Interpolation with 8-point stnecil to a 1D
! array of real values, ax(lx:mx), with a target point at (x) in
! the index unit coordinates. If the target is outside the region of the
! grid that allows for centered 8-point interpolation,
! the logical Failure Flag, FF, is returned .true. and no interpolation
! is attempted; otherwise, for a normal condition, the interpolated values
! is returned as "a".
! Real 8*8 array, rp8, is the constant coefficients matrix precomputed in
! subroutine getrp8
!===========================================================================
use pkind, only: spi,dp
implicit none
integer(spi),                   intent(in ):: lx,mx
real(dp),dimension(0:7,0:7),    intent(in ):: rp8
real(dp),dimension(lx:mx),      intent(in ):: ax
real(dp),                       intent(in ):: x
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp),dimension(-3:4):: xp,wx
real(dp)                :: rx
integer(spi)            :: i,ix
!============================================================================
ix=floor(x)
ff=ix<lx+3 .or. ix>=mx-3
if(ff)return
rx=x-ix
call pown(7,rx,xp)
wx=matmul(xp,rp8)
a=dot_product(wx,ax(ix-3:ix+4))
end subroutine abswi8_1
!============================================================================
subroutine abswi8_2(lx,mx,ly,my,rp8,axy,x,y, a,ff)!                  [abswi8]
!============================================================================
! Apply the B-Spline-Weighted Interpolation with 8-point stnecil to a 2D
! array of real values, axy(lx:mx,ly:my), with a target point at (x,y) in
! the index unit coordinates. If the target is outside the region of the
! grid that allows for centered 8-point interpolation in both directions,
! the logical Failure Flag, FF, is returned .true. and no interpolation
! is attempted; otherwise, for a normal condition, the interpolated values
! is returned as "a".
! Real 8*8 array, rp8, is the constant coefficients matrix precomputed in
! subroutine getrp8
!===========================================================================
use pkind, only: spi,dp
implicit none
integer(spi),                   intent(in ):: lx,mx,ly,my
real(dp),dimension(0:7,0:7),    intent(in ):: rp8
real(dp),dimension(lx:mx,ly:my),intent(in ):: axy
real(dp),                       intent(in ):: x,y
real(dp),                       intent(out):: a
logical,                        intent(out):: ff
!----------------------------------------------------------------------------
real(dp),dimension(-3:4):: ax,xp,yp,wx,wy
real(dp)                :: rx,ry
integer(spi)            :: i,ix,iy
!============================================================================
ix=floor(x); iy=floor(y)
ff=ix<lx+3 .or. ix>=mx-3 .or. iy<ly+3 .or. iy>=my-3
if(ff)return
rx=x-ix; ry=y-iy
call pown(7,rx,xp); call pown(7,ry,yp)
wx=matmul(xp,rp8); wy=matmul(yp,rp8)
do i=-3,4; ax(i)=dot_product(wy,axy(ix+i,iy-3:iy+4)); enddo
a=dot_product(wx,ax)
end subroutine abswi8_2

end module pbswi

!#

