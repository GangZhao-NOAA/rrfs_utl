program regrid_esg2rll_iplib
!================================================================================
 use omp_lib
 use netcdf
 use pkind, only: dp, sp, dpi, spi
 use pietc, only: dtor,rtod
 use pesg
 use pbswi, only: abswi2,abswi4,abswi6,abswi8,getrp4,getrp6,getrp8
 use ip_mod                               ! ip lib
 use grib_mod                             ! g2 lib(grib)
 use bacio_module                         ! prerequisite for g2 lib (grib)

 implicit none

 character(len=100)      :: input_file
 character(len=100)      :: input_file_nc
 character(len=100)      :: latlon_file_nc
 character(len=100)      :: output_file_nc
 character(len=100)      :: outfile_chk
 integer                 :: ncid, varid
 integer                 :: varid_glon, varid_glat
 integer                 :: varid_rlon, varid_rlat
 integer                 :: varid_orig, varid_ipol
 integer                 :: varid_diff
 integer                 :: dimid_x, dimid_y, dimid_time
 character(len=20)       :: varname_nc

 real                    :: fill
 real                    :: diff, maxdiffx, maxdiffy
 integer                 :: nret, iopt
 character(len=4)        :: cgrid_i
 integer                 :: badpts
 integer                 :: ii, jj, nn
 integer                 :: iii, jjj
 integer                 :: n_bitmap

 integer                 :: iunit, iret, lugi
 integer                 :: iunito_chk

! integer                 :: ip, ipopt(20)     ! parameters for interpolation (iplib)
 integer                 :: interp_opt

 integer                 :: j, jdisc, jpdtn, jgdtn, k
 integer                 :: jids(200), jgdt(200), jpdt(200)
!integer                 :: idim_output, jdim_output
 integer, allocatable    :: igdtmpl_i(:)
 integer                 :: igdtlen_i
 integer                 :: igdtnum_i

 integer                 :: idim_input, jdim_input ! dimension size pre-set in namelist or hardwired in code
 integer                 :: mi
 integer                 :: imdl_i, jmdl_i         ! dimension size read in grib2 data file
 integer                 :: npts_i
 logical                 :: unpack
 integer, allocatable    :: ibi(:)
!integer                 :: ibi
!logical*1, allocatable  :: input_bitmap(:,:)      ! 2D Data in 1D slice
 logical*1, allocatable  :: input_bitmap(:)        ! 2D Data in 1D slice
!real, allocatable       :: input_data(:,:)        ! 2D Data in 1D slice
 real, allocatable       :: input_data(:)          ! 2D Data in 1D slice

 type(gribfield)         :: gfld_input

!--- used in gdswzd
 real, allocatable :: xpts_i(:,:), ypts_i(:,:)
 real, allocatable :: glat_i(:,:), glon_i(:,:)
 real, allocatable :: crot_i(:,:), srot_i(:,:)
 real, allocatable :: xlon_i(:,:), xlat_i(:,:)
 real, allocatable :: ylon_i(:,:), ylat_i(:,:), area_i(:,:)

 real(dp), allocatable :: xpts_i2(:,:), ypts_i2(:,:)  ! x/y in ESG x/y coordinates

 real(dp), allocatable :: rotlat(:), rotlon(:)     ! rotated lat and lon

!real(dp) :: rp4(4,4)
 real(dp) :: rp6(6,6)
 real(dp) :: rp8(8,8)

! integer                 :: mo
! integer                 :: no
 integer                 :: km
 integer(spi)            :: ipts_o, jpts_o         ! dimension size read in netcdf data file
 integer                 :: npts_o

! integer                 :: igdtlen_o
! integer                 :: igdtnum_o
! integer, allocatable    :: igdtmpl_o(:)

! integer, allocatable    :: ibo(:)
!!integer                 :: ibo
!!real, allocatable       :: output_data(:,:)           ! 2D Data in 1D slice
! real, allocatable       :: output_data(:)             ! 2D Data in 1D slice
!!logical*1, allocatable  :: output_bitmap(:,:)         ! 2D Data in 1D slice
! logical*1, allocatable  :: output_bitmap(:)           ! 2D Data in 1D slice
!!real, allocatable       :: output_glat(:,:)           ! 2D Data in 1D slice
! real, allocatable       :: output_glat(:)             ! 2D Data in 1D slice
!!real, allocatable       :: output_glon(:,:)           ! 2D Data in 1D slice
! real, allocatable       :: output_glon(:)             ! 2D Data in 1D slice

 real, allocatable       :: glon_o(:,:), glat_o(:,:)
 real(dp), allocatable       :: slmask_o(:,:)
 real(dp), allocatable       :: howv_o(:,:)                 ! data on ESG grid
 real(dp), allocatable       :: howv_o1(:,:)                ! data on Rotated Lat-Lon grid
 real(dp), allocatable       :: howv_o2(:,:)                ! data on Rotated Lat-Lon grid
 real(dp), allocatable       :: howv_diff(:,:)              ! data on Rotated Lat-Lon grid
 real, allocatable       :: xpts_o2(:,:), ypts_o2(:,:)  ! x/y in ESG x/y coordinates
 real, allocatable       :: diff_xy(:,:)

!integer                 :: xdimid, ydimid, timedimid

!integer                 :: n_center

 logical           :: l_chk_gdswzd     ! if call gdswzd with iopt=-1 to check
 logical           :: l_chk_bitmap     ! check bitmap with False and the value of data
 logical           :: l_chk_undefi     ! check undefined value

! Used for gnomonic transformation (Jim Purser)
 real(dp) :: A, Kappa, pazi, plon, plat
 real(dp) :: delx, dely
 real(dp) :: dlon, dlat
 logical :: ff
!logical :: misha
 real(dp), dimension(2) :: xm
!real(dp), parameter          :: zero=0_dp
!real(dp), parameter          :: one=1_dp
 real(dp), parameter          :: two=2_dp
 real(dp), parameter          :: lam=0.8
 real(dp)                     :: m_arcx,m_arcy,q
 real(dp)                     :: m_delx,m_dely
!real(dp)                     :: delxre,delyre
 real(dp)                     :: arcx,arcy

 integer                      :: nxh,nyh, nx,ny, nxm,nym, lx,ly


 interface

   subroutine baopenr(iunit, input_file, iret)
     integer,          intent(in   ) :: iunit
    character(len=*), intent(in   ) :: input_file
     integer,          intent(  out) :: iret
   end subroutine baopenr

   subroutine baclose(iunit, iret)
     integer,          intent(in   ) :: iunit
     integer,          intent(  out) :: iret
   end subroutine baclose

   subroutine getgb2(lugb, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
                     unpack, k, gfld, iret)
     use grib_mod
     integer,               intent(in   ) :: lugb, lugi, j, jdisc, jpdtn, jgdtn
     integer, dimension(:), intent(in   ) :: jids(*), jpdt(*), jgdt(*)
     logical,               intent(in   ) :: unpack
     integer,               intent(  out) :: k, iret
     type(gribfield),       intent(  out) :: gfld
   end subroutine getgb2

   subroutine gf_free(gfld)
     use grib_mod
     type(gribfield),  intent(in   ) :: gfld
   end subroutine gf_free

 end interface

!---------------------------------------------------------------------------
! 0. Initialisation (dimensions, etc.)

! 1. Set up the grids (input/output)
! 1.1 Rotated Latlon (RLL) grid (output grid)
!     read grib2 file for grid template and grid info
!     read in background data (2D field) for output grid 
!     (to cover the grids which are not covered by the input grid data)

!     set up the geo latlon of RLL grid

! 1.2 Extended Schimdt Gnomonic (ESG) grid (input grid)
!     read netcdf file
!     read grid latlon of ESG grid
!     read ESG grid information from global attributes of netcdf file
!     call bestesg_geo to set up the ESG grid parameters
!     read the input data (2D field)


! 2. Set up x-y coordinates based on ESG grid
!    (the origin of x-y coordinates is the center of ESG project grid)
! 2.1 convert geo latlon of ESG to x-y in ESG-based x-y coordinates

! 2.2 convert geo latlon of RLL to x-y in ESG-based x-y coordinates
!
! 3. Interpolation from input grid (ESG) to oupput grid (RLL)




!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  l_chk_gdswzd = .true.
  l_chk_bitmap = .false.
  l_chk_undefi = .false.

  interp_opt = 2     ! 2: abswi2;
                     ! 4: abswi4;
                     ! 6: abswi6;
                     ! 8: abswi8;
!---------------------------------------------------------------------------
! open the grib 2 file containing data to be interpolated.
! for this example, there are only one data record (HTSGW).
!---------------------------------------------------------------------------

 iunit=9
 input_file="./input.data.grib2"
 call baopenr(iunit, input_file, iret)
 if (iret /= 0) then
   write(6,*) 'return from baopenr: ',iret
   stop 'Error: baopenr failed.'
 end if

!---------------------------------------------------------------------------
! prep for call to g2 library to degrib data. the data are on a 
! rotated-lat/lon grid with i/j dimension of 360/181. 
!---------------------------------------------------------------------------

 idim_input = 4881 ! the i/j dimensions of input grid
 jdim_input = 2961
 mi         = idim_input * jdim_input   ! total number of pts, input grid

 jdisc   = -1         ! search for any discipline
 jpdtn   = -1         ! search for any product definition template number

!jgdtn   =  0         ! search for grid definition template number 0 - regular lat/lon grid
 jgdtn   =  1         ! search for grid definition template number 1 - rotated lat/lon grid

 jids    = -9999      ! array of values in identification section, set to wildcard
 jgdt    = -9999      ! array of values in grid definition template 3.m
!jgdt(8) = idim_input ! search for grid with i/j of 4881/2961
!jgdt(9) = jdim_input
 jpdt    = -9999      ! array of values in product definition template 4.n
 unpack  = .true.     ! unpack data
 lugi    = 0          ! no index file
!lugi    = iunit

 nullify(gfld_input%idsect)
 nullify(gfld_input%local)
 nullify(gfld_input%list_opt)
 nullify(gfld_input%igdtmpl)  ! holds the grid definition template information
 nullify(gfld_input%ipdtmpl)
 nullify(gfld_input%coord_list)
 nullify(gfld_input%idrtmpl)
 nullify(gfld_input%bmap)     ! holds the bitmap
 nullify(gfld_input%fld)      ! holds the data

!---------------------------------------------------------------------------
! degrib the data.  non-zero "iret" indicates a problem during degrib.
!---------------------------------------------------------------------------

 km = 1                  ! number of records to interpolate

 do j = 0, (km-1)    ! number of records to skip

   call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld_input, iret)

   if (iret /= 0) then
     write(6,*) 'return from sub getgb2: ', iret
     stop 'Error: getb2 failed.'
   else
     write(6,*) ' checking the data in array gfld_input read by getgb2: '
     write(6,*) ' version: ',    gfld_input%version
     write(6,*) ' discipline: ', gfld_input%discipline
     write(6,*) ' idsectlen: ',  gfld_input%idsectlen
     write(6,*) ' idsect(:):, ', gfld_input%idsect(1:gfld_input%idsectlen)
     write(6,*) ' locallen: ',   gfld_input%locallen
     if (gfld_input%locallen > 0) &
       write(6,*) ' local(:): ',   gfld_input%local(1:gfld_input%locallen)
     write(6,*) ' ifldnum: ',    gfld_input%ifldnum
     write(6,*) ' griddef: ',    gfld_input%griddef
     write(6,*) ' ngrdpts: ',    gfld_input%ngrdpts
     write(6,*) ' numoct_opt: ', gfld_input%numoct_opt
     write(6,*) ' interp_opt: ', gfld_input%interp_opt
     write(6,*) ' num_opt: ',    gfld_input%num_opt
     write(6,*) ' igdtnum: ',    gfld_input%igdtnum
     write(6,*) ' igdtlen: ',    gfld_input%igdtlen
     write(6,*) ' igdtmpl(:): ', gfld_input%igdtmpl(1:gfld_input%igdtlen)
     write(6,*) ' ipdtnum: ',    gfld_input%ipdtnum
     write(6,*) ' ipdtlen: ',    gfld_input%ipdtlen
     write(6,*) ' ipdtmpl(:): ', gfld_input%ipdtmpl(1:gfld_input%ipdtlen)
     write(6,*) ' num_coord: ',  gfld_input%num_coord
     if (gfld_input%num_coord > 0) &
       write(6,*) ' coord_list(:) ', gfld_input%coord_list(:)
     write(6,*) ' ndpts: ',      gfld_input%ndpts
     write(6,*) ' idrtnum: ',    gfld_input%idrtnum
     write(6,*) ' idrtlen: ',    gfld_input%idrtlen
     write(6,*) ' idrtmpl(:): ', gfld_input%idrtmpl(1:gfld_input%idrtlen)
     write(6,*) ' unpacked: ',   gfld_input%unpacked
     write(6,*) ' expanded: ',   gfld_input%expanded
     write(6,*) ' ibmap: ',      gfld_input%ibmap
     write(6,*) ' length of bmap(:) ',     size(gfld_input%bmap)
     write(6,*) ' fld(1): ',       gfld_input%fld(1)
!    n_center = nint(gfld_input%ngrdpts/2.0)
     write(6,*) ' fld(ngrdpts/2): ', gfld_input%fld(nint(gfld_input%ngrdpts / 2.0))
   end if

   imdl_i = gfld_input%igdtmpl(8)               ! x-dim size of rotated latlon grid
   jmdl_i = gfld_input%igdtmpl(9)               ! y-dim size of rotated latlon grid
   npts_i = imdl_i * jmdl_i

!--- input grid template information
   igdtnum_i = jgdtn
   write(cgrid_i,'(I4.4)') igdtnum_i
   if (igdtnum_i == 0) then
     write(6,*)"input model grid is defined on regular lat-lon grid"
   else if (igdtnum_i == 1) then
     write(6,*)"input model grid is defined on rotated lat-lon grid"
   else
     write(6,*)"input model grid is defined on the grid with grid number = ",igdtnum_i
   end if
   igdtlen_i = gfld_input%igdtlen
   allocate(igdtmpl_i(igdtlen_i))
   igdtmpl_i(:) = gfld_input%igdtmpl(:)

!---------------------------------------------------------------------------
! does input data have a bitmap?
!---------------------------------------------------------------------------
   allocate(ibi(km))
!  allocate(input_bitmap(npts_i, km))
   allocate(input_bitmap(npts_i))
!  allocate(input_data(npts_i, km))
   allocate(input_data(npts_i))

   if (gfld_input%ibmap==0) then  ! input data has bitmap
     ibi                   = 1        ! tell ipolates to use bitmap
     input_bitmap(:)       = gfld_input%bmap
   else                           ! no bitmap, data everywhere
     ibi                   = 0        ! tell ipolates there is no bitmap
     input_bitmap(:)       = .true.
   endif

   input_data(:)           = gfld_input%fld  ! the input data field

   if (l_chk_bitmap .or. l_chk_undefi) then

     n_bitmap = 0
     do nn = 1, npts_i

       if ( .not. input_bitmap(nn) .and. l_chk_bitmap ) then
!      if (       input_bitmap(nn) .and. l_chk_bitmap) then
         n_bitmap = n_bitmap + 1
         if ( mod(n_bitmap, 1000) == 0 ) then
           write(6,*) '(input) nn  bitmap  data: ', nn, input_bitmap(nn), input_data(nn), n_bitmap
         end if
       end if

       if ( abs(input_data(nn)) .ge. 50.0 .and. l_chk_undefi ) then
           write(6,*) '(input) nn  bitmap  abs(data)>=50.0): ', nn, input_bitmap(nn), input_data(nn)
       end if

     end do
   end if
 
!---------------------------------------------------------------------------
! first, call gdswzd to calculate lat/lon for each grid point of the input grid
! (rotated-latlon grid).
!---------------------------------------------------------------------------

   fill = -9999.

   allocate (xpts_i(imdl_i,jmdl_i),ypts_i(imdl_i,jmdl_i))
   allocate (glat_i(imdl_i,jmdl_i),glon_i(imdl_i,jmdl_i))
   allocate (crot_i(imdl_i,jmdl_i),srot_i(imdl_i,jmdl_i))
   allocate (xlon_i(imdl_i,jmdl_i),xlat_i(imdl_i,jmdl_i))
   allocate (ylon_i(imdl_i,jmdl_i),ylat_i(imdl_i,jmdl_i))
   allocate (area_i(imdl_i,jmdl_i))

   xpts_i = fill
   ypts_i = fill
   glat_i = fill
   glon_i = fill
   crot_i = fill
   srot_i = fill
   xlon_i = fill
   xlat_i = fill
   ylon_i = fill
   ylat_i = fill
   area_i = fill

   iopt = 0

   write(6,*) ' call gdswzd with iopt=', iopt
   call gdswzd(igdtnum_i, igdtmpl_i, igdtlen_i, iopt, npts_i, fill, xpts_i, ypts_i, glon_i, glat_i, &
               nret, crot_i, srot_i, xlon_i, xlat_i, ylon_i, ylat_i, area_i)

   if (nret /= npts_i) then
     write(6,*)'ERROR. WRONG NUMBER OF POINTS RETURNED FROM gdswzs (iopt=0): ',nret,npts_i
     stop 33
   else
     write(6,*)' successfully done with gdswzs (iopt=0): ',nret,npts_i
   endif

   print*,'LAT/LON POINT(1,1):   ',glat_i(1,1),glon_i(1,1)
   print*,'LAT/LON POINT(1,JM):  ',glat_i(1,jmdl_i),glon_i(1,jmdl_i)
   print*,'LAT/LON POINT(IM,1):  ',glat_i(imdl_i,1),glon_i(imdl_i,1)
   print*,'LAT/LON POINT(IM,JM): ',glat_i(imdl_i,jmdl_i),glon_i(imdl_i,jmdl_i)
   print*,'xpts/ypts POINT( 1, 1):  ',xpts_i(1,1),           ypts_i(1,1)
   print*,'xpts/ypts POINT( 1,JM):  ',xpts_i(1,jmdl_i),      ypts_i(1,jmdl_i)
   print*,'xpts/ypts POINT(IM, 1):  ',xpts_i(imdl_i,1),      ypts_i(imdl_i,1)
   print*,'xpts/ypts POINT(IM,JM):  ',xpts_i(imdl_i,jmdl_i), ypts_i(imdl_i,jmdl_i)

   outfile_chk = "./grid" // trim(cgrid_i) // ".iopt0.bin"
   iunito_chk = 101
   open (iunito_chk, file=trim(outfile_chk), access='direct', recl=imdl_i*jmdl_i*4)
   write(iunito_chk, rec=1) real(glat_i,4)
   write(iunito_chk, rec=2) real(glon_i,4)
   write(iunito_chk, rec=3) real(xpts_i,4)
   write(iunito_chk, rec=4) real(ypts_i,4)
   write(iunito_chk, rec=5) real(crot_i,4)
   write(iunito_chk, rec=6) real(srot_i,4)
   write(iunito_chk, rec=7) real(xlon_i,4)
   write(iunito_chk, rec=8) real(xlat_i,4)
   write(iunito_chk, rec=9) real(ylon_i,4)
   write(iunito_chk, rec=10) real(ylat_i,4)
   write(iunito_chk, rec=11) real(area_i,4)
   close (iunito_chk)

   if (l_chk_gdswzd) then
!--- the first call to gdswzd computed the lat/lon at each point.  now,
!--- given that lat/lon, compute the i/j map coordinate indices.  it should be reversable.

     iopt = -1
     xpts_i=fill
     ypts_i=fill

     write(6,*) ' call gdswzd with iopt=', iopt
     call gdswzd(igdtnum_i, igdtmpl_i, igdtlen_i, iopt, npts_i, fill, xpts_i, ypts_i, glon_i, glat_i, &
                 nret, crot_i, srot_i, xlon_i, xlat_i, ylon_i, ylat_i, area_i)

     if (nret /= npts_i) then
       write(6,*)'ERROR. WRONG NUMBER OF POINTS RETURNED FROM gdswzs (iopt=-1): ',nret,npts_i
       stop 34
     else
       write(6,*)' successfully done with gdswzs (iopt=-1): ',nret,npts_i
     endif

     outfile_chk = "./grid" // trim(cgrid_i) // ".ioptm1.bin"
     iunito_chk = 101
     open (iunito_chk, file=trim(outfile_chk), access='direct', recl=imdl_i*jmdl_i*4)
     write(iunito_chk, rec=1) real(glat_i,4)
     write(iunito_chk, rec=2) real(glon_i,4)
     write(iunito_chk, rec=3) real(xpts_i,4)
     write(iunito_chk, rec=4) real(ypts_i,4)
     write(iunito_chk, rec=5) real(crot_i,4)
     write(iunito_chk, rec=6) real(srot_i,4)
     write(iunito_chk, rec=7) real(xlon_i,4)
     write(iunito_chk, rec=8) real(xlat_i,4)
     write(iunito_chk, rec=9) real(ylon_i,4)
     write(iunito_chk, rec=10) real(ylat_i,4)
     write(iunito_chk, rec=11) real(area_i,4)
     close (iunito_chk)

!------------------------------------------------------------------------------
!    did the second call to gdswzd work?
!
!    note: the gdswzdcb routine works on a grid that is tilted 45 degrees for
!          Arakawa-E grid (eg, Eta model), so the internal i/j's do not match 
!          the normal convention.  account for this.
!          (see example in
!          "EMC_NCEPLIBS-ip/tests/reg_tests/gdswzd/sorc/gdswzd_driver.f90")
!          The following checkness could not work for E-grid.
!------------------------------------------------------------------------------

     maxdiffx = -99999.
     maxdiffy = -99999.

     badpts=0
     do jj = 1, jmdl_i            ! conventional grid-index i of grid points
     do ii = 1, imdl_i            ! conventional grid-index j of grid points
       diff = abs(float(ii)-xpts_i(ii,jj))
       maxdiffx = max(maxdiffx, diff)
       if ( diff > .01) then
         print*,'BAD X POINT: ',ii,jj,xpts_i(ii,jj),ypts_i(ii,jj)
         badpts=badpts+1
       endif 
       diff = abs(float(jj)-ypts_i(ii,jj))
       maxdiffy = max(maxdiffy, diff)
       if ( diff  > .01) then
         print*,'BAD Y POINT: ',ii,jj,xpts_i(ii,jj),ypts_i(ii,jj)
         badpts=badpts+1
       endif 
     enddo
     enddo

     if (badpts > 0) &
       write(6,'(1X,A,I6)') "NUMBER OF BAD POINTS(mismatch from two calls of gdswzdcb with iopt=0 & -1) : ", badpts
     write(6,*) 'MAX DIFFERENCES IN X/Y CALCULATIONS: ', maxdiffx, maxdiffy

   end if    ! l_chk_gdswzd == .true.
!
!--- read in the lat/lon of output grid from netcdf file
!
   ipts_o = 3950              ! x-dim size of ESG grid
   jpts_o = 2700              ! y-dim size of ESG grid
   npts_o = ipts_o * jpts_o
   write(6,*) 'output grid dimensions -- ipts_o:', ipts_o, ' jpts_o:', jpts_o
   allocate(glat_o(ipts_o, jpts_o))
   allocate(glon_o(ipts_o, jpts_o))
   allocate(slmask_o(ipts_o, jpts_o))
   allocate(howv_o(ipts_o, jpts_o))

   latlon_file_nc = "./latlon.data.nc"
   call check( nf90_open(latlon_file_nc, nf90_nowrite, ncid) )

   varname_nc = "geolon"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, glon_o ))
   write(6,*) 'read in ', trim(adjustl(varname_nc)),'=',glon_o(1,1),glon_o(1,jpts_o),glon_o(ipts_o,jpts_o),glon_o(ipts_o,1)

   varname_nc = "geolat"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, glat_o ))
   write(6,*) 'read in ', trim(adjustl(varname_nc)),'=',glat_o(1,1),glat_o(1,jpts_o),glat_o(ipts_o,jpts_o),glat_o(ipts_o,1)

   call check( nf90_close(ncid) )

   input_file_nc = "./input.data.nc"
   call check( nf90_open(input_file_nc, nf90_nowrite, ncid) )

   varname_nc = "slmsk"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, slmask_o ))
   write(6,*) 'read in ', trim(adjustl(varname_nc)),'=',slmask_o(1,1),slmask_o(1,jpts_o),slmask_o(ipts_o,jpts_o),slmask_o(ipts_o,1)

   varname_nc = "howv"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, howv_o ))
   write(6,*) 'read in ', trim(adjustl(varname_nc)),'=',howv_o(1,1),howv_o(1,jpts_o),howv_o(ipts_o,jpts_o),howv_o(ipts_o,1)

   call check( nf90_close(ncid) )

!-------------------------------------------------------------------------------!

   lx = -3962
   ly = -2712
   delx = 0.0129656309291837_dp      ! half grid/cell size, unit of degree
   delx = delx * dtor                ! degree to radian
   dely = 0.0132456552576884_dp
   dely = dely * dtor
   plat = 55.0_dp
   plon = 247.5_dp     ! -112.5_dp;  
   pazi = 0.0_dp

   nxh=-lx
   nyh=-ly
   nx=2*nxh
   ny=2*nyh
   nxm=nx-1
   nym=ny-1

   arcx=delx*nxh
   arcy=dely*nyh
   write(6,'(1X,A,7(1x,D))')' before calling bestesg_geo: A,Kappa,delx,dely,plat,plon,pazi=',A,Kappa,delx,dely,plat,plon,pazi
   call bestesg_geo(lam,arcx,arcy, a,kappa,m_arcx,m_arcy,q,ff)
   m_delx=m_arcx/nxh ! Map-space grid steps in x
   m_dely=m_arcy/nyh ! Map-space grid steps in y
   delx=m_delx
   dely=m_dely
   write(6,'(1X,A,7(1x,D))')' after  calling bestesg_geo: A,Kappa,delx,dely,plat,plon,pazi=',A,Kappa,delx,dely,plat,plon,pazi
   write(6,'(1X,A,7(1x,D))')' after  calling bestesg_geo: m_arcx, m_arcy, q, ff =',m_arcx, m_arcy, q, ff
!
   A=0.183131392268429_dp
   Kappa=-0.265835885178773_dp
   delx=0.0129656309291837_dp    ! in deg., not rad. (*6370*pi/180=1.4415 km, half grid/cell size)
   delx = delx * dtor            ! deg. to radian
   dely=0.0132456552576884_dp
   dely = dely * dtor
   plat=55.0_dp             ! center lat of gnomonic grid
   plon=-112.5_dp           ! center lon of gnomonic grid  ! -112.5_dp = 247.5_dp
   pazi=0.0_dp
   m_arcx = delx * nxh
   m_arcy = dely * nyh
   write(6,'(1X,A,7(1x,D))')' using the pre-defined     : A,Kappa,delx,dely,plat,plon,pazi=',A,Kappa,delx,dely,plat,plon,pazi
   write(6,'(1X,A,7(1x,D))')'       the pre-defined     : m_arcx, m_arcy        =',m_arcx, m_arcy
!-------------------------------------------------------------------------------!
!--- convert lat/lon of rotated latlon grid points to x/y in the ESG grid x-y
!    coordinates (with A/Kappa/plat/plon/pazi/delx/dely of ESG)
   allocate (xpts_i2(imdl_i,jmdl_i),ypts_i2(imdl_i,jmdl_i))
   xpts_i2 = fill
   ypts_i2 = fill
   do jj=1,jmdl_i
   do ii=1,imdl_i
     dlat=glat_i(ii,jj)
     dlon=glon_i(ii,jj)
     call gtoxm_ak_dd_g(A,Kappa,plat,plon,pazi,two*delx,two*dely,dlat,dlon,xm,ff) !  multiply delx/dely by 2.0 to get values on compuational grid
!    call gtoxm_ak_dd_g(A,Kappa,plat,plon,pazi,    delx,    dely,dlat,dlon,xm,ff) !  if not multiply delx/dely by 2.0 to get values on compuational grid
     xpts_i2(ii,jj)=xm(1)
     ypts_i2(ii,jj)=xm(2)
     xpts_i2(ii,jj) = (ipts_o+1)/two + xpts_i2(ii,jj)  ! move origin of X/Y coordinate from center
                                                       ! of ESG grid to lower-left corner of ESG grid
     ypts_i2(ii,jj) = (jpts_o+1)/two + ypts_i2(ii,jj)
   enddo
   enddo
   print*,'Rotated LatLon grid dimensions ipts,jpts=',imdl_i,jmdl_i
   print*,'max/min xpts of RLL in ESG x-y coordinates =',maxval(xpts_i2),minval(xpts_i2)
   print*,'max/min ypts of RLL in ESG x-y coordinates =',maxval(ypts_i2),minval(ypts_i2)
   print*,'xpts/ypts POINT( 1, 1):  ',xpts_i2(1,1),           ypts_i2(1,1)
   print*,'xpts/ypts POINT( 1,JM):  ',xpts_i2(1,jmdl_i),      ypts_i2(1,jmdl_i)
   print*,'xpts/ypts POINT(IM, 1):  ',xpts_i2(imdl_i,1),      ypts_i2(imdl_i,1)
   print*,'xpts/ypts POINT(IM,JM):  ',xpts_i2(imdl_i,jmdl_i), ypts_i2(imdl_i,jmdl_i)

!--- convert lat/lon of ESG grid points to x/y in the ESG grid x-y
!    coordinates (with A/Kappa/plat/plon/pazi/delx/dely of ESG)
   allocate(xpts_o2(ipts_o, jpts_o), ypts_o2(ipts_o, jpts_o))
   allocate(diff_xy(ipts_o, jpts_o))
   xpts_o2 = fill
   ypts_o2 = fill
   do jj=1,jpts_o
   do ii=1,ipts_o
     dlat=glat_o(ii,jj)
     dlon=glon_o(ii,jj)
!    routine from Jim Purser.  Transform lat/lon to x,y for the ESG regional grid.
     call gtoxm_ak_dd_g(A,Kappa,plat,plon,pazi,two*delx,two*dely,dlat,dlon,xm,ff)  !    multiply delx/dely by 2.0 to get values on compuational grid
     xpts_o2(ii,jj)=xm(1) ! origin at map projection center
     ypts_o2(ii,jj)=xm(2) ! origin at map projection center
     xpts_o2(ii,jj) = (ipts_o+1)/two + xpts_o2(ii,jj)  ! convert to lower left hand corner of grid
     ypts_o2(ii,jj) = (jpts_o+1)/two + ypts_o2(ii,jj)  ! convert to lower left hand corner of grid
     diff_xy(ii,jj)=sqrt((xpts_o2(ii,jj)-float(ii))**2 + (ypts_o2(ii,jj)-float(jj))**2)
   enddo
   enddo
   print*,'ESG grid dimensions ipts,jpts=',ipts_o,jpts_o
   print*,'max/min difference of x-y coordiante ESG grid =',maxval(diff_xy),minval(diff_xy)
   print*,'max/min xpts ESG in ESG x-y coordiantes =',maxval(xpts_o2),minval(xpts_o2)
   print*,'max/min ypts ESG in ESG x-y coordiantes =',maxval(ypts_o2),minval(ypts_o2)

!--- interpolation from ESG to RLL
   allocate(howv_o1(imdl_i, jmdl_i))
   allocate(howv_o2(imdl_i, jmdl_i))
   allocate(howv_diff(imdl_i, jmdl_i))
   nn = 0
   do jj = 1, jmdl_i
   do ii = 1, imdl_i
     nn = nn + 1
     howv_o1(ii,jj) = input_data(nn)                   ! intialized with firstguess
     howv_o2(ii,jj) = input_data(nn)                   ! intialized with firstguess
   end do
   end do
   
   if(interp_opt.eq.2)then
     print*,'Using BSWI-2, interp_opt= ', interp_opt
     nn = 0
     do jj=1,jmdl_i
     do ii=1,imdl_i
       iii = dint(xpts_i2(ii,jj))
       jjj = dint(ypts_i2(ii,jj))
       nn = nn + 1
!      howv_o2(ii,jj) = input_data(nn)                   ! intialized with firstguess
       if(iii .gt. 0 .and. jjj .gt. 0 .and. iii .le. (ipts_o - 1) .and. jjj .le. (jpts_o - 1) )then
         call abswi2(1,ipts_o,1,jpts_o,howv_o,xpts_i2(ii,jj),ypts_i2(ii,jj),howv_o2(ii,jj),ff)
       end if
!      howv_diff(ii,jj) = howv_o2(ii,jj) - howv_o1(ii,jj)
!!      if ( ibi(1) == 1 .and. .not. input_bitmap(nn) ) howv_o2(ii,jj) = 0.0
       if (                   .not. input_bitmap(nn) ) then
!!        howv_o2(ii,jj) = 0.0
         howv_o2(ii,jj) = howv_o1(ii,jj)
       end if
       howv_diff(ii,jj) = howv_o2(ii,jj) - howv_o1(ii,jj)
     enddo
     enddo
   elseif(interp_opt.eq.6)then
     print*,'Using BSWI-6, interp_opt= ', interp_opt
     nn = 0
     call getrp6(rp6)
     do jj=1,jmdl_i
     do ii=1,imdl_i
       nn = nn + 1
!      howv_o2(ii,jj) = input_data(nn)                   ! intialized with firstguess
       call abswi6(1,ipts_o,1,jpts_o,rp6,howv_o,xpts_i2(ii,jj),ypts_i2(ii,jj),howv_o2(ii,jj),ff)
!      howv_diff(ii,jj) = howv_o2(ii,jj) - howv_o1(ii,jj)
!      if ( ibi(1) == 1 .and. .not. input_bitmap(nn) ) howv_o2(ii,jj) = 0.0
       if (                   .not. input_bitmap(nn) ) then
         howv_o2(ii,jj) = 0.0
       end if
       howv_diff(ii,jj) = howv_o2(ii,jj) - howv_o1(ii,jj)
     enddo
     enddo
   elseif(interp_opt.eq.8)then
     print*,'Using BSWI-8, interp_opt= ', interp_opt
     nn = 0
     call getrp8(rp8)
     do jj=1,jmdl_i
     do ii=1,imdl_i
       nn = nn + 1
!      howv_o2(ii,jj) = input_data(nn)                   ! intialized with firstguess
       call abswi8(1,ipts_o,1,jpts_o,rp8,howv_o,xpts_i2(ii,jj),ypts_i2(ii,jj),howv_o2(ii,jj),ff)
!      howv_diff(ii,jj) = howv_o2(ii,jj) - howv_o1(ii,jj)
!      if ( ibi(1) == 1 .and. .not. input_bitmap(nn) ) howv_o2(ii,jj) = 0.0
       if (                   .not. input_bitmap(nn) ) then
         howv_o2(ii,jj) = 0.0
       end if
       howv_diff(ii,jj) = howv_o2(ii,jj) - howv_o1(ii,jj)
     enddo
     enddo
   else
     print*, "Unknown Interp_opt=", interp_opt, " stop!"
     stop(11)
   endif

   write(6,'(1X,A,2(1X,F))') "max/min diff of howv on RLL grid : ", maxval(howv_diff), minval(howv_diff)

!--- output the data (howv_01/howv_o2/howv_diff) on rotated latlon grid  to a netcdf file
!- Create the netcdf file
   output_file_nc = "howv_rll.nc"
   call check(nf90_create(trim(output_file_nc), nf90_netcdf4, ncid))   

!- Define the dimensions
   call check(nf90_def_dim(ncid,    "X",         imdl_i,    dimid_x))
   call check(nf90_def_dim(ncid,    "Y",         jmdl_i,    dimid_y))
   call check(nf90_def_dim(ncid, "Time", nf90_unlimited, dimid_time))

!- Define true lat/lon variables
!  (no rotated lat/lon defined here, which should be included and will be added later.)
   call check(nf90_def_var(ncid, "rotlon", nf90_double, (/dimid_x /),         varid_rlon))
   call check(nf90_def_var(ncid, "rotlat", nf90_double, (/dimid_y/),          varid_rlat))
   call check(nf90_def_var(ncid, "geolon", nf90_double, (/dimid_x, dimid_y/), varid_glon))
   call check(nf90_def_var(ncid, "geolat", nf90_double, (/dimid_x, dimid_y/), varid_glat))

!- Define data variables
   call check(nf90_def_var(ncid, "howv_orig", nf90_double, (/dimid_x, dimid_y, dimid_time/), varid_orig))
   call check(nf90_def_var(ncid, "howv_ipol", nf90_double, (/dimid_x, dimid_y, dimid_time/), varid_ipol))
   call check(nf90_def_var(ncid, "howv_diff", nf90_double, (/dimid_x, dimid_y, dimid_time/), varid_diff))

!- Add the attributes
   call check(nf90_put_att(ncid, nf90_global, 'description', 'wave height data fistguess/analysis/increment'))
   call check(nf90_put_att(ncid, nf90_global, 'note', 'Rotated Lat/Lon Grid'))
   call check(nf90_put_att(ncid, nf90_global, 'Projection', 'Rotated Lat/Lon Grid'))
   call check(nf90_put_att(ncid, nf90_global, 'lat_ll_rll', -37.000000))
   call check(nf90_put_att(ncid, nf90_global, 'lat_ur_rll',  37.000000))
   call check(nf90_put_att(ncid, nf90_global, 'dlat_rll',     0.025000))
   call check(nf90_put_att(ncid, nf90_global, 'lon_ll_rll', 299.000000))
   call check(nf90_put_att(ncid, nf90_global, 'lon_ur_rll',  61.000000))
   call check(nf90_put_att(ncid, nf90_global, 'dlon_rll',     0.025000))
   call check(nf90_put_att(ncid, nf90_global, 'latitude_center', real(plat)))
   call check(nf90_put_att(ncid, nf90_global, 'longitude_center', real(plon)))
   call check(nf90_put_att(ncid, nf90_global, 'latitude_south_pole_rotated', -35.0))
   call check(nf90_put_att(ncid, nf90_global, 'longitude_south_pole_rotated', 247.0))
   call check(nf90_put_att(ncid, nf90_global, 'azimuth_rotated', 0.0))
   call check(nf90_put_att(ncid, varid_glon,  'description', 'geographical longitude'))
   call check(nf90_put_att(ncid, varid_glon,  'units', 'degree_east'))
   call check(nf90_put_att(ncid, varid_glat,  'description', 'geographical latitude'))
   call check(nf90_put_att(ncid, varid_glat,  'units', 'degree_north'))
   call check(nf90_put_att(ncid, varid_rlon,  'description', 'rotated longitude'))
   call check(nf90_put_att(ncid, varid_rlon,  'units', 'degree_east'))
   call check(nf90_put_att(ncid, varid_rlat,  'description', 'rotated latitude'))
   call check(nf90_put_att(ncid, varid_rlat,  'units', 'degree_north'))
   call check(nf90_put_att(ncid, varid_orig,  'units', 'meters'))
   call check(nf90_put_att(ncid, varid_orig,  'description', 'Significant Wave height (firstguess)'))
   call check(nf90_put_att(ncid, varid_ipol,  'units', 'meters'))
   call check(nf90_put_att(ncid, varid_ipol,  'description', 'Significant Wave height (analysis)'))
   call check(nf90_put_att(ncid, varid_diff,  'units', 'meters'))
   call check(nf90_put_att(ncid, varid_diff,  'description', 'Significant Wave height (increment)'))

!- End definition of variables
   call check(nf90_enddef(ncid))

!- Write the data
   call check(nf90_put_var(ncid, varid_glon, real(glon_i,8)))
   call check(nf90_put_var(ncid, varid_glat, real(glat_i,8)))

   allocate(rotlon(imdl_i))
   allocate(rotlat(jmdl_i))
   do iii = 1, imdl_i
     rotlon(iii) = 299.0_dp + 0.025_dp*(iii-1)
     if (rotlon(iii) >= 360.0_dp ) rotlon(iii) = rotlon(iii) - 360.0_dp
     if (rotlon(iii) <    0.0_dp ) rotlon(iii) = rotlon(iii) + 360.0_dp
   end do
   do jjj = 1, jmdl_i
     rotlat(jjj) = -37.0_dp + 0.025_dp*(jjj-1)
   end do
   call check(nf90_put_var(ncid, varid_rlon, rotlon))
   call check(nf90_put_var(ncid, varid_rlat, rotlat))

   call check(nf90_put_var(ncid, varid_orig, howv_o1))
   call check(nf90_put_var(ncid, varid_ipol, howv_o2))
   call check(nf90_put_var(ncid, varid_diff, howv_diff))

!- Close the dataset
   call check(nf90_close(ncid))

!--- clean the memory
   deallocate(input_bitmap, input_data)

   deallocate(igdtmpl_i)
   deallocate(xpts_i, ypts_i)
   deallocate(glat_i, glon_i)
   deallocate(crot_i, srot_i)
   deallocate(xlon_i, xlat_i)
   deallocate(ylon_i, ylat_i)
   deallocate(area_i)
!  deallocate (ibi)

   deallocate(glat_o, glon_o)
   deallocate(howv_o1)
   deallocate(howv_o2)
   deallocate(howv_diff)
   deallocate(slmask_o)
!  deallocate(output_glat)
!  deallocate(output_glon)
!  deallocate(output_data)
!  deallocate(output_bitmap)
!  deallocate (ibo)

 enddo

!--- close grib file
 call baclose(iunit, iret)

 if (iret /= 0) stop 'Error: baclose failed.'

!--- free the memory usage by gfld
 call gf_free(gfld_input)

!================================================================================
      contains
!-----------------------------------------------------------------------
!
      subroutine check(status)
      integer,intent(in) :: status
!
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
      end subroutine check
!
!-----------------------------------------------------------------------
!

end program regrid_esg2rll_iplib
