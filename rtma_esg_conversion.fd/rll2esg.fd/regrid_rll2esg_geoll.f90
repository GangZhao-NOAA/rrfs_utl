program regrid_rll2esg_iplib
!================================================================================
 use omp_lib
 use netcdf
 use pkind, only:dp, sp, dpi, spi
 use ip_mod                               ! ip lib
 use grib_mod                             ! g2 lib(grib)
 use bacio_module                         ! prerequisite for g2 lib (grib)

 implicit none

 character(len=100)      :: input_file
 character(len=100)      :: input_file_nc
 character(len=100)      :: latlon_file_nc
 character(len=100)      :: outfile_chk
 integer                 :: ncid, varid
 character(len=20)       :: varname_nc

 real                    :: fill
 real                    :: diff, maxdiffx, maxdiffy
 integer                 :: nret, iopt
 character(len=4)        :: cgrid_i
 integer                 :: badpts
 integer                 :: ii, jj, nn
 integer                 :: n_bitmap

 integer                 :: iunit, iret, lugi
 integer                 :: iunito_chk

 integer                 :: ip, ipopt(20)     ! parameters for interpolation

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

 integer                 :: mo, no
 integer                 :: km
 integer                 :: ipts_o, jpts_o         ! dimension size read in netcdf data file
 integer                 :: npts_o

 integer                 :: igdtlen_o
 integer                 :: igdtnum_o
 integer, allocatable    :: igdtmpl_o(:)

 integer, allocatable    :: ibo(:)
!integer                 :: ibo
!real, allocatable       :: output_data(:,:)           ! 2D Data in 1D slice
 real, allocatable       :: output_data(:)             ! 2D Data in 1D slice
!logical*1, allocatable  :: output_bitmap(:,:)         ! 2D Data in 1D slice
 logical*1, allocatable  :: output_bitmap(:)           ! 2D Data in 1D slice
!real, allocatable       :: output_glat(:,:)           ! 2D Data in 1D slice
 real, allocatable       :: output_glat(:)             ! 2D Data in 1D slice
!real, allocatable       :: output_glon(:,:)           ! 2D Data in 1D slice
 real, allocatable       :: output_glon(:)             ! 2D Data in 1D slice

 real, allocatable       :: glon_o(:,:), glat_o(:,:)
 real, allocatable       :: slmask_o(:,:)
 real, allocatable       :: data_o(:,:)

 integer                 :: xdimid, ydimid, timedimid
 integer                 :: howv_varid

!integer                 :: n_center

 logical           :: l_chk_gdswzd     ! if call gdswzd with iopt=-1 to check
 logical           :: l_chk_bitmap     ! check bitmap with False and the value of data
 logical           :: l_chk_undefi     ! check undefined value

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

  l_chk_gdswzd = .true.
  l_chk_bitmap = .false.
  l_chk_undefi = .true.

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

   imdl_i = gfld_input%igdtmpl(8)
   jmdl_i = gfld_input%igdtmpl(9)
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
     write(6,*) 'There are bitmap data associated with the data.'
     ibi                   = 1        ! tell ipolates to use bitmap
     input_bitmap(:)       = gfld_input%bmap
   else                           ! no bitmap, data everywhere
     write(6,*) 'There is NO bitmap data associated with the data.'
     ibi                   = 0        ! tell ipolates there is no bitmap
     input_bitmap(:)       = .true.
   endif

   input_data(:)           = gfld_input%fld  ! the input data field

   write(6,*) 'reset input_bitmap to be true everywhere.'
   input_bitmap(:)       = .true.

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
   print*,'xpts/ypts POINT(1,1):  ',xpts_i(1,1),           ypts_i(1,1)
   print*,'xpts/yptsPOINT(1,JM):  ',xpts_i(1,jmdl_i),      ypts_i(1,jmdl_i)
   print*,'xpts/yptsPOINT(IM,1):  ',xpts_i(imdl_i,1),      ypts_i(imdl_i,1)
   print*,'xpts/yptsPOINT(IM,JM): ',xpts_i(imdl_i,jmdl_i), ypts_i(imdl_i,jmdl_i)

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
   ipts_o = 3950
   jpts_o = 2700
   npts_o = ipts_o * jpts_o
   write(6,*) 'output grid dimensions -- ipts_o:', ipts_o, ' jpts_o:', jpts_o
   allocate(glat_o(ipts_o, jpts_o))
   allocate(glon_o(ipts_o, jpts_o))
   allocate(data_o(ipts_o, jpts_o))
   allocate(slmask_o(ipts_o, jpts_o))

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
   call check( nf90_open(input_file_nc, nf90_write, ncid) )

   varname_nc = "slmsk"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, slmask_o ))
   write(6,*) 'read in ', trim(adjustl(varname_nc)),'=',slmask_o(1,1),slmask_o(1,jpts_o),slmask_o(ipts_o,jpts_o),slmask_o(ipts_o,1)

!---------------------------------------------------------------------------
!  the output "grid" is ESG grid, treated as a series of random station points.
!  in this case, set the grid definition template number of a negative number.
!  the grid definition template array information is not used, so set to a flag value.
!---------------------------------------------------------------------------
   igdtnum_o = -1 
   igdtlen_o =  1
   allocate(igdtmpl_o(igdtlen_o))
   igdtmpl_o =  -9999

   mo = npts_o
   no = mo

   allocate (ibo(km))
   allocate (output_glat(mo))
   allocate (output_glon(mo))
   allocate (output_data(mo))
   allocate (output_bitmap(mo))

   nn = 0
   do jj = 1, jpts_o
   do ii = 1, ipts_o
     nn = nn + 1
     output_glat(nn) = glat_o(ii,jj)
     output_glon(nn) = glon_o(ii,jj)
   end do
   end do

!---------------------------------------------------------------------------
!  setup arguments for ipolates (scalar interpolation) call.
!---------------------------------------------------------------------------

   ip       = 0                         ! bilinear interpolation
   ipopt(:) = 0                         ! options for bilinear
   ipopt(1) = 75                        ! set minimum mask to 75%

!---------------------------------------------------------------------------
!  call ipolates to interpolate scalar data.  non-zero "iret" indicates
!  a problem.
!---------------------------------------------------------------------------

   write(6,*) ' call ipolates_grib2 ... '
!  call ipolates_grib2_single_field(ip, ipopt, igdtnum_i, igdtmpl_i, igdtlen_i, &
!  call ipolates(ip, ipopt, igdtnum_i, igdtmpl_i, igdtlen_i, &
   call ipolates_grib2(ip, ipopt, igdtnum_i, igdtmpl_i, igdtlen_i, &
                 igdtnum_o, igdtmpl_o, igdtlen_o, &
                 mi, mo, km, ibi, input_bitmap, input_data, no, output_glat, &
                 output_glon, ibo, output_bitmap, output_data, iret)

   if (iret /= 0) then
     write(6,*) ' ipolates_grib2 failed.'
     stop
   else
     write(6,*) ' ipolates_grib2 is done successfully.'
   end if

   nn = 0
   do jj = 1, jpts_o
   do ii = 1, ipts_o
       nn = nn + 1
       data_o(ii,jj) = output_data(nn)
!      if ( .not. output_bitmap(nn) .or. slmask_o(ii,jj) > 0.1 ) then
!        data_o(ii,jj) = -0.01
!        data_o(ii,jj) =  0.0
!      end if
   end do
   end do

!--- inquire dimension id
   call check( nf90_inq_dimid(ncid, "xaxis_1",    xdimid) )   
   call check( nf90_inq_dimid(ncid, "yaxis_1",    ydimid) )   
   call check( nf90_inq_dimid(ncid, "Time",    timedimid) )

!--- define new variable
   call check( nf90_redef(ncid) )
!  call check( nf90_def_var(ncid, "howv", nf90_double, &    ! rrfs sfc_data.nc changes data type to NF90_FLOAT sometime in 2023
   call check( nf90_def_var(ncid, "howv", nf90_float,  &
                      (/ xdimid, ydimid, timedimid /), howv_varid )  )
   call check( nf90_enddef(ncid) )

!--- output data to new variable
   call check( nf90_put_var(ncid, howv_varid, data_o))
   write(6,*) ' create new variable howv and append it to input netcdf file '

   if (l_chk_bitmap .or. l_chk_undefi) then

     nn = 0
     n_bitmap = 0
     do jj = 1, jpts_o
     do ii = 1, ipts_o

       nn = nn + 1
       data_o(ii,jj) = output_data(nn)

       if ( .not. output_bitmap(nn) .and. l_chk_bitmap ) then
!      if (       output_bitmap(nn) .and. l_chk_bitmap ) then
         n_bitmap = n_bitmap + 1
         if ( mod(n_bitmap, 1000) == 0 ) then
           write(6,'(1x,A,I,L,F,F,F,I)') '(output) nn  bitmap slmask  data: ', nn, output_bitmap(nn), slmask_o(ii,jj), output_data(nn), data_o(ii,jj), n_bitmap
         end if
       end if

       if ( abs(output_data(nn)) .ge. 50.0 .and. l_chk_undefi ) then
           write(6,*) '(output) nn  bitmap  abs(data)>=50.0): ', nn, output_bitmap(nn), output_data(nn)
       end if

     end do
     end do

   end if

   call check( nf90_close(ncid) )

!--- clean the memory
   deallocate(input_bitmap, input_data)

   deallocate(igdtmpl_i)
   deallocate(xpts_i, ypts_i)
   deallocate(glat_i, glon_i)
   deallocate(crot_i, srot_i)
   deallocate(xlon_i, xlat_i)
   deallocate(ylon_i, ylat_i)
   deallocate(area_i)
   deallocate (ibi)

   deallocate(glat_o, glon_o)
   deallocate(output_glat)
   deallocate(output_glon)
   deallocate(output_data)
   deallocate(output_bitmap)
   deallocate (ibo)

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

end program regrid_rll2esg_iplib
