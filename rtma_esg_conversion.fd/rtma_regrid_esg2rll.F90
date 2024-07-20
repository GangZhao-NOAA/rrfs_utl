program rtma_regrid_esg2rll
!================================================================================
  use omp_lib
  use netcdf                               ! netcdf lib
  use pkind, only: dp, sp, dpi, spi        ! Jim Purser's lib
  use pietc, only: dtor,rtod               ! Jim Purser's lib
  use gdswzd_mod, only: gdswzd             ! ip lib (for conversion from earth to grid coor or vice versa)
!-- required when using iplib v4.0 or higher
#ifndef IP_V3
! use ip_mod                               ! ip lib (for interpolation)
! use ipolates_mod                         ! ip lib (for interpolation)
! use ip_grid_descriptor_mod
! use ip_grids_mod
! use ip_grid_mod
! use ip_grid_factory_mod
#endif
  use grib_mod                             ! g2 lib(grib)
  use bacio_module                         ! prerequisite for g2 lib (grib)
  use pbswi, only: abswi2
! use pbswi, only: abswi4, abswi6, abswi8
  use mod_rtma_regrid, only: rotated_gridopts, variable_options, esg_gridopts
  use mod_rtma_regrid, only: set_esg_gridopts, set_variable_options,                  &
                             check_varopts_grb2, set_time4data,                       &
                             check_grbmsg, set_rllgridopts,                           &
                             set_bitmap_grb2, check_data_1d_with_bitmap,              &
                             ll_to_xy_esg
#ifdef IP_V3
  use mod_rtma_regrid, only: gdt2gds_rll
#endif

  implicit none

  real(dp),     parameter :: undef_real = -9999.00_dp
  integer,      parameter :: undef_int  = -9999

  type(gribfield)         :: gfld_input
  type(rotated_gridopts)  :: rll_opts
  type(variable_options)  :: var_opts
  type(esg_gridopts)      :: esg_opts

  character(len=100)      :: input_data_rll_file_grb2
  character(len=100)      :: input_data_esg_file_nc
  character(len=100)      :: input_data_esg_file_fgs_nc
  character(len=100)      :: esg_grid_spec_file_nc
  character(len=100)      :: output_data_rll_file_nc
  character(len=20)       :: varname_nc
  character(len=20)       :: varname_input
 
  integer                 :: ii, jj
  integer                 :: iii, jjj
! integer                 :: nn, nnn
 
  integer                 :: iunit, iret, lugi


  integer                 :: adate(5)               ! year/month/day/hour/minute
  character(len=12)       :: cdate                  !yyyymmddhhmn

! integer                 :: ip, ipopt(20)     ! parameters for interpolation
 
  integer                 :: j, jdisc, jpdtn, jgdtn, k
  integer                 :: jids(200), jgdt(200), jpdt(200)
  integer                 :: km

  integer                 :: igdtlen_o
  integer                 :: igdtnum_o
  integer, allocatable    :: igdtmpl_o(:)
! integer                 :: mo
  integer                 :: imdl_o, jmdl_o         ! dimension size read in grib2 data file
  integer                 :: npts_o
  logical                 :: unpack
  integer, allocatable    :: ibo(:)
  logical*1, allocatable  :: output_bitmap(:)        ! 2D Data in 1D array
  real(dp), allocatable   :: output_data(:)          ! 2D Data in 1D array
! real(dp), allocatable   :: output_glat(:)          ! 2D Data in 1D array
! real(dp), allocatable   :: output_glon(:)          ! 2D Data in 1D array

  real(dp)                :: fill
  integer                 :: nret, iopt

#ifdef IP_V3
  integer                 :: igdt_grb2(5) 
  integer                 :: idefnum                 ! The number of entries in array ideflist, i.e. number
                                                     ! of rows (or columns) for which optional grid points are defined.
  integer,  allocatable   :: ideflist(:)             ! integer array containing the number of grid points contained in
                                                     ! each row (or column). To handle the irregular grid stuff.
  integer                 :: kgds_grb1_rll_o(200)
  integer                 :: igrid_grb1
  real(dp), allocatable   :: glat1d_o(:), glon1d_o(:)
  real(dp), allocatable   :: xpts1d_o(:), ypts1d_o(:)
#else
! type(grib2_descriptor)  :: desc_grb2
! class(ip_grid), allocatable :: ip_grid_rll_o
#endif
  real(dp), allocatable   :: glat_o(:,:), glon_o(:,:)
  real(dp), allocatable   :: xpts_o(:,:), ypts_o(:,:)

! real(dp), allocatable   :: slmask_o(:,:)
  real(dp), allocatable   :: data_o(:,:)
  real(dp), allocatable   :: data_fgs_o(:,:)
  real(dp), allocatable   :: data_tmp_o(:,:)
  real(dp), allocatable   :: diff_xy(:,:)
  logical*1, allocatable  :: output_bitmap_2d(:,:)
  real(dp), allocatable   :: rotlon(:), rotlat(:)

! integer                 :: mi, ni
  integer                 :: ipts_i, jpts_i         ! dimension size read in netcdf data file
  integer                 :: ipt2_i, jpt2_i
  integer                 :: npts_i
  real(dp), allocatable   :: xpts_i(:,:), ypts_i(:,:)
  real(dp), allocatable   :: glat_i(:,:), glon_i(:,:)
  real(dp), allocatable   :: slmask_i(:,:)
  real(dp), allocatable   :: data_i(:,:)
  real(dp), allocatable   :: data_fgs_i(:,:)

  integer                 :: ncid, varid
  integer                 :: xdimid, ydimid, timedimid                        ! for read in netcdf
  integer                 :: xt_dimid, yt_dimid
  integer                 :: dimid_x, dimid_y, dimid_time                     ! for write out netcdf
  integer                 :: varid_data
  integer                 :: varid_rlon, varid_rlat
  integer                 :: varid_glon, varid_glat

!---
  logical*1      :: l_clean_bitmap    ! if true, then set bitmap = true everywhere
  logical*1      :: verbose
  logical*1      :: l_increment_intrp ! true : interpolation with increment (when regrdding from esg to rll)
                                               ! false: interpolation with full variable
  integer        :: interp_opt        ! 2: BSWI-2 interpolation scheme (only available scheme for now)

  logical        :: ff

  character(100) :: fname_nml
  logical        :: f_exist 
  integer        :: lunin_nml

  namelist/setup/varname_input, verbose, l_clean_bitmap, l_increment_intrp, interp_opt

!-----------------------------------------------------------------------
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
! 0. reading the namelist
!---------------------------------------------------------------------------
!--- initialising the namelist variables first
  l_clean_bitmap    = .false.    ! do not reset bitmap to be true in whole domain (default)
  verbose           = .false.
  varname_input     = ''
  l_increment_intrp = .false. ! regridding with full variable (default)
  interp_opt        = 2

  f_exist = .false.
  lunin_nml = 10
!-- reading namelist
  fname_nml = 'esg2rll_namelist'
  inquire(file=trim(adjustl(fname_nml)),exist=f_exist)
  if (f_exist) then
    write(6,*) 'reading from namelist: ', trim(adjustl(fname_nml))
    open(lunin_nml, file=trim(adjustl(fname_nml)), form='formatted', status='old')
    read(lunin_nml, setup)
    close(lunin_nml)
    write(6,*) "checking the setup info in namelist:"
    write(6,setup)
  else
    write(6,'(1x,2A)')  &
         "Abort...,  failed to find the required namelist file: ", trim(adjustl(fname_nml))
    stop(99)
  end if
  if (trim(adjustl(varname_input)) == '') then
    write(6,*) "Abort...,  must set varname_input in namelist (e.g., varname_input='howv')" 
    stop(1)
  else
    var_opts%varname = trim(adjustl(varname_input))
    call set_variable_options(var_opts, iret)
    if ( iret /= 0 ) then
      write(6,'(1x,3A)') 'This program cannot process this variable : ',       &
            trim(adjustl(var_opts%varname)), ', task is ABORTED !!!!'
      stop(2)
    end if
  end if

!---------------------------------------------------------------------------
! 1. Read the output grid (RLL) info from grib2 file (on rotated grid)
!---------------------------------------------------------------------------
!   1.1 opening the grib 2 file containing data to be interpolated.
!       Note: for this example, there are only one data record 
!             (HTSGW or GUST) in grib2 file.
!---------------------------------------------------------------------------
 iunit=9
 input_data_rll_file_grb2="./input_data_rll.grib2"
 call baopenr(iunit, input_data_rll_file_grb2, iret)
 if (iret /= 0) then
   write(6,*) 'return from baopenr: ',iret
   stop 'Error: baopenr failed.'
 end if

!---------------------------------------------------------------------------
!   1.2 preparing for call to g2 library to degrib data. 
!       Note: the data are assumed to be on a rotated-lat/lon grid 
!             with i/j dimension of 360/181. 
!---------------------------------------------------------------------------
 jdisc   = -1         ! search for any discipline
 jpdtn   = -1         ! search for any product definition template number
 jgdtn   = -1         ! search for grid definition template number
                      ! 0 - regular lat/lon grid is expected.
                      ! 1 - rotated lat/lon grid is expected.
 jids    = -9999      ! array of values in identification section, set to wildcard
 jgdt    = -9999      ! array of values in grid definition template 3.m
 jpdt    = -9999      ! array of values in product definition template 4.n
 unpack  = .true.     ! unpack data
 lugi    = 0          ! no index file (if using index file, set lugi = iunit)

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
!   1.3 degrib the data.  non-zero "iret" indicates a problem during degrib.
!---------------------------------------------------------------------------
 km = 1              ! number of records to interpolate (in this example, only one record)

 do j = 0, (km-1)

   call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld_input, iret)

   if (iret /= 0) then
     write(6,'(1x,A,I4,A)') 'return from sub getgb2: ', iret, ' --> Error: getb2 failed.'
     stop(1)
   end if
   if ( verbose ) call check_grbmsg(gfld_input)

!--- check up if the variable in grib2 file matches the requried variable
   call check_varopts_grb2(gfld_input,var_opts,iret)
   if ( iret /= 0 ) then
     write(6,'(1x,3A)') 'Warning: Required variable name [',                   &
           trim(adjustl(var_opts%varname)),                                    &
           '] does not match the variable in grib2 file. Task is ABORTED ... '
     stop(3)
   end if

!--- date/time info of input date
   call set_time4data(gfld_input, adate, cdate)

!--- input grid template information
!
   imdl_o = gfld_input%igdtmpl(8)  ! x-dimension of output grid
   jmdl_o = gfld_input%igdtmpl(9)  ! y-dimension of output grid
   npts_o = imdl_o * jmdl_o        ! total dimension of output grid
   write(6,'(1x,A,3(1x,I8))') &
        'dimension of output grib2 grid (Nx, Ny, Nx*Ny) = ', imdl_o, jmdl_o, npts_o

!  jgdtn was set to be -1 so getgb2 would search for grid template number.
!  So after getgb2, jgdtn is still -1, but the true grid template number 
!  is gfld_input%igdtnum.
   igdtnum_o = gfld_input%igdtnum
   if (igdtnum_o == 1) then
     write(6,'(1x,A)') "output model grid is defined on rotated lat-lon grid, as expected."
     call set_rllgridopts(gfld_input, rll_opts)
   else
     write(6,'(1x,A,I12.12)') &
       "output model grid is defined on the grid with grid template number = ",igdtnum_o
     write(6,'(1x,A)') ' However, a Rotated Lat-Lon Grid is expected. So Abort this running ...'
     stop(-1)
   end if

   igdtlen_o = gfld_input%igdtlen
   allocate(igdtmpl_o(igdtlen_o))
   igdtmpl_o(:) = gfld_input%igdtmpl(:)      ! grid template (used by grib2)

!---------------------------------------------------------------------------
!   1.4 checking up with output model data (and its bitmap, etc.) 
!     Note:
!       These output data might be used to fill the area where the regridded
!       data could not cover (e.g., area undefined for ESG grid)
!---------------------------------------------------------------------------
!--- output data field decoded in grib2 file
   allocate(output_data(npts_o))
   output_data(:)           = gfld_input%fld  ! the output data field

!---  does grib2 data have a bitmap?
   allocate(ibo(1))
   allocate(output_bitmap(npts_o))
   l_clean_bitmap = .False.                  ! do not re-set bitmap to be true everywhere
   call set_bitmap_grb2(gfld_input,npts_o,l_clean_bitmap,ibo,output_bitmap)

!---  checking if existing any abnormal data values and counting data with false bitmap
   write(6,*)'----------------------------------------------------------'
   write(6,*)'checking the output data read from grib2 fie:'
   call check_data_1d_with_bitmap(var_opts,npts_o,output_data,ibo,output_bitmap)

!---------------------------------------------------------------------------
!   1.5 calculate lat/lon of input grid points (rotated-latlon grid here)
!       note: 
!            to use gdswzd, this type of grid (given igdtnum, igdtmpl)
!            should be recognizable by gdswzd.
!---------------------------------------------------------------------------
   iopt = 0                         ! option used in gdswzd:
                                    ! 0: calculating grid(i/j) and earth coords (lat/lon) of all grid points
                                    ! 1: calculating earth coords (lat/lon) of selected grid coordinates
                                    !-1: calculating grid coordinates of selected earth coords (lat/lon)
#ifdef IP_V3
!-- when using iplib v3.0 or lower
   allocate (xpts1d_o(npts_o),ypts1d_o(npts_o))
   allocate (glat1d_o(npts_o),glon1d_o(npts_o))
   allocate (xpts_o(imdl_o,jmdl_o),ypts_o(imdl_o,jmdl_o))
   allocate (glat_o(imdl_o,jmdl_o),glon_o(imdl_o,jmdl_o))
   fill = -9999.0_dp
   xpts1d_o = fill ; ypts1d_o = fill ;  ! Grid x & y point coords
   glat1d_o = fill ; glon1d_o = fill ;  ! Earth lat & lon in degree
!-------------------------------------------------------------------------------------!
!    As required by IP lib v3.x or older, gdszwd needs an array gds(200)              !
!    to save the grid information and parameters, and that array gds(200)             !
!    follows GRIB1 GDS info. So need to convert grid informaton from                  !
!    GRIB2 Grid Description Section (and its Grid Definition Template)                !
!    to GRIB1 GDS info (similarly as decoded by w3fi63)                               !
!-------------------------------------------------------------------------------------!
   igdt_grb2(1) = 0         ! Source of Grid Definition (0: then specified in Code Table 3.1)
   igdt_grb2(2) = npts_o    ! Number of Data Points in the defined grid
   igdt_grb2(3) = 0         ! Number of octets needed for each additional grid definition (if 0: using regular grid)
   igdt_grb2(4) = 0         ! Interpetation of list of for optional points definition (Table 3.11)
   igdt_grb2(5) = igdtnum_o ! GrdDefTmplt number(Table 3.1): 1--> rotated latlon grid (Template 3.1)
   idefnum     = 0          ! no irregular grid stuff
   allocate(ideflist(1))
   ideflist(:) = 0
   call gdt2gds_rll(igdt_grb2, igdtmpl_o, idefnum, ideflist, kgds_grb1_rll_o, igrid_grb1, iret)

!  lrot = 0                             ! return Vector Rotations (if 1)
!  lmap = 0                             ! return Map Jacobians    (if 1)
   call gdswzd(kgds_grb1_rll_o, iopt, npts_o, fill,                             &
!              lrot, lmap,                                                      &  ! not sure if needed for ip lib v3.x
               xpts1d_o, ypts1d_o, glon1d_o, glat1d_o, nret) 
   if (nret /= npts_o) then
     write(6,'(1x,2(A,1x,I8))') &
           'ERROR: Checking --> NUMBER OF VALID POINTS RETURNED FROM GDSWZS (iopt=0) ',             &
           nret, ' DOES NOT MATCH TOTAL NUMBER of GRID POINTS ', npts_o
     stop(4)
   endif
!--- reshaping 1-D output array to 2-D array
   xpts_o = reshape(xpts1d_o, (/imdl_o, jmdl_o/), order=(/1,2/))   
   ypts_o = reshape(ypts1d_o, (/imdl_o, jmdl_o/), order=(/1,2/))   
   glat_o = reshape(glat1d_o, (/imdl_o, jmdl_o/), order=(/1,2/))   
   glon_o = reshape(glon1d_o, (/imdl_o, jmdl_o/), order=(/1,2/))   
   deallocate(xpts1d_o, ypts1d_o)
   deallocate(glat1d_o, glon1d_o)
#else
!-- when using iplib v4.0 or higher
   allocate (xpts_o(imdl_o,jmdl_o),ypts_o(imdl_o,jmdl_o))
   allocate (glat_o(imdl_o,jmdl_o),glon_o(imdl_o,jmdl_o))
!--- creating the grid info from grib2 template info in grib2 file
!  allocate (xpts1d_o(npts_o),ypts1d_o(npts_o))
!  allocate (glat1d_o(npts_o),glon1d_o(npts_o))
!  desc_grb2 = init_descriptor(igdtnum_o, igdtlen_o, igdtmpl_o)
!  call init_grid(ip_grid_rll_o, desc_grb2)
!  call gdswzd(ip_grid_rll_o, iopt, npts_o, fill, xpts1d_o, ypts1d_o,                 &
!              glon1d_o, glat1d_o, nret) 
   call gdswzd(igdtnum_o, igdtmpl_o, igdtlen_o, iopt, npts_o, fill, xpts_o, ypts_o,                 &
               glon_o, glat_o, nret) 
!              crot_o, srot_o, xlon_o, xlat_o, ylon_o, ylat_o, area_o)  !<-- optional arguments
   if (nret /= npts_o) then
     write(6,'(1x,2(A,1x,I8))') &
           'ERROR: Checking --> NUMBER OF VALID POINTS RETURNED FROM GDSWZS (iopt=0) ',             &
           nret, ' DOES NOT MATCH TOTAL NUMBER of GRID POINTS ', npts_o
     stop(4)
   endif
#endif
   if ( verbose ) then
     write(6,'(1x,A,4(1x,A1,F8.3,1x,F8.3,A1))')                                                     &
               'LAT/LON   at RLL domain corners (1,1), (1,JM), (IM,1) & (IM,JM): ',                 &
               '(',glat_o(1,1),glon_o(1,1),')',                                                     &
               '(',glat_o(1,jmdl_o),glon_o(1,jmdl_o), ')',                                          &
               '(',glat_o(imdl_o,1),glon_o(imdl_o,1), ')',                                          &
               '(',glat_o(imdl_o,jmdl_o),glon_o(imdl_o,jmdl_o),')'
     write(6,'(1x,A,4(1x,A1,F8.3,1x,F8.3,A1))')                                                     &
               'XPTS/YPTS at RLL domain corners (1,1), (1,JM), (IM,1) & (IM,JM): ',                 &
               '(',xpts_o(1,1),ypts_o(1,1),')',                                                     &
               '(',xpts_o(1,jmdl_o),ypts_o(1,jmdl_o), ')',                                          &
               '(',xpts_o(imdl_o,1),ypts_o(imdl_o,1), ')',                                          &
               '(',xpts_o(imdl_o,jmdl_o),ypts_o(imdl_o,jmdl_o),')'
   end if

!---------------------------------------------------------------------------
! 2. Read information of the output ESG grid in fv3_grid_specification file (netcdf)
!---------------------------------------------------------------------------
   write(6,'(1x,A)')'=================================================================='
   esg_grid_spec_file_nc = "./fv3_grid_spec_esg.nc"
   call check( nf90_open(esg_grid_spec_file_nc, nf90_nowrite, ncid) )
!--- inquire dimension id and the dimension size
   call check( nf90_inq_dimid(ncid, "grid_xt", xt_dimid) )             ! cell center
   call check( nf90_inquire_dimension(ncid, xt_dimid, len = ipts_i) )   
   call check( nf90_inq_dimid(ncid, "grid_yt", yt_dimid) )             ! cell center
   call check( nf90_inquire_dimension(ncid, yt_dimid, len = jpts_i) )   

   npts_i = ipts_i * jpts_i
   write(6,'(3(1x,A,I8))') 'intput ESG grid dimensions -- ipts_i= ', ipts_i,    &
        ' jpts_i= ', jpts_i, ' npts_i= ', npts_i

   allocate(glon_i(ipts_i, jpts_i))
   allocate(glat_i(ipts_i, jpts_i))

   varname_nc = "grid_lont"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, glon_i ))
   write(6,'(1x,3A,4(1x,F12.6))') 'read-in ', trim(adjustl(varname_nc)), ' at 4 corners(ll->lu->ur->lr) =', &
        glon_i(1,1),glon_i(1,jpts_i),glon_i(ipts_i,jpts_i),glon_i(ipts_i,1)

   varname_nc = "grid_latt"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, glat_i ))
   write(6,'(1x,3A,4(1x,F12.6))') 'read-in ', trim(adjustl(varname_nc)), ' at 4 corners(ll->lu->ur->lr) =', &
        glat_i(1,1),glat_i(1,jpts_i),glat_i(ipts_i,jpts_i),glat_i(ipts_i,1)

   call check( nf90_close(ncid) )

!---------------------------------------------------------------------------
! 3. Read input data on ESG grid (netcdf file)
!---------------------------------------------------------------------------
   input_data_esg_file_nc = "./input_data_esg.nc"
   call check( nf90_open(input_data_esg_file_nc, nf90_nowrite, ncid) )
!--- inquire dimension id of output data file
   call check( nf90_inq_dimid(ncid, "xaxis_1",    xdimid) )   
   call check( nf90_inquire_dimension(ncid, xdimid, len = ipt2_i) )   
   call check( nf90_inq_dimid(ncid, "yaxis_1",    ydimid) )   
   call check( nf90_inquire_dimension(ncid, ydimid, len = jpt2_i) )   
   call check( nf90_inq_dimid(ncid, "Time",    timedimid) )
!--- check if dimensions of input data file (netcdf) match the dimensions of fv3_grid_spec file
   if ( ipt2_i /= ipts_i .or. jpt2_i /= jpts_i ) then
     write(6,'(1x,3A)') 'WARNING --> dimensions of input data file ', input_data_esg_file_nc, &
          ' do NOT match dimensions of fv3_grid_spec file. <-- Warning'
     write(6,'(1x,A,2(1x,I8))') 'dimension of input data    file : ', ipt2_i, jpt2_i
     write(6,'(1x,A,2(1x,I8))') 'dimension of fv3_grid_spec file : ', ipts_i, jpts_i
     write(6,*) ' Check the dimenesions above. Now ABORT the task ...'
     stop(5)
   end if
!--- read sea-land mask in input data file (ESG grid)
   allocate(slmask_i(ipts_i, jpts_i))
   varname_nc = "slmsk"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, slmask_i ))
   write(6,'(1x,3A,5(1x,F12.6))') 'read-in ', trim(adjustl(varname_nc)),       &
        ' at 4 corners & center =',slmask_i(1,1),slmask_i(1,jpts_i),           &
        slmask_i(ipts_i,jpts_i),slmask_i(ipts_i,1),slmask_i(ipts_i/2,jpts_i/2)

!--- read the interested data on ESG grid
   allocate(data_i(ipts_i, jpts_i))
   varname_nc = var_opts%var_ncf
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, data_i ))
   write(6,'(1x,3A,5(1x,F12.6))') 'read-in full data ', trim(adjustl(varname_nc)), &
        ' (on ESG grid at 4 corners and center)= ', data_i(1,1),data_i(1,jpts_i),  &
        data_i(ipts_i,jpts_i),data_i(ipts_i,1),data_i(ipts_i/2,jpts_i/2)

   call check( nf90_close(ncid) )
   write(6,'(1X,A,3(1X,F15.6))') "max/min/mean full variable data on ESG input grid : ",    &
         maxval(data_i), minval(data_i), sum(data_i)/npts_i

!---- reading the firstguess on ESG grid if regridding the increments, not the full variable
   if ( l_increment_intrp ) then   ! if l_increment_intrp is true, then read the firstguess data

!--- read the firstguess data on ESG grid
     input_data_esg_file_fgs_nc = "./input_data_esg_fgs.nc"
     call check( nf90_open(input_data_esg_file_fgs_nc, nf90_nowrite, ncid) )
!--- inquire dimension id of output data file
     call check( nf90_inq_dimid(ncid, "xaxis_1",    xdimid) )   
     call check( nf90_inquire_dimension(ncid, xdimid, len = ipt2_i) )   
     call check( nf90_inq_dimid(ncid, "yaxis_1",    ydimid) )   
     call check( nf90_inquire_dimension(ncid, ydimid, len = jpt2_i) )   
     call check( nf90_inq_dimid(ncid, "Time",    timedimid) )
!--- check if dimensions of firstguess data file (netcdf) match the dimensions of fv3_grid_spec file
     if ( ipt2_i /= ipts_i .or. jpt2_i /= jpts_i ) then
       write(6,'(1x,3A)') &
             'Dimensions of firstguess file do NOT match dimensions of fv3_grid_spec file.', &
             'dimension of input data    file : ', ipt2_i, jpt2_i,                           &
             'dimension of fv3_grid_spec file : ', ipts_i, jpts_i, '   Abort ...'
       stop(5)
     end if
!--- read the firstguess data on ESG grid
     allocate(data_fgs_i(ipts_i, jpts_i))
     varname_nc = var_opts%var_ncf
     call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
     call check( nf90_get_var(ncid, varid, data_fgs_i ))
     write(6,'(1x,3A,5(1x,F12.6))') 'read-in full firstguess data ', trim(adjustl(varname_nc)), &
          ' (on ESG grid at 4 corners and center)= ', data_fgs_i(1,1),data_fgs_i(1,jpts_i),     &
          data_fgs_i(ipts_i,jpts_i),data_fgs_i(ipts_i,1),data_fgs_i(ipts_i/2,jpts_i/2)
     call check( nf90_close(ncid) )
     write(6,'(1X,A,3(1X,F15.6))') "max/min full firstguess data on ESG input grid : ", &
           maxval(data_fgs_i), minval(data_fgs_i), sum(data_fgs_i)/npts_i

!---- calculate the increments -->  subtracting firstguess from analysis
     data_i(:,:) = data_i(:,:) - data_fgs_i(:,:)
     deallocate(data_fgs_i)
     write(6,'(1X,A,2(1X,F15.6))') "max/min increments on input ESG grid : ",         &
           maxval(data_i), minval(data_i)

   end if ! if ( l_increment_intrp )

!---------------------------------------------------------------------------
! 3. Transfer input data from ESG grid to Rotated Grid (RLL) 
!---------------------------------------------------------------------------
!   Note:
!        the input "grid" is on ESG grid, is approximately treated as 
!        structured grid, and set up x-y ESG coordinates, then use 
!        2-D interpolation. 
!---------------------------------------------------------------------------
!   3.1 set up the ESG grid parameters
!---------------------------------------------------------------------------
   call set_esg_gridopts(esg_opts)
   if ( verbose) then
     write(6,'(1X,A,7(1x,D))')'Pre-defined ESG parameters: A,Kappa,delx,dely,plat,plon,pazi=', &
           esg_opts%A, esg_opts%Kappa, esg_opts%delx, esg_opts%dely, esg_opts%plat,            &
           esg_opts%plon, esg_opts%pazi
   end if

!---------------------------------------------------------------------------
!   3.2.1 convert lat/lon of ESG grid points (input grid) 
!         to its x-y ESG coordinates
!---------------------------------------------------------------------------
   allocate (xpts_i(ipts_i,jpts_i),ypts_i(ipts_i,jpts_i))
   xpts_i = fill; ypts_i = fill;
   call ll_to_xy_esg(ipts_i, jpts_i, ipts_i, jpts_i, esg_opts, glat_i, glon_i, xpts_i, ypts_i)
   if ( verbose ) then
     write(6,'(1x,A,2(1x,I8))') 'Dimensions of ESG grid: ', ipts_i, jpts_i
     write(6,'(1x,A,4(/,A,2(1x,F12.4)))') &
           'X/Y of ESG grid corner points in ESG x/y coords: ',                 &
          'lower-left :',xpts_i(1,1),           ypts_i(1,1),                    &
          'upper-left :',xpts_i(1,jpts_i),      ypts_i(1,jpts_i),               &
          'lower-right:',xpts_i(ipts_i,1),      ypts_i(ipts_i,1),               &
          'upper-right:',xpts_i(ipts_i,jpts_i), ypts_i(ipts_i,jpts_i)
!--- to check if there is differences between the grid indices and 
!       the computed x/y coordinates for ESG grid itseld
     allocate(diff_xy(ipts_i,jpts_i))
     do jj = 1, jpts_i
       do ii = 1, ipts_i
         diff_xy(ii,jj) = sqrt( (xpts_i(ii,jj)-float(ii))**2 + (ypts_i(ii,jj)-float(jj))**2 )
       end do
     end do
     write(6,'(1x,A,2(1x,F8.4))') &
          'Max/Min differences of X-Y ESG coords (computed vs. index)= ',       &
          maxval(diff_xy),minval(diff_xy)
     if ( maxval(diff_xy) > 0.01_dp ) then
       write(6,'(1x,A)') &
         'The computed X/Y coordinates are very different to its grid indices. Stop. Please Check ...'
       stop(10)
     end if
     deallocate(diff_xy)
   end if
!---------------------------------------------------------------------------
!   3.2.2 convert lat/lon of rotated-latlon (RLL) grid points (output grid)
!         to x-y ESG coordinates
!---------------------------------------------------------------------------
   if (allocated(xpts_o)) deallocate(xpts_o)
   if (allocated(ypts_o)) deallocate(ypts_o)
   allocate (xpts_o(imdl_o,jmdl_o),ypts_o(imdl_o,jmdl_o))
   xpts_o = fill; ypts_o = fill;
   call ll_to_xy_esg(ipts_i, jpts_i, imdl_o, jmdl_o, esg_opts, glat_o, glon_o, xpts_o, ypts_o)
   if ( verbose ) then
     write(6,'(1x,2(A,2(1x,I8)))') 'Dimensions of RLL grid: ', imdl_o, jmdl_o
     write(6,'(1x,2(A,2(1x,I8)))') 'Dimensions of ESG grid: ', ipts_i, jpts_i
     write(6,'(1x,A,2(1x,F13.3))') 'max/min X of RLL in ESG x-y coordinates: ', &
                                    maxval(xpts_o),minval(xpts_o)
     write(6,'(1x,A,2(1x,F13.3))') 'max/min Y of RLL in ESG x-y coordinates: ', &
                                    maxval(ypts_o),minval(ypts_o)
     write(6,'(1x,A,4(/,A,2(1x,F12.4)))') &
           'X/Y of RLL grid corner points in ESG x/y coords: : ',               &
          'lower-left :',xpts_o(1,1),           ypts_o(1,1),                    &
          'upper-left :',xpts_o(1,jmdl_o),      ypts_o(1,jmdl_o),               &
          'lower-right:',xpts_o(imdl_o,1),      ypts_o(imdl_o,1),               &
          'upper-right:',xpts_o(imdl_o,jmdl_o), ypts_o(imdl_o,jmdl_o)
   end if
!---------------------------------------------------------------------------
!   3.3 interpolation from ESG to RLL (in ESG X-Y coordinates)
!---------------------------------------------------------------------------
   allocate(data_o(imdl_o, jmdl_o))
   allocate(data_tmp_o(imdl_o, jmdl_o))
   allocate(data_fgs_o(imdl_o, jmdl_o))
   allocate(output_bitmap_2d(imdl_o, jmdl_o))
!--- reshaping 1-D output array to 2-D array
   data_fgs_o = reshape(output_data, (/imdl_o, jmdl_o/), order=(/1,2/))   
   output_bitmap_2d = reshape(output_bitmap, (/imdl_o, jmdl_o/), order=(/1,2/))   
   data_tmp_o = 0.0_dp
   data_o = data_fgs_o  ! intialization by filling with original/firstguess data on RLL grid from grib2 file 

   if(interp_opt.eq.2)then

     print*,'Using BSWI-2, interp_opt= ', interp_opt

!    nn = 0
     do jj=1,jmdl_o
       do ii=1,imdl_o

!        nn = nn + 1
!        data_fgs_o(ii,jj) = output_data(nn)
!        output_bitmap_2d(ii,jj) = output_bitmap(nn)

         iii = dint(xpts_o(ii,jj))
         jjj = dint(ypts_o(ii,jj))

         call abswi2(1,ipts_i,1,jpts_i,data_i,xpts_o(ii,jj),ypts_o(ii,jj),data_tmp_o(ii,jj),ff)

!---     only use regridded data at that grid point if it is inside ESG grid domain
         if(iii .ge. 1 .and. jjj .ge. 1 .and. iii .le. (ipts_i - 1) .and. jjj .le. (jpts_i - 1) )then
           if ( l_increment_intrp ) then
             data_o(ii,jj) = data_fgs_o(ii,jj) + data_tmp_o(ii,jj)    ! regridding with increment data_i
           else
             data_o(ii,jj) = data_tmp_o(ii,jj)                        ! regridding with full variable data_i
           end if
         else
           data_o(ii,jj) = data_fgs_o(ii,jj)   ! if outside ESG grid domain, using orig/fgs value in grib2 file
         end if
!---     calculate the regridded analysis increment on output (RLL) grid
         data_tmp_o(ii,jj) = data_o(ii,jj) - data_fgs_o(ii,jj)

       enddo
     enddo
!--- check the regridded data on the output grid (RLL)
     write(6,'(1X,A,2(1X,F15.6))') "max/min regridded increments      on output RLL grid : ", &
           maxval(data_tmp_o), minval(data_tmp_o)
     write(6,'(1X,A,2(1X,F15.6))') "max/min orig/firstguess full data on output RLL grid : ", &
           maxval(data_fgs_o), minval(data_fgs_o)
     write(6,'(1X,A,2(1X,F15.6))') "max/min regridded full data       on output RLL grid : ", &
           maxval(data_o), minval(data_o)

   else

     write(6,'(1x,A,I4,A)') "Unknown Interp_opt=", interp_opt,                  &
           " Current code only accept interp_opt = 2. Task is abnormally terminated!"
     stop(7)

   endif

!  check the regridded data on the output grid (RLL)
   write(6,'(1x,A)')'=================================================================================='
   write(6,'(1X,A,2(1X,F15.6))') "max/min orig/firstguess full data on output RLL grid : ", &
         maxval(data_fgs_o), minval(data_fgs_o)
   write(6,*) "  ====> before set contraint to regridded data <===  "
   write(6,'(1X,A,2(1X,F15.6))') "max/min regridded increments on output RLL gridi: ", &
         maxval(data_tmp_o), minval(data_tmp_o)
   write(6,'(1X,A,2(1X,F15.6))') "max/min regridded full data on output RLL  grid : ", &
         maxval(data_o), minval(data_o)

!--- applying the specific contraints to the regridded variables
!---   non-negative feature
   if ( var_opts%varname == "howv" ) then
     where(data_o   .lt. 0.0  ) data_o = 0.0_dp       ! wave height >=0
   else if ( var_opts%varname == "gust" ) then
     where(data_o   .lt. 0.0  ) data_o = 0.0_dp       ! wind gust >=0
   end if
   data_tmp_o = data_o - data_fgs_o                   ! re-compute the analysis increments
!---  using bitmap (from original grib2 data on output RLL grid) to mask out invalid grid point
!       Note:
!            for significant wave height (howv), this action masks out data over land;
!            for 10-m wind gust(gust), this action masks out data outside ESG domain.
!  where( .not. output_bitmap_2d ) data_o = data_fgs_o
   where( .not. output_bitmap_2d )
     data_o     = undef_real
     data_tmp_o = undef_real
   end where
   write(6,*) "  ====> after set contraint to regridded data <===  "
   write(6,'(1X,A,2(1X,F15.6))') "max/min regridded increments on output RLL grid : ", &
         maxval(data_tmp_o, MASK=data_tmp_o .ne. undef_real),                          &
         minval(data_tmp_o, MASK=data_tmp_o .ne. undef_real)
   write(6,'(1X,A,2(1X,F15.6))') "max/min regridded full data on output RLL  grid : ", &
         maxval(data_o, MASK=data_o .ne. undef_real),                                  &
         minval(data_o, MASK=data_o .ne. undef_real)
   write(6,'(1x,A)')'=================================================================================='

!---------------------------------------------------------------------------
! 4. Output the regrided data on RLL grid to a new file (netcdf)
!---------------------------------------------------------------------------
!--- Create the netcdf file
   output_data_rll_file_nc = "./output_data_rll.nc"
   call check(nf90_create(trim(output_data_rll_file_nc), nf90_netcdf4, ncid))   

!  call check( nf90_redef(ncid) )
!- Define the dimensions
   call check(nf90_def_dim(ncid,    "X",         imdl_o,    dimid_x))
   call check(nf90_def_dim(ncid,    "Y",         jmdl_o,    dimid_y))
   call check(nf90_def_dim(ncid, "Time", nf90_unlimited, dimid_time))

!- Define rotated lat/lon variables
   call check(nf90_def_var(ncid, "rotlon", nf90_double, (/dimid_x /),         varid_rlon))
   call check(nf90_def_var(ncid, "rotlat", nf90_double, (/dimid_y/),          varid_rlat))
!- Define true earth lat/lon variables
   call check(nf90_def_var(ncid, "geolon", nf90_double, (/dimid_x, dimid_y/), varid_glon))
   call check(nf90_def_var(ncid, "geolat", nf90_double, (/dimid_x, dimid_y/), varid_glat))

!- Define data variables
   varname_nc=trim(adjustl(var_opts%var_ncf))
   call check(nf90_def_var(ncid, trim(adjustl(varname_nc)), nf90_double,        &
              (/dimid_x, dimid_y, dimid_time/), varid_data,                     &
              contiguous=.false.,                                               &
              chunksizes=(/imdl_o, jmdl_o, 1/),                                 &
              shuffle = .true., fletcher32 = .true.,                            &
              endianness = nf90_endian_little) ) 
   call check( nf90_def_var_fill(ncid, varid_data, 0, undef_real) )  ! set FillValue

!- Add the attributes
   call check(nf90_put_att(ncid, nf90_global, 'description', 'RRFS-3DRTMA 2-D Field on NA-3km domain'))
   call check(nf90_put_att(ncid, nf90_global, 'note', 'Rotated Lat/Lon Grid'))
   call check(nf90_put_att(ncid, nf90_global, 'Projection', 'Rotated Lat/Lon Grid'))
   call check(nf90_put_att(ncid, nf90_global, 'lon_ll_rll', rll_opts%llcnr(1)))
   call check(nf90_put_att(ncid, nf90_global, 'lat_ll_rll', rll_opts%llcnr(2)))
   call check(nf90_put_att(ncid, nf90_global, 'lon_ur_rll', rll_opts%urcnr(1)))
   call check(nf90_put_att(ncid, nf90_global, 'lat_ur_rll', rll_opts%urcnr(2)))
   call check(nf90_put_att(ncid, nf90_global, 'dlon_rll',   rll_opts%dlon))
   call check(nf90_put_att(ncid, nf90_global, 'dlat_rll',   rll_opts%dlat))
   call check(nf90_put_att(ncid, nf90_global, 'latitude_domain_center', rll_opts%ctr_lat))
   call check(nf90_put_att(ncid, nf90_global, 'longitude_domain_center', rll_opts%ctr_lon))
   call check(nf90_put_att(ncid, nf90_global, 'latitude_south_pole_rotated', rll_opts%sp_lat))
   call check(nf90_put_att(ncid, nf90_global, 'longitude_south_pole_rotated', rll_opts%sp_lon))
   call check(nf90_put_att(ncid, nf90_global, 'azimuth_rotated', 0.0))
   call check(nf90_put_att(ncid, varid_glon,  'description', 'Earth geographical longitude'))
   call check(nf90_put_att(ncid, varid_glon,  'units', 'degree_east'))
   call check(nf90_put_att(ncid, varid_glat,  'description', 'Earth geographical latitude'))
   call check(nf90_put_att(ncid, varid_glat,  'units', 'degree_north'))
   call check(nf90_put_att(ncid, varid_rlon,  'description', 'rotated longitude'))
   call check(nf90_put_att(ncid, varid_rlon,  'units', 'degree_east'))
   call check(nf90_put_att(ncid, varid_rlat,  'description', 'rotated latitude'))
   call check(nf90_put_att(ncid, varid_rlat,  'units', 'degree_north'))
   call check(nf90_put_att(ncid, varid_data,  'units', var_opts%units))
   call check(nf90_put_att(ncid, varid_data,  'description', var_opts%description))

!--- End definition of variables
   call check( nf90_enddef(ncid) )

!--- Write the data to new netcdf file
   call check(nf90_put_var(ncid, varid_glon, real(glon_o,8)))
   call check(nf90_put_var(ncid, varid_glat, real(glat_o,8)))

   allocate(rotlon(imdl_o))
   allocate(rotlat(jmdl_o))
   do iii = 1, imdl_o
     rotlon(iii) = rll_opts%llcnr(1) + rll_opts%dlon*(iii-1)
     if (rotlon(iii) >= 360.0_dp ) rotlon(iii) = rotlon(iii) - 360.0_dp
     if (rotlon(iii) <    0.0_dp ) rotlon(iii) = rotlon(iii) + 360.0_dp
   end do
   do jjj = 1, jmdl_o
     rotlat(jjj) = rll_opts%llcnr(2) + rll_opts%dlat*(jjj-1)
   end do
   call check(nf90_put_var(ncid, varid_rlon, rotlon))
   call check(nf90_put_var(ncid, varid_rlat, rotlat))

   call check(nf90_put_var(ncid, varid_data, data_o))

!- Close the dataset
   call check( nf90_close(ncid) )

!------------------------------------------------------------------------
! 5. Finalize
!-----------------------------------------------------------------------!
!--- clean the memory
   deallocate(igdtmpl_o)
   deallocate(output_data)
   deallocate(output_bitmap)
   deallocate(ibo)
   deallocate(xpts_o, ypts_o)
   deallocate(glat_o, glon_o)
   deallocate(output_bitmap_2d)
   deallocate(data_o)
   deallocate(data_fgs_o)
   deallocate(data_tmp_o)
   deallocate(rotlon, rotlat)

   deallocate(xpts_i, ypts_i)
   deallocate(glat_i, glon_i)
   deallocate(slmask_i)
   deallocate(data_i)
 enddo

!--- close grib file
 call baclose(iunit, iret)
 if (iret /= 0) stop 'Error: baclose failed.'

!--- free the memory usage by gfld
 call gf_free(gfld_input)

!-----------------------------------------------------------------------
      contains
!
      subroutine check(status)
      integer,intent(in) :: status
!
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped for nf90 error."
      end if
      return
      end subroutine check
!-----------------------------------------------------------------------
!
!================================================================================
end program rtma_regrid_esg2rll
