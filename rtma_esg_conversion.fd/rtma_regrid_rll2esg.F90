program rtma_regrid_rll2esg
!================================================================================
  use omp_lib
  use netcdf                               ! netcdf lib
  use pkind, only: dp, sp, dpi, spi        ! Jim Purser's lib
  use pietc, only: dtor,rtod               ! Jim Purser's lib
! use gdswzd_mod, only: gdswzd             ! ip lib (for conversion from earth to grid coor or vice versa)
!-- required when using iplib v4.0 or higher
#ifndef IP_V3
! use ip_mod                               ! ip lib v4 or above (for interpolation)
  use ipolates_mod                         ! ip lib v4 or above (for interpolation)
#endif
  use grib_mod                             ! g2 lib(grib)
  use bacio_module                         ! prerequisite for g2 lib (grib)
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

  character(len=100)      :: input_data_rll_file_grb2
  character(len=100)      :: esg_grid_spec_file_nc
  character(len=100)      :: output_data_esg_file_nc
  character(len=20)       :: varname_nc
  character(len=20)       :: varname_input
 
 
  integer                 :: iunit, iret, lugi

  integer                 :: ip, ipopt(20)     ! parameters for interpolation
 
  integer                 :: j, jdisc, jpdtn, jgdtn, k
  integer                 :: jids(200), jgdt(200), jpdt(200)
  integer, allocatable    :: igdtmpl_i(:)
  integer                 :: igdtlen_i
  integer                 :: igdtnum_i

#ifdef IP_V3
  integer                 :: igdt_grb2(5) 
  integer                 :: idefnum                 ! The number of entries in array ideflist, i.e. number
                                                     ! of rows (or columns) for which optional grid points are defined.
  integer,  allocatable   :: ideflist(:)             ! integer array containing the number of grid points contained in
                                                     ! each row (or column). To handle the irregular grid stuff.
  integer                 :: kgds_rll_i(200)         ! grib1 GDS for rotated_latlon grid
  integer                 :: kgds_esg_o(200)         ! grib1 GDS for esg grid (used in ip lib v3.x and older)
  integer                 :: igrid_rll_i             ! ncep re-defined grib1 grid number
#else
  integer                 :: igdtlen_o
  integer                 :: igdtnum_o
  integer, allocatable    :: igdtmpl_o(:)
#endif
 
  integer                 :: mi
  integer                 :: imdl_i, jmdl_i         ! dimension size read in grib2 data file
  integer                 :: npts_i
  logical                 :: unpack
  integer, allocatable    :: ibi(:)
  logical*1, allocatable  :: input_bitmap(:,:)      ! 2D array to match ipolates_grib2
  real(dp), allocatable   :: input_data(:,:)        ! 2D array to match ipolates_grib2

  integer                 :: adate(5)               ! year/month/day/hour/minute
  character(len=12)       :: cdate                  ! yyyymmddhhmn

  integer                 :: km

  integer                 :: mo, no
  integer                 :: ipts_o, jpts_o         ! dimension size read in netcdf data file
  integer                 :: ipt2_o, jpt2_o         ! dimension size read in netcdf data file
  integer                 :: npts_o

  integer, allocatable    :: ibo(:)
  logical*1, allocatable  :: output_bitmap(:,:)      ! 2D array to match ipolates_grib2
  real(dp), allocatable   :: output_data(:,:)        ! 2D array to match ipolates_grib2
  real(dp), allocatable   :: output_glat(:)          ! 2D Data in 1D array
  real(dp), allocatable   :: output_glon(:)          ! 2D Data in 1D array

  real, allocatable       :: glon_o(:,:), glat_o(:,:)
  real, allocatable       :: slmask_o(:,:)
  real, allocatable       :: data_o(:,:)

  integer                 :: ncid, varid
  integer                 :: xdimid, ydimid, timedimid
  integer                 :: xt_dimid, yt_dimid
  integer                 :: data_varid

  logical*1               :: l_clean_bitmap   ! if true, then set bitmap = true everywhere
  logical*1               :: verbose

  character(100) :: fname_nml
  logical        :: f_exist 
  integer        :: lunin_nml

  namelist/setup/varname_input, verbose, l_clean_bitmap
#ifdef IP_V3
  external :: ipolates
#endif

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

  f_exist = .false.
  lunin_nml = 10
!-- reading namelist
  fname_nml = 'rll2esg_namelist'
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
! 1. Read the input data and the input grid info (grib2 file, on rotated grid)
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
     stop(2)
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
   imdl_i = gfld_input%igdtmpl(8)  ! x-dimension of input grid
   jmdl_i = gfld_input%igdtmpl(9)  ! y-dimension of input grid
   npts_i = imdl_i * jmdl_i        ! total dimension of input grid
   mi     = npts_i                 ! total number of pts, input grid
   write(6,'(1x,A,3(1x,I8))')                                                  &
        'dimension of input grid (Nx, Ny, Nx*Ny) = ', imdl_i, jmdl_i, npts_i

!  jgdtn was set to be -1 so getgb2 would search for grid template number.
!  So after getgb2, jgdtn is still -1, but the true grid template number 
!  is gfld_input%igdtnum.
   igdtnum_i = gfld_input%igdtnum
   if (igdtnum_i == 1) then
     write(6,'(1x,A)') "input model grid is defined on rotated lat-lon grid, as expected."
     call set_rllgridopts(gfld_input, rll_opts)
   else
     write(6,'(1x,A,I12.12)') &
       "input model grid is defined on the grid with grid template number = ",igdtnum_i
     write(6,'(1x,A)') ' However, a Rotated Lat-Lon Grid is expected. So Abort this running ...'
     stop(-1)
   end if

   igdtlen_i = gfld_input%igdtlen
   allocate(igdtmpl_i(igdtlen_i))
   igdtmpl_i(:) = gfld_input%igdtmpl(:)

!---------------------------------------------------------------------------
!   1.4 checking up with input model data (and its bitmap, etc.) 
!---------------------------------------------------------------------------
!--- input data field decoded in grib2 file
!  allocate(input_data(npts_i))
   allocate(input_data(npts_i, 1))
   input_data(:,1)         = gfld_input%fld  ! the input data field

!---  does input data have a bitmap?
   allocate(ibi(1))
   allocate(input_bitmap(npts_i, 1))
!  allocate(input_bitmap(npts_i))
   l_clean_bitmap = .False.                  ! do not re-set bitmap to be true everywhere
   call set_bitmap_grb2(gfld_input,npts_i,l_clean_bitmap,ibi,input_bitmap(:,1))

!---  checking if existing any abnormal data values and counting data with false bitmap
   write(6,*)'----------------------------------------------------------'
   write(6,*)'checking the input data read from grib2 fie:'
   write(6,*)'----------------------------------------------------------'
   call check_data_1d_with_bitmap(var_opts,npts_i,input_data(:,1),ibi,input_bitmap(:,1))
 
#ifdef IP_V3
!-------------------------------------------------------------------------------------!
!   1.5 converting GRIB2 GDT info to GIB1 GDS info for                                !
!    compatibility backwards with IP lib v3 and older                                 !
!    As required by IP lib v3.x or older, subroutine ipolates requires grib1 GDS info,!
!    (in IP lib v4 & above, ipolates works with either grib1 GDS or grib2 GDT info),  !
!    so need to convert grid informaton from GRIB2 Grid Description Section (GDS) info!
!    (including its Grid Definition Template) to GRIB1 GDS info                       !
!    (similarly as decoded by w3fi63)                                                 !
!-------------------------------------------------------------------------------------!
   igdt_grb2(1) = 0         ! Source of Grid Definition (0: then specified in Code Table 3.1)
   igdt_grb2(2) = npts_i    ! Number of Data Points in the defined grid
   igdt_grb2(3) = 0         ! Number of octets needed for each additional grid definition (if 0: using regular grid)
   igdt_grb2(4) = 0         ! Interpetation of list of for optional points definition (Table 3.11)
   igdt_grb2(5) = igdtnum_i ! GrdDefTmplt number(Table 3.1): 1--> rotated latlon grid (Template 3.1)
   idefnum      = 0         ! no irregular grid stuff
   allocate(ideflist(1))
   ideflist(:)  = 0         ! no irregular grid stuff
   call gdt2gds_rll(igdt_grb2, igdtlen_i, igdtmpl_i, kgds_rll_i, igrid_rll_i, iret)
   deallocate(ideflist)
   if ( verbose ) write(6,*) ' checking kgds calculated by gdt2gds_rll : ', kgds_rll_i
#endif

!---------------------------------------------------------------------------
! 2. Read information of the output ESG grid in fv3_grid_specification file (netcdf)
!---------------------------------------------------------------------------
   esg_grid_spec_file_nc = "./fv3_grid_spec_esg.nc"
   call check( nf90_open(esg_grid_spec_file_nc, nf90_nowrite, ncid) )
!--- inquire dimension id and the dimension size
   call check( nf90_inq_dimid(ncid, "grid_xt", xt_dimid) )             ! cell center
   call check( nf90_inquire_dimension(ncid, xt_dimid, len = ipts_o) )   
   call check( nf90_inq_dimid(ncid, "grid_yt", yt_dimid) )             ! cell center
   call check( nf90_inquire_dimension(ncid, yt_dimid, len = jpts_o) )   

   npts_o = ipts_o * jpts_o
   write(6,*)'=========================================================='
   write(6,*)'- Checking the output grid (ESG grid)                    -'
   write(6,*)'----------------------------------------------------------'
   write(6,'(3(1x,A,I8))') 'output ESG grid dimensions -- ipts_o= ',  &
         ipts_o, ' jpts_o= ', jpts_o, ' npts_o= ', npts_o

   allocate(glon_o(ipts_o, jpts_o))
   allocate(glat_o(ipts_o, jpts_o))

   varname_nc = "grid_lont"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, glon_o ))
   write(6,'(1x,3A,4(1x,F12.6))') 'read-in ', trim(adjustl(varname_nc)), ' at 4 corners =', &
        glon_o(1,1),glon_o(1,jpts_o),glon_o(ipts_o,jpts_o),glon_o(ipts_o,1)

   varname_nc = "grid_latt"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, glat_o ))
   write(6,'(1x,3A,4(1x,F12.6))') 'read-in ', trim(adjustl(varname_nc)), ' at 4 corners =', &
        glat_o(1,1),glat_o(1,jpts_o),glat_o(ipts_o,jpts_o),glat_o(ipts_o,1)

   call check( nf90_close(ncid) )

!---------------------------------------------------------------------------
! 3. Transfer input data from Rotated Grid (RLL) to ESG grid
!      with interpolation
!---------------------------------------------------------------------------
!   Note:
!        the output "grid" is ESG grid, treated as a series of random 
!        station points. In this case, set the grid definition template 
!        number of a negative number. The grid definition template array 
!        information is not used, so set to a flag value.
!---------------------------------------------------------------------------
!   3.1 setup arguments for ipolates (scalar interpolation) call
!---------------------------------------------------------------------------
   mo = npts_o
   no = mo
   allocate (ibo(1))
   allocate (output_glat(mo))
   allocate (output_glon(mo))
   allocate (output_data(mo,1))
   allocate (output_bitmap(mo,1))
!--- rehsaping 2D array to 1D array required by ipolates
   output_glon = reshape(glon_o, (/npts_o/))
   output_glat = reshape(glat_o, (/npts_o/))

   ip       = 0                         ! bilinear interpolation
   ipopt(:) = 0                         ! options for bilinear
   ipopt(1) = 75                        ! set minimum mask to 75%

!---------------------------------------------------------------------------
!   3.2 call ipolates to interpolate scalar data.
!       non-zero "iret" indicates a problem.
!    Note:
!         the output "grid" is ESG grid, treated as a series of random 
!         station points in the interpolation with IP lib.
!         In this case, set the grid definition template 
!         number of a negative number. The grid definition template array 
!         information is not used, so set to a flag value.
!---------------------------------------------------------------------------
#ifdef IP_V3
!-------------------------------------------------------------------------------------!
!    As required by IP lib v3.x or older, ipolates requires grib1 GDS info,           !
!-------------------------------------------------------------------------------------!
!--- ESG grid points are treated as random station points in the interpolation.
   kgds_esg_o(:) =  0
   kgds_esg_o(1) = -1       ! KGDSO(1)<0 IMPLIES RANDOM STATION POINTS
   write(6,*) ' call ipolates with IP lib v3.x and older) ... '
   call ipolates(ip, ipopt, kgds_rll_i, kgds_esg_o,                          &
                 mi, mo, km, ibi, input_bitmap, input_data, no, output_glat, &
                 output_glon, ibo, output_bitmap, output_data, iret)
#else
!-------------------------------------------------------------------------------------!
!    In IP lib v4 and above, ipolates works with either grib1 GDS or grib2 GDT info.  !
!-------------------------------------------------------------------------------------!
!--- ESG grid points are treated as random station points in the interpolation.
   igdtnum_o = -1           ! set the grid definition template number of a negative number
   igdtlen_o =  1
   allocate(igdtmpl_o(igdtlen_o))
   igdtmpl_o =  -9999       ! grid definition template array info is not use, set to a flag value.
   write(6,*) ' call ipolates(==>ipolates_grib2 with IP lib v4.x and above) ... '
   call ipolates(ip, ipopt, igdtnum_i, igdtmpl_i, igdtlen_i,                 &
                 igdtnum_o, igdtmpl_o, igdtlen_o,                            &
                 mi, mo, km, ibi, input_bitmap, input_data, no, output_glat, &
                 output_glon, ibo, output_bitmap, output_data, iret)
   deallocate(igdtmpl_o)
#endif
   if (iret /= 0) then
     write(6,'(1x,A,I4,A)') ' ipolates failed with returned value iret = ', iret, '.'
     stop(7)
   else
     write(6,*) ' ipolates was done successfully.'
   end if
!---  checking if existing any abnormal data values and counting data with false bitmap
   write(6,*)'----------------------------------------------------------'
   write(6,*)'checking the regridded data interpolated by subroutine ipolates_grib2:'
   call check_data_1d_with_bitmap(var_opts,npts_o,output_data,ibo,output_bitmap)

!---------------------------------------------------------------------------
! 4. Output(appending) the regrided output data to the existing ESG grid file (netcdf)
!---------------------------------------------------------------------------
   output_data_esg_file_nc = "./output_data_esg.nc"
   call check( nf90_open(output_data_esg_file_nc, nf90_write, ncid) )
!--- inquire dimension id of output data file
   call check( nf90_inq_dimid(ncid, "xaxis_1",    xdimid) )   
   call check( nf90_inquire_dimension(ncid, xdimid, len = ipt2_o) )   
   call check( nf90_inq_dimid(ncid, "yaxis_1",    ydimid) )   
   call check( nf90_inquire_dimension(ncid, ydimid, len = jpt2_o) )   
   call check( nf90_inq_dimid(ncid, "Time",    timedimid) )
!--- check if dimension size of output data file match the dimension size of fv3_grid_spec file
   if ( ipt2_o /= ipts_o .or. jpt2_o /= jpts_o ) then
     write(6,'(1x,3A)') 'WARNING --> dimensions of output data file ', output_data_esg_file_nc, &
          ' do NOT match dimensions of fv3_grid_spec file. <-- Warning'
     write(6,'(1x,A,2(1x,I8))') 'dimension of output data   file : ', ipt2_o, jpt2_o
     write(6,'(1x,A,2(1x,I8))') 'dimension of fv3_grid_spec file : ', ipts_o, jpts_o
     write(6,*) ' Check the dimenesions above. Now ABORT the task ...'
     stop(5)
   end if
!--- read sea-land mask in ESG data file
   allocate(slmask_o(ipts_o, jpts_o))
   varname_nc = "slmsk"
   call check( nf90_inq_varid(ncid, trim(adjustl(varname_nc)), varid) )
   call check( nf90_get_var(ncid, varid, slmask_o ))
   write(6,'(1x,3A,5(1x,F12.6))') 'read-in ', trim(adjustl(varname_nc)),       &
        ' at 4 corners & center =',slmask_o(1,1),slmask_o(1,jpts_o),           &
        slmask_o(ipts_o,jpts_o),slmask_o(ipts_o,1),slmask_o(ipts_o/2,jpts_o/2)

!--- reshaping 1-D output array to 2-D array
   allocate(data_o(ipts_o, jpts_o))
   data_o = reshape(output_data(:,1), (/ipts_o, jpts_o/), order=(/1,2/))

!--- applying some constraints to the regridded variables
   if ( var_opts%varname == "howv" ) then
     where(slmask_o .gt. 0.01 ) data_o = 0.0_dp       ! re-set wave height to be 0 over the land area
     where(data_o   .lt. 0.0  ) data_o = 0.0_dp       ! wave height >=0
   else if ( var_opts%varname == "gust" ) then
     where(data_o   .lt. 0.0  ) data_o = 0.0_dp       ! wind gust >=0
   end if

!--- define new variable for the regridded data and output
   call check( nf90_redef(ncid) )
   varname_nc=trim(adjustl(var_opts%varname))
!  call check( nf90_def_var(ncid, trim(adjustl(varname_nc)), nf90_double,  &
   call check( nf90_def_var(ncid, trim(adjustl(varname_nc)), nf90_float,   &
                           (/ xdimid, ydimid, timedimid /), data_varid,    &
                           contiguous=.false.,                             &
                           chunksizes=(/ipts_o, jpts_o, 1/),               &
                           shuffle = .true., fletcher32 = .true.,          &
                           endianness = nf90_endian_little) )
   call check( nf90_def_var_fill(ncid, data_varid, 1, -9999.0) )  ! set No FillValue

   call check( nf90_enddef(ncid) )
!--- output data to new variable
   call check( nf90_put_var(ncid, data_varid, data_o))
   write(6,'(3(1X,A))') 'create new variable [', trim(adjustl(varname_nc)),&
         '] and append it to input netcdf file.'

   call check( nf90_close(ncid) )

!------------------------------------------------------------------------
! 5. Finalize
!-----------------------------------------------------------------------!
!--- clean the memory
   deallocate(igdtmpl_i)
   deallocate(input_data)
   deallocate(input_bitmap)
   deallocate (ibi)
   deallocate (ibo)
   deallocate(glat_o, glon_o)
   deallocate(output_glat)
   deallocate(output_glon)
   deallocate(output_data)
   deallocate(output_bitmap)

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
end program rtma_regrid_rll2esg
