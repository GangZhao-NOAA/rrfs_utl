module mod_rtma_regrid
!--- module interface
  use pkind, only: dp, sp, dpi, spi        ! Jim Purser's lib
  use pietc, only: dtor
  use grib_mod
  use pesg,  only: gtoxm_ak_dd_g 

  implicit none

  private

  public :: rotated_gridopts
  type :: rotated_gridopts
       real(dp) :: ctr_lon              ! earth longitude of domain center (unit: deg)
       real(dp) :: ctr_lat              ! earth latitude  of domain center (unit: deg)
       real(dp) :: sp_lon               ! earth longitude of south pole (unit: deg)
       real(dp) :: sp_lat               ! earth latitude  of south pole (unit: deg)
       real(dp) :: dlon                 ! grid spacing along x-direction/rotated-longitude (deg)
       real(dp) :: dlat                 ! grid spacing along y-direction/rotated-latitude  (deg)
       real(dp) :: llcnr(2)             ! rotated lon/lat of lower-left  corner of domain (deg)
       real(dp) :: urcnr(2)             ! rotated lon/lat of upper-right corner of domain (deg)
       integer  :: nx                   ! grid dimension in x-direction
       integer  :: ny                   ! grid dimension in y-direction
  end type rotated_gridopts

  public :: variable_options
  type :: variable_options
       character(len=20)  :: varname
       character(len=10)  :: var_ncf
       character(len=10)  :: var_grb
       character(len=100) :: des_grb
       character(len=100) :: description
       character(len=10)  :: units
       real(dp)           :: lower_bound
       real(dp)           :: upper_bound
  end type variable_options

  public :: esg_gridopts
  type :: esg_gridopts
       real(dp)           :: A
       real(dp)           :: Kappa
       real(dp)           :: delx                     ! in degree, not radian (*6370*pi/180=1.4415 km, half grid/cell size)
       real(dp)           :: dely                     ! in degree, not radian
       real(dp)           :: plat                     ! center lat of gnomonic grid
       real(dp)           :: plon                     ! center lon of gnomonic grid  ! -112.5_dp = 247.5_dp
       real(dp)           :: pazi             
  end type esg_gridopts

  public :: set_esg_gridopts
  public :: set_variable_options
  public :: check_varopts_grb2
  public :: set_time4data
  public :: check_grbmsg
  public :: set_rllgridopts
  public :: set_bitmap_grb2
  public :: check_data_1d_with_bitmap
  public :: ll_to_xy_esg
#ifdef IP_V3
  public :: gdt2gds_rll
#endif

!-----------------------------------------------------------------------
!
!---- subroutines used in rtma_regrid main code
  contains
      subroutine set_esg_gridopts(esg_opts)
        type(esg_gridopts),         intent(inout) :: esg_opts
!
        esg_opts%A      = 0.183131392268429_dp
!       esg_opts%A      = 0.183035309495599_dp     ! from Ben's rrfs-workflow v0.9.1 (but with larger error ????)

        esg_opts%Kappa  =-0.265835885178773_dp
!       esg_opts%Kappa  =-0.265943041019571_dp     ! from Ben's rrfs-workflow v0.9.1

        esg_opts%delx   = 0.0129656309291837_dp    ! in degree, not radian (*6370*pi/180=1.4415 km, half grid/cell size)
!       esg_opts%delx   = 0.0129665672419664_dp    ! from Ben's rrfs-workflow v0.9.1

        esg_opts%delx   = esg_opts%delx * dtor     ! convert degree to radian

        esg_opts%dely   = 0.0132456552576884_dp
!       esg_opts%dely   = 0.0132464375726782_dp     ! from Ben's rrfs-workflow v0.9.1

        esg_opts%dely   = esg_opts%dely * dtor
        esg_opts%plat   = 55.0_dp                  ! center lat of gnomonic grid
        esg_opts%plon   =-112.5_dp                 ! center lon of gnomonic grid  ! -112.5_dp = 247.5_dp
        esg_opts%pazi   = 0.0_dp
        return
      end subroutine set_esg_gridopts
!
      subroutine set_variable_options(var_opts,iret)
        type(variable_options),  intent(inout) :: var_opts
        integer,                 intent(  out) :: iret  

        if ( trim(adjustl(var_opts%varname)) .eq. 'howv' .or.                  &
             trim(adjustl(var_opts%varname)) .eq. 'HOWV'        ) then
          var_opts%varname     = 'howv'
          var_opts%lower_bound =   0.0
          var_opts%upper_bound = 100.0
          var_opts%var_grb     = 'HTSGW'
          var_opts%des_grb     = 'Significant Height of Combined Wind Waves and Swell'
          var_opts%var_ncf     = 'howv'
          var_opts%description = 'Significant Height of Combined Wind Waves and Swell'
          var_opts%units       = 'm'
          iret = 0
        else if ( trim(adjustl(var_opts%varname)) .eq. 'gust' .or.             &
                  trim(adjustl(var_opts%varname)) .eq. 'GUST'        ) then
          var_opts%varname     = 'gust'
          var_opts%lower_bound =   0.0
          var_opts%upper_bound = 100.0
          var_opts%var_grb     = 'GUST'
          var_opts%des_grb     = 'Wind Speed (Gust) 10-m'
          var_opts%var_ncf     = 'gust'
          var_opts%description = 'Surface Wind Gust at 10 meters'
          var_opts%units       = 'm/s'
          iret = 0
        else
          var_opts%varname     = 'unknown'
          var_opts%lower_bound = -999999999.0
          var_opts%upper_bound =  999999999.0
          var_opts%var_grb     = 'UNKNOWN'
          var_opts%des_grb     = 'Unknown'
          var_opts%var_ncf     = 'unknown'
          var_opts%description = 'Unknown'
          var_opts%units       = 'unknown'
          iret = -1
        end if

        if ( iret .eq. 0 ) write(6,'(1x,3A,2(F12.6,A))')                       &
          ' input variable name is ', trim(adjustl(var_opts%varname)),         &
          ' its pre-set range is [', var_opts%lower_bound, ' ~ ', var_opts%upper_bound, '].'
        return
      end subroutine set_variable_options
!
      subroutine check_varopts_grb2(gfld,var_opts,iret)
        use grib_mod
        type(gribfield),   intent(in   ) :: gfld
        type(variable_options),       intent(in   ) :: var_opts
        integer,           intent(  out) :: iret

        if ( gfld%ipdtnum     .eq.  0 .and.           & ! 0: Anl or fcst in a horizontal layer, see Grib2 Code Table 4.0
             gfld%discipline  .eq.  0 .and.           & ! Discipline 0 : Meteorological Products, see Table 4.1
             gfld%ipdtmpl(1)  .eq.  2 .and.           & ! Cateory 2:  Momentum, See Table 4.1 
             gfld%ipdtmpl(2)  .eq. 22 .and.           & ! Number 22:  GUST, see Table 4.2-0-2
             gfld%ipdtmpl(10) .eq.  1 .and.           & ! Prod.Templt. 4.0 Octet. 23 --> 1: Ground or Water Surface, see Table 4.5
             var_opts%varname .eq. 'gust' ) then
          iret = 0                  ! specified variable name matches the variable in grib2 file
        else if ( gfld%ipdtnum    .eq.  0 .and.      & ! 0: Anl or fcst in a horizontal layer, see Grib2 Code Table 4.0
             gfld%discipline  .eq. 10 .and.           & ! Discipline 10 : Oceanographic Products, see Table 4.1
             gfld%ipdtmpl(1)  .eq.  0 .and.           & ! Cateory 0:  Waves, See Table 4.1 
             gfld%ipdtmpl(2)  .eq.  3 .and.           & ! Number  3:  HTSGW, see Table 4.2-10-0
             gfld%ipdtmpl(10) .eq.  1 .and.           & ! Prod.Templt. 4.0 Octet. 23 --> 1: Ground or Water Surface, see Table 4.5
             var_opts%varname .eq. 'howv' ) then
          iret = 0                  ! specified variable name matches the variable in grib2 file
        else
          iret = -1                 ! variable name mis-matches the variable in grib2 file
        end if
        return
      end subroutine check_varopts_grb2
!
      subroutine set_time4data(gfld,adate,cdate)
        type(gribfield),   intent(in   ) :: gfld
        integer,           intent(  out) :: adate(5)               ! year/month/day/hour/minute
        character(len=12), intent(  out) :: cdate                  !yyyymmddhhmn
  
        integer                 :: itt
  
        adate(1)=gfld%idsect( 6)         ! Year(4digits)
        adate(2)=gfld%idsect( 7)         ! Month
        adate(3)=gfld%idsect( 8)         ! Day
        adate(4)=gfld%idsect( 9)         ! Hour
        adate(5)=gfld%idsect(10)         ! Minute
!         furtherly adjusting the time by info in Product template
        if ( gfld%ipdtmpl(8) .le. 4) then
          itt=5-gfld%ipdtmpl(8)
          adate(itt)=adate(itt) + gfld%ipdtmpl(9)
          write(6,'(1x,A,I4,1x,A,I4)') 'date time is adjusted by ', gfld%ipdtmpl(9), &
                'with unit indicator-->', gfld%ipdtmpl(9)
        else
          write(6,'(1x,A,1x,I4,1x,A)') ' * * * Warning: checking the Indicator of unit of time range: ', &
                gfld%ipdtmpl(8), ' No further adjustment to time. * * * * * * '
        end if
        write(cdate( 1: 4),'(I4.4)') adate( 1)
        write(cdate( 5: 6),'(I2.2)') adate( 2)
        write(cdate( 7: 8),'(I2.2)') adate( 3)
        write(cdate( 9:10),'(I2.2)') adate( 4)
        write(cdate(11:12),'(I2.2)') adate( 5)
        write(6,'(1x,A,1x,A4,4(A1,A2))') 'date of input model data: ',   &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),'_',cdate(9:10),':',cdate(11:12)
        return
      end subroutine set_time4data
!
      subroutine check_grbmsg(gfld)
        type(gribfield),    intent(in   ) :: gfld
        write(6,*) ' checking the data in array gfld read by getgb2: '
        write(6,*) ' version: ',    gfld%version
        write(6,*) ' discipline: ', gfld%discipline
        write(6,*) ' idsectlen: ',  gfld%idsectlen
        write(6,*) ' idsect(:):, ', gfld%idsect(1:gfld%idsectlen)
        write(6,*) ' locallen: ',   gfld%locallen
        if (gfld%locallen > 0) &
          write(6,*) ' local(:): ',   gfld%local(1:gfld%locallen)
        write(6,*) ' ifldnum: ',    gfld%ifldnum
        write(6,*) ' griddef: ',    gfld%griddef
        write(6,*) ' ngrdpts: ',    gfld%ngrdpts
        write(6,*) ' numoct_opt: ', gfld%numoct_opt
        write(6,*) ' interp_opt: ', gfld%interp_opt
        write(6,*) ' num_opt: ',    gfld%num_opt
        write(6,*) ' igdtnum: ',    gfld%igdtnum
        write(6,*) ' igdtlen: ',    gfld%igdtlen
        write(6,*) ' igdtmpl(:): ', gfld%igdtmpl(1:gfld%igdtlen)
        write(6,*) ' ipdtnum: ',    gfld%ipdtnum
        write(6,*) ' ipdtlen: ',    gfld%ipdtlen
        write(6,*) ' ipdtmpl(:): ', gfld%ipdtmpl(1:gfld%ipdtlen)
        write(6,*) ' num_coord: ',  gfld%num_coord
        if (gfld%num_coord > 0) &
          write(6,*) ' coord_list(:) ', gfld%coord_list(:)
        write(6,*) ' ndpts: ',      gfld%ndpts
        write(6,*) ' idrtnum: ',    gfld%idrtnum
        write(6,*) ' idrtlen: ',    gfld%idrtlen
        write(6,*) ' idrtmpl(:): ', gfld%idrtmpl(1:gfld%idrtlen)
        write(6,*) ' unpacked: ',   gfld%unpacked
        write(6,*) ' expanded: ',   gfld%expanded
        write(6,*) ' ibmap: ',      gfld%ibmap
        write(6,*) ' length of bmap(:) ',     size(gfld%bmap)
        write(6,*) ' fld(        1): ', gfld%fld(1)
        write(6,*) ' fld(ngrdpts/2): ', gfld%fld(nint(gfld%ngrdpts / 2.0))
        write(6,*) ' fld(ngrdpts  ): ', gfld%fld(gfld%ngrdpts)
        return
      end subroutine check_grbmsg
!
      subroutine set_rllgridopts(gfld,rll_gridopts)
        type(gribfield),         intent(in   ) :: gfld
        type(rotated_gridopts),  intent(  out) :: rll_gridopts
        rll_gridopts%sp_lon   = gfld%igdtmpl(21)/1000000.0
        rll_gridopts%sp_lat   = gfld%igdtmpl(20)/1000000.0
        if (rll_gridopts%sp_lat .le. 0.0_dp) then
          rll_gridopts%ctr_lat = 90.0_dp + rll_gridopts%sp_lat
          rll_gridopts%ctr_lon = rll_gridopts%sp_lon
        else
          rll_gridopts%ctr_lat = 90.0_dp - rll_gridopts%sp_lat
          if ( rll_gridopts%sp_lon .gt. 180.0_dp ) then
            rll_gridopts%ctr_lon = rll_gridopts%sp_lon - 180.0_dp
          else
            rll_gridopts%ctr_lon = rll_gridopts%sp_lon + 180.0_dp
          end if
        end if
        rll_gridopts%dlon     = gfld%igdtmpl(17)/1000000.0    ! Di -- i direction increment
        rll_gridopts%dlat     = gfld%igdtmpl(18)/1000000.0    ! Dj -- j direction increment
        rll_gridopts%llcnr(1) = gfld%igdtmpl(13)/1000000.0
        rll_gridopts%llcnr(2) = gfld%igdtmpl(12)/1000000.0
        rll_gridopts%urcnr(1) = gfld%igdtmpl(16)/1000000.0
        rll_gridopts%urcnr(2) = gfld%igdtmpl(15)/1000000.0
        rll_gridopts%nx       = gfld%igdtmpl( 8)              ! grid dimension in x/i-direcion
        rll_gridopts%ny       = gfld%igdtmpl( 9)              ! grid dimension in y/j-direcion
        write(6,'(1x, A, 2(1x,F8.3))') ' checking rotated grid parameters: south pole lon/lat: ',  &
              rll_gridopts%sp_lon, rll_gridopts%sp_lat
        write(6,'(1x, A, 2(1x,F8.3))') ' checking rotated grid parameters: domain center lon/lat: ',  &
              rll_gridopts%ctr_lon, rll_gridopts%ctr_lat
        write(6,'(1x, A, 2(1x,F8.3))') ' checking rotated grid parameters: grid-spacing dlon/dlat: ',  &
              rll_gridopts%dlon, rll_gridopts%dlat
        write(6,'(1x, A, 2(1x,F8.3))') ' checking rotated grid parameters: lower-left  corner lon/lat: ',  &
              rll_gridopts%llcnr(1), rll_gridopts%llcnr(2)
        write(6,'(1x, A, 2(1x,F8.3))') ' checking rotated grid parameters: upper-right corner lon/lat: ',  &
              rll_gridopts%urcnr(1), rll_gridopts%urcnr(2)
        write(6,'(1x, A, 2(1x,  I8))') ' checking rotated grid parameters: grid dimension in x/y-direction: ',  &
              rll_gridopts%nx, rll_gridopts%ny
        return
      end subroutine set_rllgridopts
!
      subroutine set_bitmap_grb2(gfld,npts,l_clean_bitmap,ibi,input_bitmap)
        type(gribfield),              intent(in   ) :: gfld
        integer,                      intent(in   ) :: npts
        logical*1,                    intent(in   ) :: l_clean_bitmap
        integer,   dimension(1),      intent(  out) :: ibi(1)
        logical*1, dimension(npts),   intent(  out) :: input_bitmap

        if (gfld%ibmap==0) then  ! input data has bitmap
          write(6,*) 'There are bitmap data associated with the data.'
          ibi                   = 1        ! tell ipolates to use bitmap
          input_bitmap(:)       = gfld%bmap
        else                           ! no bitmap, data everywhere
          write(6,*) 'There is NO bitmap data associated with the data.'
          ibi                   = 0        ! tell ipolates there is no bitmap
          input_bitmap(:)       = .true.
        endif

        if ( l_clean_bitmap ) then
          write(6,*) ' Warning --> reset input_bitmap to be true everywhere.'
          ibi                   = 0
          input_bitmap(:)       = .true.
        end if

        return
      end subroutine set_bitmap_grb2
!
      subroutine check_data_1d_with_bitmap(var_opts,npts,input_data,ibi,input_bitmap)
        type(variable_options),       intent(in   ) :: var_opts
        integer,                      intent(in   ) :: npts
        real(dp),  dimension(npts),   intent(in   ) :: input_data          ! 2D Data in 1D slice
        integer,   dimension(1),      intent(in   ) :: ibi(1)
        logical*1, dimension(npts),   intent(in   ) :: input_bitmap

        integer,   dimension(npts) :: index_bitmap
        integer    :: nn, n_bitmap, n_valid
        real(dp)   :: sum_valid

        if ( ibi(1) == 0 ) write(6,*) &
          ' bitmap associated with this data is NOT used or set as TRUE everywhere.'

        n_bitmap = 0
        n_valid = 0
        sum_valid = 0.0_dp
        index_bitmap(:) = -1
        do nn = 1, npts
!---    note: 
!            input_bitmap is  true --> field data at this grid point is valid.
!            input_bitmap is false --> field data at this grid point is invalid.
          if ( .not. input_bitmap(nn) ) then
            n_bitmap = n_bitmap + 1
            index_bitmap(n_bitmap) = nn
          else
            sum_valid = sum_valid + input_data(nn)
            n_valid = n_valid + 1
          end if

          if ( input_data(nn) .lt. var_opts%lower_bound .or.                       &
               input_data(nn) .gt. var_opts%upper_bound      ) then
              write(6,'(1x,A,A,1x,I8,1x,L2,1x,F12.5)')                            &
                'checking input data which is out of range --> ',                 &
                'nn bitmap data_value): ', nn, input_bitmap(nn), input_data(nn)
          end if
        end do

        write(6,'(1x, A, 1x, I8)') &
              ' total number of points with false bitmap : ', n_bitmap
        write(6,'(1x,A,3(1x,L2,1x,F12.5))')                                                             &
              ' checking the bitmap and data values of the first, middle and last invalid data : ',     &
               input_bitmap(index_bitmap(1)), input_data(index_bitmap(1)),                               &
               input_bitmap(index_bitmap(n_bitmap/2)), input_data(index_bitmap(n_bitmap/2)),             &
               input_bitmap(index_bitmap(n_bitmap)), input_data(index_bitmap(n_bitmap))
        write(6,'(1x, A, 1x, I8, 3(1x,F12.5))')                                    &
              'stats of data (masked by bitmap    ) -- size max min ave: ', n_valid,        &
              maxval(input_data, MASK=(input_bitmap)),                             &
              minval(input_data, MASK=(input_bitmap)), sum_valid/n_valid
        write(6,'(1x, A, 1x, I8, 3(1x,F12.5))')                                        &
              'stats of data (not-masked by bitmap) -- size max min ave: ', size(input_data),     &
              maxval(input_data), minval(input_data), sum(input_data)/size(input_data)  

        return
      end subroutine check_data_1d_with_bitmap
!
      subroutine ll_to_xy_esg(nx_esg, ny_esg, nx_rll, ny_rll, esg_opts, lats, lons, x, y)
        integer,                                 intent(in   ) :: nx_esg, ny_esg   ! domain size of esg grid domain (to define the ESG X/Y coordinates)
        integer,                                 intent(in   ) :: nx_rll, ny_rll   ! domain size of another grid domain (output grid)
        type(esg_gridopts),                      intent(in   ) :: esg_opts
        real(dp),    dimension(nx_rll, ny_rll),  intent(in   ) :: lats, lons
        real(dp),    dimension(nx_rll, ny_rll),  intent(inout) :: x, y

        integer    :: ii, jj
        real(dp)   :: dlat, dlon
        real(dp), dimension(2) :: xm
        logical :: ff
        real(dp), parameter          :: two=2_dp

        real(dp)   :: A
        real(dp)   :: Kappa
        real(dp)   :: delx                     ! in degree, not radian (*6370*pi/180=1.4415 km, half grid/cell size)
        real(dp)   :: dely                     ! in degree, not radian
        real(dp)   :: plat                     ! center lat of gnomonic grid
        real(dp)   :: plon                     ! center lon of gnomonic grid  ! -112.5_dp = 247.5_dp
        real(dp)   :: pazi             

        A     = esg_opts%A
        Kappa = esg_opts%Kappa
        plat  = esg_opts%plat
        plon  = esg_opts%plon
        pazi  = esg_opts%pazi
        delx  = esg_opts%delx
        dely  = esg_opts%dely
        do jj=1,ny_rll
          do ii=1,nx_rll
            dlat=lats(ii,jj)
            dlon=lons(ii,jj)
            call gtoxm_ak_dd_g(A,Kappa,plat,plon,pazi,two*delx,two*dely,dlat,dlon,xm,ff) !  multiply delx/dely by 2.0 to get values on compuational grid
            x(ii,jj)=xm(1)
            y(ii,jj)=xm(2)
            x(ii,jj) = (real(nx_esg)-1.0)/2.0 + x(ii,jj) + 1.0  ! Relocate the origin of X/Y coordinate from center
                                                                ! of ESG grid domain to its lower-left corner.
            y(ii,jj) = (real(ny_esg)-1.0)/2.0 + y(ii,jj) + 1.0  ! Plus (1.0, 1.0) is because the coordinates of
                                                                ! lower-left corner is (1.0, 1.0)
          enddo
        enddo

        return
      end subroutine ll_to_xy_esg
!
#ifdef IP_V3
      subroutine gdt2gds_rll(igdt, igdtlen, igdtmpl, kgds, igrid, iret)
!---    Purpose:
!         Covnert rotated latlon grid information from a GRIB2 grid tempalte info
!         to GRIB1 GDS info.  The code is based on subroutine gdt2gds in g2 lib,
!         which does not process igdt(5)=1, also is refered to subrotuine init_grib1
!         and init_grib2 in module ip_rot_equid_cylind_grid_mod of ip lib.
!         Actually, the lat/lon of center grid point in grib1 grid and grib2 grid
!         are based on init_grib1 and init_grib2. The indices used in gdt2gds.F90
!         seem to be wrong.
!         see: https://www.nco.ncep.noaa.gov/pmb/docs/on388/
!              for GDS, see
!                https://www.nco.ncep.noaa.gov/pmb/docs/on388/section2.html
!              for detailed GDS info, see
!                https://www.nco.ncep.noaa.gov/pmb/docs/on388/tabled.html
!
!        Out:
!         kgds: GRIB1 GDS as described in [NCEPLIBS-w3emc w3fi63() function]
!               (https://noaa-emc.github.io/NCEPLIBS-w3emc/w3fi63_8f.html).
!         igrid: NCEP predefined GRIB1 grid number. Set to 255, if not an NCEP grid.
!         iret Error return value:  0: No error.
!                                   1: Unrecognized GRIB2 GDT number.
!
        implicit none
  
        integer, intent(in   ) :: igdtlen
        integer, intent(in   ) :: igdt(5), igdtmpl(igdtlen)
        integer, intent(  out) :: kgds(200)
        integer, intent(  out) :: igrid, iret

        integer :: kgds72(200), kgds71(200), idum(200), jdum(200)
        integer :: ierr, j

        integer    :: iopt
        integer    :: iscale, iscale_gb2, iscale_gb1
        real(dp)   :: lon_sp_rll, lat_sp_rll
        real(dp)   :: rlon,  rlat            ! earth   lat/lon
        real(dp)   :: rlonr, rlatr           ! rotated lat/lon

        external :: w3fi71, r63w72

        iret = 1
        idum = 0
        kgds(1:200) = 0
        if (igdt(5) .eq. 1) then             ! grid number = 1  (in grib2) for Rotated Lat / Lon grid
          kgds( 1) = 205                     ! grid number =205 (in grib1) for Arakawa Staggerred for Non-E Stagger grid
          kgds( 2) = igdtmpl(8)              ! Ni (IM in init_grib)
          kgds( 3) = igdtmpl(9)              ! Nj (JM in init_grib)
          iscale = igdtmpl(10) * igdtmpl(11)
          if ( iscale == 0 ) then
            iscale_gb2 = 10**6
            iscale_gb1 = 10**3
          else
            write(6,'(1x,A,1x,I6.6)') &
                 'gdt2gds_rll::Abort ==> due to Not recognized iscale from grib2 GDT info: iscale= ', iscale
            stop(11)
          end if
          lat_sp_rll = real(igdtmpl(20), dp)/real(iscale_gb2, dp)  ! latitude  of rotated south pole
          lon_sp_rll = real(igdtmpl(21), dp)/real(iscale_gb2, dp)  ! longitude of rotated south pole
!---      first grid point: converting rotated lat/lon (in grib2 gdt) to regular earth lat/lon (for grib1 gds)
          iopt = -1                          ! rotated lat/lon --> regular earth lat/lon
          rlatr = real(igdtmpl(12), dp)/real(iscale_gb2, dp)  ! latitude  of 1st grid point (rotated value in grib2 GDT)
          rlonr = real(igdtmpl(13), dp)/real(iscale_gb2, dp)  ! longitude of 1st grid point (rotated value in grib2 GDT)
          call rll_trans_iplib(lon_sp_rll, lat_sp_rll, iopt, rlonr, rlatr, rlon, rlat)
          write(6,'(1x,2(A,2(1x,F18.9)))') 'gdt2gds_rll::rotated lat lon : ',   &
                rlatr, rlonr, ' (iplib) ==> earth lat lon : ', rlat,  rlon
          kgds( 4) = nint(rlat * real(iscale_gb1, dp)) ! Lat of 1st grid point (earth lat/lon coordinate; RLAT1 in init_grib)
          if ( rlon < 0.0 ) rlon = rlon + 360.0_dp
          kgds( 5) = nint(rlon * real(iscale_gb1, dp)) ! Lon of 1st grid point (earth lat/lon coordinate; RLON1 in init_grib)

          kgds( 6) = 0                       ! resolution and component flags: IROT in init_grib)
!         if (igdtmpl(1)==2) kgds(6) = 64
!         if (btest(igdtmpl(14), 4).OR.btest(igdtmpl(14), 5)) kgds(6) = kgds(6) + 128
!         if (btest(igdtmpl(14), 3)) kgds(6) = kgds(6) + 8
          kgds( 6) = igdtmpl(14)              ! resolution and component flags: IROT in init_grib)

          kgds( 7) = (lat_sp_rll + 90.0_dp) * real(iscale_gb1, dp) ! Earth Latitude  of rotated center point (rotated south pole lat + 90.0;  RLAT0 in init_grib)
          kgds( 8) = lon_sp_rll * real(iscale_gb1, dp)          ! Earth Longitude of rotated center point (=rotated south pole lon;  RLON0 in init_grib)

          kgds( 9) = real(igdtmpl(17), dp)/real(iscale_gb1, dp)   ! Di: x-increment, DLONS in init_grib
          kgds(10) = real(igdtmpl(18), dp)/real(iscale_gb1, dp)   ! Dj: y-increment, DLATS in init_grib

          kgds(11) = igdtmpl(19)              ! Scanning mode (nscan in init_grib)

!---      last grid point: converting rotated lat/lon (in grib2 gdt) to regular earth lat/lon (for grib1 gds)
          iopt = -1                          ! rotated lat/lon --> regular earth lat/lon
          rlatr = real(igdtmpl(15), dp)/real(iscale_gb2, dp)  ! latitude  of last grid point (rotated value in grib2 GDT)
          rlonr = real(igdtmpl(16), dp)/real(iscale_gb2, dp)  ! longitude of last grid point (rotated value in grib2 GDT)
          call rll_trans_iplib(lon_sp_rll, lat_sp_rll, iopt, rlonr, rlatr, rlon, rlat)
          write(6,'(1x,2(A,2(1x,F18.9)))') 'gds2gds_rll: rotated lat lon : ',   &
                rlatr, rlonr, ' (iplib) ==> earth lat lon : ', rlat,  rlon
          kgds(12) = nint(rlat * real(iscale_gb1, dp)) ! Lat of last grid point (earth lat/lon coordinate) (RLAT2 in init_grib)
          if ( rlon < 0.0 ) rlon = rlon + 360.0_dp
          kgds(13) = nint(rlon * real(iscale_gb1, dp)) ! Lon of last grid point (earth lat/lon coordinate) (RLON2 in init_grib)

          kgds(14) = 0
          kgds(15) = 0
          kgds(16) = 0
          kgds(17) = 0
          kgds(18) = 0
          kgds(19) = 0
          kgds(20) = 255
          kgds(21) = 0
          kgds(22) = 0
          iret = 0
        else
          write(6,'(1x, A, I5.5)') 'gdt2gds_rll: Unrecognized GRIB2 GDT = 3.', igdt(5)
          iret = 1
          kgds(1:22) = 0
          return
        endif
!
!       Can we determine NCEP grid number ?
!
        igrid = 255
        do j = 254, 1, -1
          !do j = 225, 225
          kgds71 = 0
          kgds72 = 0
          call w3fi71(j, kgds71, ierr)
          if (ierr.ne.0) cycle
          ! convert W to E for longitudes
          if (kgds71(3) .eq. 0) then    ! lat / lon
            if (kgds71(7) .lt. 0) kgds71(7) = 360000 + kgds71(7)
            if (kgds71(10) .lt. 0) kgds71(10) = 360000 + kgds71(10)
          elseif (kgds71(3) .eq. 1) then    ! mercator
            if (kgds71(7) .lt. 0) kgds71(7) = 360000 + kgds71(7)
            if (kgds71(10) .lt. 0) kgds71(10) = 360000 + kgds71(10)
          elseif (kgds71(3) .eq. 3) then     ! lambert conformal
            if (kgds71(7) .lt. 0) kgds71(7) = 360000 + kgds71(7)
            if (kgds71(9) .lt. 0) kgds71(9) = 360000 + kgds71(9)
            if (kgds71(18) .lt. 0) kgds71(18) = 360000 + kgds71(18)
          elseif (kgds71(3) .eq. 4) then     ! Guassian lat / lon
            if (kgds71(7) .lt. 0) kgds71(7) = 360000 + kgds71(7)
            if (kgds71(10) .lt. 0) kgds71(10) = 360000 + kgds71(10)
          elseif (kgds71(3) .eq. 5) then     ! polar stereographic
            if (kgds71(7) .lt. 0) kgds71(7) = 360000 + kgds71(7)
            if (kgds71(9) .lt. 0) kgds71(9) = 360000 + kgds71(9)
          endif
          call r63w72(idum, kgds, jdum, kgds72)
          if (kgds72(3) .eq. 3) kgds72(14) = 0    ! lambert conformal fix
          if (kgds72(3) .eq. 1) kgds72(15:18) = 0    ! mercator fix
          if (kgds72(3) .eq. 5) kgds72(14:18) = 0    ! polar str fix
          !           print *, ' kgds71(', j, ') =  ',  kgds71(1:30)
          !           print *, ' kgds72        =  ',  kgds72(1:30)
          if (all(kgds71 .eq. kgds72) ) then
             igrid = j
             exit
          endif
        enddo
        write(6,'(1x,A,I6)') 'sub gdt2gds_rll:: igrid = ', igrid

        return
      end subroutine gdt2gds_rll
!-----------------------------------------------------------------------
      subroutine rll_trans_iplib(lon_sp_rll, lat_sp_rll, iopt, lon_in, lat_in, lon_out, lat_out)
!       purpose:
!         conversion between the earth latitude/longitude and the rotated latitude/longitude.
!           iopt =  1: converting earth lat/lon to rotated lat/lon
!           iopt = -1: converting rotated lat/lon to earth lat/lon
!       notes:
!         the algorithm used to do the conversion between rotated latlon and regular earth latlon
!         is based on the code in SUBROUTINE GDSWZD_ROT_EQUID_CYLIND of ip_rot_equid_cylind_grid_mod.F90
!         in IP lib verson 4.3.0
!
        implicit none
!---- parameters
        real(dp),    parameter :: PI = dacos(-1.0_dp)
        real(dp),    parameter :: D2R = PI/180.0_dp    ! degree ==> radian
        real(dp),    parameter :: R2D = 180.0_dp/PI    ! radian ==> degree

        real(dp),    intent(in   ) :: lon_sp_rll    ! earth longitude of south pole after rotated (un
        real(dp),    intent(in   ) :: lat_sp_rll    ! latitude  of input (unit: deg)
        real(dp),    intent(in   ) :: lon_in        ! longitude of input (unit: deg)
        real(dp),    intent(in   ) :: lat_in        ! latitude  of input (unit: deg)
        integer,     intent(in   ) :: iopt          ! option to control the directon of transform
                                                    !   = 1 : regular lat/lon to rotated lat/lon
                                                    !   =-1 : rotated lat/lon to regular lat/lon

        real(dp),    intent(  out) :: lon_out       ! longitude of output (unit: deg)
        real(dp),    intent(  out) :: lat_out       ! latitude  of output (unit: deg)

!---- local variables
        REAL(DP) :: RLAT0, RLON0   ! latitude/longitude of center point
        REAL(DP) :: CLAT0, SLAT0   ! SIN/COS of latitude of center point
        REAL(DP) :: RLAT,  RLON    ! earth latitude/longitude
        REAL(DP) :: SLAT,  CLAT,  CLON
        REAL(DP) :: RLATR, RLONR   ! earth latitude/longitude
        REAL(DP) :: SLATR, CLATR, CLONR
        REAL(DP) :: HS
!-----------------------------------------------------------------------------!
!---- center point
        RLAT0=lat_sp_rll+90.0_dp        ! center point latitude =  south pole lat + 90
        RLON0=lon_sp_rll                ! center point longitude = south pole lon
        IF ( RLON0 < 0.0_dp ) RLON0=RLON0+360.0_dp
        CLAT0=COS(RLAT0*D2R)
        SLAT0=SIN(RLAT0*D2R)

        IF ( IOPT == 1 ) THEN           ! IOPT=1: Earth lat/lon ==> Rotated lat/lon
            RLAT=lat_in
            RLON=lon_in
            IF ( RLON < 0.0_dp ) RLON=RLON+360.0_dp
            HS=SIGN(1._dp,MOD(RLON-RLON0+180._dp+3600._dp,360._dp)-180._dp)
            CLON=COS((RLON-RLON0)*D2R)
            SLAT=SIN(RLAT*D2R)
            CLAT=COS(RLAT*D2R)
            SLATR=CLAT0*SLAT-SLAT0*CLAT*CLON
            IF(SLATR.LE.-1) THEN
                CLATR=0._dp
                RLONR=0.
                RLATR=-90.
            ELSEIF(SLATR.GE.1) THEN
                CLATR=0._dp
                RLONR=0.
                RLATR=90.
            ELSE
                CLATR=SQRT(1-SLATR**2)
                CLONR=(CLAT0*CLAT*CLON+SLAT0*SLAT)/CLATR
                CLONR=MIN(MAX(CLONR,-1._dp),1._dp)
                RLONR=HS*R2D*ACOS(CLONR)
                RLATR=R2D*ASIN(SLATR)
            ENDIF
            lat_out=RLATR
            lon_out=RLONR
        ELSEIF ( IOPT == -1 ) THEN      ! IOPT=-1: Rotated lat/lon ==> Earth lat/lon
            RLATR=lat_in
            RLONR=lon_in
            IF(RLONR > 180.0_dp) RLONR=RLONR-360.0_dp    ! in range (-180.0, 180.0)
            IF(RLONR <= 0._dp) THEN
               HS=-1.0_dp
            ELSE
               HS=1.0_dp
            ENDIF
            CLONR=DCOS(RLONR*D2R)
            SLATR=DSIN(RLATR*D2R)
            CLATR=DCOS(RLATR*D2R)
            SLAT=CLAT0*SLATR+SLAT0*CLATR*CLONR
            IF(SLAT.LE.-1._DP) THEN
                CLAT=0._DP
                CLON=DCOS(RLON0*D2R)
                RLON=0._DP
                RLAT=-90._DP
            ELSEIF(SLAT.GE.1) THEN
                CLAT=0._DP
                CLON=DCOS(RLON0*D2R)
                RLON=0._DP
                RLAT=90._DP
            ELSE
                CLAT=SQRT(1._DP-SLAT**2)
                CLON=(CLAT0*CLATR*CLONR-SLAT0*SLATR)/CLAT
                CLON=MIN(MAX(CLON,-1._dp),1._dp)
                RLON=REAL(MOD(RLON0+HS*R2D*DACOS(CLON)+3600._DP,360._dp))
                RLAT=REAL(R2D*DASIN(SLAT))
            ENDIF
            lat_out=RLAT
            lon_out=RLON
         ELSE
             WRITE(6,*) ' unrecognized option for opt, which must be either for regular to rotated (iopt=1) or vice versa (iopt=-1) '
             STOP 999
         ENDIF

         RETURN
!-----------------------------------------------------------------------
      end subroutine rll_trans_iplib
#endif
!-----------------------------------------------------------------------
!
!================================================================================
end module mod_rtma_regrid
