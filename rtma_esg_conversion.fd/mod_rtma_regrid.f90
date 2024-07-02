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
!-----------------------------------------------------------------------
!
!================================================================================
end module mod_rtma_regrid
