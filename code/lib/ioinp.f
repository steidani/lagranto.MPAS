c     ------------------------------------------------------------
c     Open input file
c     ------------------------------------------------------------

      subroutine input_open (fid,filename)

      use netcdf

      implicit none

c     Declaration of subroutine parameters
      integer      fid              ! File identifier
      character*80 filename         ! Filename

c     Declaration of auxiliary variables
      integer      ierr

      ierr = NF90_OPEN(TRIM(filename),nf90_nowrite, fid)
      IF ( ierr /= nf90_NoErr ) PRINT *,NF90_STRERROR(ierr)

      end


c     ------------------------------------------------------------
c     Read information about the grid and wind
c     ------------------------------------------------------------
      
      subroutine input_grid 
     >                   (fid,fieldname,datestr,f3,p3,ps,
     >                    mdv,xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,mode) 

c     mode = 1: grid parameters only
c     mode = 2: f3, p3 and ps (field + grid)

      use netcdf

      implicit none

c     Declaration of subroutine parameters 
      integer      fid                 ! File identifier
      real         xmin,xmax,ymin,ymax ! Domain size
      real         dx,dy               ! Horizontal resolution
      integer      nx,ny,nz            ! Grid dimensions
      real         p3(nx,ny,nz)        ! Staggered levels
      real         f3(nx,ny,nz)        ! Meteorological field
      real         ps(nx,ny)           ! Surface pressure / topography
      character*80 datestr             ! Date/time to be read from file
      character*80 fieldname           ! Variable from which to take grid info
      integer      mode                ! Reading mode
      real         mdv                 ! Missing data value

c     Numerical epsilon
      real         eps
      parameter    (eps=0.0001)

c     Flags for grid trafos
      integer      vertical_swap
      integer      closear

c     Auxiliary varaibles
      integer      ierr,stat
      integer      i,j,k
      character*80 varname
      character*80 newname
      real         lon(2000),lat(2000)
      real         lev(2000)
      integer      vardim(4)
      integer      varid,lonid,latid,levid,timid
      integer      nxf,nyf,nzf
      integer      dimids (1000)
      character*80 dimname(1000)
      real,allocatable, dimension (:,:,:) :: f3_tmp,p3_tmp
      real,allocatable, dimension (:,:  ) :: f2_tmp,ps_tmp
      real         delta,dlon
      integer      ndim
      integer      itime
      character*80 hour_req

c     --- Convert date string (YYYYMMDD_HH) into netCDF time info (HH) --------------------

      read(datestr(9:10),  *) hour_req

c     --- Read grid info -----------------------------------------------

c     lon (longitude)
      ierr = nf90_inq_dimid(fid,'lon', lonid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_dimension(fid, lonid, len = vardim(1))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_INQ_VARID(fid,'lon',varid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_GET_VAR(fid,varid,lon(1:vardim(1)))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)

c     lat (latitude)
      ierr = nf90_inq_dimid(fid,'lat', latid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_dimension(fid, latid, len = vardim(2))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_INQ_VARID(fid,'lat',varid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_GET_VAR(fid,varid,lat(1:vardim(2)))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
        
c     lev (vertical level)
      ierr = nf90_inq_dimid(fid,'lev', levid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_dimension(fid, levid, len = vardim(3))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_INQ_VARID(fid,'lev',varid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_GET_VAR(fid,varid,lev(1:vardim(3)))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)

c     time
      ierr = nf90_inq_dimid(fid,'time', timid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_dimension(fid, timid, len = vardim(4))
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)

c     --- Allocate memory for fields -----------------------------------

      allocate(p3_tmp(vardim(1),vardim(2),vardim(3)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3_tmp***'
      allocate(f3_tmp(vardim(1),vardim(2),vardim(3)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array f3_tmp***'
      allocate(f2_tmp(vardim(1),vardim(2)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array f2_tmp***'
      allocate(ps_tmp(vardim(1),vardim(2)),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ps_tmp***'

c     --- Decide on the grid structure ---------------------------------      
      
c     Check whether field must be closed 
      delta = lon(vardim(1)) - lon(1) - 360.
      dlon  = ( lon(vardim(1)) - lon(1) ) / real( vardim(1) - 1 )
      if (abs(delta+dlon).lt.eps) then
          closear = 1
      else
          closear = 0
      endif   
       
c     Decide whether fields need to be vertically swaped (index 1 = near surface)   
      vertical_swap = 0
      if ( lev(1).lt.lev(vardim(3)) ) then
        vertical_swap = 1
      endif

c     --- For mode 1, return grid parameters; nothing else to do -------

      if ( mode.eq.1 ) then
      
        nx     = vardim(1)
        ny     = vardim(2)
        nz     = vardim(3)
        xmin   = lon(1)
        xmax   = lon(nx)
        ymin   = lat(1)
        ymax   = lat(ny)
        dx     = (xmax-xmin)/real(nx-1)
        dy     = (ymax-ymin)/real(ny-1)
        mdv    = -999.
        
        if ( closear.eq.1 ) then
            nx   = nx + 1
            xmax = xmax + dx
        endif
        
        goto 100
        
      endif

c     --- Get time index -----------------------------------------------
      !!! Dani: read only wanted time step?
      !!! Better: compare directly datestr with time(i): create datestring in file with CDO?
      itime = -1
      
      if (hour_req.eq.'00') then
	itime = 1
	goto 10
      elseif (hour_req.eq.'06') then
	itime = 2
	goto 10
      elseif (hour_req.eq.'12') then
	itime = 3
	goto 10
      elseif (hour_req.eq.'18') then
	itime = 4
	goto 10
      endif

      if ( itime.eq.-1) then
        print*,' ERROR: time not found on netCDF file'
        print*,datestr
        stop
      endif
 10   continue

c     --- Get grid (p3+ps) ---------------------------------------------

c     Read 3D grid (z grid: height (in m) of pressure levels)
      varname = 'height'
      ierr = NF90_INQ_VARID(fid,varname,varid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = NF90_GET_VAR(fid,varid,p3_tmp,
     >           start = (/ 1,1,1,itime /),
     >           count = (/ vardim(1),vardim(2),vardim(3),1 /) )
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      
c     Set the topography - lowest level of p3_tmp
      do i=1,nx
        do j=1,ny
           ps_tmp(i,j) = p3_tmp(i,j,1)
        enddo
      enddo
    
c     --- Read meteorological wind field (f3) -------------------------------

c     Set variable name on MPAS netCDF
      if ( fieldname.eq.'U' ) then
         varname = 'uzonal'
      elseif ( fieldname.eq.'V' ) then
         varname = 'umeridional'
      elseif ( fieldname.eq.'OMEGA' ) then
         varname = 'w'
      else
         varname = fieldname
      endif

c     Read meta-info for variable (dimensions)
      ierr = NF90_INQ_VARID(fid,varname,varid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_variable(fid, varid, ndims  = ndim)
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      ierr = nf90_inquire_variable(fid, varid,
     >                                   dimids = dimids(1:ndim))
      IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      do i=1,ndim
           ierr = nf90_inquire_dimension(fid, dimids(i),
     >                                           name = dimname(i) )
           IF(ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      enddo
      
c     Check whether dimensions are OK
      if ( (ndim.ne.3).and.(ndim.ne.4).or. 
     >     ( dimname(1).ne.'lon'  ).or.       
     >     ( dimname(2).ne.'lat'  ).or. 
     >     ( dimname(3).ne.'lev'  ).and.(ndim.eq.4).or.
     >     ( dimname(3).ne.'time' ).and.(ndim.eq.3).or.
     >     ( dimname(4).ne.'time' ).and.(ndim.eq.4) ) 
     >then
        print*,' ERROR: dimension mismatch on netCDF file...'
        stop
      endif
      
c     Read data from file - in case of 2D, blow it up to 3D     
      if ( ndim.eq.3 ) then
        ierr = NF90_GET_VAR(fid,varid,f2_tmp,
     >           start = (/ 1,1,itime /),
     >           count = (/ vardim(1),vardim(2),1 /) )
        IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
        do i=1,vardim(1)
          do j=1,vardim(2)
            do k=1,vardim(3)
               f3_tmp(i,j,k) = f2_tmp(i,j)
            enddo
          enddo
        enddo
      else
        ierr = NF90_GET_VAR(fid,varid,f3_tmp,
     >           start = (/ 1,1,1,itime /),
     >           count = (/ vardim(1),vardim(2),vardim(3),1 /) )
        IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
      endif
      
c     ---- Set the final output fields ---------------------------------

c     Copy field to output arrays
      do i=1,vardim(1)
        do j=1,vardim(2)        
            do k=1,vardim(3)
              if ( vertical_swap.eq.1 ) then
                  p3(i,j,k) = p3_tmp(i,j,vardim(3)-k+1)
                  f3(i,j,k) = f3_tmp(i,j,vardim(3)-k+1)
               else
                  p3(i,j,k) = p3_tmp(i,j,k)
                  f3(i,j,k) = f3_tmp(i,j,k)
               endif
            enddo
            ps(i,j) = ps_tmp(i,j)
        enddo
      enddo
      
c     Close arrays
      if ( closear.eq.1 ) then
        do j=1,vardim(2)
          do k=1,vardim(3)
             p3(vardim(1)+1,j,k) = p3(1,j,k)
             f3(vardim(1)+1,j,k) = f3(1,j,k)
          enddo
          ps(vardim(1)+1,j) = ps(1,j)
        enddo
      endif
       
c     --- Deallocate memory --------------------------------------------      
      
 100  continue      
      deallocate(p3_tmp,stat=stat)
      deallocate(f3_tmp,stat=stat)
      deallocate(f2_tmp,stat=stat)
      deallocate(ps_tmp,stat=stat)
      
      end

c     ------------------------------------------------------------
c     Close input file
c     ------------------------------------------------------------

      subroutine input_close(fid)

c     Close the input file with file identifier <fid>.

      use netcdf

      implicit none

c     Declaration of subroutine parameters
      integer fid

c     Auxiliary variables
      integer ierr

      ierr = NF90_CLOSE(fid)
      IF( ierr /= nf90_NoErr) PRINT *,NF90_STRERROR(ierr)
 
      end

