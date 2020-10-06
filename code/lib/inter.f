c      *************************************************************
c      * This package provides a general interpolaton routine      *
c      *************************************************************	

c     The main interface routines are:
c         get_index3,4 : get the grid indices for interpolation
c         int_index3,4 : interpolate to the grid position

c     -------------------------------------------------------------
c     Get index in grid space for interpolation
c     -------------------------------------------------------------

      subroutine get_index4 (rid,rjd,rkd,xpo,ypo,ppo,rtp,
     >                       vert0,vert1,surf0,surf1,mode,
     >                       nx,ny,nz,lonw,lats,dlon,dlat,misdat)

c     Purpose:
c        This subroutine determines the indices (rid,rjd,rkd) in grid 
c        space for a point in physical space (xpo,ypo,ppo). The 
c        horizontal grid is specified by the south-west point (lats,lonw)
c        and the grid spacing (dlat,dlon). The vertical grid is given
c        by <vert(n1,n2,n3)>. The lower boundary (typicall surface 
c        pressure) is given by <surf(n1,n2)>.
c     Arguments:
c        rid,rjd,rkd  real  output   grid location to be interpolated to
c        xpo,ypo,ppo  real  input    physical coordinates
c        rtp          real  input    relative time position (0=beginning, 1=end)
c        n1,n2,n3     int   input    grid dimensions in x-, y- and p-direction
c        lats,lonw    real  input    south and west boundary of grid space
c        vert         real  input    vertical coordinate grid
c        surf         real  input    lower boundary (surface pressure)
c        mode         int   input    direction of vertical axis (p=1,th=-1)
c                                        1: linear, 1 -> nz (th)
c                                        2: linear, nz -> 1 (pv)
c                                        3: binary (p)

      implicit none

c     Declartion of function parameters
      integer   nx,ny,nz
      real      xpo,ypo,ppo,rtp
      real      vert0(nx*ny*nz),vert1(nx*ny*nz)
      real      surf0(nx*ny)   ,surf1(nx*ny*nz)
      real      rid,rjd,rkd
      real      dlat,dlon,lats,lonw
      real      misdat
      integer   mode

c     Set numerical parameters
      real      eps
      parameter (eps=1.e-8)

c     Auxiliary variables
      real      rid0,rjd0,rkd0,rid1,rjd1,rkd1

c     Externals 
      real      int_time
      external  int_time

c     Get the inidices
      if (abs(rtp).lt.eps) then
         call  get_index3 (rid,rjd,rkd,xpo,ypo,ppo,mode,
     >                     vert0,surf0,nx,ny,nz,lonw,lats,dlon,dlat)
      elseif (abs(rtp-1.).lt.eps) then
         call  get_index3 (rid,rjd,rkd,xpo,ypo,ppo,mode,
     >                     vert1,surf1,nx,ny,nz,lonw,lats,dlon,dlat)
      else
         call  get_index3 (rid0,rjd0,rkd0,xpo,ypo,ppo,mode,
     >                     vert0,surf0,nx,ny,nz,lonw,lats,dlon,dlat)
         call  get_index3 (rid1,rjd1,rkd1,xpo,ypo,ppo,mode,
     >                     vert1,surf1,nx,ny,nz,lonw,lats,dlon,dlat)
         rid = int_time (rid0,rid1,rtp,misdat)
         rjd = int_time (rjd0,rjd1,rtp,misdat)
         rkd = int_time (rkd0,rkd1,rtp,misdat)

      endif

      end

c     -------------------------------------------------------------
c     Interpolate to an arbitrary position in grid space and time
c     -------------------------------------------------------------

      real function int_index4 (ar0,ar1,n1,n2,n3,rid,rjd,rkd,rtp,misdat)

c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid including a linear time-interpolation. 
c     Arguments:
c        rid,rjd,rkd  real  output   grid location to be interpolated to
c        xpo,ypo,ppo  real  input    physical coordinates
c        n1,n2,n3     int   input    grid dimensions in x-, y- and p-direction
c        lats,lonw    real  input    south and west boundary of grid space
c        vert         real  input    vertical coordinate grid
c        surf         real  input    lower boundary (surface pressure)

      implicit none

c     Declartion of function parameters
      integer   n1,n2,n3
      real      ar0(n1*n2*n3),ar1(n1*n2*n3)
      real      rid,rjd,rkd
      real      rtp
      real      misdat

c     Set numerical parameters
      real      eps
      parameter (eps=1.e-8)

c     Externals  
      real      int_index3,int_time
      external  int_index3,int_time

c     Auxiliary variables
      real      val0,val1,val

c     Do the 3d-interpolation
      if (abs(rtp).lt.eps) then
         val = int_index3 (ar0,n1,n2,n3,rid,rjd,rkd,misdat)
      elseif (abs(rtp-1.).lt.eps) then
         val = int_index3 (ar1,n1,n2,n3,rid,rjd,rkd,misdat)
      else
         val0 = int_index3 (ar0,n1,n2,n3,rid,rjd,rkd,misdat)
         val1 = int_index3 (ar1,n1,n2,n3,rid,rjd,rkd,misdat)
         val  = int_time (val0,val1,rtp,misdat)
      endif

c     Return value
      int_index4 = val

      return
      end


c     -------------------------------------------------------------
c     Interpolate to an arbitrary position in grid space
c     -------------------------------------------------------------

      real function int_index3 (ar,n1,n2,n3,rid,rjd,rkd,misdat)

c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The interpolation includes the 
c        testing of the missing data flag 'misdat'. If one dimension
c        is 1, a 2d-interpolation is performed; if two dimensions
c        are 1, it is a 1d-interpolation; if all three dimensions are
c        1, no interpolation is performed and the input value is
c        returned.
c     Arguments:
c        ar        real  input   input data array
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)

      implicit none

c     Declartion of function parameters 
      integer   n1,n2,n3
      real      ar(n1*n2*n3)
      real      rid,rjd,rkd
      real      misdat

c     Set numerical parameters
      real      eps
      parameter (eps=1.e-8)

c     Local variables
      integer   i,j,k,ip1,jp1,kp1
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k
      real      ri,rj,rk
      real      val000,val001,val010,val011,val100,val101,val110,val111
      real      frc000,frc001,frc010,frc011,frc100,frc101,frc110,frc111
      real      frc
      real      mdv
      real      val

c     Elementary test for dimensions
      if ( (n1.lt.1).or.(n2.lt.1).or.(n3.lt.1) ) then
         print*,'Invalid grid dimensions ',n1,n2,n3
         stop
      endif

c     Activate or inactive the missing data check (quick and dirty)
      if (misdat.ne.0.) then
         mdv = misdat
      else
         mdv = 257.22725394015
      endif

c     Bring the indices into the grid space
      ri = amax1(1.,amin1(float(n1),rid))
      rj = amax1(1.,amin1(float(n2),rjd))
      rk = amax1(1.,amin1(float(n3),rkd))

c     Get the index of the west-south-bottom corner of the box
      i   = min0(int(ri),n1-1)
      ip1 = i+1
      j   = min0(int(rj),n2-1)
      jp1 = j+1
      k   = min0(int(rk),n3-1)
      kp1 = k+1

c     Special handling for 2d arrays
      if (n3.eq.1) then
         k=1
         kp1=1
      endif

c     Get location relative to grid box
      if ( i.ne.ip1 ) then
         frac0i = ri-float(i)
         frac1i = 1.-frac0i
      else
         frac0i = 0.
         frac1i = 1.
      endif
      if ( j.ne.jp1 ) then
         frac0j = rj-float(j)
         frac1j = 1.-frac0j
      else
         frac0j = 0.
         frac1j = 1.
      endif
      if ( k.ne.kp1 ) then
         frac0k = rk-float(k)
         frac1k = 1.-frac0k
      else
         frac0k = 0.
         frac1k = 1.
      endif

c     On a grid point - take the grid point value 
      if ( ( abs(frac0i).lt.eps ).and.
     >     ( abs(frac0j).lt.eps ).and.
     >     ( abs(frac0k).lt.eps ) ) then
         
         val = ar( i + n1*(j -1) + n1*n2*(k -1) )
         goto 100
         
      endif

c     Init the fractions
      frc000 = frac1i * frac1j * frac1k
      frc001 = frac0i * frac1j * frac1k
      frc010 = frac1i * frac0j * frac1k
      frc011 = frac0i * frac0j * frac1k
      frc100 = frac1i * frac1j * frac0k
      frc101 = frac0i * frac1j * frac0k
      frc110 = frac1i * frac0j * frac0k
      frc111 = frac0i * frac0j * frac0k

c     Init the values
      val000 = ar( i   + n1*(j  -1) + n1*n2*(k  -1) )
      val001 = ar( ip1 + n1*(j  -1) + n1*n2*(k  -1) )
      val010 = ar( i   + n1*(jp1-1) + n1*n2*(k  -1) )
      val011 = ar( ip1 + n1*(jp1-1) + n1*n2*(k  -1) )
      val100 = ar( i   + n1*(j  -1) + n1*n2*(kp1-1) )
      val101 = ar( ip1 + n1*(j  -1) + n1*n2*(kp1-1) )
      val110 = ar( i   + n1*(jp1-1) + n1*n2*(kp1-1) )
      val111 = ar( ip1 + n1*(jp1-1) + n1*n2*(kp1-1) )

c     Handle missing data
      if ( abs(val000-mdv).lt.eps ) frc000 = 0.
      if ( abs(val001-mdv).lt.eps ) frc001 = 0.
      if ( abs(val010-mdv).lt.eps ) frc010 = 0.
      if ( abs(val011-mdv).lt.eps ) frc011 = 0.
      if ( abs(val100-mdv).lt.eps ) frc100 = 0.
      if ( abs(val101-mdv).lt.eps ) frc101 = 0.
      if ( abs(val110-mdv).lt.eps ) frc110 = 0.
      if ( abs(val111-mdv).lt.eps ) frc111 = 0.

c     Build the final value
      frc = frc000 + frc001 + frc010 + frc011 + 
     >      frc100 + frc101 + frc110 + frc111   
      if ( frc.gt.0. ) then
         val = 1./frc * ( frc000 * val000 + frc001 * val001 +
     >                    frc010 * val010 + frc011 * val011 +
     >                    frc100 * val100 + frc101 * val101 +
     >                    frc110 * val110 + frc111 * val111 )
      else
         val = misdat
      endif

c     Return the value 
 100  continue

      int_index3 = val

      end


c     -------------------------------------------------------------
c     Time interpolation
c     -------------------------------------------------------------

      real function int_time (val0,val1,reltpos,misdat)

c     Purpose:
c        This subroutine interpolates linearly in time between two
c        values.
c     Arguments:
c        val0      real  input   value at time 0
c        val1      real  input   value at time 1
c        reltpos   real  input   relative time (between 0 and 1)
c        misdat    real  input   missing data flag (on if misdat<>0)

      implicit none

c     Declaration of parameters
      real      val0
      real      val1
      real      reltpos
      real      misdat

c     Numerical epsilon
      real      eps
      parameter (eps=1.e-8)

c     Local variables
      real      val
      real      mdv

c     Activate or inactive the missing data check (quick and dirty)
      if (misdat.ne.0.) then
         mdv = misdat
      else
         mdv = 257.22725394015
      endif

c     Do the linear interpolation
      if ( abs(reltpos).lt.eps ) then
         val = val0
      elseif ( abs(reltpos-1.).lt.eps ) then
         val = val1
      elseif ( (abs(val0-mdv).gt.eps).and.
     >         (abs(val1-mdv).gt.eps) ) then
         val = (1.-reltpos)*val0+reltpos*val1
      else
         val = mdv
      endif

c     Return value
      int_time = val

      end


c     -------------------------------------------------------------
c     Get the position of a physical point in grid space
c     -------------------------------------------------------------

      subroutine get_index3 (rid,rjd,rkd,xpo,ypo,ppo,mode,
     >                       vert,surf,nx,ny,nz,lonw,lats,dlon,dlat)

c     Purpose:
c        This subroutine determines the indices (rid,rjd,rkd) in grid 
c        space for a point in physical space (xpo,ypo,ppo). The 
c        horizontal grid is specified by the south-west point (lats,lonw)
c        and the grid spacing (dlat,dlon). The vertical grid is given
c        by <vert(n1,n2,n3)>. The lower boundary (typicall surface 
c        pressure) is given by <surf(n1,n2)>.
c     Arguments:
c        rid,rjd,rkd  real  output   grid location to be interpolated to
c        xpo,ypo,ppo  real  input    physical coordinates
c        n1,n2,n3     int   input    grid dimensions in x-, y- and p-direction
c        lats,lonw    real  input    south and west boundary of grid space
c        vert         real  input    vertical coordinate grid
c        surf         real  input    lower boundary (surface pressure)
c        mode         int   input    direction of vertical axis 
c                                        1: linear, 1 -> nz (th)
c                                        2: linear, nz -> 1 (pv)
c                                        3: binary (z)
c                                        4: binary (p)
 

      implicit none

c     Declartion of function parameters
      integer   nx,ny,nz
      real      vert(nx*ny*nz)
      real      surf(nx*ny)
      real      rid,rjd,rkd
      real      xpo,ypo,ppo
      real      dlat,dlon,lats,lonw
      integer   mode

c     Numerical epsilon
      real      eps
      parameter (eps=1.e-8)

c     Local variables
      integer   i,j,k
      real      ppo0,ppo1,ppom,psur
      integer   i0,im,i1
      
c     Externals 
      real      int_index3
      external  int_index3

c     Get the horizontal grid indices
      rid=(xpo-lonw)/dlon+1.
      rjd=(ypo-lats)/dlat+1.

c     Two-dimensional interpolation on horizontal plane: return
      if ( nz.eq.1 ) then
         rkd = 1.
         goto 100
      endif
         
c     Lowest-level interpolation: return
      if ( abs(ppo).lt.eps ) then
         rkd = 1.
         goto 100
      endif

c     Get the pressure at the lowest level and at the surface 
      ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real(1),0.)
      psur = int_index3(surf,nx,ny, 1,rid,rjd,real(1),0.)

c     The point is between the surface and the lowest level: return
      if ( (ppo.ge.ppo0).and.(ppo.le.psur).or.
     >     (ppo.le.ppo0).and.(ppo.ge.psur) )
     >then 
         psur = int_index3(surf,nx,ny, 1,rid,rjd,real(1),0.)
c         rkd  = (psur-ppo)/(psur-ppo0)
          rkd = 1.
         goto 100
      endif

c     Full-level search (TH): linear ascending scanning through all levels
      if ( mode.eq.1 ) then
        
         ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real(1),0.)
         rkd=0
         do i=1,nz-1
            ppo1 = int_index3(vert,nx,ny,nz,rid,rjd,real(i+1),0.)
            if ( (ppo0.lt.ppo).and.(ppo1.ge.ppo) ) then
               rkd=real(i)+(ppo0-ppo)/(ppo0-ppo1)
               goto 100
            endif
            ppo0 = ppo1
         enddo

c     Full-level search (PV): linear descending scanning through all levels
      elseif ( mode.eq.2 ) then

         ppo1 = int_index3(vert,nx,ny,nz,rid,rjd,real(nz),0.)
         rkd=0
         do i=nz-1,1,-1
            ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real(i),0.)
            if ( (ppo1.gt.ppo).and.(ppo0.le.ppo) ) then
               rkd=real(i)+(ppo0-ppo)/(ppo0-ppo1)
               goto 100
            endif
            ppo1 = ppo0
         enddo

c     Full-level search (Z):  binary search
      elseif ( mode.eq.3 ) then

         rkd  = 0
         i0   = 1
         i1   = nz
         ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real( 1),0.)
         ppo1 = int_index3(vert,nx,ny,nz,rid,rjd,real(nz),0.)
         
         do while ( i1.gt.(i0+1) )
            im   = (i0+i1)/2
            ppom = int_index3(vert,nx,ny,nz,rid,rjd,real(im),0.)
            if (ppom.gt.ppo) then
               i1   = im
               ppo1 = ppom
            else
               i0   = im
               ppo0 = ppom
            endif

         enddo
            
         rkd=real(i0)+(ppo0-ppo)/(ppo0-ppo1)

c     Full-level search (P):  binary search
      elseif ( mode.eq.4 ) then

         rkd  = 0
         i0   = 1
         i1   = nz
         ppo0 = int_index3(vert,nx,ny,nz,rid,rjd,real( 1),0.)
         ppo1 = int_index3(vert,nx,ny,nz,rid,rjd,real(nz),0.)
         
         do while ( i1.gt.(i0+1) )
            im   = (i0+i1)/2
            ppom = int_index3(vert,nx,ny,nz,rid,rjd,real(im),0.)
            if (ppom.lt.ppo) then
               i1   = im
               ppo1 = ppom
            else
               i0   = im
               ppo0 = ppom
            endif

         enddo
            
         rkd=real(i0)+(ppo0-ppo)/(ppo0-ppo1)

      endif

c     Exit point for subroutine
 100  continue

      end

