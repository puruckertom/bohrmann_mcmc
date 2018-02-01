! Array Visualizer code added...search on 'fagl'

C   Code type: FORTRAN subroutine

C   Compiler: Fortran77

C   Programmed by: C. Unkrich

C   Date: 4/95

C   Revised: 7/95 mods for microtopography by R.E. Smith
C!! slight revision in calculating depth to send to infilt routine 2/12/02

C    Laminar flow added 8/2005 by C. Unkrich - laminar flow is used up to a
C    transition depth if a 'LAM' tag & value appear in the input block

C   Description:

C     Simulates overland flow using a four point implicit finite-difference
C     approximation to the kinematic wave equation written for one-dimensional
C     unsteady flow on a sloping plane. Rainfall is modeled as a spatially
C     uniform time, intensity step function. The infiltration function is
C     coupled to the routing equations to provide a realistic interaction
C     between rainfall, runoff, and infiltration. Modifed for micro-topograpy,
C     using hydraulic radius rather than mean depth. An algorithm has been
C     added which monitors the wave celerity and adjusts the computational
C     time increment to maintain numerical accuracy.

C     For output, volumetric rainfall rate from the upstream area for each
C     time step is added to the volumetric rate over the plane and stored.
C     Dividing by contributing area gives the integrated rainfall rate over
C     the area represented by the plane and its upstream contributors.
C------------------------------------------------------------------------------

C   Arguments:

C     dtm(4)        real        smallest dt to satisfy the Courant condition
C                               during each quarter of the simulation.
C
C     xnu           real        kinematic viscosity (m**2/s or ft**2/s)
C -----------------------------------------------------------------------------
C  for other variables see module definition file
C------------------------------------------------------------------------------

C   Input Block (input file labels in parenthesis):

C     id (ID)          int      identifier,

C     xleng (L)        real     length (m or ft),

C     width (W)        real     width (m or ft),

C     slope (SL)       real     slope,

C     rs (MA or CH)    real     Manning or Chezy discharge coefficient,

C     x (X)            real     x coordinate,

C     y (Y)            real     y coordinate,

C     dintr (IN)       real     interception depth (mm or in),

C     cover (CA)       real     canopy cover fraction,

C     amr (RE)         real     average microtopographic relief (mm or in),

C     ams (SPA)        real     average microtopographic spacing (m or ft),

C     fmin (KS)        real     saturated hydraulic conductivity (mm or in/hr),
C------------------------------------------------------------------------------

C   Subroutines/functions:

C     interp        (rain.for)
C     infilt        (infilt.for)
C     infil0             "
C     iter          (miscel.for)
C     sed0          (kinsed.for)
C     kinsed             "
C     sedfin             "
C     errh2         (plane.for)
C     erpaub             "
C     qfunc              "
C     wprf               "
C     zfu                "
C     dqdaf              "
C     progss        (K2RunW.for)
C     getr4         (reader.for)
C     getstr             "
C     new           (miscel.for)
C     old                "
C     get                "
C     store              "
C     stora              "
C     geta               "
C     qwrt          (writer.for)
C     errxit        (K2RunW.for)
C------------------------------------------------------------------------------


      subroutine plane (dtm, xnu)
     & !qbal, delt, limit, clen, units, sed, diag, cour,                 

      use runpars
      use elpars
      use multip
      use itpars

C**      USE AVDEF
C                                                                     arguments
C------------------------------------------------------------------------------

      real,dimension(4) ::  dtm
      real              ::  xnu

!      logical sed, diag, cour

!      integer limit, units
!      dimension qbal(10)

C                                                               local variables
C------------------------------------------------------------------------------

      logical :: up, unif, laminar

      integer :: nd, ierr, nk, j, jp1, ms, upq, outq(2), i, k, idup,
     &        inflow(12), nu, upr, outr, ngage, nzro=0!, nwon=1 !, id

      character(LEN=18) :: msg
      character(LEN=6) :: trace, source
C** new units character  22/4/03
      character(LEN=2) :: chl

      common /plane1/ dt, dx, j, jp1, qlav, alpha(20), p(20), h1(20), 
     &                h2(20), q1(20), q2(20), h2c, ams, wpc, c1, c2, bw, 
     &                qb, omega
C
      real,dimension(10,2) :: qup
      real,dimension(20,2) :: ql
      real,dimension(20) ::  topw, slp, fav, qt, wid, afac, h1a
      real,dimension(1000) :: t, z, r, zn
      real, parameter :: dtmin = 0.
      real :: zro = 0.
      real :: htrans, qtrans

C**      INTEGER            :: fagl_delay =  900 ! milliseconds, must be < 1000
C**      INTEGER            :: fagl_start =    0 ! minutes, when to open AV window
C**      INTEGER            :: fagl_limit =  400 ! minutes, when to close AV window
C**      INTEGER            :: fagl_id    =   -1 ! element to display (0 = all, -1 = none)
C**      INTEGER            :: fagl_vals(8)
C**      INTEGER            :: fagl_s0
C**      INTEGER            :: fagl_ds
C**      CHARACTER(LEN=100) :: fagl_str
C**      LOGICAL            :: fagl_flag

C------------------------------------------------------------------------------

      external errh2

      external erpaub

      data topw /20*1./
      data qt /20*0./
      data qup /20*0./
C      data msg /'       (PLANE)'/,trace/'PLANE '/

      omega = 0.8        !  0.8 to match chan
      tol = 1.e-7
      ltyp = 0
      pof = conv*3600.
C                                                          get input parameters
C------------------------------------------------------------------------------

C                                                                    IDENTIFIER
      call geti4 ('ID', 0, id, ierr)

      if (ierr .eq. 3) then

        call errxit ('PLANE', 'invalid identifier (ID)')

      else if (ierr .ge. 1) then

        call errxit ('PLANE', 'missing identifier (ID)')

      end if

      if(DIAG) write(99,'(" Diagnostics for plane ",i3)')ID

      msg(1:16) = typname(0)
      write (msg(15:18), '(i4)') id
C
      do ms = 1, 5

        if (msg(ms:ms) .ne. ' ') go to 1

      end do

1     if(nele .gt. 5) call progss (msg(ms:18))

C                                                                        LENGTH
      call getr4 ('L', 0, xleng, ierr)

      if (ierr .gt. 0) call errxit (msg(ms:14), 'length (L) not found')

      nk = max1 ((15. * xleng / clen), 5.)
C                                                       number of spatial nodes
      if (nk .gt. 15) nk = 15
C                                                            distance increment
      dx = xleng / (float(nk) - 1.)

C** 22/4/03  dx check
      dxlim = dxlmt
      chl = 'm.'
      if(units .ge. 2) then
       dxlim = dxlmt/0.3
       chl = 'ft'
      end if
      if(dx .ge. dxlim) then
        write(files(3), 1990) id, dx, chl!  , dxlim
 1990 format(/"  Plane ",i3,":  based on length and parameter CLEN,"/
     &"  the numerical increment of ",f6.1,1x,a2," is too large for",/
     &"  realistic numerical solution of the flow equation, and may "/
     &"  give  misleading results."/)
      end if
      
      If(diag) then
        write(99,393) xleng,dlab,nk,dx,dlab
 393    format(' Plane is ',f8.1,1x,a2,' long with ',I2,' nodes and DX o
     &f',f8.4,1x,a2/)
      end if
 !     
C                                                                         WIDTH
      call getr4 ('W', 0, width, ierr)

      if (ierr .gt. 0) then
C                                                          area (sq.m or sq.ft)
        call getr4 ('AR', 0, area, ierr)

        if (ierr .gt. 0) then

          if (units .eq. 1) then
C                                                               area (hectares)
            call getr4 ('HA', 0, area, ierr)

            if (ierr .gt. 0) call errxit
     &                       (msg(ms:14), 'width (W) or area not found')

            area = area * 10000.

          else
C                                                                  area (acres)
            call getr4 ('AC', 0, area, ierr)

            if (ierr .gt. 0) call errxit
     &                       (msg(ms:14), 'width (W) or area not found')

            area = area * 43560.

          end if
        end if

        width = area / xleng

      else

        area = xleng * width

      end if

C                                                                         SLOPE
      call getr4 ('SL', 0, slope, ierr)

      if (ierr .gt. 0) call errxit (msg(ms:14), 'slope (SL) not found')

C                                                                 LAMINAR FLOW?
      call getr4 ('LAM', 0, rlam, ierr)

      if (ierr .eq. 0) then

        if (xnu .eq. 0.) call errxit (msg(ms:14),
     &    'laminar flow option requires temperature in global block')

        if (units .eq. 1) then
          grav = 9.81
        else
          grav = 32.2
        end if

!       DARCY-WEISBACH LAMINAR FLOW

        alam = slope * 8. * grav / (xnu * rlam)
        plam = 3.

        laminar = .TRUE.

      else

        laminar = .FALSE.

      end if
C                                                                      MANNING?
      call getr4 ('MA', 0, rs, ierr)

      if (ierr .eq. 0) then

        pturb = 1.66667
        fac = 1.
        if (units .eq. 2) fac = 1.49
        aturb = fac * sqrt (slope) / (rs * rmn)

	  if (laminar) then

!         TRANSITON POINT FOR MANNING-LAMINAR

          htrans = (fac*xnu*rlam/(8.*grav*rs*slope**0.5))**0.75
          qtrans = (fac*slope**0.5*htrans**(5./3.))/rs

	  end if

      else
C                                                                        CHEZY?
        call getr4 ('CH', 0, rs, ierr)

        if (ierr .eq. 0) then

          pturb = 1.5
          aturb = rs * rmn * sqrt (slope)

	    if (laminar) then

!           TRANSITION POINTS FOR CHEZY-LAMINAR

            htrans = (rs*xnu*rlam/(8.*grav * slope**0.5))**(2./3.)
            qtrans = rs*slope**0.5*htrans**1.5

	    end if

        else

          call errxit
     &    (msg(ms:14), 'hydraulic roughness (MA or CH) not found')

        end if
      end if
C                                                                      RAINGAGE
      ngage = 0

      call geti4 ('RA', 0, ngage, ierr)

      if (ierr .gt. 0) then
C                                                          Get X, Y coordinates
        call getr4 ('X', 0, x, ierr)
        if (ierr .gt. 0) x = 0.

        call getr4 ('Y', 0, y, ierr)
        if (ierr .gt. 0) y = 0.

      end if
C                                                                  INTERCEPTION
      call getr4 ('IN', 0, dintr, ierr)

      if (ierr .gt. 0) dintr = 0.

      if (dintr .gt. 0.) then

        dintr = rmin * dintr
C                                                                  CANOPY COVER
        call getr4 ('CA', 0, cover, ierr)
C!! change .eq. to .le.
        if (ierr .gt. 0 .or. cover .le. 0.)
     &  call errxit (msg(ms:14),'canopy cover fraction (CA) not found')
C!! check for out of bounds:
        if(cover .gt. 1.0) then
          call errxit (msg(ms:14),' canopy cover fraction cannot be more
     & than 1.0')
        end if

      end if
C                                                                        RELIEF
      amr = 0.
      ams = 0.

      call getr4 ('RE', 0, amr, ierr)

      if (ierr .eq. 0) then

	  if (laminar)  call errxit (msg(ms:14),
     &    'cannot specify microtopography with laminar flow option')
C                                                                       SPACING
        call getr4 ('SPA', 0, ams, ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), 'average microtopographic spacing (SPA) not found')

        amr = amr / conv

      end if

      up = .false.
      coarea = 0.
      nu = 0
      upr = 0
C                                                     upstream plane identifier
      call geti4 ('UP', 0, idup, ierr)

      if (ierr .eq. 0) then

        up = .true.
        nu = 1

      end if
C                                              get indices to storage locations
C                                              for rainfall, inflow and outflow
C------------------------------------------------------------------------------
      call new (id, 'RA', outr)

      if (up) then

        call old (idup, 'RA', upr, ierr)
        call old (idup, 'RO', upq, ierr)
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), 'upstream inflow not found')
        inflow(1) = idup

      end if

      call new (id, 'RO', outq(1))
C                                                      update contributing area
C------------------------------------------------------------------------------
      coarea = 0.
      if (up) call geta (idup, coarea)
      coarea = coarea + area
      call stora (id, coarea)

C                                              interpolate rainfall intensities
C------------------------------------------------------------------------------

      call interp (x, y, t, z, nd, sat, ngage)

      d1 = 0.
      t1 = 0.
      di = 0.
      dinact = 0.
      dintr = dintr / conv   ! convert to m.
      tfinl = float (limit - 1) * delt   !this is in seconds
      itlim = limit
      i = 0
      ndnew = nd
C!!
      zn(1) = 0.
      depth = z(nd)
      z(nd+1) = z(nd)

      do j = 1, nd

        i = i + 1

        if (i .eq. ndnew) go to 10
C                                              convert depth z into intensity r
        t2 = t(i+1)
        dt = t2 - t1
        d2 = z(i+1)
        dd = d2 - d1

        if (dintr .gt. d1) then
C                                            reduce intensity by cover fraction
          if (d2 .gt. dintr) then
C                                               reduction over partial interval
            r(i) = dd * (1. - cover) / dt
C!!  intermediate time
            adt = dt*(dintr - d1)/dd
C                                                             insert breakpoint
            r(i+1) = dd / dt
C                                                reset rain array indices
            do k = nd, i+1, -1
              t(k+1) = t(k)
              z(k+1) = z(k)
            end do
C!!  new 10/03  find array zn which is the z array less interception
            z(i+1) = z(i) + dd*adt/dt
C            
            zn(i+1) = zn(i) + r(i)*adt
            t(i+1) = t(i) + adt 
            i = i + 1
            zn(i+1) = zn(i) + r(i)*(dt-adt)
            ndnew = nd + 1
            di = dintr

          else

            r(i) = dd * (1. - cover) / dt
            di = d2
C!!  new 10/03
            zn(i+1) = zn(i) + r(i)*dt

          end if

        else
C                                                               no interception
          r(i) = dd / dt
          zn(i+1) = zn(i) + r(i)*dt    !! also new 10/03

        end if

        if (tfinl .gt. t1 .and. tfinl .le. t2) then
C                                                  compute final rainfall depth
          dtfin = (t2 - tfinl) / dt
          depth = d2 - dtfin * dd

          if (dintr .gt. 0.) dinact = di * cover

        end if

        t1 = t2
        d1 = d2

      end do

10    continue
      if(ndnew .lt. 200) then
        t(ndnew+1) = tfinl + 1.
        r(ndnew+1) = 0.
      end if
      
C                                                       final depth intercepted
      dintr = dinact
C                                                 compute integrated rain rates
C------------------------------------------------------------------------------
      qrmax = 0.
      k = 1
      qrf = r(1) * area
C!!  new 9/03
      lmn = limit - 1
      zlr = zn(1)
C!!
      do i = 1, lmn

C!! replaced 9/03:
C        if (time .gt. t(k+1)) then
C          k = k + 1
C          if(k .gt. ndnew) then   
C            stop ' rainfall undefined at long time '
CC            r(k) = 0.
CC            ndnew = k
C          end if
C          qrf = r(k) * area
C        end if
C!! new code for qrf to account for rainsteps less than DELT:
        timen = delt * float (i)   ! seconds
        do while(t(k) .lt. timen)
          k = k+1
        end do
C interpolate cumulative depth at time: 
C!!  use zn rather than z to account for interception  (10/03)
        km = k-1
        z2 = zn(km) + (zn(k) - zn(km))*(timen - t(km))/(t(k) - t(km))
        qrf = (z2 - zlr)/DELT * area
        zlr = z2
!! end new code
        qrain = qrf

        if (upr .gt. 0) then
          call get (upr, i, qr)
          qrain = qrain + qr
        end if

        call store (outr, i, qrain)
        if (qrain .gt. qrmax) qrmax = qrain

      end do
C!!  add value for i=limit (maybe unnecessary)
      call store (outr, limit, zro)

      do j = 1, 15
        qbal(j) = 0.
      end do
C                                                 plane are in m2 or ft2
      qbal(1) = area                  
C                                                               rainfall volume
      qbal(2) = area * depth
C                                                            volume intercepted
      qbal(8) = area * dintr

C                                     read / initiallize infiltration variables
C------------------------------------------------------------------------------

      call infil0 (id, 1, nk, msg(ms:14), diag, sat)

      if (sed) then
C                                                              ... and sediment
C------------------------------------------------------------------------------

        do j = 1, nk

          slp(j) = slope
          wid(j) = width

        end do

        call sed0 ( 1, id, msg(ms:14), inflow, nu, nzro,
     &             nk, omega, dx, slp, width)

      end if

C                        compute geometric variables related to microtopography
C------------------------------------------------------------------------------

      bf = 0.1

      if (amr .gt. 0.) then         ! relief height = amr
C
        if(ams .le. 0.) then
          stop ' cannot have zero spacing for relief '
        end if
        bm = .5 * (1. - bf) * ams
        ss = amr / bm
        cpz = sqrt(1. + ss * ss)
        bw = bf * ams
        ac = amr * (bm + bw)
        wpc = bw + cpz * 2. * bm ! critical wetted perim for one rill
        rmc = ac / wpc
        h2c = ac / ams
        uqc = alpha(1) * rmc**(p(1) - 1.) * ac
        wz = bf
        fnw = ams / width
        c1 = 2. / ss
        c2 = c1 * cpz

      else

        h2c = 0.
        uqc = 0.
        bw = 1.
        wpc = 1.
        ams = 1.
        wz = 1.

      end if
C                                                   initiallize local variables
C------------------------------------------------------------------------------
C                                                                    Q=0 at t=0
      do j = 1, nk
        h1(j) = 0.
        h2(j) = 0.
        q1(j) = 0.
        q2(j) = 0.
        ql(j,1) = 0.
        fav(j) = r(1)
        topw(j) = wz
      end do

	if (laminar) then
        do j = 1, nk
          p(j)     = plam
          alpha(j) = alam
          end do
	else
        do j = 1, nk
          p(j)     = pturb
          alpha(j) = aturb
        end do
	end if

      time = 0.
      xc = 0.
      k = 1
      vin = 0.
      unif = .true.
      nc = 1
      sumrf = 0.

      p1 = limit / 4
      p2 = limit / 2
      p3 = 3 * limit / 4
      qub1 = 0.                           ! default
C!!  new  RES      
      if (up) then  
        call get (upq, 1, qub1)
      end if
      vin = 0.5 * qub1
      q1s = qub1 / width
C estimate first a2u
      a2u = (q1s * ams / alpha(1) * wpc**(p(1) - 1.))**(1. / p(1))
C      h1(1) = a2u/ams
C!! end new      
C!! new location so that upstream record is not written over
      call store (outq(1), 1, zro)
      qout1 = 0.
      qp = 0.
      tp = 0.
      qsp = 0.
      tsp = 0.
      wsi = 0.
      wsmax = 0.
      vcm = 0.
      afac = 1.0
      itot = 0
      smql = 0.
C!! added balance component for upstream info.  RES 03/12/02
      smup = 0.

C**      fagl_flag = .FALSE.

C                                 start fixed time step loop (i, delt, step)...
C------------------------------------------------------------------------------

      do i = 2, limit

        itm = i
        qub2 = 0.

        dvf_delt = 0.
        rf_delt = 0.
C                                                                        inflow
        if (up) call get (upq, i, qub2)

        dqub = qub2 - qub1
        qup(1,1) = qub2
        quba = 0.5*(qub2 + qub1)
        step = float (i - 1) * delt
        dtp = delt

C**        IF( fagl_id .NE. -1 )THEN
C**          IF( .NOT. fagl_flag .AND. INT(step/60.) .LT. fagl_limit) THEN
C**            IF( fagl_id .EQ. 0 .OR. fagl_id .EQ. id )THEN
C**              IF( INT( step/60. ) .GE. fagl_start )THEN
C**                CALL faglStartWatch( h2(1:nk), istat )
C**                CALL faglShow( h2(1:nk), istat )
C**                fagl_flag = .TRUE.
C**              END IF
C**            END IF
C**          END IF
C**        END IF
C                                start variable time step loop (k, dt, time)...
C------------------------------------------------------------------------------
C!!  change limiting error band:  (could be even larger)
20      if (time .gt. (t(k) - 1.e-2)) then
C                                                   last time = last breakpoint
C                                              next time = next fixed time step
          rf = r(k)
          dtp = step - time
          k = k + 1

        end if
C                                            time elapsed since previous "step"
        dstep = time - (step - delt)

        if (time .lt. t(k) .and. time + dtp .gt. t(k)) then
C                                                                   next time =
C                                                               next breakpoint
          dtp = t(k) - time
          time = t(k)

        else

          time = step

        end if
C                                         skip this time step if it's too small
	  if (dtp .lt. DTMIN) then
            if (time .ne. step) then
              go to 20
	    else
	      cycle
            end if
	  end if

        if (cour) then
C                                                     compute next dt using max
C                                                         celerity from last dt
          nc = nint (vcm * dtp / dx)
          if (nc .eq. 0) nc = 1
          dt = dtp / float (nc)

          if (dt .lt. DTMIN) then
            nc = int (dtp / DTMIN)
            if (nc .eq. 0) nc = 1
            dt = dtp / float (nc)
          end if

        else

          dt = dtp

        end if

        dte = 0.
        sumrf = sumrf  + rf*dt
C
        do j = 2, nk
          jm = j-1
C          afac(j) = 1.
          afac(j) = 1.
          qui = q1(jm)
          if(j .eq. 2) qui = 0.5*(qub1+qub2)/width
C          flm = qui/dx/width
          hav = 0.5*(h1(j) + h1(jm))
C!!  21/12/02
          h1a(j) = hav
C** following removed 21/12/02 
C          If(h1(j) .le. 0.) then
C            h1a(j) = flm*dt
C          else
C            h1a(j) = (h1(jm)*2. + h1(j))/3.
C          End If
C!!          if(h1a(j) .lt. flm*dt) h1a(j) = flm*dt          ! new 12/02
C**  end of modifications
          if(hav .gt. 0.) then
            afac(j) = 1.0
            if( amr .gt. 0.) then
C!!new formulation                     reduce infiltration during recession
              afac(j) = 0.5*(topw(j) + topw(jm))
C              if (hav .lt. 0.5*amr) afac(j) = sqrt(2.* hav / amr)
            end if
C
          end if
C               temporarily set ql equal to inflow rate only for passing to INFILT
C!!          ql(j,1) = rf  !  changed 21/12/02
          ql(j,1) = qui/dx
        end do
C                                                 get average infiltration rate
C!!                               afac returns to represent coeff. for qlat
C!!    replace dt with dtp/ 2/18/03
        call infilt ( afac, rf, h1a, ql, dtp, fav, 1)
C
        do j = 2, nk
C                                               compute net lateral inflow rate
            diff = rf - fav(j)
            reldif = abs(diff)
            if(reldif .lt. 1.e-10) then      ! redefine ql() array
              ql(j,1) = 0.
            else
              ql(j,1) = diff
            end if
C
        end do
C
C        if(diag .and. i .lt. 50) then
C           write(99,919) 'smr ', rf*po, sumrf
C           write(99,919)' h1',(h1(j),j=1,nk)
C           write(99,919)' ql',(ql(j,1),j=1,nk)
 919  format(2x,a3,(7g12.4))
C        end if
C                                                 start Courant loop (dt, m)...
C------------------------------------------------------------------------------

50      do m = 1, nc

          vcm = 0.
C                                           time elapsed since beginning of dtp
          dte = dte + dt
C                                                time elapsed since last "step"

          deltr = dstep + dte
          qub = qub1 + dqub * deltr / delt
C                                                            interpolate inflow
          if (qub .gt. 0.) then

            q2(1) = qub / width        ! units L^2/T
            unif = .false.
C                                               compute upstream boundary depth
C------------------------------------------------------------------------------
            if (amr .lt. 1.E-8) then
C                                                            no microtopography
              if (laminar) then

                if (q2(1) .le. qtrans) then
	            p(1)     = plam
	            alpha(1) = alam
 	          else
	            p(1)     = pturb
	            alpha(1) = aturb
	          end if
	        end if

              h2(1) = (q2(1) / alpha(1))**(1. / p(1))
C                                                                      celerity
              vc = p(1) * alpha(1) * h2(1)**(p(1) - 1.)

            else
C                                                               microtopography
              qb = qub * fnw

              if (qb .ge. uqc) then

                a2u = (qb / alpha(1) * wpc**(p(1) - 1.))**(1. / p(1))

              else

                if(qub .gt. 5. * qub1) then

                  a2u = ((qb / alpha(1)) * (2. * c2 * c2 / c1)**
     &                  (.5 * (p(1) - 1.)))**(1. / (.5 * (p(1) + 1.)))

                end if
                source = 'PLanUB'
                call iter (erpaub, a2u, 0., 1000., ierr, source)

                if (ierr .gt. 0) call errxit
     &                           (msg(ms:14), 'no convergence (erpaub)')

              end if

              h2(1) = a2u / ams
              zf = zfu (bw, a2u, c1)
              wpr2 = wprf (bw, zf, c1, c2)
              trace(3:6) = 'UPBd'
              vc = dqdaf (alpha(1), a2u, zf, wpr2, p(1), c2, trace)

C                                                              compute celerity
            end if

          else

            q2(1) = 0.
            h2(1) = 0.
            vc = 0.

          end if

          if (vc .gt. vcm) vcm = vc
C                                                               zone A solution
C------------------------------------------------------------------------------
          if (unif .and. ql(2,1) .gt. 0.) then

!*          8/2005 added for laminar flow option

            if (laminar) then

              if (h2(2) .le. htrans) then
	          p(2)     = plam
	          alpha(2) = alam
	        else
	          p(2)     = pturb
	          alpha(2) = aturb
	        end if
	      end if

            h2(2) = h1(2) + ql(2,1) * dt
C                                                              compute celerity
            if (h2(2) .ge. h2c) then

              a2u = h2(2) * ams
              vc = p(2) * alpha(2) * (a2u / wpc)**(p(2) - 1.)
              wpr2 = wpc
              q2(2) = alpha(2) * h2(2)**p(2) / (wpr2 / ams)**(p(2) - 1.)

             else
C
              if(h2(2) .gt. 0.) then
                a2u = h2(2) * ams
                zf = zfu (bw, a2u, c1)
                wpr2 = wprf (bw, zf, c1, c2)
                q2(2) = qfunc (alpha(2), a2u, wpr2, p(2)) / ams
                trace(3:6) = 'zonA'
                vc = dqdaf (alpha(2), a2u, zf, wpr2, p(2), c2, trace)
              else
                h2(2) = 0.
              end if
C
            end if
            if(h2(2) .le. 0.) exit

            xc = xc + vc * dt
            if (vc .gt. vcm) vcm = vc

            if (xc .lt. dx) then
C                                                                  solution is
C                                                                     in zone A
!            if(diag ) write(99,222) vc,xc,a2u,bw,zf,wpr2
! 222  Format(" geom",8g13.4)

              do j = 3, nk
                h2(j)    = h2(2)
                q2(j)    = q2(2)
	          p(j)     = p(2)
	          alpha(j) = alpha(2)
              end do

              if (laminar) then
                do j = 3, nk
	            p(j)     = p(2)
	            alpha(j) = alpha(2)
                end do
	        end if
C                                                 bypass iterative calculations
              go to 70

            else
C                                                     solution is out of zone A
              unif = .false.

              go to 55

            end if

          end if
C only gets to this point if ql() negative
          if (amr .gt. 0.) then
C                  under microtopography, reduce infiltration during recession 
C------------------------------------------------------------------------------
C            if (h2(1) .eq. 0. .and. h1(2) .gt. 0. .and.
C     &                                    ql(2) .lt. 0.) ql(1) = 0.

            do j = 2, nk
                ql(j,1) = ql(j,1) * afac(j)
            end do

          end if
C                                                         start spatial loop...
C------------------------------------------------------------------------------

55        do j = 1, nk - 1

            jp1 = j + 1
            htest = h2(j) + h1(jp1) + h1(j)
            qlav = ql(jp1,1) !0.5 * (ql(j) + ql(jp1))
C                                                               test for runoff
            if (htest .le. 1.e-8 .and. qlav .le. 1.e-11) then

              h2(jp1) = 0.
              q2(jp1) = 0.
C!!
              ql(j,1) = 0.
              vc = 0.

            else
C                                                iterative solution for h2(j+1)
              h2jp1 = (h2(j) + h1(jp1)) / 2.

!*            8/2005 added for laminar flow option

              if (laminar) then

                if ((2.* h1(jp1) + h2(j)) / 3. .le. htrans) then
	            p(jp1)     = plam
	            alpha(jp1) = alam
	          else
	            p(jp1)     = pturb
	            alpha(jp1) = aturb
	          end if
	        end if

              source = 'PlanH2'
              write(source(5:6),'(i2)') jp1
              hmin = min(-h2jp1,-tol)
C              hmin = -10.*tol
              call iter (errh2, h2jp1, hmin, 10., ierr, source)

              if (ierr .gt. 0) then
                if(diag) write(99,*) jp1, hmin,h2jp1,ql(j,1),ql(jp1,1)
                close(99)
                call errxit(msg(ms:14), 'no convergence (errh2)')
              else if (j .eq. 2) then
                itot = itot - ierr
              end if
              if(j .eq. 1) write(99,'(i3," itrs")') ierr

              if (h2jp1 .lt. 1.E-7) h2jp1 = 0.        !!!!!!!!!!
              h2(jp1) = h2jp1
C                                             au2 = area of flow per unit width
              a2u = h2jp1 * ams

              if (h2jp1 .ge. h2c) then

                wpr = wpc
                vc = p(jp1) * alpha(jp1) *
     &               (h2jp1 * ams / wpc)**(p(jp1) - 1.)

              else

                zf = zfu (bw, a2u, c1)
                wpr = wprf (bw, zf, c1, c2)
C                                                              compute celerity
                trace(3:6) = 'H2IT'
                vc = dqdaf (alpha(jp1), a2u, zf, wpr, p(jp1), c2, trace)

              end if
            end if
C                                           vcm is max. celerity for current dt
            if (vc .gt. vcm) vcm = vc

          end do
C                                                       ... end of spatial loop
C------------------------------------------------------------------------------

C                                                       check Courant condition
C------------------------------------------------------------------------------

C**      IF( fagl_flag )THEN
C**        CALL DATE_AND_TIME( VALUES=fagl_vals )
C**        fagl_s0 = fagl_vals(8)
C**        DO
C**          CALL DATE_AND_TIME( VALUES=fagl_vals )
C**          fagl_ds = fagl_vals(8) - fagl_s0
C**          IF( fagl_ds < 0 ) fagl_ds = 1000 + fagl_ds
C**          IF( fagl_ds >= fagl_delay )EXIT
C**        END DO

C**        WRITE(fagl_str,"('h2: t=',f8.2,' dt=',f6.2)")
C**     +  (step+dte)/60.,dt/60.
C**        CALL faglUpdate ( h2(1:nk), istat )
C**        CALL faglName( h2(1:nk), fagl_str, istat )
C**      END IF

          if (vcm * dt / dx .gt. 1. .and. dt .gt. DTMIN) then
C                                                                     reduce dt
            dtprev = dt
            dtr = dtp - dte + dt
            nc = int (2. * vcm * dtr / dx)
            if (nc .eq. 0) nc = 1
            dt = dtr / float (nc)

            if (dt .lt. DTMIN) then
              nc = int (dtr / DTMIN)
              if (nc .eq. 0) nc = 1
              dt = dtr / float (nc)
            end if

            if (cour) then
C                                                recompute with new subinterval
              dte = dte - dtprev
              go to 50

            else
C                                             find minimum dt over each quarter
              if (i .le. p1) then
                if (dt .lt. dtm(1)) dtm(1) = dt
              else if (i .le. p2) then
                if (dt .lt. dtm(2)) dtm(2) = dt
              else if (i .le. p3) then
                if (dt .lt. dtm(3)) dtm(3) = dt
              else
                if (dt .lt. dtm(4)) dtm(4) = dt
              end if
C                                                               don't recompute
              nc = 1
              dt = dtp

            end if
          end if

70        continue
C
C!! added for balance component RES  03/12/02:
          smup = smup + 0.5*(q1(1)+q2(1))*dt
C
          hsum = 0.5*(h2(1) + h2(nk))
C
          do j = 1, nk
C!!           moved up from sed option below              calculate top widths
            if (h2(j) .lt. h2c) then

              a2u = h2(j) * ams
              z2 = zfu (bw, a2u, c1)
              topw(j) = z2 / ams
            else
              topw(j) = 1.

            end if
C                                                            infiltrated volume
            if(j .ge. 2) then
              if(j .lt. nk) hsum = hsum + h2(j)

	        if (h2(j) .eq. 0.) then

                dvf = 0.5 * ((h2(j-1) - h1(j) - h1(j-1)) / dt -
     &                       (q2(j-1) + q1(j-1) - q1(j)) / dx) + rf

	        else
                dvf = rf -  ql(j,1)    !)0.5 * (ql(j-1) +
	        end if

!              if (h1(j-1) .le. 1.e-6) dvf = amin1 (rf, dvf)

              qbal(6) = qbal(6) + dt * dvf
              smql = smql +  ql(j,1) * dt * dx

              dvf_delt = dvf_delt + dt * dvf
              rf_delt = rf_delt + dt * rf
            end if

            h1(j) = h2(j)
            q1(j) = q2(j)
          end do
C
          storh = hsum*dx
          smin = qbal(6)*dx
          qbal(10) = qbal(10) + dt * 0.5 * (q1(nk) + q2(nk))  !outflow summation
C!! added smup to balance calcs  RES  03/12/02
          sterr =  smql + smup - qbal(10) - storh   !sumrf*xleng - smin

          if(diag) then
C            write(99,'(" iter2: ",i4, g12.4)') itot, h1r
C            write(99,'("  ql",8g12.4)') (ql(j,1),j=1,nk)
            write(99,'("  h2",8g12.4)') (h2(j),j=1,nk)
            write(99,'("  q2",8g12.4)') (q2(j),j=1,nk)
           if(amr .gt. 0.)write(99,'(" tw ",8g12.4)') (topw(j),j=1,nk)
  101 format("    step  t(min)  sums(L^2): RF       Infil        Excess 
     &       OutQ         Stor         Bal")
            write(99,101)
            write(99,'(" ba:",i4,8g13.5)') i-1,time/60.,sumrf*xleng,
     &  smin, smql, qbal(10), storh, sterr
C
          end if
C                                                     ending time of current dt
          timec = time - dtp + dte

          if (q2(nk) .gt. qp) then
C                                                                peak discharge
            qp = q2(nk)
            tp = timec

          end if

          qpw = qp * width
C          
C                                                                      sediment
C------------------------------------------------------------------------------
          if (sed) then

C                                                                        kinsed
            call kinsed ( i, dt, deltr, ql, h2, q2,
     &                   topw, rf, qup, qt, timec)

          end if
C                                                                outflow volume
        end do
C                                                        ...end of Courant loop
C------------------------------------------------------------------------------

        if (time .ne. step) go to 20
C                                             ...end of variable time step loop
C------------------------------------------------------------------------------

C                                                               store discharge
        qout2 = width * q2(nk)

        call store (outq(1), i, qout2)
C!!  revised RES 3/12/02                                                                 inflow volume
        if(i .lt. limit) vin = vin + qub2

        qub1 = qub2
C  
        if(nele .le. 5) call progss(msg(ms:18))
C

C**        IF( fagl_flag )THEN
C**          IF( INT( step/60. ) .GE. fagl_limit )THEN
C**            CALL faglEndWatch( h2(1:nk), istat )
C**            CALL faglClose( h2(1:nk), istat )
C**            fagl_flag = .FALSE.
C**          END IF
C**        END IF

        dvf_delt = dvf_delt*1000.*3600./delt/FLOAT(nk-1)
        WRITE( 88, '(F10.2,",",F10.2,",",F10.2)' )
     &  step/60., dvf_delt, rf_delt*1000.*3600./delt/FLOAT(nk-1)

!**
!        WRITE(20,'(F8.0,",",F12.4,",",L1)') step/60.,q2(nk)*width,unif 
!        WRITE(20,'(15(F12.8))') (h2(j), j = 1, nk)
!**
      end do
C                                               ...end of fixed time step loop
C------------------------------------------------------------------------------

C                                                  finish sediment computations
C      if(i .ge. limit) stop ' before sedfin '
      if (sed) call sedfin

C                                                                 inflow volume
      vin = (delt / 2.) * (2. * vin + qub2)
C                                                                outflow volume
      qbal(10) = qbal(10) * width

      do j = 2, nk - 1
C                                                                storage volume
        qbal(9) = qbal(9) + h2(j)

      end do

      qbal(9) = (dx / 2.) * (h2(1) + 2. * qbal(9) + h2(nk)) * width

C                                                            infiltrated volume
      qbal(6) = max (qbal(6) * width * dx, 0.)

C                                                         pass values to writer

      call qwrt (id, 1, msg(ms:18), qpw, tp, vin, coarea,
     &            outq, 0, outr, qrmax)

      return

      end

C------------------------------------------------------------------------------


      subroutine errh2 (h2jp1, ferr, dferr)


C   Computes the residual ferr and derivative dferr at h2jp1 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for flow on a sloping plane.  Includes microtopog. geometry

      integer j, jp1
      character(LEN=6) :: trace
      common /plane1/ dt, dx, j, jp1, qlav, alpha(20), p(20), h1(20), 
     &                h2(20), q1(20), q2(20), h2c, ams, wpc, c1, c2, bw, 
     &                qb, omega
      data trace/'errh2 '/

      dhdt = .5 * (h2jp1 - h1(jp1) + h2(j) - h1(j)) / dt


      if (h2jp1 .ge. h2c .or. h2c .le. 0.) then

        a2u = h2jp1 * ams
        a2pos = max(0., a2u)
        h2p = max(0.,h2jp1)
        q2(jp1) = alpha(jp1) * h2p **p(jp1) / (wpc / ams)**(p(jp1) - 1.)  !m^2/s

        dferr = .5 / dt + omega * alpha(jp1) * p(jp1) *
     &          (a2pos / wpc)**(p(jp1) - 1.) / dx

      else

        a2jp1 = h2jp1 * ams
        zfjp1 = zfu (bw, a2jp1, c1)   ! shape parameter
        wp2jp1 = wprf(bw, zfjp1, c1, c2)  ! wetted perim
        q2(jp1) = qfunc (alpha(jp1), a2jp1, wp2jp1, p(jp1)) / ams

        dferr = .5 / dt + omega * dqdaf (alpha(jp1), a2jp1, zfjp1,
     &                                   wp2jp1, p(jp1), c2, trace)
     &                                   / dx / ams

      end if
      
      dqdx = (omega * (q2(jp1) - q2(j)) +
     &       (1. - omega) * (q1(jp1) - q1(j))) / dx

      ferr = dhdt + dqdx - qlav
C      if(jp1 .eq. 2) then
C        write(99,88)  h2jp1, q2(jp1), dhdt, dqdx, h2(j), h1(j), h1(jp1)
   88 format(2x,7g13.4)
C      end if

      return

      end

C------------------------------------------------------------------------------


      subroutine erpaub (aub, ferr, dferr)


C   Computes the residual and derivative of the normal discharge function
C   to obtain the upstream area for a given upstream flow for trapezoidal
C   microtopography

      integer j, jp1
      character*6 trace
      common /plane1/ dt, dx, j, jp1, qlav, alpha(20), p(20), h1(20), 
     &                h2(20), q1(20), q2(20), h2c, ams, wpc, c1, c2, bw, 
     &                qb, omega
      data trace/'erpaub'/

      zpd = zfu (bw, aub, c1)
      wpr = wprf (bw, zpd, c1, c2)

      ferr = alpha(1) * aub**p(1) / wpr**(p(1) - 1.) - qb

      dferr = dqdaf (alpha(1), aub, zpd, wpr, p(1), c2, trace)

      return

      end

C------------------------------------------------------------------------------


      real function qfunc (al, a, wp, p)

C   Computes discharge as a function of area 'a' and wetted perimeter 'wp'

      real al, a, wp, p

      if(a .gt. 0.) then
        qfunc = al * a * (a / wp)**(p - 1.)
      else
        qfunc = 0.
      end if

      return

      end

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      real function wprf (w, z, c1, c2)

C   Computes the wetted perimeter of a trazpezoid of bottom width 'w'

      real w, z, c1, c2

      wprf = w + c2 * (z - w) / c1

      return

      end



      real function zfu (w, a, c1)

C   Computes geometric shape parameter 'z'

      real w, a, c1

      if( a .ge. 0.) then
        zfu = sqrt (w * w + a * 2. * c1)
      else
        zfu = w
      end if

      return

      end



      real function dqdaf (al, a2, d, wpr, p, cz, sourc)

C   Computes derivative of 'q' with resp to 'a2'

      real al, a2, d, wpr, p, cz
      character*6  sourc
C
      if(a2 .le. 0.) then
        dqdaf = 0.
      else
        dqdaf = al *(a2 / wpr)**( p - 1.) * (p -
     &        ((p - 1.) * a2 * cz / (d * wpr)))
      end if
      return

      end
