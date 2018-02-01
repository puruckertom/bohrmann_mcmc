C   Code type: FORTRAN subroutine

C   Compiler: Microsoft FORTRAN v5.1

C   Programmed by: C. Unkrich

C   Date: 2/96; 2/02

C   Description:

C     Simulates open channel flow using a four-point implicit finite-difference
C     approximation to the kinematic wave equation written for one-dimensional
C     unsteady flow in a compound trapezoidal channel, where the cross section
C     may include an "overbank" section which floods when the main channel is
C     overtopped. Hydraulic and sediment inputs include concentrated upstream
C     inflow, uniformly distributed lateral inflow, and rainfall. Infiltration
C     is coupled to the routing equations to provide a realistic interaction
C     between flow conditions and infiltration. A baseflow component can be
C     specified as a discharge value at the downstream end of the channel,
C     whereupon baseflow is linearized between the upstream and downstream
C     values, and infiltration parameters are ignored.

C     An algorithm has been added which monitors the wave celerity and, if
C     enabled,  adjusts the computational time increment to maintain numerical
C     accuracy.

C     For output, volumetric rainfall rate from upstream areas for each time
C     step are added and stored.  Dividing by contributing area gives the
C     integrated rainfall rate over the area contributing to the channel.
C------------------------------------------------------------------------------

C   Arguments:

C     qbal(10)      real        array with element volume balance components;
C                               qbal(1)  = area              (sq.m/sq.ft),
C                               qbal(2)  = rainfall          (cu.m/cu.ft),
C                               qbal(3)  = baseflow          (cu.m/cu.ft),
C                               qbal(6)  = overbank infil.       "       ,
C                               qbal(7)  = chan. infiltration(cu.m/cu.ft),
C                               qbal(9)  = storage           (cu.m/cu.ft),
C                               qbal(10) = outflow           (cu.m/cu.ft),
C                               qbal(11) = overbank area     (sq.m/sq.ft),
c                               qbal(12) = overbank rainfall (cu.m/cu.ft),
C                               qbal(13) = chan. bottom infil      "
C
C     diag          log         .true. = write diagnostic info to unit 99,

C     filen         int         unit number of input file (for reader),

C     dtm(4)        real        smallest dt to satisfy the Courant condition
C                               during each quarter of the simulation,
C------------------------------------------------------------------------------
C
C  Global parameters:
C     CHRain        logical     simulate rain on a channel width at depth DE

C     delt          real        user-specified time step (sec),

C     limit         int         number of time steps,

C     clen          real        characteristic length (m or ft),

C     cour          log         .true.  adjust time step for Courant condition,
C                               .false. do not adjust (but report),

C     units         int         1 = metric, 2 = english,

C     sed           log         .true. = route sediment,
C ---------------------------------------------------------------------------

C  Input Block (input file labels in parenthesis):

C     id (ID)        int*4         identifier (up to 6 digits),

C     xleng (LE)     real          length (m or ft),

C     idup (UP)      int*4         identifiers of up to 10 upstream
C                                  contributing elements,

C     idlat (LA)     int*4         identifier(s) of 1 or 2 lateral elements
C                                  (for a compound channel, the 1st element
C                                  listed is assumed to contribute to the main
C                                  channel, and the second element if present
C                                  to the overbank section),

C     wool (WO)      log*1         switch to apply "Woolhiser" effective
C                                  wetted perimeter function to infiltration
C                                  loss: "YES" = .true., "NO" = .false.,
C                                  default is .false.,

C     wco (WC)       real          coefficient for Woolhiser function to over-
C                                  ride default value of 0.15,

C     qbf (QB)       real          baseflow discharge at end of channel (cu. m
C                                  or cu. ft),

C     Each of the following hydraulic variables may have two values, one for
C     the upstream x-section and one for the downstream x-section, specified as
C     a list or column with two values, with the upstream value being first in
C     the list or on top of the column.


C     width (WI)        real       bottom width (m or ft) where:

C     slope (SL)        real       slope along flow path, 

C     rs (MAN or CH)    real       Manning or Chezy discharge coefficient,

C     z1 (SS1)          real       bank side slope, ratio of vertical to horizontal,

C     z2 (SS2)          real       bank side slope (for a compound section,
C                                  refers to the overbank side,

C     depth (DE)        real       representative channel depth, m or ft, at which 
C                                  width RWIDTH for computation of rain volume is used.  
C                                  default zero, so that rwidth = width
C------------------------------------------------------------------------------

C  Input Block -- OVERBANK:

C     thresh (THR)      real       threshold depth to overbank, m or ft,

C     z3 (SS3)          real       overbank side slope, V/H,

C     zb) (LS)          real       overbank lateral bottom slope,

C     width (WI)        real       bottom width (m or ft) where:

C     slope (SL)        real       slope in flow direction, 

C     rs (MA or CH)     real       Manning or Chezy discharge coefficient,
C------------------------------------------------------------------------------

C  Subroutines/functions:

C     infilt        (infilt.for)
C     infil0            "
C     iter          (miscel.for)
C     sed0          (kinsed.for)
C     kinsed             "
C     sedfin             "
C     progss        (K2RunW.for)
C     reader        (reader.for)
C     getr4              "
C     getstr             "
C     geti4              "
C     new           (miscel.for)
C     old                "
C     get                "
C     store              "
C     geta               "
C     stora              "
C     qwrt          (writer.for)
C     errxit        (errxit.for)
C------------------------------------------------------------------------------

      subroutine channl ( dtm)

      use runpars        ! global parameters
      use elpars
      use multip
      use itpars
C                                                                     arguments
C------------------------------------------------------------------------------
C
      integer  filen

      dimension  dtm(4) 

C                                                               local variables
C------------------------------------------------------------------------------

      logical cmpd, wool, unif, bflow, man, oblat, minflg ! rain,

      integer :: ierr, m, nk, j, jp1, nd, k, nu, nl, nlu, ir, ms,
     &        leng, upq,  i, inflow, idup, upr, idlat, outr  
            !,id, ltyp, idoblat, nchan

      character msg*18, c*1, string*10, source*6
      character(LEN=1), dimension(2) :: ul = (/'M','f'/)
C** new chl
      character(LEN=1) :: chl
      real :: zeron=0. , vone = 1.0                                          ! 10/02

      common /chan/ a1(20,2), a2(20,2), q1(20,2), q2(20,2), fav(20,2),
     &              wpr1(20,2), wpr2(20,2), co1(20,2), co2(20,2),
     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20,2),
     &              ewpr2(20,2), a0(20), wpr0(20), tw0(20), hb(20),
     &              ahb(20), wprhb(20), twhb(20), co1hb(20), co2hb(20),
     &              qlat(2), j, jp1, omega, qub, dt, dx, cmpd, dqbf,
     &              width(20,2), alpha(20,2), thresh(20), p, wool, wco,
     &              rf(20,2), qltr(20), afov(20), rwidth(20,2)

      integer,dimension(2) :: outq, latq, latr 
      real,dimension(1000) :: t, z, r, q, v, h
      real,dimension(2) :: dqu, qu, qu1, qu2, z1, z2, z3, qlat1, qlat2,
     &                     dthr, dqlat, depth, zb, quba
C
      dimension  qt(20), topw(20,2), slope(20,2), rs(2,2), inflow(12), 
     &           upq(10,2), qup(10,2), upr(10), ql(20,2), !, rwidth(20,2)
     &           idup(10), idlat(2), afac(20), hr1(20,2)
     
      real, parameter :: dtmin = 1.0
C                                                                   common data
C------------------------------------------------------------------------------

      external errhub, errch2 , errha

C      real :: omega = 0.8

C      data msg /'       (CHANNEL)'/


C                                                    get basic input parameters
C------------------------------------------------------------------------------
      filen = files(1)     ! unit number for input parameters
      omega = 0.8
      tol = 1.e-7          ! new 10.01 : using this as lower (neg) iter limit allows 
C                            iteration to find a zero result if appropriate
C                                                                    IDENTIFIER
      call geti4 ('ID', 0, id, ierr)

      if (ierr .eq. 3) then

        call errxit ('CHANNEL', 'invalid identifier (ID)')

      else if (ierr .ge. 1) then

        call errxit ('CHANNEL', 'missing identifier (ID)')

      end if

      if(DIAG) then
        jd = 2
        write(99,'(" Diagnostics for channel ",i3)') ID
      else
        jd = 100
      end if
C
      msg(1:16) = typname(1)
      write (msg(15:18), '(i4)') id

      do ms = 1, 5

        if (msg(ms:ms) .ne. ' ') go to 1

      end do

1     continue
      if(nele .gt. 5) call progss (msg(ms:18))

      nu = 0
C                                                upstream element identifier(s)
      do i = 1, 12
        inflow(i) = 0
      end do
C
      do k = 1, 10

        call geti4 ('UP', k, idup(k), ierr)

        if (ierr .eq. 3) call errxit
     &  (msg(ms:18), 'invalid upstream identifier (UP)')

        if (ierr .eq. 0) then
C                                                        upstream inflow (main)
          nu = nu + 1
          inflow(nu) = idup(nu)  ! new location

        else

          go to 11

        end if
      end do

11    nl = 0  ! this is the number of lateral inputs
      nm = 0  ! this is  "     "        "        "   that go into the main channel
      nlu = 0

      do k = 1, 2
C                                                    up to 2 lateral inflow(s)
        call geti4 ('LA', k, idlat(k), ierr)

        if (ierr .eq. 3) call errxit
     &  (msg(ms:18), 'invalid lateral element identifier (LA)')

        if (ierr .eq. 0) then

          nl = nl + 1
          nm = nm + 1
          inflow(nl+10) = idlat(nl)
          nlu = 1

        else
 
          go to 12

        end if
      end do
C
12    continue

      call getr4 ('LE', 0, xleng, ierr)
C                                                                       LENGTH
      if (ierr .gt. 0)
     &call errxit (msg(ms:18),'length (LE) not found')
C                                                       compute number of nodes
      nk = nint (amax1 ((15. * xleng / clen), 5.))
      if (nk .gt. 15) nk = 15
      nkm = nk - 1
      rnk = float (nkm) 
      itlim = limit
C                                                     compute spatial increment
      dx = xleng / rnk
C** 22/4/03  dx check      
      dxlim = dxlmt
      chl = 'm.'
      if(units .ge. 2) then
       dxlim = 3.33*dxlmt
       chl = 'ft'
      end if
      if(dx .ge. dxlim) then
        write(files(3), 1990) id, dx, chl!  , dxlim
 1990 format(/"  Channel ",i3,":  based on length and parameter CLEN,"/
     &"  the numerical increment of ",f6.1,1x,a2," is too large for",/
     &"  realistic numerical solution of the flow equation, and may "/
     &"  give  misleading results."/)
      end if
 !     
      if(diag) write(99,'(" length div. into ",i2," segments of ",f8.2,
     &1x,a2," each.")') nkm, dx, dlab
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
C                                                   apply Woolhiser's function?
      call getstr ('WO', 0, c, m, ierr)

      if (ierr .eq. 0 .and. (c .eq. 'y' .or. c .eq. 'Y')) then

        wool = .true.

        call getr4 ('WC', 0, val, ierr)
C                                                                   coefficient
        if (ierr .eq. 0) then

          wco = val

        else
C                                                               default is 0.15
          wco = 0.15

        end if
C                                                         conversion for metric
        if (units .eq. 1) wco = wco * 0.5520869

      else

        wool = .false.

      end if
C                                          baseflow discharge at end of channel
      call getr4 ('QB', 0, qbf, ierr)
C*2/5/2004 Distinguish between baseflow = 0 and no baseflow -- set 'bflow' here
      if (ierr .gt. 0) then
        bflow = .false.
        qbf = 0.
	else
        bflow = .true.
      end if

      cmpd = .false.
      nchan = 1
      ltyp = 1

      call getstr ('TY', 0, c, leng, ierr)
C                                                                         TYPE
      if (ierr .eq. 0 .and. (c .eq. 'C' .or. c .eq. 'c')) then

        cmpd = .true.
        nchan = 2
        ltyp = 4
        msg(1:14) = typname(ltyp)
        write (msg(15:18), '(i4)') id

      end if
      
      depth(1) = 0.   ! default
      depth(2) = 0.
C
      if(chrain) then
      
        call getr4( 'DE', 1, dep, ierr)                               !   DEPTH
C   this is the depth at which effective width for rain calculations is found    
        if(ierr .eq. 0) then
          depth(1) = max(0.,dep)
          call getr4( 'DE', 2, dep, ierr)
          if(ierr .eq. 0) then
            depth(2) = max(0.,dep)
          else
            depth(2) = depth(1)
          end if
        end if
      
      end if

C                       loop to read upstream, downstream section parameters...
C------------------------------------------------------------------------------

      do j = 1, 2

        call getr4 ('WI', j, width(j,1), ierr)

        if (ierr .gt. 0) then
C                                                                       WIDTH
          if (j .eq. 1) then

            call errxit (msg(ms:18), 'width (WI) not found')

          else
C                                                     set downstream = upstream
            width(2,1) = width(1,1)

          end if
        end if
C                                                                       SLOPE
        call getr4 ('SL', j, slope(j,1), ierr)

        if (ierr .gt. 0) then

          if (j .eq. 1) then

            call errxit
     &      (msg(ms:18), 'main section slope (SL) not found')

          else

            slope(2,1) = slope(1,1)

          end if
        end if
C                                                                       MANNING?
        call getr4 ('MA', j, rs(j,1), ierr)

        if (ierr .eq. 0) then

          p = 1.66667

          if (units .eq. 1) then
            fac = 1.
          else
            fac = 1.49
          end if

          man = .true.

        else
C                                                                        CHEZY?
          call getr4 ('CH', j, rs(j,1), ierr)

          if (ierr .eq. 0) then

            p   = 1.5
            fac = 1.
            man = .false.

          else if (j .eq. 1) then

            call errxit
     &      (msg(ms:18), 'hydraulic roughness (MA or CH) not found')

          else

            rs(2,1) = rs(1,1)

          end if
        end if
C                                                    apply roughness multiplier
        rs(j,1) = rs(j,1) * rmn
C                                                                side slope SS1
        call getr4 ('SS1', j, z1(j), ierr)

        if (ierr .gt. 0) then

          if (j .eq. 1) then

            call errxit
     &      (msg(ms:18), 'main section side slope (SS1) not found')

          else

            z1(2) = z1(1)

          end if
        end if

        call getr4 ('SS2', j, z2(j), ierr)
C                                                                side slope SS2
        if (ierr .gt. 0) then

          if (j .eq. 1) then

            call errxit
     &      (msg(ms:18), 'main section side slope (SS2) not found')

          else

            z2(2) = z2(1)

          end if
        end if

      end do
C                                            ...end of upstream/downstream loop
C------------------------------------------------------------------------------

3     continue

C                                              interpolate hydraulic parameters
C------------------------------------------------------------------------------

      if (width(1,1) .eq. 0.) width(1,1) = .01
      width(nk,1) = width(2,1)
      if (width(nk,1) .eq. 0.) width(nk,1) = .01
      dwidth = (width(nk,1) - width(1,1)) / rnk
C      rwidth(nk,1) = rwidth(2,1)
      if (man) then
        alpha(1,1) = fac * sqrt (slope(1,1)) / rs(1,1)
        alpha(nk,1) = fac * sqrt (slope(2,1)) / rs(2,1)
      else
        alpha(1,1) = fac * sqrt (slope(1,1)) * rs(1,1)
        alpha(nk,1) = fac * sqrt (slope(2,1)) * rs(2,1)
      end if
      dslope = (slope(2,1) - slope(1,1)) / rnk
      slope(nk,1) = slope(2,1)
      drs = (rs(2,1) - rs(1,1)) / rnk
      co1(1,1) = 1. / z1(1) + 1. / z2(1)
      co1(nk,1) = 1. / z1(2) + 1. / z2(2)
!      dco1 = (co1(nk,1) - co1(1,1)) / rnk
      dz1 = (z1(2) - z1(1)) / rnk
      dz2 = (z2(2) - z2(1)) / rnk
      co2(1,1) = sqrt (1. + 1. / z1(1)**2) + sqrt (1. + 1. / z2(1)**2)
      co2(nk,1) = sqrt (1. + 1. / z1(2)**2) + sqrt (1. + 1. / z2(2)**2)
      drwid = 0.
      if(chrain) then
        rwidth(1,1) = width(1,1) + depth(1)*co1(1,1)
        rwidth(nk,1) = width(nk,1) + depth(2)*co1(nk,1)
        drwid = (rwidth(nk,1) - rwidth(1,1)) / rnk
      else
        rwidth(1,1) = 0.
        rwidth(nk,1) = 0.
      end if

      do j = 2, nk - 1

        rj = float (j - 1)
        width(j,1) = width(1,1) + dwidth * rj
        rwidth(j,1) = rwidth(1,1) + drwid * rj
        slope(j,1) = slope(1,1) + dslope * rj
        if (man) then
          alpha(j,1) = fac * sqrt (slope(j,1)) / (rs(1,1) + drs * rj)
        else
          alpha(j,1) = fac * sqrt (slope(j,1)) * (rs(1,1) + drs * rj)
        end if
        z1j = z1(1) + dz1 * rj
        z2j = z2(1) + dz2 * rj
C*2/3/2004 co1 is nonlinear with respect to z1 and z2
!        co1(j,1) = co1(1,1) + dco1 * rj
        co1(j,1) = 1. / z1j + 1. / z2j
        co2(j,1) = sqrt (1. + 1. / z1j**2) + sqrt (1. + 1. / z2j**2)

      end do
C                                                          initialize variables
C------------------------------------------------------------------------------

      do j = 1, nk

        h2(j) = 0.
        wpr2(j,1) = width(j,1)
        wpr1(j,1) = width(j,1)        
        wpr2(j,2) = 0.                ! default setting with no overbank
        wpr1(j,2) = 0.
        width(j,2) = 0.                 ! new 12/12/01
        ewpr1(j,2) = width(j,2)
        rwidth(j,2) = 0.
        afov(j) = 0.
!       
        ewpr1(j,1) =  width(j,1)
        thresh(j) = 100.

        if (wool) ewpr1(j,1) = 0.05 * width(j,1)  ! starting value

        do k = 1, 2
          h1(j,k) = 0.
          a1(j,k) = 0.
          a2(j,k) = 0.
          q1(j,k) = 0.
          q2(j,k) = 0.
          fav(j,k) = 0.
          topw(j,k) = 0.
          rf(j,k) = 0.
          ql(j,k) = 0.
        end do

      end do
C                                                    initialize volumes to zero
      do j = 1, 15
        qbal(j) = 0.
      end do

      do k = 1, 2
        qlat(k) = 0.
        qlat1(k) = 0.
        qlat2(k) = 0.
        dqlat(k) = 0.
        qu1(k) = 0.
        qu2(k) = 0.
      end do

      qtot1 = 0.
      qp = 0.
      tp = 0.
      vin = 0.
      xc = 0.
      unif = .true.

      p1 = limit / 4
      p2 = limit / 2
      p3 = 3 * limit / 4

      iretn = 0      ! count number of courant resets of dt
      nc = 1
      coarea = 0.
      qrain = 0.
C      tfinl = float (limit - 1) * delt   !this is in seconds
      barea = xleng * 0.5*(width(1,1) + width(nk,1)) ! this is channel bottom area
      call interp (x, y, t, z, nd, sat, ngage)

      if (chrain) then
C                                                                      rainfall
C------------------------------------------------------------------------------
        area = xleng * 0.5 * (rwidth(1,1) + rwidth(nk,1)) 
        rarea = area
        qbal(1) =  area
        coarea = coarea + area

       do ir = 2, nd 

           r(ir-1) = (z(ir) - z(ir-1)) / (t(ir) - t(ir-1))
C           write(77,'(i3,2f10.4)')ir, z(ir), t(ir)/60.
        end do
C                                                               rainfall volume
C        r(nd+1) = r(nd)               ! to insure tfin period is covered
C**        qbal(2) = area * z(nd)   !this assumes run goes to end of rain record
C
      else
 
        rarea = 0.
        qbal(1) = xleng * (width(1,1) + width(nk,1)) /2.
        nd = 2
        t(1) = 0.
        r(1) = 0.
        t(2) = float (limit) * delt
        r(2) = 0.
C**1/22/03 SAT will be set to -1 in interp() if not specified in rainfall file
C        sat = -1.              ! default no init. sat set yet
        
      end if

      call new (id, 'RA', outr)
      call new (id, 'RO', outq(1))
      if (cmpd) call new (id, 'ro', outq(2))

C                                     read / initialize infiltration variables
C------------------------------------------------------------------------------

      call infil0 (id, 1, nk, msg(ms:18), diag, sat)
C                                                               ...and sediment
      if (sed) call sed0 ( 1, id, msg(ms:18), inflow, nu, nm,
     &                    nk, omega, dx, slope, vone)

C
      if (cmpd) then
C                                                              compound channel
C------------------------------------------------------------------------------

C                                                              read input block
        call reader (filen, string, ierr)

        if (ierr .gt. 0 .or. string(1:8) .ne. 'OVERBANK')
     &  call errxit (msg(ms:18), 'overbank input block not found')
     
        oblat = .false.

        if (nl .lt. 2) then
C                                              look for possible lateral inflow
C                                              associated with overbank section
          nl = nl + 1

          call geti4 ('LA', 1, idlat(nl), ierr)

          if (ierr .eq. 3) call errxit
     &    (msg(ms:18), 'invalid lateral element identifier (LA)')

          if (ierr .eq. 0) then
            oblat = .true.       ! this lateral input assigned to overbank
            if(nl .gt. 2) call errxit
     &    (msg(ms:18), 'too many total lateral elements (LA)')
            nlu = 2  
          else
            nl = nl - 1   ! no input to overbank found - reset NL

          end if
        
        end if
C                       loop to read upstream, downstream section parameters...
C------------------------------------------------------------------------------

        do j = 1, 2

          call getr4 ('WI', j, width(j,2), ierr)
C                                                                overbank width
          if (ierr .gt. 0) then

            if (j .eq. 1) then

              call errxit (msg(ms:18), 'overbank width (WI) not found')

            else
C                                                     set downstream = upstream
              width(2,2) = width(1,2)

            end if
          end if

          if(chrain) then
            rwidth(1,2) = width(1,2)
            rwidth(2,2) = width(2,2)
          end if
C                                                                        SLOPE
          call getr4 ('SL', j, slope(j,2), ierr)

          if (ierr .gt. 0) then

            if (j .eq. 1) then
              slope(1,2) = slope(1,1)
C              call errxit (msg(ms:16), 'overbank slope (SL) not found')

            else

              slope(2,2) = slope(nk,1)

            end if
C
          end if
C                                                                       MANNING?
          call getr4 ('MA', j, rs(j,2), ierr)

          if (ierr .eq. 0) then

            if (.not. man) call errxit
     &      (msg(ms:18), 'must also use Chezy for overbank')

          else
C                                                                       CHEZY
            call getr4 ('CH', j, rs(j,2), ierr)

            if (ierr .eq. 0) then

              if (man) call errxit
     &        (msg(ms:18), 'must also use Manning for overbank')

            else if (j .eq. 1) then

              call errxit
     &        (msg(ms:18), 'overbank hydraulic roughness not found')

            else

              rs(2,2) = rs(1,2)

            end if
          end if
C                                                                side slope SS3
          call getr4 ('SS3', j, z3(j), ierr)

          if (ierr .gt. 0) then

            if (j .eq. 1) then

              call errxit
     &        (msg(ms:18), 'overbank side slope (SS3) not found')

            else

              z3(2) = z3(1)

            end if
          end if

          call getr4 ('SSL', j, zb(j), ierr)
C                                                          lateral bottom slope
          if (ierr .gt. 0) then

            if (j .eq. 1) then

              call errxit
     &        (msg(ms:18),'missing overbank lateral bottom slope (SSL)')

            else

              zb(2) = zb(1)

            end if
          end if

          if (zb(j) .eq. 0.) zb(j) = 0.0001
C                                                               threshold depth
          call getr4 ('THR', j, dthr(j), ierr)

          if (ierr .gt. 0) then

            if (j .eq. 1) then

              call errxit
     &        (msg(ms:18), 'threshold depth (THR) not found')

            else

              dthr(2) = dthr(1)

            end if
          end if
        end do
C                                              interpolate hydraulic parameters
C------------------------------------------------------------------------------

        if (width(1,2) .eq. 0.) width(1,2) = .01
        width(nk,2) = width(2,2)
        if (width(nk,2) .eq. 0.) width(nk,2) = .01
        dwidth = (width(nk,2) - width(1,2)) / rnk
        rwidth(nk,2) = rwidth(2,2)
        drwid = (rwidth(nk,2) - rwidth(1,2)) / rnk
        if (man) then
          alpha(1,2) = fac * sqrt (slope(1,2)) / rs(1,2)
          alpha(nk,2) = fac * sqrt (slope(2,2)) / rs(2,2)
        else
          alpha(1,2) = fac * sqrt (slope(1,2)) * rs(1,2)
          alpha(nk,2) = fac * sqrt (slope(2,2)) * rs(2,2)
        end if
        dslope = (slope(2,2) - slope(1,2)) / rnk
        slope(nk,2) = slope(2,2)
        drs = (rs(2,2) - rs(1,2)) / rnk
        thresh(1) = dthr(1)
        thresh(nk) = dthr(2)
        dthresh = (dthr(2) - dthr(1)) / rnk
        co1(1,2) = 1. / z3(1)
        co1(nk,2) = 1. / z3(2)
        co1a(1) = 1. / z1(1)
        co1a(nk) = 1. / z1(2)
        co2(1,2) = sqrt (1. + 1. / z3(1)**2)
        co2(nk,2) = sqrt (1. + 1. / z3(2)**2)
        dz3 = (z3(2) - z3(1)) / rnk
        co2a(1) = sqrt (1. + 1. / z1(1)**2)
        co2a(nk) = sqrt (1. + 1. / z1(2)**2)
        a0(1) = thresh(1) * width(1,1) + 0.5 *
     &          thresh(1) * thresh(1) * co1(1,1)
        a0(nk) = thresh(nk) * width(nk,1) + 0.5 *
     &           thresh(nk) * thresh(nk) * co1(nk,1)
        da0 = (a0(nk) - a0(1)) / rnk
        wpr0(1) = width(1,1) + thresh(1) * co2(1,1)
        wpr0(nk) = width(nk,1) + thresh(nk) * co2(nk,1)
        tw0(1) = width(1,1) + thresh(1) * co1(1,1)
        tw0(nk) = width(nk,1) + thresh(nk) * co1(nk,1)
        dtw0 = (tw0(nk) - tw0(1)) / rnk
        dzb = (zb(2) - zb(1)) / rnk
        hb(1) = width(1,2) * zb(1)
        hb(nk) = width(nk,2) * zb(2)
        dhb = (hb(nk) - hb(1)) / rnk
        co1hb(1) = 1. / zb(1)
        co1hb(nk) = 1. / zb(2)
        co2hb(1) = sqrt (1. + 1. / zb(1)**2)
        co2hb(nk) = sqrt (1. + 1. / zb(2)**2)
        ahb(1) = 0.5 * (hb(1) * hb(1)) / zb(1)
        ahb(nk) = 0.5 * (hb(nk) * hb(nk)) / zb(2)
        wprhb(1) = hb(1) * co2hb(1)
        wprhb(nk) = hb(nk) * co2hb(nk)
        twhb(1) = hb(1) / zb(1)
        twhb(nk) = hb(nk) / zb(2)
        dtwhb = (twhb(nk) - twhb(1)) / rnk

        do j = 2, nk - 1

          rj = float (j - 1)
          z1j = z1(1) + dz1 * rj
          z3j = z3(1) + dz3 * rj
          width(j,2) = width(1,2) + dwidth * rj
          rwidth(j,2) = rwidth(1,2) + drwid * rj
          slope(j,2) = slope(1,2) + dslope * rj
          if (man) then
            alpha(j,2) = fac * sqrt (slope(j,2)) / (rs(1,2) + drs * rj)
          else
            alpha(j,2) = fac * sqrt (slope(j,2)) * (rs(1,2) + drs * rj)
          end if
          thresh(j) = thresh(1) + dthresh * rj
          co1(j,2) = 1. / z3j
          co1a(j) = 1. / z1j
          co2(j,2) = sqrt (1. + 1. / z3j**2)
          co2a(j) = sqrt (1. + 1. / z1j**2)
          wpr0(j) = width(j,1) + thresh(j) * co2(j,1)
          a0(j) = a0(1) + da0 * rj
          tw0(j) = tw0(1) + dtw0 * rj
          zbj = zb(1) + dzb * rj
          hb(j) = hb(1) + dhb * rj
          co1hb(j) = 1. / zbj
          ahb(j) = 0.5 * (hb(j) * hb(j)) / zbj
          co2hb(j) = sqrt (1. + 1. / zbj**2)
          wprhb(j) = hb(j) * co2hb(j)
          twhb(j) = twhb(1) + dtwhb * rj

        end do
C
        if(chrain) then
          xarea = xleng* 0.5*(rwidth(1,2) + rwidth(nk,2))
          rarea = rarea + xarea
          coarea = coarea + xarea
          qbal(11) = xarea
C          qbal(1) = qbal(1) + xarea
          qbal(12) = xarea*z(nd)
C          qbal(2) = qbal(2) + xarea*z(nd)
        else
          qbal(1) = qbal(1) + xleng * 0.5*(width(1,2) + width(nk,2))
        end if

C                             read / initialize overbank infiltration variables
C------------------------------------------------------------------------------

        call infil0 (id, 2, nk, msg(ms:18),  diag, sat)
C                                                               ...and sediment
        if (sed) call sed0 ( 2, id, msg(ms:18), inflow, nu, nl,
     &                       nk, omega, dx, slope, vone)


      end if       !  end of compound section read
C                                              get indices to storage locations
C                                              for rainfall, inflow and outflow
C------------------------------------------------------------------------------

      do k = 1, nu

        upr(k) = 0

        call old (idup(k), 'RA', upr(k), ierr)

        call old (idup(k), 'RO', upq(k,1), ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:18), 'upstream inflow not found')

C moved        inflow(k) = idup(k)

        call old (idup(k), 'ro', upq(k,2), ierr)

        if (ierr .eq. 0) inflow(k) = -idup(k)

      end do

      do k = 1, nl

        latr(k) = 0

        call old (idlat(k), 'RA', latr(k), ierr)

        call old (idlat(k), 'RO', latq(k), ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:18), 'lateral inflow not found')

C moved        inflow(10+k) = idlat(k)

      end do
      
C                                                      update contributing area
C------------------------------------------------------------------------------
      do k = 1, nu

        call geta (idup(k), area)
        coarea = coarea + area

      end do

      do k = 1, nl

        call geta (idlat(k), area)
        coarea = coarea + area

      end do
C                                                       store contributing area
      call stora (id, coarea)

C                                                 compute integrated rain rates
C------------------------------------------------------------------------------
      qrmax = 0.
      k = 1
      qrf = r(1) * rarea
C!!  new 9/03 RES
      lmn = limit - 1
      zlr = z(1)
C!!              11/02
      sumrch = 0.

c11 replaced 9/03:
C      do i = 1, limit
C        time = delt * float (i - 1)
C        if (time .gt. t(k+1) .and. (k+1) .le. nd) then
C          k = k + 1
C          qrf = r(k) * rarea
C        end if
C!!  now new code:
      do i = 1, lmn
        timen = delt*float(i)   ! seconds
        do while (t(k) .lt. timen)
          k = k+1
        end do
        km = k-1
        znr = z(km) + (z(k) - z(km))*(timen - t(km))/(t(k) - t(km))
        qrf = (znr - zlr)/delt * rarea
        zlr = znr
C !! end new code  
C!!  corrected 12/02.  can't use z for cases where run stops before rain
        if(i .gt. 1) sumrch = sumrch + qrf*delt

        qrain = qrf

        do j = 1, nu
          if (upr(j) .gt. 0) then
            call get (upr(j), i, qru)
            qrain = qrain + qru
          end if
        end do

        do j = 1, nl
          if (latr(j) .gt. 0) then
            call get (latr(j), i, qr)
            qrain = qrain + qr
          end if
        end do
C
        call store (outr, i, qrain)
        if (qrain .gt. qrmax) qrmax = qrain

      end do
C!!
      qbal(2) = sumrch
C                              obtain initial upstream inflow into main channel
C------------------------------------------------------------------------------

      do k = 1, nu

        call get (upq(k,1), 1, qup(k,1))
        qu1(1) = qu1(1) + qup(k,1)
C!! this should be done here and not just in the baseflow segment below
        if(qu1(1) .gt. 0.) then
          qub = qu1(1)
          source = 'ChnUB1'
C*2/3/2004 index j needs to be set to one (see changes to errhub)
	    j = 1
          call iter (errhub, h1(1,1), 0., 100., ierr, source)
          if (ierr .gt. 0) call errxit
     &      (msg(ms:18), 'no convergence (errhub)')

          if (cmpd .and. h1(1,1) .gt. thresh(j)) call errxit
     &      (msg(ms:18),  'incoming baseflow > threshold')

	    unif = .false.
        end if

      end do

      dqbf = 0.
C**  baseflow and upstream inflow should not be treated equally here:
C*2/5/2004 'bflow' is now set when baseflow is read from the parameter file
      if (bflow) then   ! .or. qu1(1) .gt. 0.  
C                                                      baseflow in main channel
C------------------------------------------------------------------------------
C!! bflow should only be set if there is baseflow, not upstream input
C        bflow = .true.
        qbal(3) = (qbf - qu1(1)) * delt * float (limit - 1)

        if (qbal(3) .lt. 0.) then
C                                            decreasing baseflow = infiltration
          qbal(7) = - qbal(3)
          qbal(3) = 0.

        end if

C         compute initial baseflow profile --  dqbf is (assumed linear) increase/
C         decrease in incremental baseflow between specified values qu1 and qbf

        dqbf = (qbf - qu1(1)) / float (nk-1)
        vsb = 0.
C                                                         use errhub to compute
C                                                                    flow depth
        a1(1,1) = a2(1,1)
        wpr1(1,1) = wpr2(1,1)
        q1(1,1) = qub
        if (sed) topw(1,1) = width(1,1) + h1(1,1) * co1(1,1)

        do j = 2, nk

          qub = qu1(1) + dqbf * float (j-1)

          if (qub .gt. 0.) then

            h1(j,1) = h1(j-1,1)
            source = 'ChnUB1'
            call iter (errhub, h1(j,1), 0., 100., ierr, source)
            if (ierr .gt. 0) call errxit
     &      (msg(ms:18), 'no convergence (errhub)')
            if (cmpd .and. h1(j,1) .gt. thresh(j)) call errxit
     &      (msg(ms:18),  'baseflow > threshold')
            a1(j,1) = a2(j,1)
            wpr1(j,1) = wpr2(j,1)
            q1(j,1) = qub
            if (sed) topw(j,1) = width(j,1) + h1(j,1) * co1(j,1)

          end if

          vsb = vsb + 0.5 * (a1(j,1) + a1(j-1,1))

        end do
C*10/10/2003 Compute baseflow contribution to lateral inflow
        dqbf = (qbf - qu1(1)) / xleng
C                                   vsb is storage in baseflow profile at t = 0
        vsb = vsb * dx
C**        vin = 0.5 * (q1(1,1) + qbf)
C                                                 compute initial concentration
        if (sed) call sedbf (1, a1, q1, topw)

!C*2/2/2004 compute initial baseflow volume
        qbstor = 0.
        do j = 2, nk - 1
          qbstor = qbstor + a1(j,1)
        end do
        qbstor = dx * (a1(1,1) + 2. * qbstor + a1(nk,1)) / 2.
        if (qbstor .lt. 0.) qbstor = 0.

!C*2/2/2004 don't allow zone A
	  unif = .false.
      end if
C                                                    outflow at t=0 is baseflow
      call store (outq(1), 1, qbf)

      if (cmpd) call store (outq(2), 1, zeron)

      time = 0.
C!! changed
      qu(1) = qu1(1)                                    ! 12/05  res
      ir = 1
      if (chrain) d1 = z(1)
      d2 = 0.
      rfi = 0.
      dt = delt
      itot = 0
      itrt = 0
      qtest = 0.
      sov = 0.
      swr = 0.
      swq = 0.
      tfr = 0.
      tfq = 0.
      vtn = 0.           ! running balance components for diag. out
      vup = 0.
      vlat = 0.
      smrn = 0.
C
      if(diag) write(99,'(" opts: compound=",L3,", basefl=",L3,", rain="
     &,L3)') cmpd, bflow, chrain
!**
      IF( nprint .EQ. 5)
     &  WRITE( xfile, '("T(min), j, Q(j), V(j), H(j), j = 1, nk")' )
!**
C                                       start fixed time step loop (i, delt)...
C------------------------------------------------------------------------------
      do i = 2, limit

        itm = i
        time = float (i - 2) * delt
C                                            time at START of current time step
        if (chrain) then
C                                                         compute rainfall rate
C                                                          as average over delt
 23       continue
           do while((time + delt) .gt. t(ir))
C          if (time + delt .gt. t(ir)) then
            ir = ir + 1
           end do
C            go to 23
C          end if

          d2 = z(ir-1) +  r(ir-1) * (time + delt - t(ir-1))
          rfi = (d2 - d1) / delt
          if (rfi .lt. 0.) rfi = 0.
          qtest = rfi
          d1 = d2

C**9/8/2004**following section moved past infilt() calls
!          do k = 1, nchan
!            do j = 1, nk
!              widu = min(ewpr1(j,k),rwidth(j,k))
!              if(k .eq. 2) widu = rwidth(j,2)
!              rf(j,k) = rfi * widu      ! effective rain input into flowing water
!            end do
!          end do
          smrn = smrn + rfi*0.5*(rwidth(1,1) + rwidth(nk,1))*xleng*delt
C
        end if
C
        qu2(1) = 0.
        qu2(2) = 0.
C                                                            upstream inflow(s)
        do k = 1, nu

          call get (upq(k,1), i, qup(k,1))

          qu2(1) = qu2(1) + qup(k,1)
          qup(k,2) = 0.

          if (inflow(k) .lt. 0) then
C                                                     overbank section upstream
            call get (upq(k,2), i, qup(k,2))

            if (cmpd) then

              qu2(2) = qu2(2) + qup(k,2)

            else

              qu2(1) = qu2(1) + qup(k,2)

            end if

          end if

        end do

        qtest = qtest + qu1(1) + qu2(1) + qbf
        if (cmpd) qtest = qtest + qu1(2) + qu2(2)
        vup = vup + 0.5*(qu2(1)+qu1(1))*delt
C                                                                lateral inflow
        qlat2(1) = 0.
        qlat2(2) = 0.
        do k = 1, nl

          call get (latq(k), i, qlat2k)

          qtest = qtest + qlat2k
          
          if(k .gt. nm  .and. oblat) then
            qlat2(2) = qlat2k/xleng          ! lateral flow onto overbank
          else
            qlat2(2) = 0.
            if(k .ge. 2) qlat1r = qlat2(1)      ! remember this for sed calcs
            qlat2(1) = qlat2(1) + qlat2k/xleng  !lateral flows lumped
          end if
C          qlat2(k) = qlat2(k) / xleng

        end do

        qt(1) = qu(1) - q2(1,1)
C**
        do j = 1, nk
          qt(j) = 0.
          qltr(j) = 0.  ! initialize
        end do
C                                                                 inflow volume
        vin = vin + qu2(1) + (qlat2(1) + qlat2(2)) * xleng
        vlat = vlat + 0.5*(qlat2(1) + qlat1(1))*xleng*delt  ! good for simple chans only
C
        if (cmpd) vin = vin + qu2(2)
C new segment
        source = 'ChnUBn'
        hulim = -tol
        qub = qu2(1)
        h2n = h2(1)
        if(h2n .le. 0.) then
          h2n = (qub/width(1,1)/alpha(1,1))**(1./p) ! init. estimate
        end if
C*2/3/2004 index j needs to be set to one (see changes to errhub)
	  j = 1
        call iter (errhub, h2n, hulim, 100., ierr, source)
        if(h2n .lt. 0.) h2n = 0.

        if (ierr .gt. 0) then
          call errxit (msg(ms:18), 'no convergence (errhubn)')
        end if
C                                      prepare to linearly  interpolate inflows
C------------------------------------------------------------------------------

        do k = 1, nchan

          dqu(k) = (qu2(k) - qu1(k)) / delt
C          quba(k) = (qu2(k) + qu1(k))/2.

        end do

        do k = 1, nlu      ! loop for each lateral side (after consolidation)

          dqlat(k) = (qlat2(k) - qlat1(k)) / delt

        end do

        if (.not. bflow) then
C                                   compute infiltration rate for both sections
C------------------------------------------------------------------------------

          if (cmpd) then
C                                                       overbank section
C            hr1(2) = h1(1,2)
            do j = 2, nk       ! new start at 2

C!!              ql(j,2) = qlat(2) / width(j,2) + rfi
              ql(j,2) = 0.5 * (qlat1(2) + qlat2(2)) / width(j,2) /dx !+ rfi 
               !lat flow onto overbank goes in as distributed, units of L/T
C
              if (h1(j,1) .gt. thresh(j)) then
                afov(j) = 1.0
C                if(j .gt. 1) then
                h1(j,2) = h1(j,1) - thresh(j)
                if(h1(j,2) .lt. hb(j)) afov(j) = ewpr1(j,2)/width(j,2)
C                  ewpr2(j,2) = afov(j)*width(j,2)
              else
                h1(j,2) = 0.
                afov(j) = 1.0
                ewpr2(j,2) = 0.
              end if

            end do
C!!
            call infilt ( afov, rfi, h1, ql, delt, fav, 2) ! for overbank only
            
            do j = 2, nk
C       ! xs from overbank into main
              qltr(j) = (ql(j,2) - fav(j,2)) * width(j,2) * afov(j)
              if(qltr(j) .lt. 0.) then
                qltr(j) = 0.
              else
                afov(j) = 1.0
              end if
            end do
C** patch
            qltr(1) = qltr(2)

          else  ! simple channel has already lumped lateral inflows

            do j = 1, nk  ! in absence of overbank:  (units L*2/t)
C              ql(j,2) = 0.
              qltr(j) = 0.  !0.5*(qlat1(2) + qlat2(2))
            end do

          end if
C                                                                  main section
C          ql(1,1) = rfi
          do j = 2, nk
            jm = j-1
            hm = (h1(j,1) + h1(jm,1))*.5
            if (cmpd .and. hm .gt. thresh(j)) then
!              ql(j,1) = qlat(1) / width(j,1) + rfi
              qnet = 0.
C
            else
C   lateral inflow into main chan includes excess from overbank
              qnet = max(qltr(j), 0.)/ewpr1(j,1)
!        
            end if
!      
C**            ql(j,1) = rfi
            afac(j) = 0.01
            if(ewpr1(j,1) .gt. 1.e-4) afac(j) = ewpr1(j,1)/width(j,1)
            afac(j) = min(afac(j),1.0)
C
            qui = q1(jm,1)
            if(j .eq. 2) then
              qui = 0.5*(qu1(1)+qu2(1))  !upper bound inflo
            end if
            flm = qui/dx/ewpr1(j-1,1)             ! equiv. vert. flux
C!!  replaced 18/12/02       if(j .eq. 2 .and. qui .gt. 0.) then
C**              dcorr = 0.5*(h2n - h1(1,1))/delt  ! correction for input into storage
               
C              if(diag) write(99,'(" f,h:",5g12.4)') flm*delt, hm, qui,
C     &   h2n, dcorr*delt
C              hm = (flm - dcorr)*delt
C            else
C              hm = 0.5*(hm + flm*delt)
C**            end if
C             if(flm .lt. hm/delt ) hm = max(flm*delt,0.)  !.and. h1(j+1,1) .le. 0.
C
C**            hr1(j,1) = hm + (0.5*(qlat1(1)+qlat2(1))/ewpr1(j,1)
C!! new 21/12/02  put inflows in equivalent rain, not equivalent depths   
            ql(j,1) = flm + (qnet + 
     &              0.5*(qlat1(1)+qlat2(1))/ewpr1(j,1))/dx  !+ rfi
            hr1(j,1) = hm
C  this is effective surface depth to be divided by delt in INFILT
            if(diag.and. j .eq. jd) write(99,'("ql ",i3, 6g12.4)') j, 
     &  qlat(1), qnet, ewpr1(j,1), flm, hm, h2n
C  
          end do

C!!          hr1(1,1) = 0.      ! not needed
          call infilt ( afac, rfi, hr1, ql, delt, fav, 1)

        end if
C
C
C**9/8/2004**following section copied from above and changed to accomodate impervious channels
        if (chrain) then
          do k = 1, nchan
            do j = 1, nk
	        if( k .eq. 1 .and. fav(j,k) .gt. 0. )then
                widu = min(ewpr1(j,k),rwidth(j,k))
	        else
                widu = rwidth(j,k)
	        end if
              rf(j,k) = rfi * widu      ! effective rain input into flowing water
            end do
          end do
        end if

        dte = 0.
        minflg = .false.
C                                                 start Courant loop (dt, m)...
C------------------------------------------------------------------------------

50      do m = 1, nc

          vcm = 0.

          if (qtest .le. 0.) then
C                                                                       no flow
C------------------------------------------------------------------------------

            do j = 1, nk

              h2(j) = 0.

              do k = 1, nchan
                a2(j,k) = 0.
                q2(j,k) = 0.
                fav(j,k) = rfi
                wpr2(j,k) = width(j,k)
                topw(j,k) = 0.
                if(wool) then
                  ewpr2(j,k) = .05 * width(j,k)
                else
                  ewpr2(j,k) = width(j,k)
                end if
              end do

            end do

            qtot2 = 0.      ! qtot1?
            vstor = 0.
            if(diag) write(99,'("  bypass...")')
C                                                                bypass routing
            go to 80

          end if
C                                               dte = time from last fixed step
C                                                     to end of current dt
          dte = dte + dt
          qub = 0.
C                                            interpolate upstream inflows
          do k = 1, nchan
            qu(k) = qu1(k) + dqu(k) * dte
            qub =  qub + qu(k)
          end do
C                                             lateral inflow is average over dt
          do k = 1, nlu
            qlat(k) = qlat1(k) + dqlat(k) * (dte - 0.5 * dt)
          end do

          if (qub .le. 0.) then
C                                                   upstream boundary depth = 0
C------------------------------------------------------------------------------

            a2(1,1) = 0.
            wpr2(1,1) = width(1,1)
            q2(1,1) = 0.
            a2(1,1) = 0.
            wpr2(1,2) = 0.
            q2(1,2) = 0.
            h2(1) = 0.
            qt(1) = 0.

            if (wool) then

!              ewpr2(1) = 0.05 * width(1,1)
              ewpr2(1,1) = 0.05 * width(1,1)   ! starting value

            else

!              ewpr2(1) = width(1,1)
              ewpr2(1,1) = width(1,1)

            end if

          else
C                                               compute upstream boundary depth
C------------------------------------------------------------------------------
            source = 'ChnUB2'
            hulim = -tol
            if(h2(1) .le. 0.) then
              h2(1) = (qub/width(1,1)/alpha(1,1))**(1./p) ! init. estimate
            end if
C*2/3/2004 index j needs to be set to one (see changes to errhub)
	      j = 1
            call iter (errhub, h2(1), hulim, 100., ierr, source)
            if(h2(1) .lt. 0.) h2(1) = 0.

            if (ierr .gt. 0) then
              call errxit (msg(ms:18), 'no convergence (errhub)')
            else !  if (j .eq. 2) then
              itot = itot - ierr
            end if

            unif = .false.
!            ewpr2(1) = wpr2(1,1)
            ewpr2(1,1) = wpr2(1,1)

            if (.not. cmpd .or. h2(1) .le. thresh(1)) then
C                                                              flow confined to
C                                                                  main channel
!              if (wool) ewpr2(1) = ewprf (h2(1), width(1,1),
              if (wool) ewpr2(1,1) = ewprf (h2(1), width(1,1),
     &                                      wpr2(1,1), wco)

              dadh1 = dadhf (width(1,1), h2(1), co1(1,1))
              dqdh1 = dqdhf (alpha(1,1), wpr2(1,1), p,
     &                         a2(1,1), dadh1, co2(1,1))
              dadh2 = 0.
              dqdh2 = 0.
              qt(1) = 0.

            else

              hh = h2(1) - thresh(1)
              dadh1 = dadhf (tw0(1), h2(1), co1a(1))
              dqdh1 = dqdhf (alpha(1,1), wpr2(1,1), p,
     &                         a2(1,1), dadh1, co2a(1))

              ewpr2(1,2) = wpr2(1,2)*afov(1)                         ! new 10/01
C              if (wool) ewpr2(1,2) = ewprf (hh, width(1,2),
C     &                                      wpr2(1,2), wco)  ! new 10/01
              if (hh .le. hb(1)) then

                dadh2 = dadhf (0., hh, co1hb(1))
                dqdh2 = dqdhf (alpha(1,2), wpr2(1,2), p,
     &                           a2(1,2), dadh2, co2hb(1))
              else

                hhh = hh - hb(1)
                dadh2 = dadhf (twhb(1), hhh, co1(1,2))
                dqdh2 = dqdhf (alpha(1,2), wpr2(1,2), p,
     &                           a2(1,2), dadh2, co2(1,2))

              end if

              qt(1) = qu(1) - q2(1,1)

            end if
C                                                          estimate celerity at
C                                                                upper boundary
            vc2 = 0.                                  ! prev undefined.  1-/02
            vc1 = 0.
            if (dadh1 .ne. 0.) vc1 = dqdh1 / dadh1
            if (dadh2 .ne. 0.) vc2 = dqdh2 / dadh2
            vc = amax1 (vc1, vc2)
            if (vc .gt. vcm) vcm = vc                ! used for courant check

          end if
                             ! end of upper bound calculations
C 
          qlt2 = (rfi - fav(2,1))*ewpr1(2,1) + qlat(1) + qltr(2)
                                                   ! new - net input at node 2
          if (unif .and. qlt2 .gt. 0.) then   ! 
C                                                               zone A solution
C------------------------------------------------------------------------------
            h2a = h1(2,1)
            j = 2
            source = 'ChanHA'
            halim = -20.*tol
            call iter (errha, h2a, halim, 10000., ierr, source)
            if(h2a .lt. 0.) h2a = 0.

            h2(2) = h2a
            dadh = dadhf (width(2,1), h2a, co1(2,1))
            dqdh = dqdhf (alpha(2,1), wpr2(2,1), p,
     &                    a2(2,1), dadh, co2(2,1))
C                                                              compute celerity
            vc = dqdh / dadh
            xc = xc + dt * vc
            if (vc .gt. vcm) vcm = vc

            if (xc .ge. dx) then
C                                                       begin kinematic routing
              unif = .false.

            else

              q2(2,1) = qf (a2(2,1), wpr2(2,1), alpha(2,1), p)

              do j = 3, nk
                h2a = h1(j,1)
                source = 'ChanH2'
                call iter (errha, h2a, -tol, 10000., ierr, source)
                if(h2a .lt. 0.) h2a = 0.

                h2(j) = h2a
                q2(j,1) = qf (a2(j,1), wpr2(j,1), alpha(j,1), p)

              end do
C                                                      bypass kinematic routing
          if(diag) write(99,'(" zoneA ",8g13.5)') 
     &  ewpr1(2,1),rf(2,1),fav(2,1),qlat(1),qlat(2)
              go to 70

            end if
          end if                              ! end zone A solution
C------------------------------------------------------------------------------
C                                                start spatial loop (j, jp1)...

          do j = 1, nk - 1

            jp1 = j + 1
C            h2jp1 = h2(jp1)
            h2jp1 = (h2(j) + h1(jp1,1)) / 2.
C            if(h1(jp1) .le. 0.) h2jp1 = 0.
            source = 'ChanCH'
            hmin = min(-20.*tol,0.)
C  route flow for general case
            call iter (errch2, h2jp1, hmin, 100., ierr, source)

            if (ierr .gt. 0) then
              if(diag)  write(99,'(I2,6g13.4)') 
     &               j,h1(j,1),h1(jp1,1),h2(j),fav(j,1)
              call errxit(msg(ms:18), 'no convergence (errch2)')
            else if (j .eq. jd) then  !  
              itrt = itrt - ierr
              if(diag .and. j .eq. jd) 
     &     write(99,177) j,h1(j,1),h2(j),h1(jp1,1),h2jp1,fav(jp1,1)
 177  format(" J    h1(j,1)      h2(j,1)       h1(jp1,1)      h2(jp)    
     &   fav(jp1)"/I2,6g13.4)
            end if

            if(h1(j,1) .gt. 0. .and. h2jp1 .lt. 0.) then
              if(fav(j,1) .gt. 0.) fav(jp1,1) = fav(jp1,1)+0.5*h2jp1/dt
            end if
            if( h2jp1 .le. 1.e-5 .and. h2(j) .le. 1.e-5) then
              h2(jp1) = 0.
              a2(jp1,1) = 0.
              q2(jp1,1) = 0.
            else
              h2(jp1) = max(h2jp1, 0.)
            end if

            if (.not. cmpd .or. h2jp1 .le. thresh(jp1)) then

              dadh1 = dadhf (width(jp1,1), h2jp1, co1(jp1,1))
              dqdh1 = dqdhf (alpha(jp1,1), wpr2(jp1,1), p,
     &                         a2(jp1,1), dadh1, co2(jp1,1))
              dadh2 = 0.
              dqdh2 = 0.
C                                                     flow in main channel only
              if (cmpd) then

                a2(jp1,2) = 0.
                q2(jp1,2) = 0.
                wpr2(jp1,2) = 0.

              end if

            else
C                                                      flow in overbank section
              hh = h2jp1 - thresh(jp1)

              dadh1 = dadhf (tw0(jp1), h2jp1, co1a(jp1))
              dqdh1 = dqdhf (alpha(jp1,1), wpr2(jp1,1), p,
     &                         a2(jp1,1), dadh1, co2a(jp1))

              if (hh .le. hb(jp1)) then

                dadh2 = dadhf (0., hh, co1hb(jp1))
                dqdh2 = dqdhf (alpha(jp1,2), wpr2(jp1,2), p,
     &                           a2(jp1,2), dadh2, co2hb(jp1))

              else

                hhh = hh - hb(jp1)
                dadh2 = dadhf (twhb(jp1), hhh, co1(jp1,2))
                dqdh2 = dqdhf (alpha(jp1,2), wpr2(jp1,2), p,
     &                           a2(jp1,2), dadh2, co2(jp1,2))

              end if
C                                            compute lateral transfer discharge
C------------------------------------------------------------------------------

60            aout = fav(jp1,1) * 0.25 * (ewpr2(jp1,1) + ewpr1(jp1,1)
     &               + ewpr2(j,1) + ewpr1(j,1))

              dadt = 0.5 * (a2(jp1,1) - a1(jp1,1) +
     &                      a2(j,1) - a1(j,1)) / dt

              dqdx = (omega * (q2(jp1,1) - q2(j,1)) +
     &               (1. - omega) * (q1(jp1,1) - q1(j,1))) / dx

              qterms = qlat(1) + rf(jp1,1) + dqbf  ! + qlat(2)  missing 12/12

              qt(jp1) = qterms - (dadt + dqdx + aout)   ! this is used in SED calcs

            end if
C                                                             estimate celerity
C------------------------------------------------------------------------------
            vc2 = 0.
            vc1 = 0. 
            if (dadh1 .ne. 0.) vc1 = dqdh1 / dadh1
            if (dadh2 .ne. 0.) vc2 = dqdh2 / dadh2
            vc = amax1 (vc1, vc2)
C                                           vcm is max. celerity for current dt
            if (vc .gt. vcm) vcm = vc

          end do
C                                               ...end of spatial loop (j, jp1)
C------------------------------------------------------------------------------

C                                                       check Courant condition
C------------------------------------------------------------------------------

          if (vcm * dt / dx .gt. 1. .and. 
     &             (dt .gt. DTMIN) .and. (.not. minflg)) then
C                                                                compute new dt
            dtprev = dt
            dtr = delt - dte + dt
            nc = int (2. * vcm * dtr / dx)
            if (nc .eq. 0) nc = 1
            dt = dtr / float (nc)

            if (dt .lt. DTMIN) then
              nc = int (dtr / DTMIN)
              if (nc .eq. 0) nc = 1
              dt = dtr / float (nc)
              minflg = .true.               ! don't try to subdivide further
                                            ! for the reaminder of this delt
            end if

            if (cour) then
C                                                 recompute profile with new dt
              iretn = iretn + 1
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
              dt = delt

            end if
          end if

70        continue

          if (.not. bflow) then
C                                                    compute infiltrated volume
C------------------------------------------------------------------------------
C
            vstor = 0.
            do j = 2, nk
              jm = j-1
C  for diagnos. bal:
              vstor = vstor + 0.5*dx*(a2(j,1) + a2(jm,1))
C
              avew = 0.25*(ewpr2(jm,1) + ewpr2(j,1) + ewpr1(j,1) 
     &              + ewpr1(jm,1))
              avrw = 0.5*(rwidth(j,1) + rwidth(jm,1))
              dwfr = max((avrw - avew),0.)
C**9/8/2004**IF block added for impervious channel
	        if( fav(j,1) .gt. 0. )then
                dfir = rfi*dwfr
	        else
	          dfir = 0.
	        end if
              dwf  = avew
              dfif = fav(j,1)*dwf
              test = dt*(dfif + dfir)
C     &             
              test2 = 0.
C
              if(ewpr2(j,1) .gt. width(j,1)) then
                dwb = dt*(0.5*fav(jm,1)*(width(jm,1) + width(j,1))
     &    + max((avrw - 0.5*(width(j,1)+width(jm,1))),0.)*rfi)
              else
                dwb = test
              end if
C
              if(cmpd) then
                avew2 = 0.25*(ewpr1(j,2) + ewpr1(jm,2) + ewpr2(j,2)
     &                  + ewpr2(jm,2))
                avrw2 = 0.5*(rwidth(j,2) + rwidth(jm,2))
C                test2 = dt*( 0.5 * (fav(j-1,2) * ewpr2(j-1,2)
C     &      + fav(j,2)*ewpr2(j,2)) + (rwidth(j,2) - ewpr2(j,2))*rfi)
C**9/8/2004**IF block added for impervious channel
	          if( fav(j,2) .gt. 0. )then
                  dwfr = dwfr + avrw2*(1.-afov(j))
                  dfir = dfir + rfi*dwfr
	           else
	             dwfr = 0.
	           end if
                 dwf = dwf + avew2
                dfif = dfif + fav(j,2)*avew2
                test2 = dt*(fav(j,2)*avew2 + avrw2*rfi)
              end if  
C
     &                  
              if (h2(j)  .le. 0.) then   !+ h2(j-1)
C                                      dry: no depth at node j
                sup =    0.5 *( a1(j,1) + a1(j-1,1) +
     &             dt *(rf(j,1) + rf(j-1,1))) + dt * (qlat(1) + qltr(j))
     &             + rfi*(rwidth(j,1) - 0.5*(ewpr1(j,1)+ewpr2(j,1)))
                if(sup .lt. test) test = sup
              end if
C
              if(cmpd .and. (h2(j) .le. thresh(j) .and. h2(j-1) .le.
     &             thresh(j-1))) then
C                                    overbank has no flow depth at j,j-1
                sup2 = 0.5 * (a1(j,2) + a1(j-1,2)) +
     &                     dt * (rfi*rwidth(j,2) +  qlat(2))
                if(sup2 .lt. test2) test2 = sup2
              end if
C
C              dvf = test + test2
C
              tfr = tfr + dfir*dt !sum infil on rain only surface
              swr = swr + dwfr*dt
              tfq = tfq + dfif*dt !sum infil under flows
              swq = swq + dwf*dt
              vtn = vtn + dt

              qbal(7) = qbal(7) + test 
              if(cmpd) qbal(6) = qbal(6) + test2
              qbal(13) = qbal(13) + dwb
              sov = sov + test2*dx
C
            end do
          end if
C------------------------------------------------------------------------------

80        continue
          timec = time + dte
C                                                     ending time of current dt

          qtot2 = q2(nk,1)
          if (cmpd) qtot2 = qtot2 + q2(nk,2)
C                                                          peak total flow rate
          if (qtot2 .gt. qp) then

            qp = qtot2
            tp = timec

          end if

          qbal(10) = qbal(10) + 0.5 * (qtot1 + qtot2) * dt
          qtot1 = qtot2
C - - - - - - - - -
            if(diag) then
              tmi = (timec)/60.                     ! time at end of time step
              write(99,'("  Diagnostic arrays:")')
              write(99,992)'h2  ',(h2(j),j=1,nk)
              write(99,992)'q2  ',(q2(j,1),j=1,nk)
C              write(99,992)'wp2 ',(ewpr2(j,2),j=1,nk)
C              write(99,992)'fav:', (fav(j,1),j=1,4)
C              write(99,992)'ewpr',(ewpr2(j,1),j=1,4)
C              write(99,'(" iter: ",2i5, g12.4)') itot, itrt
              tbal = vup + vlat + smrn - vstor - qbal(7)*dx - qbal(10)
              write(99,993)jd,ul(units),ul(units),ul(units),tmi,rf(jd,1)
     &         ,fav(jd,1), q2(jd,1), ql(jd,1), qlat(1) !, sov
             write(99,991)jd, itrt, qbal(7)*dx, vup, vstor, tbal,qlat(1)
     &          ,qlat(2),smrn
  993  format(1x,I2," t(min)      rf(",a1,"/s)     f(",a1,"/s)      q2("     
     &,a1,"^3)     ql         qlat"/(8g12.4))
  992  format(1x,a4,(8g12.4))
  991  format(' j itr    Qb7          SumIn        stor         bal     
     &  QL1         QL2'/1x,i2,i4,(8g12.4))
            end if
C - - - - - - - - - - - -
C                                                                      sediment
C------------------------------------------------------------------------------
          if (sed) then   ! prepare variables for call to kinsed:

            do j = 1, nk
C
              if(nm .le. 1) then
                ql(j,1) = qlat2(1)
                ql(j,2) = qltr(j)
              else    ! reseparate ql for use in kinsed:
                ql(j,1) = qlat1r
                ql(j,2) = qlat2(1) - qlat1r
              end if
C                                                             compute top width
              if (.not. cmpd .or. h2(j) .le. thresh(j)) then

                topw(j,1) = width(j,1) + h2(j) * co1(j,1)
                topw(j,2) = 0.

              else

                hh = h2(j) - thresh(j)
                topw(j,1) = tw0(j) + hh * co1a(j)
C
                if (hh .le. hb(j)) then

                  topw(j,2) = hh * co1hb(j)

                else

                  hhh = hh - hb(j)
                  topw(j,2) = twhb(j) + hhh * co1(j,2)

                end if
              end if
            end do
C                                                                        kinsed
            call kinsed ( i, dt, dte, ql, a2, q2,
     &                   topw, rfi, qup, qt, timec)

          end if
C                                      reassign variables for next time substep
C------------------------------------------------------------------------------

          qtest = 0.

          do j = 1, nk

            h1(j,1) = h2(j)

            if (cmpd) then

              h1(j,2) = h2(j) - thresh(j)
              if (h1(j,2) .lt. 0.) h1(j,2) = 0.

            end if

            do k = 1, nchan

              a1(j,k) = a2(j,k)
              wpr1(j,k) = wpr2(j,k)
              q1(j,k) = q2(j,k)
              qtest = qtest + q2(j,k)
              ewpr1(j,k) = ewpr2(j,k) ! <-

            end do
          end do
        end do
C                                                ...end of Courant loop (dt, m)
C------------------------------------------------------------------------------

        if (cour) then
        
          minflg = .false.
C                                        next dt based on previous max celerity
          nc = nint (vcm * delt / dx)
          if (nc .eq. 0) nc = 1
          dt = delt / float (nc)
          dte = 0.

          if (dt .lt. DTMIN) then
            nc = int (delt / DTMIN)
            if (nc .eq. 0) nc = 1
            dt = delt / float (nc)
            minflg = .true.            ! don't subdivide further
          end if

        end if
C                                                                 store outflow
        call store (outq(1), i, q2(nk,1))
        if (cmpd) call store (outq(2), i, q2(nk,2))
C                                                               store rain rate
        do k = 1, nchan

          qu1(k) = qu2(k)

        end do

        do k = 1, nlu

          qlat1(k) = qlat2(k)

        end do

!**
        IF( nprint .EQ. 5)THEN
          DO j = 1, nk
            IF( a2(j,1) .GT. 0. )THEN
              v2 = q2(j,1)/a2(j,1)
            ELSE
              v2 = 0.
            END IF
            WRITE( xfile,'(F10.2,",",I3,",",G13.6,",",G13.6,",",G13.6)')
     &      time/60., j, q2(j,1), v2, h2(j)
          END DO
          IF( i .LE. 1000 )THEN
            q(i) = q2(nk,1)
            v(i) = v2
            h(i) = h2(nk)
          END IF
        END IF  
!**

C
        if(nele .le. 5) call progss (msg(ms:18))
C 
      end do
C                                                ...end of fixed time step loop
      if(diag .and. cour) then
        write(99,'(/I4," resets of dt by courant criterea. ")') iretn
      end if
C
      if(vtn .gt. 0.) then
       avwr = swr/vtn
       avwq = swq/vtn
      else
       avwr = 0.
       avwq = 0.
      end if
C
      avwi = 0.
      if(qbal(7) .gt. 0.)
     &  avwi = (avwr*tfr + avwq*tfq)/qbal(7)  ! weighted mean width of infil 
      avai = avwi*xleng
C      dina = (tfq+tfr)*dx/(avwi*xleng)*conv
C      write(77,'(" avwid (m): ",3f7.3,5g12.4)')
C     &          avai, avwi, avwq, barea, tfq, qbal(1), qbal(7), qbal(13)

C------------------------------------------------------------------------------

C                                                  finish sediment computations
      if (sed) call sedfin
C                                                                 inflow volume
      vin = qu2(1) + (qlat2(1) + qlat2(2)) * xleng + 2. * vin
      if (cmpd) vin = vin + qu2(2)
C*2/3/2004 Add baseflow to inflow volume
      vin = (vin / 2.) * delt + qbal(3)
C                                                                outflow volume
C*2/2/2004 double-counting baseflow
C      qbal(10) = qbal(10) + delt * float (limit -1) * qbf
C                                                            infiltrated volume
      if (.not. bflow) then
        qbal(7) = max (qbal(7) * dx, 0.)
        if(cmpd) qbal(6) = max (qbal(6) * dx, 0.)
      end if
      qbal(13) = max (qbal(13) * dx, 0.)

      do j = 2, nk - 1
C                                                                storage volume
        qbal(9) = qbal(9) + a2(j,1)
        if (cmpd) qbal(9) = qbal(9) + a2(j,2)

      end do
C*2/2/2004 subtract initial baseflow profile
      qbal(9) = dx * (a2(1,1) + 2. * qbal(9) + a2(nk,1)) / 2. - qbstor
      if (qbal(9) .lt. 0.) qbal(9) = 0.
C                                                         pass values to writer
C------------------------------------------------------------------------------

      call qwrt (id, nchan, msg(ms:18), qp, tp, vin,
     &           coarea, outq, ltyp, outr, qrmax)

      IF( nprint .EQ. 5)THEN
        WRITE( xfile, '("T(min), Q(nk), V(nk), H(nk)")' )
        DO i = 2, MIN( limit, 1000 )
          time = FLOAT( i - 2 ) * delt
          WRITE( xfile,'(F10.2,",",G13.6,",",G13.6,",",G13.6)')
     &    time/60., q(i), v(i), h(i)
        END DO
        CLOSE( xfile ) 
      END IF

      return

      end

C------------------------------------------------------------------------------


      subroutine errhub (hub, ferr, dferr)


C  Computes the residual and derivative of the depth hub in the uniform
C  discharge equation for a trapezoidal open channel.
C*2/3/2004 Replace first dimension of a2, wpr2, etc. with index 'j' so that
C          depth can be computed for any spatial node.  This is necessary for
C          computing the initial baseflow profile when width(1,1) ne width(nk,1)

      logical cmpd, wool

      integer j, jp1

      common /chan/ a1(20,2), a2(20,2), q1(20,2), q2(20,2), fav(20,2),
     &              wpr1(20,2), wpr2(20,2), co1(20,2), co2(20,2),
!     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20),
!     &              ewpr2(20), a0(20), wpr0(20), tw0(20), hb(20),
     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20,2),
     &              ewpr2(20,2), a0(20), wpr0(20), tw0(20), hb(20),
     &              ahb(20), wprhb(20), twhb(20), co1hb(20), co2hb(20),
     &              qlat(2), j, jp1, omega, qub, dt, dx, cmpd, dqbf,
     &              width(20,2), alpha(20,2), thresh(20), p, wool, wco,
     &              rf(20,2), qltr(20), afov(20), rwidth(20,2)


      if (.not. cmpd .or. hub .le. thresh(1)) then

        a2(j,1) = areaf (0., hub, width(j,1), co1(j,1))
        wpr2(j,1) = cwprf (width(j,1), hub, co2(j,1))
        q2(j,1) = qf (a2(j,1), wpr2(j,1), alpha(j,1), p)
        dadh = dadhf (width(j,1), hub, co1(j,1))
        dq2adh = dqdhf (alpha(j,1), wpr2(j,1), p,
     &                    a2(j,1), dadh, co2(j,1))
        dq2bdh = 0

        if (cmpd) then

          a2(j,2) = 0.
          wpr2(j,2) = 0.
          q2(j,2) = 0.

        end if

      else

        hh = hub - thresh(j)
        a2(j,1) = areaf (a0(j), hh, tw0(j), co1a(j))
        wpr2(j,1) = cwprf (wpr0(j), hh, co2a(j))
        q2(j,1) = qf (a2(j,1), wpr2(j,1), alpha(j,1), p)
        dadh = dadhf (tw0(j), hh, co1a(j))
        dq2adh = dqdhf (alpha(j,1), wpr2(j,1), p,
     &                    a2(j,1), dadh, co2a(j))

        if (hh .le. hb(j)) then

          wpr2(j,2) = cwprf (0., hh, co2hb(j))
          a2(j,2) = areaf (0., hh, 0., co1hb(j))
          dadh = dadhf (0., hh, co1hb(j))
          dq2bdh = dqdhf (alpha(j,2), wpr2(j,2), p,
     &                      a2(j,2), dadh, co2hb(j))

        else

          hhh = hh - hb(j)
          wpr2(j,2) = cwprf (wprhb(j), hhh, co2(j,2))
          a2(j,2) = areaf (ahb(j), hhh, twhb(j), co1(j,2))
          dadh = dadhf (twhb(j), hhh, co1(j,2))
          dq2bdh = dqdhf (alpha(j,2), wpr2(j,2), p,
     &                      a2(j,2), dadh, co2(j,2))

        end if

        q2(j,2) = qf (a2(j,2), wpr2(j,2), alpha(j,2), p)

      end if

      ferr = q2(j,1) - qub
      if (cmpd) ferr = ferr + q2(j,2)
      dferr = dq2adh + dq2bdh

      return

      end

C------------------------------------------------------------------------------


      subroutine errch2 (h2jp1, ferr, dferr)


C  Computes the residual and derivative of the depth h2jp1 in a finite
C  difference of kinematic flow equation for a trapezoidal open channel.


      logical cmpd, wool

      integer j, jp1

      common /chan/ a1(20,2), a2(20,2), q1(20,2), q2(20,2), fav(20,2),
     &              wpr1(20,2), wpr2(20,2), co1(20,2), co2(20,2),
!     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20),
!     &              ewpr2(20), a0(20), wpr0(20), tw0(20), hb(20),
     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20,2),
     &              ewpr2(20,2), a0(20), wpr0(20), tw0(20), hb(20),
     &              ahb(20), wprhb(20), twhb(20), co1hb(20), co2hb(20),
     &              qlat(2), j, jp1, omega, qub, dt, dx, cmpd, dqbf,
     &              width(20,2), alpha(20,2), thresh(20), p, wool, wco,
     &              rf(20,2), qltr(20), afov(20), rwidth(20,2)

! New declaration

      real dewpr(2)

      if (.not. cmpd .or.  h2jp1 .le. thresh(jp1)) then

        a2(jp1,1) = areaf (0., h2jp1, width(jp1,1), co1(jp1,1))
        wpr2(jp1,1) = cwprf (width(jp1,1), h2jp1, co2(jp1,1))
        q2(jp1,1) = qf (a2(jp1,1), wpr2(jp1,1), alpha(jp1,1), p)
        da2adh = dadhf (width(jp1,1), h2jp1, co1(jp1,1))
        dq2adh = dqdhf (alpha(jp1,1), wpr2(jp1,1), p,
     &                    a2(jp1,1), da2adh, co2(jp1,1))
        da2bdh = 0.
        dq2bdh = 0.

        if (cmpd) then

          a2(jp1,2) = 0.
          wpr2(jp1,2) = 0.
          q2(jp1,2) = 0.
          ewpr2(jp1,2) = 0.
          dewpr(2) = 0.

        end if

        if (wool) then

!          ewpr2(jp1) = ewprf (h2jp1, width(jp1,1), wpr2(jp1,1), wco)
          ewpr2(jp1,1) = ewprf (h2jp1, width(jp1,1), wpr2(jp1,1), wco)

! New function

          dewpr(1) = dewprf (h2jp1, width(jp1,1), ewpr2(jp1,1),
     &                       co2(jp1,1), wco)

!          dewpr = (width(jp1,1) + 2. * h2jp1 * co2(jp1,1))
!     &            / (wco * sqrt(width(jp1,1)))

        else

!          ewpr2(jp1) = wpr2(jp1,1)
          ewpr2(jp1,1) = wpr2(jp1,1)
          dewpr(1) = co2(jp1,1)

        end if

      else       ! compound and h .gt. thresh:

        hh = h2jp1 - thresh(jp1)
        wpr2(jp1,1) = cwprf (wpr0(jp1), hh, co2a(jp1))
        a2(jp1,1) = areaf (a0(jp1), hh, tw0(jp1), co1a(jp1))
        q2(jp1,1) = qf (a2(jp1,1), wpr2(jp1,1), alpha(jp1,1), p)
        da2adh = dadhf (tw0(jp1), hh, co1a(jp1))
        dq2adh = dqdhf (alpha(jp1,1), wpr2(jp1,1),
     &           p, a2(jp1,1), da2adh, co2a(jp1))

        if (hh .le. hb(jp1)) then

          wpr2(jp1,2) = cwprf (0., hh, co2hb(jp1))
          a2(jp1,2) = areaf (0., hh, 0., co1hb(jp1))
          da2bdh = dadhf (0., hh, co1hb(jp1))
          dq2bdh = dqdhf (alpha(jp1,2), wpr2(jp1,2),
     &             p, a2(jp1,2), da2bdh, co2hb(jp1))

          if (wool) then
            ewpr2(jp1,2) = ewprf(hh, width(jp1,2), wpr2(jp1,2), wco)
            dewpr(2) = dewprf(hhh, width(jp1,2), ewpr2(jp1,2),
     &                         co1hb(jp1), wco)
          else
            ewpr2(jp1,2) = wpr2(jp1,2)
            dewpr(2) =  co1hb(jp1)
          endif

        else

          hhh = hh - hb(jp1)
          wpr2(jp1,2) = cwprf (wprhb(jp1), hhh, co2(jp1,2))
          a2(jp1,2) = areaf (ahb(jp1), hhh, twhb(jp1), co1(jp1,2))
          da2bdh = dadhf (twhb(jp1), hhh, co1(jp1,2))
          dq2bdh = dqdhf (alpha(jp1,2), wpr2(jp1,2),
     &             p, a2(jp1,2), da2bdh, co2(jp1,2))

          if (wool) then
            ewpr2(jp1,2) = ewprf(hhh, width(jp1,2), wpr2(jp1,2), wco)
            dewpr(2) = dewprf(hhh, width(jp1,2), ewpr2(jp1,2),
     &                        co2(jp1,2), wco)
          else
            ewpr2(jp1,2) = wpr2(jp1,2)
            dewpr(2) = co2(jp1,2)
          endif

        end if

        q2(jp1,2) = qf (a2(jp1,2), wpr2(jp1,2), alpha(jp1,2), p)
!        ewpr2(jp1) = wpr2(jp1,1)
        ewpr2(jp1,1) = wpr2(jp1,1)
        dewpr(1) = co2a(jp1)

      end if

      comega = 1. - omega
      avew = 0.5*(omega*(ewpr2(j,1) + ewpr2(jp1,1)) + 
     &            comega*(ewpr1(j,1) + ewpr1(jp1,1)))
C      aout = 0.5 * (fav(jp1,1) * 0.5 *(ewpr2(jp1,1) + ewpr1(jp1,1))
C     &              + fav(jp1,2) * afov(jp1)* rwidth(jp1,2) !(ewpr2(jp1,2) + ewpr1(jp1,2))
C     &            + fav(j,1) * 0.5 * (ewpr2(j,1) + ewpr1(j,1))
C     &	            + fav(j,2) * afov(j)*rwidth(j,2)) !
      aout = fav(jp1,1) * avew + fav(jp1,2)*afov(jp1)*0.5*(rwidth(j,2) 
     &                       + rwidth(jp1,2))

      dadt = 0.5 * (a2(jp1,1) - a1(jp1,1) + a2(j,1) - a1(j,1) +
     &       a2(jp1,2) - a1(jp1,2) + a2(j,2) - a1(j,2)) / dt

      dqdx = (omega * (q2(jp1,1) - q2(j,1) + q2(jp1,2) - q2(j,2)) +
     &       comega * (q1(jp1,1) - q1(j,1) + q1(jp1,2) - q1(j,2))) / dx
C
      qterms = qlat(1) + qlat(2) + rf(jp1,1) + rf(jp1,2)*afov(jp1) +dqbf  !+ rf(jp1,2)

      ferr = dadt + dqdx + aout - qterms    ! 

      dferr = 0.5 * (da2adh + da2bdh) / dt + omega * (dq2adh +
     &        dq2bdh) / dx + 0.5 * omega *(fav(jp1,1) * dewpr(1)) 
C
C      if(jp1 .eq. 2) then
C        write(99,109) jp1, h2jp1, q2(jp1,1), dadt, dqdx, q2(j,1),
C     &  q1(j,1), q1(jp1,1)
C        write(99,108),a2(jp1,1),q1(jp1,1),q2(j,1)
C        write(99,110) fav(j,1),fav(jp1,1),rf(jp1,1),ewpr2(jp1,1),qlat(1)
C      end if
  108 format(" h,q1&2",8g12.4)
  109 format(" itj ",i2,7g13.4)
C  110 format(" f&ew",8g12.4)
      return

      end

C------------------------------------------------------------------------------


      subroutine errha (hza, ferr, dferr)

C  Computes the residual and derivative for depth hza in an equation for flow
C  in zone A of the solution space (main channel only).

      logical cmpd, wool

      integer j, jp1
! ewpr
      common /chan/ a1(20,2), a2(20,2), q1(20,2), q2(20,2), fav(20,2),
     &              wpr1(20,2), wpr2(20,2), co1(20,2), co2(20,2),
!     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20),
!     &              ewpr2(20), a0(20), wpr0(20), tw0(20), hb(20),
     &              co1a(20), co2a(20), h1(20,2), h2(20), ewpr1(20,2),
     &              ewpr2(20,2), a0(20), wpr0(20), tw0(20), hb(20),
     &              ahb(20), wprhb(20), twhb(20), co1hb(20), co2hb(20),
     &              qlat(2), j, jp1, omega, qub, dt, dx, cmpd, dqbf,
     &              width(20,2), alpha(20,2), thresh(20), p, wool, wco,
     &              rf(20,2), qltr(20), afov(20), rwidth(20,2)


      a2(j,1) = areaf (0., hza, width(j,1), co1(j,1))
      dadh = dadhf (0., hza, co1(j,1))
      wpr2(j,1) = cwprf (width(j,1), hza, co2(j,1))
      dwpr2 = co2(j,1)
      dewpr2 = dwpr2

      if (wool) then
!        ewpr2(j) = ewprf (hza, width(j,1), wpr2(j,1), wco)
        ewpr2(j,1) = ewprf (hza, width(j,1), wpr2(j,1), wco)

! new function call

        dewpr2 = dewprf (hza, width(j,1), ewpr2(j,1), co2(j,1), wco)

!        denom = wco * sqrt (width(j,1))

!        if (hza .lt. denom) then
C                                       compute derivative of ewpr2 w/resp to h
C------------------------------------------------------------------------------

!          if (ewpr2(j) .eq. .05 * width(j,1)) then

!              dewpr2 = 0.

!          else

!              dewpr2 = 2. * hza * co2(j,1) / denom

!          end if

!        end if

      else

!        ewpr2(j) = wpr2(j,1)
        ewpr2(j,1) = wpr2(j,1)
      end if

       ferr = a2(j,1) - a1(j,1) - dt * ((qlat(1)  + rf(j,1)) - !+ qlat(2)
     &       0.5 * fav(j,1) * (ewpr1(j,1) + ewpr2(j,1)))
C                                                        derivative w/resp to h
      dferr = dadh + dt * fav(j,1) * 0.5 * dewpr2

      return

      end

C------------------------------------------------------------------------------


      function areaf (a0, h, tw0, co1)

      areaf = a0

      if(h .ge. 0.) then
        areaf = a0 + h * tw0 + 0.5 * h * h * co1
      end if

      return

      end



      function cwprf (wp0, h, co2)

      if(h .ge. 0.) then
        cwprf = wp0 + h * co2
      else
        cwprf = wp0
      end if

      return

      end



      function qf (ar, wp, al, p)

      if(ar .ge. 0.) then
        qf = al * ar**p / wp**(p - 1.)
      else
        qf = 0.
      end if

      return

      end



      function dadhf (tw0, h, co1)

      if(h .ge. 0.) then
        dadhf = tw0 + h * co1
      else
        dadhf = tw0
      end if

      return

      end



      function dqdhf (al, wp, p, ar, dadh, co2)
C
      if(ar .gt. 0.) then
        dqdhf = al * (wp**(p - 1.) * p * dadh * ar**(p - 1.) -
     &          ar**p * (p - 1.) * co2 * wp**(p - 2.)) /
     &          wp**(2. * p - 2.)
      else
        dqdhf = 0.
      end if

      return

      end


      function ewprf (h, w, wp, wco)

      if(h .ge. 0.) then
        ewprf = amax1 (amin1 (h / (wco * sqrt(w)), 1.) * wp, .05 * w)
      else
        ewprf = .05*w
      end if

      return

      end


!     New function: derivative of eff wetted perim (ewpr) w/resp to h

      function dewprf (h, w, ewpr, co2, wco)

      if (ewpr < .05 * w) then
        dewprf = 0.
      else
        denom = wco * sqrt(w)
        if (h .lt. denom) then
          hu = MAX(h,0.)
          dewprf = 2. * hu * co2 / denom
        else
          dewprf = co2
        end if
      end if

      return

      end
