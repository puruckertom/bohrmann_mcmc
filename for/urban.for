C   Code type: FORTRAN subroutine

C   Compiler: Microsoft FORTRAN v5.1

C   Programmed by: C. Unkrich

C   Date: 2/97

C   Description:

C     Represents a composite urban element consisting of up to six overland
C     flow areas representing various combinations of pervious and impervious
C     surfaces contributing to a paved, crowned street.  The element can
C     recieve upstream inflow (into the sreet) but not lateral inflow.  The
C     relative proportions of the six overland flow areas are specified as
C     fractions of the total element area.  All six areas do not need to be
C     present, except for connecting areas if the corresponding indirectly
C     connected area is specified.  The element is modeled as rectangular:

C      +--------------+--------------+--------------+--------------+
C      |              |              |  Indirectly  |              |
C      |              |              |  Connected   |  Indirectly  |
C      |  Directly    |  Directly    |  Impervious  |  Connected   |
C      |  Connected   |  Connected   +--------------+  Pervious    |
C      |  Pervious    |  Impervious  |              |              |
C      |              |              |  Connecting  +--------------+
C      |              |              |  Pervious    |  Connecting  |
C      |              |              |              |  Impervious  |
C      |              |              |              |              |
C      =============================================================
C      |                                                           |
C      |          Street (half)                                    |
C      |                                                           |
C      +----------------------------- Street Crown ----------------+



C     For output, volumetric rainfall rate from upstream areas for each time
C     step are added and stored.  Dividing by contributing area gives the
C     integrated rainfall rate over the area contributing to the element.

C     Upstream inflow cannot originate from a compound channel.

C     Sediment production is assumed to be negligible, and sediment throughflow
C     is not simulated.
C------------------------------------------------------------------------------

C   Arguments:

C     qbal(15)      real        array with element volume balance components;
C                               qbal(1)  = area             (sq.m/sq.ft),
C                               qbal(2)  = rainfall         (cu.m/cu.ft),
C                               qbal(7)  = infiltration     (cu.m/cu.ft),
C                               qbal(8)  = interception     (cu.m/cu.ft), 
C                               qbal(9)  = storage          (cu.m/cu.ft),
C                               qbal(10) = outflow          (cu.m/cu.ft),

C     diag          log         .true. = write diagnostic info to unit 99,

C     dtm(4)        real        smallest dt to satisfy the Courant condition
C                               during each quarter of the simulation,
C------------------------------------------------------------------------------
C
C  other run variables are defined in the global modules in KINMOD.FOR file
C
C -----------------------------------------------------------------------------
C
C  Input Block (input file labels in parenthesis):

C     id (ID)            int*4   identifier (up to 6 digits),

C     idup (UP)          int*4   identifiers of up to 10 upstream
C                                contributing elements,

C     xleng (LE)         real    length of element in the downslope
C                                direction exclusive of the street (m or ft),

C     width (WI)         real    TWO widths listed in order of:
C                                (1) total width of element,
C                                (2) width of street half (m or ft),

C     slope (SL)         real    TWO slopes listed in order of:
C                                (1) average slope of all overland flow areas,
C                                (2) street slope,

C     rs (MA or CH)      real    THREE Manning or Chezy discharge coefficients
C                                listed in order of:
C                                (1) impervious surfaces,
C                                (2) pervious surface,
C                                (3) street,

C     zc (CS)            real    lateral street slope from crown,

C     fract(1) (DCP)     real    directly connected pervious
C     fract(2) (DCI)     real    directly connected impervious
C     fract(3) (ICI)     real    indirectly connected impervious
C     fract(4) (CP)      real    connecting pervious
C     fract(5) (ICP)     real    indirectly connected pervious
C     fract(6) (CI)      real    connecting impervious

C     dintr, dintri (IN) real    TWO interception depths listed in order of:
C                                (1) pervious surfaces,
C                                (2) impervious surface,
C------------------------------------------------------------------------------

C  Subroutines/functions:

C     infilt        (infilt.for)
C     infil0            "
C     iter          (miscel.for)
C     progss        (K2RunW.for)
C     reader        (reader.for)
C     getr4              "
C     getstr             "
C     geti4              "
C     new           (clerk.for)
C     old                "
C     get                "
C     store              "
C     geta               "
C     stora              "
C     qwrt          (writer.for)
C     errxit        (K2RunW.for)
C------------------------------------------------------------------------------

C      subroutine urban (qbal, diag, dtm)  !, delt, limit, clen, units, cour
C**  11/04/03
      subroutine urban (dtm)

      use multip
      use runpars
C**  
      use itpars
      use elpars    ! need to get id into this module for use in EVENT
C                                                                     arguments
C------------------------------------------------------------------------------
C** 11/04/03
C      logical diag   !, cour

      dimension  dtm(4)  !qbal(15),

C                                                               local variables
C------------------------------------------------------------------------------

      logical :: unif, infil, intcep, upflow, addq, first, lwarn

      integer ierr, m, nk, j, jp1, nd, k, nu, nl, ir, ms, upq,outc(5,2),
     &        outq(2), i, idup, upr, outr, latq, out999, up888 !units, id,  
C**
      character(LEN=2),dimension(5):: attr
      character(len=6) :: source
C** 22/4/03
      character(LEN=2) :: chl
      character(len=18) :: msg
C** 4/03  0. = zro
      real :: zro = 0.
      real, dimension(5) :: wso = (/0.,0.,0.,0.,0./)

      dimension t(1000), r(1000), z(1000), upq(10), slope(2), fav(20),
     &          upr(10), idup(10), fract(6), alph(3), ql(20), afac(20),
     &          h_avg(20), q_sup(20)  !4/03 remove wso(5), outc(5,2) attr here

C                                                                   common data

      common /vars/ a1(20), a2(20), q1(20), q2(20), wpr1(20), wpr2(20),
     &              co1, co2, h1(20), h2(20), j, jp1, omega, qub, dt,
     &              dx, alpha, p, rfi, qlav


C------------------------------------------------------------------------------

      data attr /'S1', 'S2', 'S3', 'S4', 'S5'/

      external errpl, errsub , errst, errsza

      tol = 1.e-6
      ltyp = 6
      lwarn = .true.
      do j = 1, 13
        qbal(j) = 0.
      end do

C                                                          get input parameters
C------------------------------------------------------------------------------

C                                                                    IDENTIFIER
      call geti4 ('ID', 0, id, ierr)

      if (ierr .eq. 3) then
        call errxit ('URBAN', 'invalid identifier (ID)')
      else if (ierr .ge. 1) then
        call errxit ('URBAN', 'missing identifier (ID)')
      end if

      msg(1:16) = typname(6)
      write (msg(15:18), '(i4)') id

      do ms = 1, 5
        if (msg(ms:ms) .ne. ' ') go to 1
      end do

1     continue
      if(nele .gt. 5) call progss (msg(ms:18))

C                                                                        LENGTH
      call getr4 ('L', 0, xleng, ierr)
      if (ierr .gt. 0) call errxit (msg(ms:14), 'length (L) not found')

C                                                                         WIDTH
      call getr4 ('W', 0, width, ierr)
      if (ierr .gt. 0) call errxit (msg(ms:14), 'width (W) not found')

      call getr4 ('W', 2, street, ierr)
      if (ierr .gt. 0)
     &call errxit (msg(ms:14), 'street (1/2) width not found')


C                                                                         SLOPE
      call getr4 ('SL', 1, slope(1), ierr)
      if (ierr .gt. 0) call errxit
     &                 (msg(ms:14), 'overland flow slope not found')

      call getr4 ('SL', 2, slope(2), ierr)
      if (ierr .gt. 0) call errxit
     &                 (msg(ms:14), 'street slope not found')

C                                                                      MANNING?
      call getr4 ('MA', 1, rs, ierr)    ! impervious surfaces

      if (ierr .eq. 0) then

        p = 1.66667
        fac = 1.
        if (units .eq. 2) fac = 1.49
        alph(1) = fac * sqrt (slope(1)) / (rs * rmn)  ! impervious alph

        call getr4 ('MA', 2, rs, ierr)     ! pervious surfaces
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), '3 roughness values must be specified')
        alph(2) = fac * sqrt (slope(1)) / (rs * rmn)

        call getr4 ('MA', 3, rs, ierr)     ! street
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), '3 roughness values must be specified')
        alph(3) = fac * sqrt (slope(2)) / (rs * rmn)

      else
C                                                                        CHEZY?
        p = 1.5

        call getr4 ('CH', 1, rs, ierr)
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), 'no roughness values found')
        alph(1) = rs * rmn * sqrt (slope(1))

        call getr4 ('CH', 2, rs, ierr)
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), '3 roughness values must be specified')
        alph(2) = rs * rmn * sqrt (slope(1))

        call getr4 ('CH', 3, rs, ierr)
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), '3 roughness values must be specified')
        alph(3) = rs * rmn * sqrt (slope(2))

      end if
C                                                                  INTERCEPTION
      call getr4 ('IN', 1, dintr, ierr)

      if (ierr .gt. 0) dintr = 0.
      if (dintr .gt. 0.) then

        dintr = rmin * dintr / conv
C                                                                  SURFACE COVER
        call getr4 ('CA', 0, cover, ierr)

        if (ierr .eq. 0 .and. cover .le. 0.) then
        call errxit (msg(ms:14),'canopy cover fraction (CA) must be >0')
        else if(ierr .gt. 0) then
          cover = 1.0
        end if
C
      end if

      call getr4 ('IN', 2, dintri, ierr)

      if (ierr .gt. 0) dintri = 0.

      dintri = rmin * dintri / conv
C                                                               GAGE
      ngage = 0

      call geti4 ('GA', 0, ngage, ierr)

      if (ierr .gt. 0) then
C                                                          Get X, Y coordinates
        call getr4 ('X', 0, x, ierr)
        if (ierr .gt. 0) x = 0.

        call getr4 ('Y', 0, y, ierr)
        if (ierr .gt. 0) y = 0.

      end if
C                                                                CROWN SLOPE ZC
      call getr4 ('CS', 0, zc, ierr)
      if (ierr .gt. 0)
     &call errxit (msg(ms:14), 'road crown slope (CS) not found')
C                                                                     FRACTIONS
      do i = 1, 6
        fract(i) = 0.
      end do
C                                                   directly connected pervious
      call getr4 ('DCP', 0, fract(1), ierr)
C                                                 directly connected impervious
      call getr4 ('DCI', 0, fract(2), ierr)
C                                               indirectly connected impervious
      call getr4 ('ICI', 0, fract(3), ierr)
C                                                           connecting pervious
      call getr4 ('CP', 0, fract(4), ierr)
      if (fract(3) .gt. 0. .and. ierr .gt. 0)
     &  call errxit (msg(ms:14), 'indirect impervious not connected')

C                                                 indirectly connected pervious
      call getr4 ('ICP', 0, fract(5), ierr)
C                                                         connecting impervious
      call getr4 ('CI', 0, fract(6), ierr)
      if (fract(5) .gt. 0. .and. ierr .gt. 0)
     &  call errxit (msg(ms:14), 'indirect pervious not connected')

      sum = 0.
      do i = 1, 6
        sum = sum + fract(i)
      end do

      if (sum .lt. 0.99 .or. sum .gt. 1.01)
     &  call errxit (msg(ms:14), 'area fractions do not sum to one')

      nu = 0
C                                            get upstream element identifier(s)
      do k = 1, 10

        call geti4 ('UP', k, idup(k), ierr)
        if (ierr .eq. 3) call errxit
     &  (msg(ms:14), 'invalid upstream identifier (UP)')

        if (ierr .eq. 0) then
          nu = nu + 1
        else
          go to 611
        end if

      end do

611   continue
C                                               rainfall from upstream elements
C------------------------------------------------------------------------------

      do k = 1, nu
        upr(k) = 0
        call old (idup(k), 'RA', upr(k), ierr)
      end do
C                                                      update contributing area
C------------------------------------------------------------------------------
      coarea = 0.
      qbal(8) = 0.

      do k = 1, nu

        call geta (idup(k), area)
        coarea = coarea + area

      end do
C                                                                      rainfall
C------------------------------------------------------------------------------
      call interp (x, y, t, z, nd, sat, ngage)

C** temp expedient:
      z(nd+1) = z(nd)
      t(nd+1) = t(nd) + delt
C
      do ir = 2, nd + 1

         r(ir-1) = (z(ir) - z(ir-1)) / (t(ir) - t(ir-1))

      end do
C                                                 compute integrated rain rates
C------------------------------------------------------------------------------
      call new (id, 'RA', outr)

      itlim = limit
      qrmax = 0.
      k = 1
      area =  width* (xleng + street)
      qrf = r(1) * area
C!!  new code, 9/03, RES
      lmn = limit - 1
      zlr = z(1)
C!! replaced, 9/03:
C      do i = 1, limit
C        time = delt * float (i - 1)
C        if (time .gt. t(k+1)) then
C          k = k + 1
C          qrf = r(k) * area
C        end if
C!! new code 9/03
      do i = 1, lmn
        timen = delt*float(i)
        do while(t(k) .lt. timen)
          k = k+1
        end do
C!! moved 9/03:
          if(k .gt. nd) then   ! artificial extend rainfall
            t(k) = tfin*60. + 1.
            r(k) = 0.
C!!            z(k) = z(nd)
            nd   = k
          end if
C!! end moved
C!! continue new code:
C interpolate cumulative depth at time:
        km = k-1
        znr = z(km) + (z(k) - z(km))*(timen - t(km))/(t(k) - t(km))
        qrf = (znr - zlr)/DELT*area         ! mean over interval
        zlr = znr
C!! end new code
        qrain = qrf

        do j = 1, nu
          if (upr(j) .gt. 0) then
            call get (upr(j), i, qr)
            qrain = qrain + qr
          end if
        end do

        call store (outr, i, qrain)
        if (qrain .gt. qrmax) qrmax = qrain

      end do
C                                                         OVERLAND FLOW ROUTING
C------------------------------------------------------------------------------

      area = xleng * width
C                                                               outflow storage
      call new (-999, 'RO', out999)
C** 4/03  0. = zro
      call store (out999, 1, zro)

      first = .true.
      omega = 0.7
C** dx check !
      dxlim = dxlmt
      chl = 'm.'
      if(units .ge. 2) then
        dxlim = dxlmt/0.3
        chl = 'ft'
      end if
C                                   ... start of cycle thru 6 possible segments
      do n = 1, 6

        if (fract(n) .eq. 0.) cycle
     
        if (n .eq. 1) then
C                                                   directly connected pervious
          wid = area * fract(1) / xleng
          xl = xleng
          alpha = alph(2)
          infil = .true.
          upflow = .false.
          outq(1) = out999
          dint = dintr
          covr = cover
          if (first) then
            addq = .false.
            first = .false.
          else
            addq = .true.
          end if

        else if (n .eq. 2) then
C                                                 directly connected impervious
          wid = area * fract(2) / xleng
          xl = xleng
          alpha = alph(1)
          infil = .false.
          upflow = .false.
          dint = dintri
          covr = 1.0
          outq(1) = out999
          if (first) then
            addq = .false.
            first = .false.
          else
            addq = .true.
          end if

        else if (n .eq. 3) then
C                                               indirectly connected impervious

          wid = area * (fract(3) + fract(4)) / xleng
          xl = area * fract(3) / wid
          alpha = alph(1)
          infil = .false.
          dint = dintri
          covr = 1.0
          upflow = .false.
          call new (-888, 'RO', outq(1))
C** 4/03  0. = zro
          call store (outq(1), 1, zro)
          addq = .false.

        else if (n .eq. 4) then
C                                                           connecting pervious
          xl = xleng - xl
          alpha = alph(2)
          infil = .true.
          upflow = .true.
          dint = dintr
          covr = cover
          outq(1) = out999
          if (first) then
            addq = .false.
            first = .false.
          else
            addq = .true.
          end if

        else if (n .eq. 5) then
C                                                 indirectly connected pervious

          wid = area * (fract(5) + fract(6)) / xleng
          xl = area * fract(5) / wid
          alpha = alph(2)
          infil = .true.
          upflow = .false.
          dint = dintr
          covr = cover
          call new (-888, 'RO', outq(1))
C** 4/03  0. = zro
          call store (outq(1), 1, zro)
          addq = .false.

        else if (n .eq. 6) then
C                                                         connecting impervious
          xl = xleng - xl
          alpha = alph(1)
          infil = .false.
          upflow = .true.
          outq(1) = out999
          dint = dintri
          covr = 1.0
          if (first) then
            addq = .false.
            first = .false.
          else
            addq = .true.
          end if
        end if

        if (upflow) call old (-888, 'RO', up888, ierr)

        nk = max1 ((15. * xl / clen), 5.)
        if (nk .gt. 15) nk = 15
        dx = xl / (float(nk) - 1.)
C** 22/4/03
        if(n .le. 2 .and. dx .ge. dxlim .and. lwarn) then
          lwarn = .false.
          write(files(3), 1990) id, dx, chl, n!  , dxlim
 1990 format(/" Urban Element",i3,":  based on its length and parameter 
     &CLEN, the"/"  numerical increment of ",f6.1,1x,a2," for fraction "
     &,i1," is too large for realistic",/"  solution of the flow equatio
     &n, and may give  misleading results."/)
        end if
 !     
C                                                   initialize local variables
C------------------------------------------------------------------------------
        do j = 1, nk
          h1(j) = 0.
          h2(j) = 0.
          q1(j) = 0.
          q2(j) = 0.
          ql(j) = 0.
          fav(j) = r(1)
          afac(j) = 1.0
        end do

        time = 0.
        xc = 0.
        k = 1
        unif = .true.
        nc = 1
        depth = 0.
        dinact = 0.

        p1 = limit / 4
        p2 = limit / 2
        p3 = 3 * limit / 4

        qub1 = 0.
        vcm = 0.
        dfir = dint
        if(diag) write(99,'(" c,d",2f12.6)') covr, dint
C                                     read / initiallize infiltration variables
C------------------------------------------------------------------------------

        if(infil) call infil0 (id, 1, nk, msg(ms:14), diag, sat) ! units, 


C                                 start fixed time step loop (i, delt, step)...
C------------------------------------------------------------------------------

        do i = 2, limit
          itm = i
          qub2 = 0.
C                                                                        inflow
          if (upflow) call get (up888, i, qub2)

          dqub = qub2 - qub1
          qup = qub2
          step = float (i - 1) * delt
          dtp = delt

C                                start variable time step loop (k, dt, time)...
C------------------------------------------------------------------------------

20        if (time .gt. (t(k)- 1.e-6)) then
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

          depth = depth + rf * dtp

	    if(dinact .lt. dint) then
	      if(depth .gt. dint) then
	        rfu = (rf * dtp - (dint - dinact) * cover) / dtp
	        dinact = dint
	      else
	        dinact = depth
	        rfu = rf * (1. - cover)
	      end if
          else
            rfu = rf
          end if

          if (cour) then
C                                                     compute next dt using max
C                                                         celerity from last dt
            nc = nint (vcm * dtp / dx)
            if (nc .eq. 0) nc = 1
            dt = dtp / float (nc)

          else

            dt = dtp

          end if
C
          if (infil) then
C
            do j = 2, nk
C** 4/03  revise call parameter definition:
              q_sup(j) = q2(j-1)/dx
              h_avg(j) = 0.5*(h1(j) + h1(j-1))
            end do
C                                                      get average infiltration
C                                                                          rate
            call infilt ( afac, rfu, h_avg, q_sup, dtp, fav, 1)
C                                                           compute net lateral
C                                                                   inflow rate
            do j = 1, nk
              diff = rfu - fav(j)
              reldif = abs(diff)
              if(reldif .lt. 1.e-10) then          ! redefine array ql:
                ql(j) = 0.
              else
                ql(j) = diff
              end if
C
            end do
          else
C                                set lateral inflow rate equal to rainfall rate
            do j = 2,nk
              ql(j) = rfu
            end do
C
          end if

          dte = 0.
C                                                 start Courant loop (dt, m)...
C------------------------------------------------------------------------------

50        do m = 1, nc
C                                           time elapsed since beginning of dtp
            vcm = 0.
C                                                time elapsed since last "step"
            dte = dte + dt
C                                                            interpolate inflow
            deltr = dstep + dte
            qub = qub1 + dqub * deltr / delt

            if (qub .gt. 0.) then
C                                               compute upstream boundary depth
C------------------------------------------------------------------------------

              q2(1) = qub / wid
              unif = .false.
              h2(1) = (q2(1) / alpha)**(1. / p)
              vc = p * alpha * h2(1)**(p - 1.)

            else

              q2(1) = 0.
              h2(1) = 0.
              vc = 0.

            end if

            if (vc .gt. vcm) vcm = vc
C                                                               zone A solution
C------------------------------------------------------------------------------
            if (unif .and. ql(2) .gt. 0.) then

              h2(2) = h1(2) + ql(2) * dt
C                                                              compute celerity
              vc = p * alpha * h2(2)**(p - 1.)
              xc = xc + vc * dt
              if (vc .gt. vcm) vcm = vc

              if (xc .lt. dx) then
C                                                         solution is in zone A
                q2(2) = alpha * h2(2)**p

                do j = 3, nk
                  h2(j) = h2(2)
                  q2(j) = q2(2)
                end do
C                                                 bypass iterative calculations
                go to 60

              else
C                                                     solution is out of zone A
                unif = .false.

              end if
            end if
C                                                         start spatial loop...
C------------------------------------------------------------------------------

            do j = 1, nk - 1

              jp1 = j + 1
              htest = h2(j) + h1(jp1) + h1(j)
              qlav = ql(jp1)  ! ql(j) is now mean lat inflow between j and j-1
C                                                               test for runoff
              if (htest .eq. 0. .and. qlav .le. 0.) then

                h2(jp1) = 0.
                q2(jp1) = 0.
                vc = 0.

              else
C                                                iterative solution for h2(j+1)
                h2jp1 = (h2(j) + h1(jp1)) / 2.
                hmin = min(-tol,-h2jp1)
                source = 'UrbPLh'
                call iter (errpl, h2jp1, hmin, 10., ierr, source)

                if (ierr .gt. 0) call errxit
     &                           (msg(ms:14), 'no convergence (errpl)')

                if (h2jp1 .lt. 1.E-8) h2jp1 = 0.
                h2(jp1) = h2jp1
C                                                              compute celerity
                vc = p * alpha * h2jp1**(p - 1.)

              end if
C                                           vcm is max. celerity for current dt
              if (vc .gt. vcm) vcm = vc

            end do
C                                                       ... end of spatial loop
C------------------------------------------------------------------------------

C                                                       check Courant condition
C------------------------------------------------------------------------------

            if (vcm * dt / dx .gt. 1.) then
C                                                                     reduce dt
              dtprev = dt
              dtr = dtp - dte + dt
              nc = int (2. * vcm * dtr / dx)
              if (nc .eq. 0) nc = 1
              dt = dtr / float (nc)

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

60          continue

            h1(1) = h2(1)
            q1(1) = q2(1)

            do j = 2, nk

              h1(j) = h2(j)
              q1(j) = q2(j)
C                                                            infiltrated volume
              dvf = rfu - ql(j)   !) 0.5 * (ql(j-1) +
              if (h1(j) .le. 0.) dvf = amin1 (rfu, dvf)
              qbal(6) = qbal(6) + dt * dvf * wid * dx

            end do

          end do
C                                                        ...end of Courant loop
C------------------------------------------------------------------------------

          if (time .ne. step) go to 20
C                                             ...end of variable time step loop
C------------------------------------------------------------------------------

C                                                               store discharge
          qout2 = 0.
          if (addq) call get (outq(1), i, qout2)
          qout2 = qout2 + wid * q2(nk)
          call store (outq(1), i, qout2)
          qub1 = qub2

          if(diag) then
            write(99,'(" n,q: ",I2,3g12.4)') n, dint, rfu, q2(nk)
C            write(99,'(" iter2: ",i4, g12.4)') itot, h1r
            write(99,'("  ql",8g12.4)') (ql(j),j=1,nk)
            write(99,'("  h2",8g12.4)') (h2(j),j=1,nk)
          end if
C
        end do
C                                                ...end of fixed time step loop
C------------------------------------------------------------------------------

        ast = 0.5*(h2(1) + h2(nk))
        do j = 2, nk-1
          ast = ast + h2(j)
        end do
        qbal(9) = qbal(9) + wid*dx*ast  ! storage in m^3

        qbal(8) = qbal(8) + dinact * area * covr * fract(n)
C                                                ...end of loop thru 6 segments
      end do
C                                                          ROUTE THROUGH STREET
C------------------------------------------------------------------------------

      alpha = alph(3)
      co1 = 1. / zc
      co2 = 1. + sqrt (1. + 1. / zc**2)
C                                                                lateral inflow
      call old (-999, 'RO', latq, ierr)

      nl = 1
C                                                               outflow storage
      call new (id, 'RO', outq(1))

C                                                               upstream inflow
      do k = 1, nu
        call old (idup(k), 'RO', upq(k), ierr)
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), 'upstream inflow not found')
      end do
C                                                       area / contibuting area
      area = width * (xleng + street)
      qbal(1) = area
      coarea = coarea + area
      call stora (id, coarea)
C                                                          initialize variables
C------------------------------------------------------------------------------

      xl = width

      nk = nint (amax1 ((15. * xl / clen), 5.))
      if (nk .gt. 15) nk = 15
      dx = xl / (float(nk) - 1.)
C**  22/4/03 dx check
      if(dx .ge. dxlim) then
        write(files(3), 1991) id, dx, chl!  , dxlim
 1991 format(/"   Urban street ",i3,":  based on its length and paramete
     &r CLEN,"/"  the numerical increment of ",f6.1,1x,a2," is too large
     & for",/"  realistic numerical solution of the flow equation, and m
     &ay "/"  give  misleading results."/)
      end if
 !     
C** following loop moved, 4/03, since NK must be redefined first:
      do j = 1, nk
        h2(j) = 0.
        wpr2(j) = 0.
        wpr1(j) = 0.
        h1(j) = 0.
        a1(j) = 0.
        a2(j) = 0.
        q1(j) = 0.
        q2(j) = 0.
      end do
      itlim = limit
      omega = 0.8
      nc = 1
      xc = 0.
      unif = .true.
      time = 0.
      dt = delt
      qtest = 0.
      d1 = 0.
      vin = 0.
      qp = 0.
      ir = 1
      qlat1 = 0.
      qu1 = 0.
C                                       start fixed time step loop (i, delt)...
C------------------------------------------------------------------------------

      do i = 2, limit
C                                            time at start of current time step
        itm = i
        time = float (i - 2) * delt
C                                                         compute rainfall rate
C                                                          as average over delt
23      if (time + delt .gt. t(ir)) then
          ir = ir + 1
          go to 23
        end if

        d2 = z(ir-1) +  r(ir-1) * (time + delt - t(ir-1))
        rfi = (d2 - d1) / delt * street
        if (rfi .lt. 0.) rfi = 0.
        qtest = qtest + rfi
        d1 = d2

C** 4/03        vin = vin + qu2
        vin = vin + qu1
C                                                                 inflow volume
        qu2 = 0.
C                                                            upstream inflow(s)
        do k = 1, nu
          call get (upq(k), i, qup)
          qu2 = qu2 + qup
        end do
C                                                                lateral inflow
        call get (latq, i, qlat2)

C                                      prepare to linearly  interpolate inflows
C------------------------------------------------------------------------------

        dqu = (qu2 - qu1) / delt
        dqlat = (qlat2 - qlat1) / delt
        dte = 0.
        qtest = qtest + qu2 + qlat2

        if (qtest .le. 0.) then
C                                                                       no flow
C------------------------------------------------------------------------------
          do j = 1, nk
            h2(j) = 0.
            a2(j) = 0.
            q2(j) = 0.
            wpr2(j) = 0.
          end do

          dt = delt
          nc = 1
          qtest = 0.
C                                                                bypass routing
          go to 100

        end if
C                                                 start Courant loop (dt, m)...
C------------------------------------------------------------------------------

80      do m = 1, nc

          vcm = 0.
C                                               dte = time from last fixed step
C                                                     to end of current dt
          dte = dte + dt
C                                                           interpolate inflows
          qub = qu1 + dqu * dte
C                                                             lateral inflow is
C                                                               average over dt
          qlav = (qlat1 + dqlat * (dte - 0.5 * dt)) / xl

          if (qub .eq. 0.) then
C                                                   upstream boundary depth = 0
C------------------------------------------------------------------------------

            a2(1) = 0.
            wpr2(1) = 0.
            q2(1) = 0.
            h2(1) = 0.
C                                               compute upstream boundary depth
C------------------------------------------------------------------------------

          else

            source = 'UrbUBs'
            call iter (errsub, h2(1), 0., 100., ierr, source)

            if (ierr .gt. 0) call errxit
     &      (msg(ms:14), 'no convergence (errsub)')

            unif = .false.

            dadh = dafct (h2(1), co1)
            dqdh = dqfct (alpha, wpr2(1), p, a2(1), dadh, co2)

C                                                          estimate celerity at
C                                                                upper boundary
            if (dadh .ne. 0.) vc = dqdh / dadh
            if (vc .gt. vcm) vcm = vc

          end if

          if (unif) then
C                                                               zone A solution
C------------------------------------------------------------------------------
            h2a = h1(2)
            j = 2
            source = 'UrbZAS'
            call iter (errsza, h2a, 0., 10000., ierr, source)

            h2(2) = h2a
            dadh = dafct (h2a, co1)
            dqdh = dqfct (alpha, wpr2(2), p, a2(2), dadh, co2)
C                                                              compute celerity
            vc = dqdh / dadh
            xc = xc + dt * vc
            if (vc .gt. vcm) vcm = vc

            if (xc .ge. dx) then
C                                                       begin kinematic routing
              unif = .false.

            else

              q2(2) = qfunct (a2(2), wpr2(2), alpha, p)

              do j = 3, nk
                h2a = h1(j)
                source = 'UrbZBS'
                call iter (errsza, h2a, 0., 10000., ierr, source)
                h2(j) = h2a
                q2(j) = qfunct (a2(j), wpr2(j), alpha, p)
              end do
C                                                      bypass kinematic routing
              go to 90

            end if
          end if
C                                                start spatial loop (j, jp1)...
C------------------------------------------------------------------------------

          do j = 1, nk - 1

            jp1 = j + 1
            h2jp1 = h2(jp1)
            source = 'UrbH2S'
            call iter (errst, h2jp1, 0., 100., ierr, source)

            if (ierr .gt. 0) call errxit
     &      (msg(ms:14), 'no convergence (errst)')

            h2(jp1) = h2jp1

            dadh1 = dafct (h2jp1, co1)
            dqdh = dqfct (alpha, wpr2(jp1), p, a2(jp1), dadh1, co2)

C                                                             estimate celerity
C------------------------------------------------------------------------------

            if (dadh1 .ne. 0.) vc = dqdh / dadh1
            if (vc .gt. vcm) vcm = vc

          end do
C                                               ...end of spatial loop (j, jp1)
C------------------------------------------------------------------------------

C                                                       check Courant condition
C------------------------------------------------------------------------------

          if (vcm * dt / dx .gt. 1.) then
C                                                                compute new dt
            dtprev = dt
            dtr = delt - dte + dt
            nc = int (2. * vcm * dtr / dx)
            if (nc .eq. 0) nc = 1
            dt = dtr / float (nc)

            if (cour) then
C                                                 recompute profile with new dt
              dte = dte - dtprev
              go to 80

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
C                                                     ending time of current dt
90        timec = time + dte
C                                                                peak flow rate
          if (q2(nk) .gt. qp) then
            qp = q2(nk)
            tp = timec
          end if
C         if(diag) write(99,'(" ch ",3g12.4)') timec/60., qlav*xl, q2(nk)

          qbal(10) = qbal(10) + 0.5 * (q1(nk) + q2(nk)) * dt


C                                      reassign variables for next time substep
C------------------------------------------------------------------------------

          qtest = 0.

          do j = 1, nk

            h1(j) = h2(j)
            a1(j) = a2(j)
            wpr1(j) = wpr2(j)
            q1(j) = q2(j)
            qtest = qtest + q2(j)

          end do

        end do
C                                                ...end of Courant loop (dt, m)
C------------------------------------------------------------------------------

        if (cour) then
C                                        next dt based on previous max celerity
          nc = nint (vcm * delt / dx)
          if (nc .eq. 0) nc = 1
          dt = delt / float (nc)
          dte = 0.

        end if
C                                                                 store outflow
        call store (outq(1), i, q2(nk))

        qu1 = qu2
        qlat1 = qlat2
        qtest = qtest + qu2 + qlat2

100     continue
C
         if(nele .le. 5) call progss (msg(ms:18))
C
      end do
C                                                ...end of fixed time step loop
C------------------------------------------------------------------------------

C                                                               rainfall volume
      qbal(2) = area * z(nd)
C                                                                 inflow volume
      vin = (qu2 + 2. * vin) / 2. * delt
C                                                                storage volume
      sare = 0.
      do j = 2, nk - 1
        sare = sare + a2(j)
      end do

      qbal(9) = qbal(9) + dx * (a2(1) + 2. * sare + a2(nk)) / 2.
      if (qbal(9) .lt. 0.) qbal(9) = 0.
C                                                         pass values to writer
C------------------------------------------------------------------------------

      call qwrt (id, 1, msg(ms:18), qp, tp, vin,
     &           coarea, outq, 6, outr, qrmax)


      if (sed) then
C                                            store zero sediment concentrations
C------------------------------------------------------------------------------
        do n = 1, nps
          call new (id, attr(n), outc(n,1))
	    do i = 1, limit
C** 4/03  0. = zro
            call store (outc(n,1), i, zro)
	    end do
        end do
C**  this is a dummy call to put in zeros
        call swrt (wso, 0., 0., 0., 0., 0., outc)
      end if

      return

      end

C------------------------------------------------------------------------------

      subroutine errpl (h2jp1, ferr, dferr)

C   Computes the residual ferr and derivative dferr at h2jp1 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for flow on a sloping plane.

      integer j, jp1
C      character(LEN=6) :: trace
      common /vars/ a1(20), a2(20), q1(20), q2(20), wpr1(20), wpr2(20),
     &              co1, co2, h1(20), h2(20), j, jp1, omega, qub, dt,
     &              dx, alpha, p, rfi, qlav


      dhdt = .5 * (h2jp1 - h1(jp1) + h2(j) - h1(j)) / dt
C
      q2(jp1) = 0.
      if(h2jp1 .ge. 0.) q2(jp1) = alpha * h2jp1**p
      dqdx = (omega * (q2(jp1) - q2(j)) +
     &       (1. - omega) * (q1(jp1) - q1(j))) / dx
C
      ferr = dhdt + dqdx - qlav
      h2pw = 0.
      if(h2jp1 .gt. 0.) h2pw = h2jp1**(p-1.)
      dferr = .5 / dt + omega * alpha * p * h2pw / dx

      return

      end

C------------------------------------------------------------------------------

      subroutine errsub (hub, ferr, dferr)

C  Computes the residual and derivative w/resp to h at hub in a uniform
C  discharge equation for the street half.

      integer j, jp1

      common /vars/ a1(20), a2(20), q1(20), q2(20), wpr1(20), wpr2(20),
     &              co1, co2, h1(20), h2(20), j, jp1, omega, qub, dt,
     &              dx, alpha, p, rfi, qlav

      a2(1) = afunct (hub,  co1)
      wpr2(1) = wfunct (hub, co2)
      q2(1) = qfunct (a2(1), wpr2(1), alpha, p)
      dadh = dafct (hub, co1)
      dq2dh = dqfct (alpha, wpr2(1), p, a2(1), dadh, co2)
      ferr = q2(1) - qub
      dferr = dq2dh

      return

      end

C------------------------------------------------------------------------------

      subroutine errst (h2jp1, ferr, dferr)

C  Computes the residual and derivative w/resp to h at h2jp1 in a finite
C  difference kinematic flow equation for the street half.

      integer j, jp1

      common /vars/ a1(20), a2(20), q1(20), q2(20), wpr1(20), wpr2(20),
     &              co1, co2, h1(20), h2(20), j, jp1, omega, qub, dt,
     &              dx, alpha, p, rfi, qlav

      a2(jp1) = afunct (h2jp1, co1)
      wpr2(jp1) = wfunct (h2jp1, co2)
      q2(jp1) = qfunct (a2(jp1), wpr2(jp1), alpha, p)
      da2dh = dafct (h2jp1, co1)
      dq2dh = dqfct (alpha, wpr2(jp1), p, a2(jp1), da2dh, co2)
      dadt = 0.5 * (a2(jp1) - a1(jp1) + a2(j) - a1(j)) / dt
      dqdx = (omega * (q2(jp1) - q2(j)) +
     &       (1. - omega) * (q1(jp1) - q1(j))) / dx
      ferr = dadt + dqdx - qlav - rfi
      dferr = 0.5 * da2dh / dt + omega * dq2dh / dx

      return

      end

C------------------------------------------------------------------------------

      subroutine errsza (hza, ferr, dferr)

C  Computes the residual and derivative w/resp to h at hza in an eqn for flow
C  in zone A of the solution space (street half, no upstream inflow).

      integer j, jp1

      common /vars/ a1(20), a2(20), q1(20), q2(20), wpr1(20), wpr2(20),
     &              co1, co2, h1(20), h2(20), j, jp1, omega, qub, dt,
     &              dx, alpha, p, rfi, qlav

      a2(j) = afunct (hza, co1)
      wpr2(j) = wfunct (hza, co2)
      ferr = a2(j) - a1(j) - dt * (qlav + rfi)
      dferr = dafct (hza, co1)

      return

      end

C------------------------------------------------------------------------------

      function afunct (h, co1)

      afunct = 0.5 * h * h * co1

      return

      end

C ----------------------------------------------------------------

      function wfunct (h, co2)

      wfunct = h * co2

      return

      end

C --------------------------------------------------------------

      function qfunct (ar, wp, al, p)

      if (ar .gt. 0. .and. wp .gt. 0.) then
        qfunct = al * ar**p / wp**(p - 1.)
      else
        qfunct = 0.
      end if

      return

      end

C -------------------------------------------------------------

      function dafct (h, co1)

      dafct = h * co1

      return

      end
C ------------------------------------------------------------

      function dqfct (al, wp, p, ar, dadh, co2)

      if (wp .gt. 0.) then
        dqfct = al * (wp**(p - 1.) * p * dadh * ar**(p - 1.) -
     &            ar**p * (p - 1.) * co2 * wp**(p - 2.)) /
     &            wp**(2. * p - 2.)
      else
        dqfct = 0.
      end if

      return

      end
