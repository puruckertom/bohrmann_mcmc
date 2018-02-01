! Code type: format FORTRAN subroutine
! Compiler: Fortran90
! Programmed by: C. Unkrich
! Date: 4/95
! Revised: 7/00
!
! Description:
!
!   Four-point implicit finite-difference method for approximating the
!   solution of the kinematic wave equation for flow in a conduit of circular
!   cross section.  The program will be terminated with a warning if inflow
!   dictates full pipe flow.  An algorithm has been added which monitors the
!   wave celerity and adjusts the computational time increment to maintain
!   numerical accuracy.
!
!   For output, volumetric rainfall rate from the upstream area for each time
!   step is added and stored.  Dividing by contributing area gives the
!   integrated rainfall rate over the area represented by the its upstream
!   contributors.
!-------------------------------------------------------------------------------
!
! Arguments:
!   dtm(4)        real        smallest dt to satisfy the Courant condition
!                (out)        during each quarter of the simulation,
!-------------------------------------------------------------------------------
!
! Input Block (input file labels in parenthesis):
!
!   id (ID)          int      identifier (1 to 9999),
!
!   xleng (L)        real     length (m or ft),
!
!   diam (D)         real     diameter (m or ft),
!
!   slope (SL)       real     slope,
!
!   rs (MA or CH)    real     Manning or Chezy discharge coefficient,
!
!   upstream(UP)     int      identifiers of input elements (up to 10)
!   input         
!-------------------------------------------------------------------------------
!  other variables defined in module definition file KINMODS
!
! Subroutines/functions:
!
!   iter          (iter.for)
!   sed0          (kinsed.for)
!   kinsed             "
!   sedfin             "
!   errth2        (pipe.for)
!   erthub             "
!   progss        (progss.for)
!   getr4         (reader.for)
!   getstr             "
!   new           (clerk.for)
!   old                "
!   get                "
!   store              "
!   stora              "
!   geta               "
!   qwrt          (writer.for)
!   errxit        (errxit.for)
!-------------------------------------------------------------------------------

      subroutine pipe( dtm )
C
      use runpars
      use elpars
      use multip
!                                                              declare arguments
!-------------------------------------------------------------------------------
      dimension  dtm(4)
!                                                        declare local variables
!-------------------------------------------------------------------------------
      integer upr, upq, outq(2), outr
      character msg*18
      character(LEN=6) :: source = 'pipe  '  
C** 22/4/03  units character
      character(LEN=2) :: chl
      dimension ql(20,2),  a2(20), q2(20), idup(10),!, q1(20) a1(20),
     &          upr(10), upq(10,2), inflow(12), slp(20), topw(20),
     &          qt(20), qup(10,2)

      external erthub, errth2

      data pi, twopi /3.1415927, 6.2831853/
!                                                                common data
!-------------------------------------------------------------------------------
      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &               j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &               qr1(20), qr2(20), thz, az, daz, dqz
C
      omega = 0.8
!                                                           get input parameters
!-------------------------------------------------------------------------------

! IDENTIFIER

      call geti4('ID', 0, id, ierr)
      if(ierr .gt. 0) call errxit('PIPE', 'missing or invalid identifier
     & (ID)')
      ltyp = 2
      rsoil = .false.   ! no infiltrating sol

		! construct identifier string as 'PIPE element  XXXX'

      msg(1:16) = typname(2)
      write (msg(15:18), '(i4)') id

		! don't include leading blanks

      do ms = 1, 5
        if(msg(ms:ms) .ne. ' ') exit
      end do

		! update computational progress displayed on standard output

      if(nele .gt. 5) call progss(msg(ms:18))

! LENGTH

      call getr4('L', 0, xleng, ierr)
      if(ierr .gt. 0) call errxit(msg(ms:18), 'length (L) not found')

		! compute number of spatial nodes: 5 <= nk <= 15

      nk = max1((15.* xleng/clen), 5.)
      if(nk .gt. 15) nk = 15

		! compute distance increment

      dx = xleng/(float(nk) - 1.)

C** 22/04/03  dx check
      dxlim = dxlmt
      chl = 'm.'
      if(units .ge. 2) then
       dxlim = dxlmt/0.3
       chl = 'ft'
      end if
      if(dx .ge. dxlim) then
        write(files(3), 1990) id, dx, chl!  , dxlim
 1990 format(/"  Pipe ",i3,":  based on its length and parameter CLEN,"/
     &"  the numerical increment of ",f6.1,1x,a2," is too large for",/
     &"  realistic numerical solution of the flow equation, and may "/
     &"  give  misleading results."/)
      end if
 !     
! DIAMETER

      call getr4('D', 0, diam, ierr)
      if(ierr .gt. 0) call errxit(msg(ms:18), 'diameter (D) not found')

! SLOPE

      call getr4('SL', 0, slope, ierr)
      if(ierr .gt. 0) call errxit(msg(ms:18), 'slope (SL) not found')

! MANNING?

      call getr4('MA', 0, rs, ierr)

      if(ierr .eq. 0) then

        power = 1.66667
        fac = 1.
        if(units .eq. 2) fac = 1.49
        alpha = fac * sqrt (slope)/(rs * rmn)

      else

! CHEZY?

        call getr4('CH', 0, rs, ierr)

        if(ierr .eq. 0) then
	  power = 1.5
          alpha = rs * rmn * sqrt(slope)
	else
	  call errxit(msg(ms:18),'hydraulic roughness (MA or CH) not fou
     &nd')
        end if

      end if

! UP - upstream element identifier(s)

      nu = 0

      do k = 1, 10

        call geti4('UP', k, idup(k), ierr)
        if(ierr .eq. 3) call errxit(msg(ms:18),'invalid upstream identif
     &ier (UP)')

        if(ierr .eq. 0) then
          nu = nu + 1
        else
          exit
        end if

      end do
!                                      get indices to storage locations
!-------------------------------------------------------------------------------

!     integrated rainfall

      call new(id, 'RA', outr)

!     loop over upstream contributing elements

      do k = 1, nu

        call old(idup(k), 'RO', upq(k,1), ierr)
        if(ierr .gt. 0) call errxit(msg(ms:18),
     &                              'upstream inflow not found')
        inflow(k) = idup(k)

!       overbank section upstream?

        call old(idup(k), 'ro', upq(k,2), ierr)
        if(ierr .eq. 0) inflow(k) = -idup(k)

!       rainfall

        upr(k) = 0
        call old(idup(k), 'RA', upr(k), ierr)

      end do

!     outflow

      call new(id, 'RO', outq(1))

!                                              initialize variables
!-------------------------------------------------------------------------------
      call store (outq(1), 1, 0.)

      do j = 1, nk
        th1(j) = 0.
        th2(j) = 0.
        ar1(j) = 0.
        a2(j) = 0.
        qr2(j) = 0.
        ar2(j) = 0.
        q2(j) = 0.
      end do
      q1nk = 0.
      itotl = 0

      do j = 1, 10
        qbal(j) = 0.
      end do

! compute limiting inflow rate (when pipe is at full capacity)

      qmax = 0.99 * alpha * pi * diam * (diam/4.)**power

      area = 0.
      qp = 0.
      tp = 0.
      vi = 0.
C      vo = 0.
C      xc = 0.
      vcm = 0.
      nc = 1
      qub1 = 0.
      qubprev = 0.

      p1 = limit / 4
      p2 = limit / 2
      p3 = 3 * limit / 4

      diamsq = diam * diam
C set up min angle parameters:
      thz = 0.
      thzu = 0.05
      az = area_f(thzu)
      daz = dadth(thzu)
      dqz = dqdth_f(thzu)
      thz = thzu
      wz = thz*diam/2.
C
      if(diag) write(99,'(" al, p : ",2g12.4)') alpha, power

!                                     read/initialize sediment variables
!-------------------------------------------------------------------------------
      if(sed) then

        do j = 1, nk
          slp(j)  = slope
	    ql(j,1) = 0.
	    ql(j,2) = 0.
	    qt(j)   = 0.
        end do

        call sed0(1, id, msg(ms:18), inflow, nu, 0, nk,
     &            omega, dx, slp, 1.)

      end if
!                                               update contributing area
!-------------------------------------------------------------------------------
      coarea = 0.

      do k = 1, nu
        call geta(idup(k), area)
        coarea = coarea + area
      end do

      call stora(id, coarea)
!                            integrate rain rates over contributing areas
!-------------------------------------------------------------------------------
      qrmax = 0.

      do i = 1, limit

        qrain = 0.

        do j = 1, nu

          if(upr(j) .gt. 0) then
            call get(upr(j), i, qr)
            qrain = qrain + qr
          end if

        end do

        call store(outr, i, qrain)

        if(qrain .gt. qrmax) qrmax = qrain

      end do
!                                start fixed time step loop (i, delt)...
!-------------------------------------------------------------------------------
      do i = 2, limit

        time = float (i - 2) * delt	!time at start of current time step

      ! obtain upstream inflow(s)

        qub2 = 0.

        do k = 1, nu

          call get(upq(k,1), i, qup(k,1))
          qub2 = qub2 + qup(k,1)

          if(inflow(k) .lt. 0) then

  ! add inflow from overbank section

            call get(upq(k,2), i, qup(k,2))
            qub2 = qub2 + qup(k,2)
          end if

        end do
			! prepare to linearly interpolate inflow

        dqub = (qub2 - qub1)/delt

        if(cour) then

! divide current time step based on previous maximum celerity

          nc = nint(vcm * delt/dx)
          if(nc .eq. 0) nc = 1
          dt = delt/float(nc)
        else
          dt = delt
        end if

        dte = 0.
!                                          start Courant loop (dt, m)...
!-------------------------------------------------------------------------------
50      do m = 1, nc

          vcm = 0.
          dte = dte + dt ! time elapsed from beginning of fixed time step

				! interpolate inflow
          qub = qub1 + dqub * dte/delt

! compute upper boundary flow angle th2(1)

          if(qub .le. 1.e-15) then

  ! no inflow

            th2(1) = thz
            ar2(1) = 0.
            qr2(1) = 0.
            topw(1) = wz

           else

          ! check for full pipe flow

            if(qub .gt. qmax)
     &  call errxit(msg(ms:18), 'free surface flow capacity exceeded')

          ! iteratively solve for th2(1)

            if(qub .gt. 10. * qubprev) then
	      th2(1) = 5. * th1(1)
	    else
	      th2(1) = th1(1)
	    endif

            if(th2(1) .eq. 0.) th2(1) = pi/10.
            source(5:6) = 'UB'
            call iter(erthub, th2(1), .0, twopi, ierr, source)          ! ITER

            if (ierr .gt. 0) then
              call errxit(msg(ms:18),' no convergence (erthub)')
            else if(ierr .lt. 0) then
              itub = -ierr
            end if
!          if(th2(1) .lt. .001) th2(1) = 0.
            ar2(1) = area_f(th2(1))
            qr2(1) = qub
            topw(1) = topw_f(th2(1))

          end if

	  qubprev = qub
!                                                   start spatial loop...
!-------------------------------------------------------------------------------
          do j = 1, nk - 1

            jp1 = j + 1

            if(th2(j) + th1(j) + th1(jp1) .le. 3.*thz) then

					  ! no flow
              th2(jp1) = thz
              ar2(jp1) = 0.
              qr2(jp1) = 0.
              topw(jp1) = wz
              vc = 0

            else

! iterative solution for th2(j+1)

              th2jp1 = amax1(0.5*(th2(j)+th1(jp1)), thz)
              source(5:6) = 'IT'
              call iter(errth2, th2jp1, .0, twopi, ierr, source)

              if (ierr .gt. 0) then
                write(99,'(2i3,4g12.4)') j, ierr, th2(j),th1(jp1),th2jp1
                call errxit(msg(ms:18),'no convergence (errth2)')
              else
               itotl = itotl - ierr
              end if
!            if(th2jp1 .lt. .011) th2jp1 = 0.
C              angl1a = (th2jp1 + th1(jp1))/2.
              th2(jp1) = th2jp1
              ar2(jp1) = area_f(th2jp1)
              qr2(jp1) = q_f(th2jp1)
              topw(jp1) = topw_f(th2jp1)

						! compute celerity
              vc = dqdth_f(th2jp1)

            end if

            if(vc .gt. vcm) vcm = vc ! vcm is maximum celerity for current dt

          end do
!                                                   ... end of spatial loop
!-------------------------------------------------------------------------------

        ! check Courant condition

          if(vcm * dt/dx .gt. 1.) then
! reduce dt
            dtprev = dt
            dtr = delt - dte + dt
            nc = int(2. * vcm * dtr/dx)
            if(nc .eq. 0) nc = 1
            dt = dtr/float(nc)

            if(cour) then

C  recompute with new subinterval
              dte = dte - dtprev
              go to 50

            else
				 ! find minimum dt over each quarter
              if(i .le. p1) then
                if(dt .lt. dtm(1)) dtm(1) = dt
              else if(i .le. p2) then
                if(dt .lt. dtm(2)) dtm(2) = dt
              else if(i .le. p3) then
                if(dt .lt. dtm(3)) dtm(3) = dt
              else
                if(dt .lt. dtm(4)) dtm(4) = dt
              end if

						! don't recompute dt
              nc = 1
              dt = delt

            end if

          end if

          timec = time + dte ! ending time of current dt

      ! check for peak flow

          if(q2(nk) .gt. qp) then
            qp = q2(nk)
            tp = timec
          end if

		! increment outflow volume (trapezoidal integration)

          qbal(10) = qbal(10) + dt * 0.5 * (q1nk + q2(nk))

        ! reset variables for next time step

          q1nk = q2(nk)
          do j = 1, nk
            qr1(j) = qr2(j)
            q2(j) = qr2(j)     ! for kinsed
            ar1(j) = ar2(j)
            a2(j) = ar2(j)     ! for kinsed
            th1(j) = th2(j)
          end do
          if(diag) then
            write(99,801)'pipe',i,itotl,qub
            write(99,800)' a2 ',(ar2(j),j=1,nk)
            write(99,800)' q2 ',(qr2(j),j=1,nk)
            write(99,800)'th2 ',(th2(j),j=1,nk)
          end if
 800  format(a4,8g12.4)
  801 format(a4,2i5,5g12.4)
        ! compute sediment concentration

          if(sed) call kinsed(i,dt,dte,ql,a2,q2,topw,0.,qup,qt,timec)

        end do
!                                                         ...end of Courant loop
!-------------------------------------------------------------------------------

      ! store current outflow

        call store(outq(1), i, q2(nk))

			! increment inflow volume

        vi = vi + qub1

      ! reset variables for next fixed time step

        qub1 = qub2
C
      if(nele .le. 5) call progss(msg(ms:18))

      end do
!                                           ...end of fixed time step loop
!-------------------------------------------------------------------------------

    ! finish sediment computations

      if(sed) call sedfin()

    ! finish inflow volume calculation

      vi = (delt/2.) * (2. * vi + qub2)

! compute volume of water remaining in the pipe (Simpson's Rule)

      vsum = 0.

      do j = 2, nk - 1
        vsum = vsum + th2(j) - sin(th2(j))
      end do

      qbal(9) = (dx/2.) * (diamsq/8.) *
     &     (th2(1) - sin(th2(1)) + 2. * vsum + th2(nk) - sin(th2(nk)))

   ! pass info to writer

      call qwrt(id, 1, msg(ms:18), qp, tp, vi, coarea,
     &  outq, 2, outr, qrmax)

      end

!-------------------------------------------------------------------------------

      subroutine erthub(thub, ferr, dferr)

! Computes the residual and derivative at angle thub of the error function
! form of the discharge equation for circular pipe flow with a free surface
! condition:

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      ferr = q_f(thub) - qub
      dferr = dqdth_f(thub)

      return
      end

!-------------------------------------------------------------------------------

      subroutine errth2(th2jp1, ferr, dferr)

! Computes the residual and derivative at angle th2jp1 of the error function
! form of the kinematic equation for circular pipe flow with a free surface
! condition

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

C
      if (th2jp1 .le. thz) then
        a2thp = 0.
        q2thp = 0.
        da2jp1 = daz
        dq2jp1 = dqz
      else
        a2thp = area_f(th2jp1)
        q2thp = q_f(th2jp1)
        da2jp1 = dadth(th2jp1)
        dq2jp1 = dqdth_f(th2jp1)
      end if
C
      ferr = (ar2(j) - ar1(j) + a2thp - ar1(jp1))/2./dt +
     &    (omega*(q2thp - qr2(j)) + (1.-omega)*(qr1(jp1) - qr1(j)))/dx
C
      dferr = da2jp1/2./dt + dq2jp1*omega/dx
C
C      write(99,103) jp1, th2(j), th1(j), th1(jp1), th2jp1, q2thp, a2thp
C      write(99,100) angl2a, angl2q, cos2a, area2q, dadx, fac2, term1 !,
C      write(99,101) term2, term3, term4, dadth2, deriv, dqdth2, d2qdtx
  100 format(" 101", 8g12.5)
  101 format(" 102", 8g12.5)
  103 format(" ths", i2,8g12.5)
       return
      end

!-------------------------------------------------------------------------------

      function area_f(theta) ! Cross sectional area

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      if(theta .gt. thz) then
        if(theta .lt. 0.01) then
          area_f = diamsq/8. * thz**3 /6. - az
        else
          area_f = diamsq * (theta - sin (theta))/8. - az
        end if
      else
        area_f = 0.
      end if

      return
      end

!----------------------------------------------------------------------------

      function dadth(theta)
C
      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      if(theta .gt. thz) then
        if(theta .gt. .01) then
          dadth = diamsq * (1 - cos (theta))/8.
        else
          dadth = diamsq * theta**2/16. 
        end if
      else
        dadth = daz
      end if

      return
      end

!----------------------------------------------------------------------------

      function wpm_f(theta)    ! Wetted perimeter

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      if(theta .ge. thz) then
        wpm_f = diam * theta/2.
      else
        wpm_f = diam * thz/2.
      end if

      return
      end

!-------------------------------------------------------------------------------

      function topw_f(theta) ! Top width

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      topw_f = amax1(diam * sin (theta/2.), 0.01)

      return
      end

!-------------------------------------------------------------------------------

      function q_f(theta)    ! Discharge

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      afu = area_f(theta)
      if(afu .gt. 0.) then
        q_f = alpha * afu**power/wpm_f(theta)**(power - 1.)
      else
        q_f = 0.
      end if

      return
      end

!-------------------------------------------------------------------------------

      function dqdth_f(theta) ! Derivative of discharge w/resp to angle theta

      common /pipe1/ th1(20), th2(20), qub, diamsq, diam, alpha, power,
     &              j, jp1, dt, dx, omega, ar1(20), ar2(20),
     &              qr1(20), qr2(20), thz, az, daz, dqz

      area = area_f(theta)

      if(area .gt. 0.) then
C      if(theta .gt. 0.05)

        fac = (2.* area/(theta * diam))**(power - 1.)
        term1 = power * dadth(theta)     !diamsq * (1.- cos(theta))/8.
        term2 = area * (power - 1.)/theta

        dqdth_f = alpha * fac * (term1 - term2)

      else             ! surogate function for small theta
        dqdth_f = 0.
C        v = diamsq*(diam/4.)**(power-1.)*alpha/8.
C        dqdth_f = .05*v*4.34*theta**(3.34)
      end if

      return
      end

