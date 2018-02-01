C   Code type: FORTRAN subroutine

C   Compiler: Microsoft FORTRAN v5.1

C   Programmed by: C. Unkrich

C   Date: 2/96

C   Description:

C     Reads a file with time and inflow data pairs, interpolates values at
C     the the user-specified time steps, and stores them. If sediment routing
C     has been requested, it will also attempt to read nps columns of sediment
C     concentration, and if that fails, will assume the incoming water carries
C     no sediment, i.e., zero concentration for all particle classes.
C------------------------------------------------------------------------------

C   Arguments:

C     qbal(15)       real       array with element volume balance components:

C                               qbal(4)  = injection       (cu.m/cu.ft),
C                               qbal(10) = outflow         (cu.m/cu.ft),

C                               in this case qbal(4) = qbal(10);

C     delt                      user-specified time step (sec),

C     limit                     number of time steps,

C     sed                       .true. = route sediment,

C     nps                       number of particle classes,

C     rho(5)                    particle densities,

C     nord(5)                   new ordering sequence of particle classes
C                               according to erodability (sed00),
C------------------------------------------------------------------------------

C  Input Block (input file labels in parenthesis):

C     id (ID)          int*4       identifier (up to 6 digits),

C     filespec (FI)    char*50     input file specification,

C     offstim (OFF)     real        time offset (min) to be added to values
C                                  read from the input file - must be positive,
C------------------------------------------------------------------------------

C  Subroutines/functions:

C     getr4         (reader.for)
C     getstr             "
C     geti4              "
C     new           (miscel.for\clerk)
C     store              "
C     stora              "
C     qwrt          (writer.for)
C     swrt          (writer.for)
C     errxit        (K2RunW.for)
C     progss        (K2RunW.for)
C------------------------------------------------------------------------------


      subroutine inject (nord) ! , delt, limit, sed, nps, rho

      use runpars
      use elpars

      logical sedmt, next   !sed, 

      integer j, ierr, n, i, k, outq(2), outc, ms, nord !, nps, limit

      character filespec*150, attr*2, test*150  !, msg*15
      character(LEN=18) :: msg

      dimension  wso(5), sd(5), sd0(5), outc(5,2),
     &          c1(5), c2(5), attr(5), dtm(4), nord(5) !, rho(5)

C      data msg /'       (INJECT)'/
      data attr /'S1', 'S2', 'S3', 'S4', 'S5'/

      zero = 0.
      qp = 0.
      tp = 0.
      qsp = 0.
      tsp = 0.
      nchan = 1
      ltyp = 5
      rsoil = .false.

      do k = 1, 4

        dtm(k) = delt

      end do

      sedmt = .FALSE.

      do k = 1, 15

        qbal(k) = 0.

      end do
C                                                                         input
C------------------------------------------------------------------------------

      call geti4 ('ID', 0, id, ierr)
C                                                                    identifier
      if (ierr .gt. 0) call errxit ('INJECT', 'id not found')

C      write (msg(1:6), '(i6.5)') id
      msg(1:16) = typname(5)
      write (msg(15:18), '(i4)') id

      do ms = 1, 5
        if (msg(ms:ms) .ne. ' ') go to 1
      end do

1     call progss (msg(ms:15))

      call getstr ('FI', 0, filespec, j, ierr)
C                                                               input file
      if (ierr .gt. 0) call errxit
     &(msg(ms:15), 'input file not specified')
C                                                                     open
      open (90, err = 100, file = filespec(1:j), status = 'old')
C!!                     proposed add
      call getr4 ('ARE', 0, areain, ierr)
      if(ierr .gt. 0.) areain = 0.
C!! end proposed add

      call getr4 ('OFF', 0, offstim, ierr)
C                                                          time offset (min)
      if (ierr .gt. 0) then

        offstim = 0.

      end if

      if (offstim .lt. 0.) call errxit (msg(ms:15), 'negative offset')

      if (sed) then
C                                         look for sediment conc. (columns 3-7)
C------------------------------------------------------------------------------

        read (90, 5, err = 110) test
5       format (a)
        num = 0
        next = .true.

        do n = 1, 150

          if (ichar (test(n:n)) .gt. 44) then

            if (next) num = num + 1
            next = .false.

          else

            next = .true.

          end if

        end do

        if (num .EQ. nps + 2) sedmt = .TRUE.

        rewind (90)

      end if
C                                                             read first record
C------------------------------------------------------------------------------
      if (sed) then

        if (sedmt) then

          read (90, *, err = 110) t1, q1, (c1(nord(n)), n = 1, nps)

          do n = 1, nps
            sd0(n) = q1 * c1(n) * rho(n)
          end do

	  else

          read (90, *, err = 110) t1, q1

          do n = 1, nps
	      c1(n)  = 0.
            sd0(n) = 0.
          end do

	  end if

      else

        read (90, *, err = 110) t1, q1

      end if

      q0 = q1

      if (t1 .ne. 0.) call errxit
     &(msg(ms:15), 'time must start at zero')

C                                                                 second record
C------------------------------------------------------------------------------
      if (offstim .gt. 0.) then

        if (q1 .ne. 0.) call errxit
     &  (msg(ms:15), 'cannot offset nonzero initial discharge')

        t2 = offstim
        q2 = 0.

        if (sed) then

          do n = 1, nps
            c1(n) = 0.
            c2(n) = 0.
            sd0(n) = 0.
          end do

        end if

      else

	  if (sed) then

          if (sedmt) then

            read (90, *, err = 110) t2, q2, (c2(nord(n)), n = 1, nps)

          else

            read (90, *, err = 110) t2, q2

            do n = 1, nps
              c2(n) = 0.
            end do

	    end if

	  else

          read (90, *, err = 110) t2, q2

        end if

      end if
C                                                               reserve storage
C------------------------------------------------------------------------------

      call new (id, 'RO', outq(1))

      call store (outq(1), 1, q1)

      vsum = 0.

      if (sed) then

        do n = 1, nps

          call new (id, attr(n), outc(n,1))

          call store (outc(n,1), 1, c1(n))

          wso(n) = 0.
          sd(n) = 0.

        end do

      end if

      dt = delt / 60.
      qout = q2         ! at time = 0.
C                                      interpolate to user-specified time steps
C------------------------------------------------------------------------------

      do i = 2, limit
C!! correct summation  RES
        if(i .lt. limit) vsum = vsum + qout

        if (sed) then

          do n = 1, nps

            wso(n) = wso(n) + sd(n)

          end do

        end if

        time = dt * float(i - 1)

40      if (time .gt. t2) then
C
          t1 = t2
          q1 = q2

          if (sed) then
          
            if (sedmt) then

              do n = 1, nps
                c1(n) = c2(n)
              end do

              read(90, *, end=60, err=110) t2, q2, (c2(nord(n)),n=1,nps)
	      
            else

              read (90, *, end = 60, err = 110) t2, q2

              do n = 1, nps
                c2(n) = 0.
              end do

            end if

          else
C                                 read input file: t(min), q(cu.m/s or cu.ft/s)
            read (90, *, end = 60, err = 110) t2, q2
          
          end if

          t2 = t2 + offstim

          go to 40

        end if
C                                                             compute discharge
        qout = q2 - (t2 - time) * (q2 - q1) / (t2 - t1)
C                                                                         store
        call store (outq(1), i, qout)

        if (qout .gt. qp) then
C                                                                peak discharge
          qp = qout
          tp = time

        end if
C                                                                      sediment
C------------------------------------------------------------------------------

        if (sed) then

          tsd = 0.

          do n = 1, nps

            cs = c2(n) - (t2 - time) * (c2(n) - c1(n)) / (t2 - t1)

            call store (outc(n,1), i, cs)
C                                                            sediment discharge
            sd(n) = cs * qout * rho(n)
C                                                      total sediment discharge
            tsd = tsd + sd(n)

          end do

          if (tsd .gt. qsp) then
C                                                                peak discharge
            qsp = tsd
            tsp = time

          end if

        end if

      end do

60    close (90)
C                                       zero discharge for remaining time steps
C------------------------------------------------------------------------------

      do k = i, limit

        call store (outq(1), k, zero)

        if (sed) then

          do  n = 1, nps

            call store (outc(n,1), k, zero)

          end do

        end if

      end do
C------------------------------------------------------------------------------

C                                                                outflow volume
      qbal(10) = delt * (q0 + 2. * vsum + qout) / 2.

      qbal(4) = qbal(10)

      tp = tp * 60.
C                                                         pass values to writer
C------------------------------------------------------------------------------

      call qwrt (id, 1, msg(ms:18), qp, tp, 0., areain, outq, 5, 0, 0.)

C                                                                     zero area
      call stora (id, 0.)

      if (sed) then

        do n = 1, nps

          wso(n) = delt * (sd0(n) + 2. * wso(n) +
     &                     sd(n) * qout * rho(n)) / 2.

        end do

        tsp = tsp * 60.
C                                                         pass values to writer
C------------------------------------------------------------------------------

        call swrt (wso, 0., 0., 0., qsp, tsp, outc)

      end if

      return

100   call errxit (msg(ms:15), 'input file not found')

110   call errxit (msg(ms:15), 'error while reading input file')

      end
