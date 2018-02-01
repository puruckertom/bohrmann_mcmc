C   Code type: FORTRAN subroutine + 1 entry point

C   Compiler: Microsoft FORTRAN v5.1

C   Programmed by: C. Unkrich

C   Date: 4/24/95

C   Description:

C     Reads time-depth/intensity pairs for up to twenty raingages;
C     entry point interp returns a set of time-depth pairs based on
C     interpolation from raingage sites surrounding the element
C     centroid. The subroutine assumes the data is depth if the
C     column or list label is 'D...' and assumes it is intensity
C     if the label is 'I...'. If only one gage is present x, y
C     coordinates are not necessary. If no rainfall data file was
C     specified, then a zero input is assumed.

C     If initial soil saturation is specified, it will also be
C     interpolated.
C------------------------------------------------------------------------------

C   Arguments:

C   -- rain --

C     filen      integer       unit number assigned to the rainfall file,

C   -- interp --

C     x0         real          x coordinate of element centroid,

C     y0         real          y coordinate of element centroid,

C     t(1000)    real          interpolated time (sec),

C     z(1000)    real          interpolated depth (m or ft),

C     ndi        integer       number of t, z data pairs,

C     sat        real          initial soil saturation (-1 = no value)
C------------------------------------------------------------------------------

C  Subroutines/functions:

C     reader     (reader.for)
C     iget4      (reader.for)
C     rget4      (reader.for)
C     errxit     (K2RunW.for)

C------------------------------------------------------------------------------


      subroutine rain (filen)  !, units, limit, delt

      use runpars
      use messagew

      logical xy

      integer m, n, i, j, k, ierr, filen, kr(100), krtemp, kk, kkp1,
     &        nrg, icase, k1, k2, k3, i1, i2, i3, ii, nd, ndi, indx,
     &        np, flag, ngage, nt2  !, units, limit

      dimension t(1000), z(1000), ti(100000), zi(100000), x(100), 
     &      y(100), a(100), b(100), d(100), nd(100), indx(100), si(100)

      character rlab*1, gageid*10, msg*20

      save nrg, ti, zi, nd, x, y, a, b, d, kr, indx, si

      tfins = tfin*60.  ! end time in seconds
      if (filen .eq. -1) then
C                                               no rainfall data -- assume zero
        nrg = 1
        nd(1) = 2
        indx(1) = 1
        ti(1) = 0.
        ti(2) = tfins
        zi(1) = 0.
        zi(2) = 0.

        return

      end if

      k = 0
      n = 0
      xy = .true.
C                                             read data sets for all rain gages
C------------------------------------------------------------------------------

      do i = 1, 100

        call reader (filen, gageid, ierr)

        if (ierr .gt. 0) then

          msg = 'rain gage '//gageid

          if (ierr .eq. 1) then
            if(i .eq. 1) then
C                                                               end of file
          write(msgstring(1),"('       Rainfile specified as:     ')")
          write(msgstring(2),"(' ',a50,' is empty.    ')") fname(2)
          write(msgstring(3),"('   Make sure file names are correct befo
     &re proceeding')")
              ibutt = 0
              iconx = 1
              nlm = 3
!              call reportw
              if(errxt) 
     & call errxit(fname(2), ' eof encountered in rainfile')
            end if
            go to 8

          else if (ierr .eq. 2) then

            call errxit (msg, 'invalid data block')

          else

            call errxit (msg, 'block too big for data buffer')

          end if

        end if

        if (i .eq. 2 .and. .not. xy)
     &  call errxit (msg, 'x or y coordinate missing or invalid')

        msg = 'rain gage '//gageid
C                                                                       x coord
        call getr4 ('X', 0, x(i), ierr)

        if (ierr .gt. 0) then

          if (i .ge. 2) then

            call errxit (msg, 'x coordinate missing or invalid')

          else

            xy = .false.

          end if

        end if
C                                                                       y coord
        call getr4 ('Y', 0, y(i), ierr)

        if (ierr .gt. 0) then

          if (i .ge. 2) then

            call errxit (msg, 'y coordinate missing or invalid')

          else

            xy = .false.

          end if

        end if
C                                                       initial soil saturation
        call getr4 ('S', 0, si(i), ierr)

        if (ierr .gt. 0) si(i) = -1.
C                                                               # of data pairs
        call geti4 ('N', 0, np, ierr)

        if (ierr .gt. 0) call errxit
     &  (msg, 'number of data pairs not specified')

        t1 = -1.
        rlab = 'D'

        call getr4 (rlab, 1, zi(n+1), ierr)

        if (ierr .gt. 0) rlab = 'I'
C                                                        index to start of next
C                                                          data set in array ti
        indx(i) = k + 1

        do j = 1, np
C                                                                      get data
          call getr4 ('T', j, t2, ierr)

          if (ierr .gt. 0) call errxit
     &    (msg, 'time value missing or invalid')

          if (t2 .le. t1) call errxit
     &    (msg, 'time not increasing')

          k = n + j

          if (k .gt. 100000) stop ' error - too much rainfall data'
C protective roundoff to nearest 0.1 min:
          ot2 = t2
          nt2 = nint(ot2*10.)
          t2 = real(nt2,4)/10.
C          write(77,'(" t, round t:", 2f10.5)') ot2, t2
C                                                                   time in sec
          ti(k) = t2 * 60.
          t1 = t2

          call getr4 (rlab, j, zi(k), ierr)

          if (ierr .gt. 0) call errxit
     &    (msg, 'depth/intensity value missing or invalid')

        end do

        if (ti(n+np) .lt. tfins) then
C                                                               add one more pt
          np = np + 1
          k = n + np
          if (k .gt. 100000) stop ' error - too much rainfall data'
          ti(k) = tfins + 1.
          zi(k) = zi(k-1)

        end if

        nd(i) = np

        if (rlab .eq. 'I') then
C                                      convert intensity (mm or in/hr) to depth
          z1 = 0.

          do j = 2, np
            k = n + j
            z2 = z1 + zi(k-1) * (ti(k) - ti(k-1)) / 3600.
            zi(k-1) = z1
            z1 = z2
          end do

          zi(k) = z2

        else
C                                              check for decreasing depth value
          do j = 2, np
            k = n + j
            if (zi(k) .lt. zi(k-1)) call errxit
     &      (msg, 'decreasing depth value')
          end do

        end if

C            subtract first depth from subsequent depths and convert to m or ft

        do j = 1, np
          k = n + j
          zi(k) = (zi(k) - zi(n+1)) / conv
        end do

        n = k

      end do

8     nrg = i - 1

      return
C------------------------------------------------------------------------------


      entry interp (x0, y0, t, z, ndi, sat, ngage)

C                                                        interpolation override
      if (ngage .gt. 0) then
        if (ngage .gt. nrg) then
          call errxit ('INTERP',
     &                 'nonexistent gage specified for override')
        else
          k1 = max(ngage,1)
          go to 20
        end if
      end if

      k1 = 1

      if (nrg .eq. 1) go to 20


      do k = 1, nrg
C                                          rank gages by distance from centroid
        kr(k) = k
        a(k) = x(k) - x0
        b(k) = y(k) - y0
        d(k) = sqrt (a(k) * a(k) + b(k) * b(k))

      end do

      do k = 1, nrg - 1

        kk = kr(1)

        do ik = 1, nrg - 1

          kkp1 = kr(ik+1)

          if (d(kkp1) .lt. d(kk)) then

            krtemp = kr(ik)
            kr(ik) = kr(ik+1)
            kr(ik+1) = krtemp

          end if

          kk = kkp1

        end do

      end do


      if (nrg .eq. 2) then

C   2 gages -- determine whether the element centroid lies within the strip
C   bounded by two (parallel) lines that pass through the gage gage locations
C   and are perpendicular to the line connecting the two points:

        k1 = kr(1)
        k2 = kr(2)
        a12 = a(k2) - a(k1)
        b12 = b(k2) - b(k1)
        vmag = sqrt (a12 * a12 + b12 * b12)
        proj = (a12 * a(k2) + b12 * b(k2)) / vmag

        if (proj .gt. 0. .and. proj .lt. vmag) then

          c7 = abs (proj) / vmag
          c8 = 1. - c7
          icase = 2

          go to 10

        else
C                                                                  use one gage
          go to 20

        end if
      end if

C   Loop thru all combinations of gages taken three at a time, such that for
C   each succeeding triple, the sum of the distances from each gage to the
C   centroid increases. For each triple, determine whether the three gages
C   surround the element centroid by checking the direction of the cross
C   products between the vector connecting a pair of gages and the vector from
C   a member of the pair to the centroid, as well as the vector connecting the
C   same member to the third gage. Opposite directions indicate the centroid
C   and the third gage lie on opposite sides of the line through the gage
C   pair. This will work only if the gages are not colinear!

      do i = 1, nrg - 2

        k1 = kr(i)

        do j = i + 1, nrg - 1

          k2 = kr(j)

          do k = j + 1, nrg

            k3 = kr(k)

C   Check for colinearity:

            test = ((y(k1) - y(k2)) / (x(k1) * y(k2) - x(k2) * y(k1))) *
     &             x(k3) + ((x(k1) - x(k2)) / (y(k1) * x(k2) - y(k2) *
     &             x(k1))) * y(k3) + 1.

            if (abs (test) .gt. 1.E-4) then

              res12 = ((y(k2) - y(k1)) * x(k3) - (x(k2) - x(k1)) * y(k3)
     &                + y(k1) * (x(k2) - x(k1)) - x(k1) *
     &                (y(k2) - y(k1))) * ((y(k2) - y(k1)) * x0 - (x(k2)
     &                - x(k1)) * y0 + y(k1) * (x(k2) - x(k1)) - x(k1) *
     &                (y(k2) - y(k1)))

              res23 = ((y(k3) - y(k2)) * x(k1) - (x(k3) - x(k2)) * y(k1)
     &                + y(k2) * (x(k3) - x(k2)) - x(k2) *
     &                (y(k3) - y(k2))) * ((y(k3) - y(k2)) * x0 - (x(k3)
     &                - x(k2)) * y0 + y(k2) * (x(k3) - x(k2)) - x(k2) *
     &                (y(k3) - y(k2)))

              res13 = ((y(k3) - y(k1)) * x(k2) - (x(k3) - x(k1)) * y(k2)
     &                + y(k1) * (x(k3) - x(k1)) - x(k1) *
     &                (y(k3) - y(k1))) * ((y(k3) - y(k1)) * x0 - (x(k3)
     &                - x(k1)) * y0 + y(k1) * (x(k3) - x(k1)) - x(k1) *
     &                (y(k3) - y(k1)))


              if (res12 .ge. 0. .and. res23 .ge. 0. .and. res13 .ge. 0.)
     &        then
C                               surrounded - compute interpolating coefficients

                c1 = (x(k3) - x(k1)) * (y(k1) - y0) - (x(k1) - x0) *
     &               (y(k3) - y(k1))

                c2 = (x(k1) - x0) * (y(k2) - y(k1)) - (x(k2) - x(k1)) *
     &               (y(k1) - y0)

                c3 = (x(k2) - x(k1)) * (y(k3) - y(k1)) - (x(k3) - x(k1))
     &               * (y(k2) - y(k1))

                c4 = c1 / c3
                c5 = c2 / c3
                c6 = 1. - c4 - c5
                icase = 3

                go to 10

              end if
            end if
          end do
        end do
      end do

C   Loop thru all combinations of gages taken two at a time, such that for each
C   succeeding pair, the sum of the distances from each gage to the centroid
C   increases. For each pair, determine whether the element centroid lies
C   within the strip bounded by two (parallel) lines that pass through the gage
C   gage locations and are perpendicular to the line connecting the two points:

      do i = 1, nrg - 1

        k1 = kr(i)

        do j = i + 1, nrg

          k2 = kr(j)

C**CU 6/2007 I don't know why res12 was calculated for this case, k3 is not even updated

C          res12 = ((y(k2) - y(k1)) * x(k3) - (x(k2) - x(k1)) * y(k3) +
C     &            y(k1) * (x(k2) - x(k1)) - x(k1) * (y(k2) - y(k1))) *
C     &            ((y(k2) - y(k1)) * x0 - (x(k2) - x(k1)) * y0 + y(k1)
C     &            * (x(k2) - x(k1)) - x(k1) * (y(k2) - y(k1)))

C          if (res12 .le. 0.) then
C                                             project vector k1->k0 onto k1->k2
            a12 = a(k2) - a(k1)
            b12 = b(k2) - b(k1)
            vmag = sqrt (a12 * a12 + b12 * b12)
            proj = (a12 * a(k2) + b12 * b(k2)) / vmag

            if (proj .gt. 0. .and. proj .lt. vmag) then

              c7 = abs (proj) / vmag
              c8 = 1. - c7
              icase = 2

              go to 10

            end if
C          end if
        end do
      end do
C                           no suitable gage pair found - use closest gage only
      k1 = kr(1)

      go to 20
C                                             boolean union of breakpoint times
C------------------------------------------------------------------------------

10    i1 = indx(k1)
      ndi = nd(k1)

      if (ndi .gt. 1000) go to 40
C                                                    put in times from 1st gage
      do k = 1, ndi

        t(k) = ti(i1+k-1)

      end do
C                                                    insert times from 2nd gage
      kk = k2
      i2 = indx(k2)
      ii = i2
      flag = 0

1109  j = 1

      do k = 0, nd(kk) - 1

        tt = ti(ii+k)

1111    if (tt .lt. t(j)) then

          if (ndi + 1 .gt. 1000) go to 40

          do m = ndi, j, -1

            t(m+1) = t(m)

          end do

          t(j) = tt
          ndi = ndi + 1

        else if (tt .eq. t(j)) then

          j = j + 1

        else

          j = j + 1

          if (j .gt. ndi) then

            if (ndi + nd(kk) - k .gt. 1000) go to 40

            do m = k, nd(kk) - 1

              ndi = ndi + 1
              t(ndi) = ti(ii+m)

            end do

            go to 1112

          end if

          go to 1111

        end if
      end do

1112  if (icase .eq. 3 .and. flag .eq. 0) then
C                                                     insert time from 3rd gage
        kk = k3
        i3 = indx(k3)
        ii = i3
        flag = 1

        go to 1109

      end if

C                                       interpolate depth for each gage at each
C                                      time, then interpolate depth at centroid
      t(1) = 0.
      z(1) = 0.

      do i = 2, ndi

        if (ti(i1+1) .lt. t(i)) i1 = i1 + 1

        d1 = zi(i1) + ((t(i) - ti(i1)) / (ti(i1+1) - ti(i1))) *
     &       (zi(i1+1) - zi(i1))

        if (ti(i2+1) .lt. t(i)) i2 = i2 + 1

        d2 = zi(i2) + ((t(i) - ti(i2)) / (ti(i2+1) - ti(i2))) *
     &       (zi(i2+1) - zi(i2))

        if (icase .eq. 2) then

          z(i) = c7 * d1 + c8 * d2
          if (z(i) .lt. 0.) z(i) = 0.

        else

          if (ti(i3+1) .lt. t(i)) i3 = i3 + 1

          d3 = zi(i3) + ((t(i) - ti(i3)) / (ti(i3+1) - ti(i3))) *
     &         (zi(i3+1) - zi(i3))

          z(i) = c6 * d1 + c4 * d2 + c5 * d3
          if(z(i) .lt. 0.) z(i) = 0.

        end if

      end do
C                                                           interpolate initial
C                                                                 soil moisture
      if (si(k1) .lt. 0. .or. si(k2) .lt. 0.) then

        sat = -1.

      else if (icase .eq. 2) then

        sat = c7 * si(k1) + c8 * si(k2)

        if (sat .lt. 0.) sat = 0.

      else if (si(k3) .lt. 0.) then

        sat = -1.

      else

        sat = c6 * si(k1) + c4 * si(k2) + c5 * si(k3)

        if (sat .lt. 0.) sat = 0.

      end if

      go to 30

20    ndi = nd(k1)

      i1 = indx(k1) - 1

      if (ndi .gt. 1000) go to 40

      do j = 1, ndi

        i = i1 + j

        t(j) = ti(i)
        z(j) = zi(i)

      end do

      sat = si(k1)

30    return

40    stop ' error (interp) -- too many breakpoints'

      end
