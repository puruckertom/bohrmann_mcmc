C   Code type: Fortran subroutine + 2 entry points

C   Compiler: Fortran77

C   Programmed by: R.E. Smith

C   Date: 11/93

C   Last modification: 4/97 by R.E. Smith

C   Description:

C     A four-point explicit finite-difference solution of the equation of
C     supply and conservation of a transported material. The solution is
C     explicit since advantage is taken of the kinematic wave solution for the
C     same time step. This routine is called from the PLANE, CHANNEL, and PIPE
C     modules. This revised version treats a sediment supply drawn from as
C     many as 5 different particle size classes, thus to deal with a watershed
C     where soils are not uniform for all elements. Revised again to handle
C     parallel processing for compound channel sections.
C------------------------------------------------------------------------------

C   Entry points:

C     sed0          reads parameters & initializes vars for current element,

C     sed00         reads & initializes global sediment parameters & vars,

C     sedbf         compute initial conc. in a baseflow profile based on
C                   transport capacity,

C     sedfin        finish sediment balance computations for current element,
C------------------------------------------------------------------------------

C   Arguments:

C     k              int          = 1 for planes,
C                                 = 2 for channels,

C     nsect          int          = 1 for main channel section & plane,
C                                 = 2 for overbank section,

C     id             int          element identifier,
C     msg            char*(*)     element id string for error messages,
C     indx          int          current time index,

C     inflow(12)     int          inflow(1:10)  = upstream element id's, where
C                                                 a neg value indicates a
C                                                 compound channel upstream,
C                                 inflow(11:12) = lateral element id's,

C     units          int          1 = metric, 2 = english,
C     qup(10,2)      real         upstream inflow rate(s) at current time,
C     dt             real         current time increment (sec),
C     dte            real         time elapsed since last user-defined time step,
C     ql(20,k)       real         lateral inflow of water ...
C                                 ... for a plane, the average rate over the time
C                                     interval, m/s or ft/s;
C                                 ... for a channel, the rate at the end of the
C                                     interval, m*m/s or ft*ft/s,
C     qil(20,k)      real         bed or soil inflow of water (cu.m/s or cu.ft/s),
C     a2(20,k)       real         area at end of interval (sq.m or sq.ft),
C     q2(20,k)       real         discharge at end of time interval,
C     topw(20,k)     real         width of flow at water surface (m or ft),
C     wid(20,k)      real         width of plane or channel,
C     rf             real         current rainfall rate (m/s or ft/s),
C     qt(20)         real         for a compound channel, the net discharge
C                                 transferred from main -> overbank section;
C                                 (-) indicates reverse flow direction,

C     time           real         current simulation time,
C     typ            int          index for type of element:
C                                 = 0 for "plane", with splash as the source of
C                                     lateral contribution,
C                                 = 1 for a channel, with no splash and
C                                     potentially some lateral inflow,
C                                 = 2 for a pipe with no lateral inflow,
C                                 = 4 for a compound channel,

C     nx             int          number of spatial nodes,
C     omeg           real         time weighting factor from finite-diff. eqn.
C                                        of calling routine,
C     delx           real         spatial increment (m or ft),
C     slp(20,k)      real         bottom slope at each spatial node.
C     units          int          1 = metric, 2 = english,
C     del            real         fixed time step (sec),
C     lim            int          number of (fixed) time steps),
C     nps            int          number of particle classes
C     dpd(5)         real         particle diameters (m or ft),
C     rho(5)         real         particle densities,
C     norder(5)      int          new ordering sequence of particle classes
C                                            according to erodability,
C     diag           log          .true. = write diagnostic info to unit 99,
C------------------------------------------------------------------------------

C   Input Block (input file labels in parenthesis):

C     temp               real     water temperature (C or F),
C     dpd(5)             real     particle diameters (mm or in),
C     rho(5)             real     particle density (g/cc),
C     pav(k) (PA)        real     erosion pavement fraction (1 = noneroding),
C     bspl (SPL)         real     rain splash coefficient (plane only),
C     bcoh(k) (CO)       real     cohesion coefficient,
C     pros(5,k) (FR)     real     fraction of total material represented by
C                                     each particle class,
C     fmin (KS)          real     Ks, mm/fr or in/hr (plane only),
C------------------------------------------------------------------------------

C   Subroutines/functions:

C     splash       (kinsed.for)
C     trnscap           "
C     vsetl             "
C     switch            "
C     getr4        (reader.for)
C     swrt         (writer.for)
C     errxit       (K2RunW.for)
C------------------------------------------------------------------------------


      subroutine kinsed ( indx, dt, dte, ql, a2, q2,
     &                    topw, rf, qup, qt, time)

      use runpars
      use elpars
      use multip
C                                                                     arguments
C------------------------------------------------------------------------------

      integer  indx, typ, nx, nsect

C      dimension inflow(12), ql(20,k), a2(20,k), q2(20,k), topw(20,k),
C     &          qup(10,2), wid(20,k), slp(20,k), qt(20), norder(5)! dp(5), rh(5),vs(5)
      real, dimension(5) :: vs
      real, dimension(20) :: qt
C      real, dimension(20,nchan) ::  a2, q2, topw, wid, slp
      real, dimension(20,nchan) ::  a2, q2, topw, slp
      real, dimension(20,2) :: ql
      real, dimension(10,2) :: qup
      integer, dimension(12) :: inflow

      character msg*(*)

C                                                               local variables
C------------------------------------------------------------------------------

      logical,save :: rpt, over

      integer :: ierr, i, im, j, jp, m, n, l, na, nup, nlat    !typl,
      integer,save :: nk, ilast, nu, nl, nm
C
      real :: zro = 0.
      real, save :: omega, comega, dx, qsp, tsp, xl, qfac, 
     &              bspl, fmin, bdep, grav, swup
      real,save,dimension(2) :: pav, bcoh
      real,dimension(2) :: qu, qlu
      real,save,dimension(5) :: vsl, pwt, wso, dwso, dwsol, wsin, dwsup, 
     &                     dwsupl
      real,save,dimension(5,2) :: csub, pros, sd, csub1, csub2, dcsub,
     &                       qcsl1, qcsl2, dqcsl, qcsl, qcslm
      real,save,dimension(20,2) :: depnet, slope, qilm, a1, q1, qil,
     &                        qll !, width
      real,save,dimension(5,20,2) :: cs, cl, cml, deprt, supl, cm
C
      integer, save, dimension(2) :: nlp
      integer, dimension(5) :: norder, nord
      integer, save, dimension(5,2) :: outc, latc
      integer, save, dimension(10,5,2) :: upc
      save :: nord
C
      character(LEN=2),dimension(5,2) :: attr = reshape((/'S1', 'S2', 
     &'S3', 'S4', 'S5', 's1', 's2', 's3', 's4', 's5'/),(/5,2/))

      integer, parameter :: plane = 0, chan = 1, pipe = 2, cmpd = 4,
     &        urban = 6

C      data qilm /40*0./, cml /200*0./

C------------------------------------------------------------------------------

C      ku = k
      if (indx .ne. ilast) then
C                                             store outgoing concentration from
C                                                   end of last fixed time step
        do n = 1, nps

          call store (outc(n,1), ilast, cs(n,nk,1))

          if (ltyp .eq. cmpd) call store (outc(n,2), ilast, cs(n,nk,2))

        end do

      end if

      dt2 = dt / 2.

      do j = 1, 2

        do n = 1, nps
          sd(n,j) = 0.
          csub(n,j) = 0.
        end do

      end do

      if (ltyp .eq. plane) then
C                                             lateral inflow is average over dt
        do i = 2, nk
          qil(i,1) = ql(i,1)
        end do
C temp fix
C        ql(1,1) = ql(2,1)

      end if

      if (nu .gt. 0) then
C                                       compute concentration at upper boundary
C------------------------------------------------------------------------------

        if (ltyp .eq. plane) then

          if (indx .ne. ilast) then
C                                                          next fixed time step
            do n = 1, nps
              csub1(n,1) = csub2(n,1)
            end do

            if (qup(1,1) .gt. 0.) then

              do n = 1, nps
                call get (upc(1,n,1), indx, csub2(n,1))
              end do

            else

              do n = 1, nps
                csub2(n,1) = 0.
              end do

            end if

            do n = 1, nps
              dcsub(n,1) = (csub2(n,1) - csub1(n,1)) / delt
            end do

          end if
C                                                                   interpolate
          do n = 1, nps
            csub(n,1) = csub1(n,1) + dcsub(n,1) * dte
          end do

        else
C                                                                channel / pipe
          if (indx .ne. ilast) then
C                                                          next fixed time step
            qu(1) = 0.
            qu(2) = 0.

            do n = 1, nps
              csub1(n,1) = csub2(n,1)
              if (ltyp .eq. cmpd) csub1(n,2) = csub2(n,1)
            end do

            do j = 1, nu

              if (qup(j,1) .gt. 0.) then

                qu(1) = qu(1) + qup(j,1)

                do n = 1, nps
                  call get (upc(j,n,1), indx, cj)
C          if(diag) write(99,'(" cj: ",2I3,g13.4)') upc(j,n,1), indx, cj
                  sd(n,1) = sd(n,1) + qup(j,1) * cj
                end do

              end if

              if (qup(j,2) .gt. 0.) then

                qu(2) = qu(2) + qup(j,2)
C                                                  inflow from overbank section
                do n = 1, nps
                  call get (upc(j,n,2), indx, cj)
                  sd(n,2) = sd(n,2) + qup(j,2) * cj
                end do

              end if

            end do

            if (ltyp .eq. cmpd) then
C                                                              compound channel
              if (qu(2) .lt. 1.E-8) then
C                                                         all upstream sediment
                if (qu(1) .gt. 0.) then
C                                                               -> main channel
                  do n = 1, nps
                    csub2(n,1) = (sd(n,1) + sd(n,2)) / qu(1)
                    csub2(n,2) = 0.
                  end do

                end if

              else
C                                             distribute sediment at confluence
                if (qt(1) .ge. 0.) then
C                                           net sediment flow: main -> overbank
                  m = 1
                  j = 2

                else
C                                                              overbank -> main
                  m = 2
                  j = 1

                end if

                absqt = abs (qt(1))
C                                                compute sediment concentration
                do n = 1, nps

                  if (qu(m) .gt. 0.) csub2(n,m) = sd(n,m) / qu(m)

                  if (qu(j) .gt. 0.) then

                    csub2(n,j) = (sd(n,j) + csub2(n,m) * absqt) /
     &                           (qu(j) + absqt)

                  else

                    csub2(n,j) = csub2(n,m)

                  end if

                end do

              end if

              do n = 1, nps
                dcsub(n,1) = (csub2(n,1) - csub1(n,1)) / delt
                dcsub(n,2) = (csub2(n,2) - csub1(n,2)) / delt
                csub(n,2) = csub1(n,2)
              end do

            else
C                                                                simple channel
              if (qu(1) .gt. 0.) then

                do n = 1, nps
                  csub2(n,1) = (sd(n,1) + sd(n,2)) / qu(1)
                end do

              end if

              do n = 1, nps
                dcsub(n,1) = (csub2(n,1) - csub1(n,1)) / delt
              end do

            end if

          end if
C                                                                   interpolate
          do n = 1, nps
            csub(n,1) = csub1(n,1) + dcsub(n,1) * dte
          end do

          if (ltyp .eq. cmpd) then

            do n = 1, nps
              csub(n,2) = csub1(n,2) + dcsub(n,2) * dte
            end do

          end if

        end if

      end if
C                                                                lateral inflow
C------------------------------------------------------------------------------

      if ((ltyp .eq. chan .or. ltyp .eq. cmpd) .and. nl .ge. 1) then
C
        qlu(1) = 0.
        qlu(2) = 0.
        do m = 1, nl
          if(m .gt. nm) then      ! create combinations for various nm, nl options
            qlu(2) = ql(2,m)
          else
            qlu(1) = qlu(1) + ql(2,m)
          end if
        end do

        j = 1
        nmx = nm
        if (over) then
          j = 2     ! this is for nl = 1 only
          nmx = 1
        end if
        m = 1
C   
        do while(m .le. nl)

          if (indx .ne. ilast) then
C                                                        next fixed time step
            do n = 1, nps
              qcsl1(n,j) = qcsl2(n,j)
              qcsl2(n,j) = 0.
            end do
C
            do while(m .le. nmx )
              do n = 1, nps
                call get (latc(n,m), indx, cj)           !? err found 2/02  1 -> m
                qcsl2(n,j) = qcsl2(n,j) + cj * ql(2,m)   ! ql(i,m) are the same for any i
              end do 
C            if(diag) write(99,'(" cj: ",i3, 2g13.4)') m, cj, ql(1,m)
              m = m + 1
            end do
            do n = 1, nps
              dqcsl(n,j) = (qcsl2(n,j) - qcsl1(n,j)) / delt
C              qcsl1(n,j) = qcsl2(n,j)                ! this seems wrong too
            end do

          else
            m = nl+1   !  to end loop
          end if
C
          do n = 1, nps
C                                           average sediment discharge
            qcsl(n,j) = qcsl1(n,j) + dqcsl(n,j) * (dte - 0.5 * dt)

          end do
C                                         average lateral inflow rate
          do i = 1, nk
            qil(i,j) = (qlu(j) + qll(i,j)) / 2.
            qll(i,j) = qlu(j)
          end do

          if(j .lt. 2) j = j + 1
C          nmx = nmx + 1

        end do

      end if

      ilast = indx
C                            BEGIN SEDIMENT ROUTING -- upper boundary condition
C------------------------------------------------------------------------------

      if (a2(1,1) .lt. 1.E-8) then
C                                                            no upstream inflow
        if (ltyp .eq. plane) then

          hb = 0.
          qcsu = splash (hb, rf, bspl, bdep)
          qiltot = rf * topw(1,1)

        else if (ltyp .ne. pipe) then

          qiltot = 0.

          do j = 1, nl
            qiltot = qiltot + qil(1,j)
          end do

        end if

        do n = 1, nps

          if (ltyp .eq. plane) then
C                                          rain splash determines concentration
            qcs1 = qcsu * pros(n,1)

          else if (ltyp .ne. pipe) then
C                                           use lateral inflow at upstream node
            qcs1 = 0.

            do j = 1, nl
              qcs1 = qcs1 + qcsl(n,j)
            end do

          end if

          denom = qiltot + vsl(n) * topw(1,1)

          if (denom .gt. 0.) then

            cs(n,1,1) = qcs1 / denom

          else

            cs(n,1,1) = 0.

          end if

          if (ltyp .eq. cmpd) cs(n,1,2) = 0.

        end do

      else                     ! positive upper bound area:
C                                        upstream flow determines concentration
        do j = 1, nchan

          do n = 1, nps

            cs(n,1,j) = csub(n,j)
            cm(n,1,j) = 0.

            if (a2(1,j) .gt. 1.E-8) then
C                                                    compute transport capacity
              v = q2(1,j) / a2(1,j)
              d = a2(1,j) / topw(1,j)
              cm(n,1,j) = trnscap (d, v, slope(1,j),
     &                             dpd(n), rho(n), grav)

            end if

          end do

        end do

      end if
C                                                 loop through spatial nodes...
C------------------------------------------------------------------------------

      do i = 2, nk

        if (qt(i) .ge. 0. .or. nchan .le. 1) then
C                                                    process main section first
          l = 1

        else
C                                                        overbank section first
          l = 2

        end if

        im = i - 1
C
        rux = amin1 (qil(i,1), 0.0)

C                                         loop for possible overbank section...
C------------------------------------------------------------------------------

        do j = 1, nchan

          sumdel = 0.
          qcsb = -.9
          qcsi = -.9       !  defaults for non-plane diagnostic print only

          if(ltyp .eq. plane) then
C                                                           compute rain splash
C                                                                     component
            hb = 0.5 * (a2(i,1) + a2(im,1))
            qcsi = splash (hb, rf, bspl, bdep)
            qcsb = splash( 0., rf, bspl, bdep)                 ! new 2/02
C                                                          erodibility index of
C                                                        largest particle class
          end if

          dcpx = cml(nlp(l),i,l) - cl(nlp(l),i,l)

          if (a2(i,l) .gt. 1.E-8) then
C                                                  velocity and hydraulic depth
            v = q2(i,l) / a2(i,l)
            d = a2(i,l) / topw(i,l)

            if (i .eq. 4) dxpr = dcpx

          end if
C                                                   loop thru particle sizes...
C------------------------------------------------------------------------------

          do n = 1, nps
C            if(pros(n,l) .lt. 0.001) cycle          ! added res 10/98

            if (a2(i,l) .gt. 1.E-8) then
C                                                            transport capacity
              cm(n,i,l) = trnscap (d, v, slope(i,l),
     &                             dpd(n), rho(n), grav)
              cp = cl(n,i,l)
              cpm = cl(n,im,l)
              pro = 0.
              upav = 1.

              if (dcpx .gt. 0. .and. n .ge. nlp(l)) then

C    largest class is erodible, so all are (potentially): determine erodability
C   factor based on availability (default erodability factor based on cohesion)

                pro = pros(n,l) * bcoh(l)

                if (deprt(n,i,l) .le. 1.E-6) then
C                                                         no deposited material
                  upav = 1. - pav(l)

                  if (ltyp .eq. pipe .or. pav(l) .ge. 0.999) upav = 0.

                else if (pav(l) .ge. .999 .and.

     &                   depnet(i,l) .gt. 1.E-5) then
C                                                          erosion of deposited
C                                                              material in pipe
                  prt = deprt(n,i,l) / depnet(i,l)
                  pro = amin1 (prt, 1.)
                  upav = 1.

                end if

              else
C                                                      at least some deposition
                if (cml(n,i,l) .lt. cp) then

                  pro = 1.
C                                                          deposition indicated
                else

                  supl(n,i,l) = 0.
C                                                  erosion during previous step
                end if

              end if
C                                           condition at beginning of time step
C------------------------------------------------------------------------------

              if (ltyp .eq. plane) then
C                                                      rain splash contribution
                spav = upav
                qcs = qcsi * pros(n,1)
C                                                 balance init. conditions
                if (a1(i,1) .lt. 1.E-8 .and. qil(i,1) .gt. 1.E-6)
     &            cp = qcs / (qil(i,1) + vsl(n))

                if (a1(im,1) .lt. 1.E-8 .and. qil(im,1) .gt. 1.E-6)
     &            cpm = qcs / (qil(im,1) + vsl(n))

                clow = qcsb /(amax1(qil(i,1),0.) + vsl(n))  !modifed 2/02
C   trial
                qcs = qcs * topw(i,1) +
     &                (1. - topw(i,1)) * clow * (qil(i,1) - rux)

              else if (ltyp .eq. chan) then
C                                                        channel lateral supply
                qcs = qcsl(n,1) + qcsl(n,2)

                if(a1(i,1) .lt. 1.E-8) cp = (qcslm(n,1) + qcslm(n,2))/
     &                                          
     &                  (qilm(i,1) + qilm(i,2) + vsl(n) * topw(i,1))

                spav = 1.
 
                clow = qcs / vsl(n) / topw(i,l)

              else if (ltyp .eq. cmpd) then
C                                                        channel lateral supply
                if (l .eq. 2) then
C                                                              overbank section
                  m = 2
C                  if (over) m = 1
                  if (a1(i,2) .lt. 1.E-8) cp = qcslm(n,m) / (qilm(i,m)
     &                                         + vsl(n) * topw(i,m))

                  if (q2(i,2) .gt. 0.) then

                    qcs = qcsl(n,m)

                  else
C                                                no flow in overbank section --
C                                                lateral inflow -> main section
                    qcs = 0.

                  end if

                else if (.not. over) then
C                                                                  main section
                  qcs = qcsl(n,1)

                  if (q2(i,2) .eq. 0.)
C                                                            all lateral inflow
C                                                               to main section
     &               qcs = qcs + qcsl(n,2)

                  if(a1(i,1) .lt. 1.E-8) cp = (qcslm(n,1) + qcslm(n,2))/
     &                                          
     &                    (qilm(i,1) + qilm(i,2) + vsl(n) * topw(i,1))

                else if (q2(i,2) .eq. 0.) then
C                                                no flow in overbank section --
C                                                 single inflow -> main section
                  qcs = qcsl(n,1)

                else

                  qcs = 0.
C                                               flow into overbank section only
                end if

                spav = 1.

                clow = qcs / vsl(n) / topw(i,l)

              else if (ltyp .eq. pipe) then
C                                                      pipe -- no lateral input
                qcs = 0.
                clow = 0.
                spav = upav

              end if
C                                            compute erosion or deposition term
C------------------------------------------------------------------------------

              cgu = vsl(n) * pro
              dfac = upav * cgu * topw(i,l)
              eterm = dfac * cm(n,i,l)

              rpt = .false.

15            dtrc = eterm + qcs * spav

              absqt = abs (qt(i))

              if (j .eq. 2) then
C                                                 compound channel -- determine
C                                           sediment transfer  between sections
                if (l .eq. 1) then

                  ct = cs(n,i,2) + cl(n,i,2)
C
                else

                  ct = cs(n,i,1) + cl(n,i,1)

                end if
C                                                        sediment is carried in
C                                                         from adjacent section
                qst = 0.5 * absqt * ct
                qto = 0.

              else
C                                                       sediment is carried out
                qto = 0.5 * absqt
                qst = -0.5 * absqt * cl(n,i,l)

              end if
C                                                         compute concentration
C------------------------------------------------------------------------------

              supc = dtrc + comega * supl(n,i,l)
              af = cp * a1(i,l) / dt
              bf = (omega * cs(n,im,l) * q2(im,l) - comega *
     &             (cp * q1(i,l) - cpm * q1(im,l))) / dx
              df = 1. / dt + omega * v / dx
              deno = a2(i,l) * df + dfac - omega * rux * topw(i,1) + qto
              cs(n,i,l) = (supc + af + bf + qst) / deno

              if (cs(n,i,l) .lt. 0.) cs(n,i,l) = 0.

C                               compute erosion or deposition for current class
C------------------------------------------------------------------------------

              eros = eterm - cs(n,i,l) * dfac

              if (pav(l) .lt. 1. .and. n .gt. nlp(l)
     &            .and. .not. rpt .and. erolim .gt. 0.) then
C                                                             check for mixture
C                                                                 limit on eros
                erop = erolim * pros(n,l) / pros(nlp(l),l)

                if (eros .gt. erop .and. eros .gt. 0.) then
C                                limit erosion to proportion of least erodable class
                  eterm = erop
C                  dfac = 0.
                  rpt = .true.

                  go to 15

                end if

              end if
C                                                           net supply rate for
C                                                          this class in sq.m/s
              sup = rux * topw(i,l) * cs(n,i,l)
              spu = eros + qcs * spav

              if (cs(n,i,l) .lt. 0.) cs(n,i,l) = clow

            else
C                                                       no flow at current node
              cm(n,i,l) = 0.
              cs(n,i,l) = 0.
              qcs = 0.
              sup = 0.
              qst = 0.
              spu = 0.
              eros = 0.
              spav = 1.
              qto = 0.

            end if
C                                         erosion / deposition depth accounting
C------------------------------------------------------------------------------

            if (n .eq. nlp(l)) erolim = eros

            if (ltyp .eq. plane) then
C                                                    interval supply in tons/m2
              dela = (comega * supl(n,i,l) +
     &                omega * sup + spu) * dt * rho(n)

            else

              dela = (comega * supl(n,i,l) + omega * sup
     &                + spu - qcs * spav) * dt * rho(n)

              if (qto .eq. 0.) then
C                                                     extra sediment carried in
                dela = dela - qst
              else
                dela = dela + 0.5 * absqt * (cl(n,i,l) + cs(n,i,l))
              end if

            end if
C                                                   erosion / deposition depth
            if (ltyp .eq. pipe .and. dela .gt. deprt(n,i,l))
     &        dela = 0.999 * deprt(n,i,l)

            deprt(n,i,l) = deprt(n,i,l) - dela
            sumdel = sumdel + dela
            supl(n,i,l) = sup

          end do
C                                             ...end loop thru particle classes
C------------------------------------------------------------------------------

          depnet(i,l) = depnet(i,l) - sumdel
          l = 3 - l

        end do
C                                                     ...end main/overbank loop
C------------------------------------------------------------------------------
      end do
C                                                           ...end spatial loop
C------------------------------------------------------------------------------

C                                            reset variables for next time step
C------------------------------------------------------------------------------
      do j = 1, nchan

        do i = 1, nk

          a1(i,j) = a2(i,j)
          q1(i,j) = q2(i,j)

          do n = 1, nps
            cl(n,i,j) = cs(n,i,j)
            cml(n,i,j) = cm(n,i,j)
          end do

        end do

      end do
       
      if(diag) then
        write(99,313) nk, qcsi, qcsb, cm(nps,nk,1) !, qcsl(nps,2)
        write(99,311)' a1m:',(a1(i,1),i=1,nk)
        write(99,311)' cm: ',(cl(nps,i,1),i=1,nk)
        if(ltyp .eq. cmpd) then
          write(99,311)' a1o:',(a1(i,2),i=1,nk)
          write(99,311)' cov:',(cl(nps,i,2),i=1,nk)
        end if
C        write(99,'(2g13.4)') depnet(5,1), depnet(5,2)
 311  format(a5/(8g11.3))
 313  format(' qcsl: ',i3,5g13.4)
      end if

      do j = 1, nl

        do i = 1, nk
          qilm(i,j) = qil(i,j)
        end do

        do n = 1, nps
          qcslm(n,j) = qcsl(n,j)
        end do

      end do
C                    compute contribution to total outflow at current time step
C------------------------------------------------------------------------------

      twso = 0.

      do n = 1, nps

        dwso(n) = qfac * q2(nk,1) * cs(n,nk,1)

        if (ltyp .eq. cmpd) dwso(n) = dwso(n) + q2(nk,2) * cs(n,nk,2)

        wso(n) = wso(n) + dt2 * (dwso(n) + dwsol(n))
        dwsol(n) = dwso(n)
        twso = twso + dwso(n) * rho(n)

      end do
C                                             check for peak sediment discharge
      if (twso .gt. qsp) then

        qsp = twso
        tsp = time

      end if
C                                                 compute contribution to total
C                                                   inflow at current time step
      do n = 1, nps

        dwsup(n) = 0.

        do j = 1, nchan

          dwsup(n) = dwsup(n) + qfac * q2(1,j) * csub(n,j)

        end do

        wsin(n) = wsin(n) + dt2 * (dwsup(n) + dwsupl(n))
        swup = swup + dwsupl(n)*dt
        dwsupl(n) = dwsup(n)
C                                                             lateral inflow is
C                                                      already averaged over dt
        do j = 1, nl

          wsin(n) = wsin(n) + dt * xl * qcsl(n,j)

        end do

      end do
      if(diag .and. ltyp .eq. chan) write(99,'(" swup,win: ",4g13.4)')
     &  q2(1,1), csub(1,1), swup, wsin(1)

      return

C------------------------------------------------------------------------------

      entry sedfin
C                                          finish sediment balance computations
C------------------------------------------------------------------------------

      do n = 1, nps
C                                                  store outgoing concentration
        call store (outc(n,1), ilast, cs(n,nk,1))

        if (ltyp .eq. cmpd) call store (outc(n,2), ilast, cl(n,nk,2))

      end do

      wsi = 0.
C                              compute weight of incoming and outgoing sediment
      do n = 1, nps

        wso(n) = wso(n) * rho(n)
        wsi = wsi + wsin(n) * rho(n)

      end do
C                                 compute weight of eroded / deposited material
      wse = 0.

      do j = 1, nchan

        wsej = 0.

        do i = 2, nk

          wsej = wsej + depnet(i,j)

        end do

        wse =  wse + wsej

      end do

      wse = dx * qfac * wse
C                                          compute weight of suspended material
      wss = 0.
C      if(nps .gt. 5 .or. nchan .gt. ku .or. nchan .le. 0) stop ' lim'
      do n = 1, nps

        do j = 1, nchan
C
          dwss1 = a1(1,j) * cs(n,1,j)

          do i = 2, nk

            dwss2 = a1(i,j) * cs(n,i,j)
            wss = wss + rho(n) * 0.5 * (dwss2 + dwss1) * dx
            dwss1 = dwss2

          end do

        end do

      end do

      wss = qfac * wss
C                                                         pass values to writer
C------------------------------------------------------------------------------

      call swrt ( wso, wsi, wse, wss, qsp, tsp, outc)


      return

C------------------------------------------------------------------------------


      entry sedbf (k, a2, q2, topw)     ! k not used


C                                       baseflow: compute initial concentration
C------------------------------------------------------------------------------

      do i = 1, nk

        d = a2(i,1) / topw(i,1)
        v = q2(i,1) / a2(i,1)

        do n = 1, nps

          cs(n,i,1) = trnscap (d, v, slope(i,1), dpd(n), rho(n), grav)

        end do

      end do

      do n = 1, nps

        call store (outc(n,1), 1, cs(n,nk,1))

      end do

      return

C------------------------------------------------------------------------------


      entry sed0 ( nsect, idv, msg,  inflow, nup, nlat,    !typ,
     &            nx, omeg, delx, slp, qf )


C      nchan = nsect
      idu = id

      if (nsect .eq. 1) then
C                                        save some arguments as local variables
C------------------------------------------------------------------------------
        swup = 0.
        nk = nx
        omega = omeg
        comega = 1. - omega
        dx = delx
        xl = dx * float (nk - 1)
        qfac = qf
        qsp = 0.
        tsp = 0.
C        typl = typ   ! use ltyp
        over = .false.
C**  22/4/03
        ilast = 1    != 0
        nm = nlat
        nl = nlat
        nu = nup

        if (ltyp .eq. chan .or. ltyp .eq. cmpd) then
C                                                                obtain lateral
C                                                               inflow block(s)
C------------------------------------------------------------------------------

          do i = 1, nm

            do n = 1, nps

              j = abs (inflow(i+10))

              call old (j, attr(n,1), latc(n,i), ierr)

              if (ierr .gt. 0) call errxit
     &        (msg, 'lateral sediment concentration not found')

            end do

          end do

          do j = 1, nchan
C                                                initialize lateral inflow vars
            do n = 1, nps
              qcsl(n,j) = 0.
              qcsl1(n,j) = 0.
              qcsl2(n,j) = 0.
              dqcsl(n,j) = 0.
            end do

            do i = 1, nk
              qil(i,j) = 0.
              qilm(i,j) = 0.
            end do

          end do

        end if
C                                               obtain upstream inflow block(s)
C------------------------------------------------------------------------------

        do i = 1, nu

          do n = 1, nps

            inv = abs (inflow(i))

            call old (inv, attr(n,1), upc(i,n,1), ierr)

            if (ierr .gt. 0) call errxit
     &      (msg, 'upstream sediment concentration not found')

            if (inflow(i) .lt. 0) then
C                                                           neg value indicates
C                                                              compound channel
              call old (inv, attr(n,2), upc(i,n,2), ierr)

              if (ierr .gt. 0) call errxit
     &        (msg, 'upstream sediment concentration not found')

            end if

          end do

        end do

C        if (ltyp .eq. plane) qfac = qf

      end if

      if (nsect .eq. 2 .and. nlat .gt. nm .and. ltyp .ne. pipe) then
C                                                      check for lateral inflow
C                                             specified in overbank input block
C------------------------------------------------------------------------------

        nl = nlat
        j = abs (inflow(nl+10))

        do n = 1, nps

          call old (j, attr(n,1), latc(n,nl), ierr)

          if (ierr .gt. 0) call errxit
     &    (msg, 'lateral sediment concentration not found')

        end do
C                                           "over" indicates the single lateral
C                                            input goes to the overbank section
        if (nl .eq. 1) over = .true.

      end if
C------------------------------------------------------------------------------

      do n = 1, nps
C                                                     obtain blocks for outflow
C                                                         concentration storage
        call new (idu, attr(n,nsect), outc(n,nsect))
C                                                            set c = 0 at t = 0
        call store (outc(n,nsect), 1, zro)

      end do

      do i = 1, nk
        slope(i,nsect) = slp(i,nsect)
C        width(i,nchan) = wid(i,nchan)
      end do

      do n = 1, nps
        pros(n,nsect) = 0.
      end do
C                                                       get sediment parameters
C------------------------------------------------------------------------------

C                                                     erosion pavement fraction
      if (ltyp .ne. pipe) then

        call getr4 ('PA', 0, pav(nsect), ierr)

        if (ierr .gt. 0) pav(nsect) = 0.

        if (pav(nsect) .lt. 1.) then

          if (ltyp .eq. plane) then
C                                                       rain splash coefficient
            call getr4 ('SPL', 0, bspl, ierr)

            if (ierr .gt. 0)
     &      call errxit (msg, 'rain splash coefficient (SP) not found')

            bspl = rmspl * bspl

            call getr4 ('KS', 1, fmin, ierr)
C                                                          Ks (for rain splash)
            if (ierr .gt. 0) fmin = 0.

            fmin = rmks * fmin / (conv * 3600.)

          end if
C                                                          cohesion coefficient
          call getr4 ('CO', 0, bcoh(nsect), ierr)

          if (ierr .gt. 0) call errxit
     &    (msg, 'soil cohesion coefficient (CO) not found')

          bcoh(nsect) = bcoh(nsect) * rmcoh
C!! add warning exit
        else if(pav(nsect) .gt. 1.00001) then
          call errxit (msg, ' surface pavement coefficient cannot be gre
     &ater than 1.0 ')
 !
        else
C!! add limit 
          pav(nsect) = min(1.0,pav(nsect))
          bspl = 0.
          fmin = 0.
          bcoh(nsect) = 0.

        end if

      end if
C                                       check that class proportions sum to one
C                                                       particle size fractions
          do n = 1, nps

            call getr4 ('FR', n, prost, ierr)

            if (ierr .gt. 0) go to 1

            pros(nord(n),nsect) = prost

          end do
  1       sum = 0.

          do n = 1, nps
            sum = sum + pros(n,nsect)
          end do

          if (abs (sum - 1.0) .gt. 0.05)
     &    call errxit (msg, 'particle size fractions do not add to 1')

          nlp(nsect) = 1
C                                                   find index of largest class
2         if (pros(nlp(nsect),nsect) .le. 0.001) then

            nlp(nsect) = nlp(nsect) + 1

            go to 2

          end if
C                                                          initialize variables
C------------------------------------------------------------------------------

      qu(nsect) = 0.

      do n = 1, nps

        sd(n,nsect) = 0.
        csub1(n,nsect) = 0.
        csub2(n,nsect) = 0.
        dcsub(n,nsect) = 0.
        dwso(n) = 0.
        dwsol(n) = 0.
        wso(n) = 0.
        wsin(n) = 0.
        dwsup(n) = 0.
        dwsupl(n) = 0.

      end do

      do i = 1, nk

        depnet(i,nsect) = 0.

        do n = 1, nps

          supl(n,i,nsect) = 0.
          deprt(n,i,nsect) = 0.
          cs(n,i,nsect) = 0.
          cl(n,i,nsect) = 0.
          cml(n,i,nsect) = 0.

        end do
      end do

      return
C------------------------------------------------------------------------------


      entry sed00 ( vs, norder) !units, del, lim, np, dp, rh,



C      delt = del
C      limit = lim
C                                                get global sediment parameters
C------------------------------------------------------------------------------

C                                                    water temperature (C or F)
      call getr4 ('T', 0, temp, ierr)

        if (ierr .gt. 0) stop ' error - temperature not specified'


      do i = 1, 5
C                                                  particle diameter (mm or in)
        call getr4 ('DI', i, dpd(i), ierr)

        if (ierr .gt. 0) go to 20
C                                                       particle density (g/cc)
        call getr4 ('DE', i, rho(i), ierr)

        if (ierr .gt. 0) stop ' error - missing particle density'

      end do

      i = 6

20    nps = i - 1

      if (nps .lt. 1) stop ' error - particle classes not defined'

      np = nps

      if (units .eq. 1) then
C                                                                        metric
        grav = 9.81
        bdep = 656.

      else
C                                                                       English
        grav = 32.2
        temp = (temp - 32.) * 5. / 9.
        bdep = 200.

      end if
C                                                  kinematic viscosity (m**2/s)
      xnu = (.0000017756 / (1. + 0.03368 *
     &       temp + 0.000221 * temp * temp))
C                                                            convert to ft**2/s
      if (units .eq. 2) xnu = xnu * 10.76496

      sumvs = 0.

      do n = 1, nps
        dpd(n) = dpd(n) / conv
        vsl(n) = vsetl (rho(n), dpd(n), xnu, grav)
        sumvs = sumvs + 1. / vsl(n)
      end do
C                                        compute relative erodibility weighting
      do n = 1, nps

        pwt(n) = 1. / vsl(n) / sumvs
        nord(n) = n

      end do

      if (nps .gt. 1) then
C                                       order classes by increasing erodibility
        j = 1
30      jp = j + 1

        if (jp .le. nps) then

          l = j

          do i = jp, nps

            if (pwt(i) .gt. 0. .and. pwt(i) .lt. pwt(l)) l = i

          end do

          if (l .ne. j) then

            call switch (rho(j), rho(l))
            call switch (vsl(j), vsl(l))
            call switch (dpd(j), dpd(l))
            call switch (pwt(j), pwt(l))
            na = nord(j)
            nord(j) = nord(l)
            nord(l) = na

          end if

          j = j + 1

          go to 30

        end if

      else

        nord(1) = 1

      end if

      do n = 1, nps
       norder(n) = nord(n)
       vs(n) = vsl(n)
      end do
C      write(77,'(" order: ",5I3)')(norder(n),n=1,nps)
      return

      end
C------------------------------------------------------------------------------


      function splash (hb, rain, bspl, bdep)

C  Estimates the rate of removal of soil material from the surface caused by
C  the direct energy of the falling rain, in terms of volume of sediment per
C  unit time per unit area, [L/t].


C         estimate splash due to rainrate and depression of splash with depth

      splash = bspl * rain * rain * exp (-bdep * hb)
      if (splash .lt. 0.) splash = 0.

      return

      end

C------------------------------------------------------------------------------

      function trnscap (h, v, slope, dpd, rho, grav)

C  Finds Sediment transport capacity concentration as a function of hydraulic
C  conditions using the Engelund-Hansen relation.  H is hydraulic depth,
C  V is local velocity, dpd is particle diam, and rho is particle density.

      if (h .gt. 0.) then
C                                                        compute shear velocity
        ustar = sqrt (grav * h * slope)
        sm = rho - 1.
        trnscap = 0.05 * v * ustar**3 / grav / grav
     &            / dpd / sm / sm / h / rho

      else

        trnscap = 0.

      end if

      return

      end
C------------------------------------------------------------------------------


      function vsetl (ss, d, vnu, grav)

C  Finds settling velocity for a round particle of diameter d and specific
C  gravity ss, in water of viscosity vnu.


      iter = 0
      ca = 24 * vnu / d
      ck = 4. / 3. * grav * (ss - 1.) * d
      cb = 3. * sqrt (vnu / d)
      cc = .34
      vs = 2.8E05 * d * d

      do i = 1, 40

        f = ca * vs + cb * vs**1.5 + cc * vs * vs - ck
        fp = ca + 1.5 * cb * sqrt (vs) + 2. * cc * vs
        dvs = f / fp

        if (abs (dvs) .gt. 0.0001 * vs) then

          vs = vs - dvs
          if (vs .le. 0.) vs = 0.5 * (vs + dvs)

        else

          vsetl = vs

          return

        end if

      end do

      write (77, 5)
5     format (' error - unable to compute settling velocity for', /,
     &        '        the given temperature and particle size', /)
      stop ' error computing Vsetl for given Temp. and size '

      end
C------------------------------------------------------------------------------


      subroutine switch (a, b)

      temp = a
      a = b
      b = temp

      return

      end
