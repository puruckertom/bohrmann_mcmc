C   Description:

C   Combines water and sediment from up to 10 upstream elements

C-------------------------------------------------------------------------------

C  Global parameters:

C     sed            log           .true. = route sediment (module runpars)

C     nps            int            number of particle classes (module runpars)

C     limit          int            number of time steps (module runpars)

C     typname        char          'ADDER' (module elpars)

C-------------------------------------------------------------------------------

C  Input Block (input file labels in parenthesis):

C     id (ID)        int*4         identifier (up to 6 digits)

C     idup (UP)      int*4         identifiers of up to 10 upstream
C                                  contributing elements

C-------------------------------------------------------------------------------

C  Subroutines/functions:

C     getr4         (reader.for)
C     getstr             "
C     geti4              "
C     new           (miscel.for)
C     old                "
C     get                "
C     store              "
C     geta               "
C     stora              "
C     errxit        (errxit.for)
C-------------------------------------------------------------------------------

      subroutine adder ()

      use runpars
      use elpars

      implicit none

      integer ierr, i, j, ms, n, nu, outr, outq 
      integer upq(10), idup(10), upr(10), upc(10,5), outc(5)

      real qrain, qu_sum, qru, qu, cu, area, coarea, qs_sum(5)

      character msg*18
      character attr*2(5)

      data attr/'S1', 'S2', 'S3', 'S4', 'S5'/


      call geti4 ('ID', 0, id, ierr)

      if (ierr .eq. 3) then

        call errxit ('ADDER', 'invalid identifier (ID)')

      else if (ierr .ge. 1) then

        call errxit ('ADDER', 'missing identifier (ID)')

      end if

      msg(1:16) = typname(1)
      write (msg(15:18), '(i4)') id

      do ms = 1, 5
        if (msg(ms:ms) .ne. ' ') go to 1
      end do

1     continue

      nu = 0

      do j = 1, 10

        call geti4 ('UP', j, idup(j), ierr)

        if (ierr .eq. 3) call errxit
     &  (msg(ms:18), 'invalid upstream identifier (UP)')

        if (ierr .eq. 0) then

          nu = nu + 1

        else

          go to 11

        end if
      end do

11    continue

      coarea = 0.

      do j = 1, nu

        upr(j) = 0

        call old (idup(j), 'RA', upr(j), ierr)

        call old (idup(j), 'RO', upq(j), ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:18), 'upstream inflow not found')

        call geta (idup(j), area)
        coarea = coarea + area

        if (sed) then

          do n = 1, nps

            call old (idup(j), attr(n), upc(j,n), ierr)

            if (ierr .gt. 0) call errxit
     &      (msg, 'upstream sediment concentration not found')

          end do

        end if

      end do

      call stora (id, coarea)

      call new (id, 'RA', outr)
      call new (id, 'RO', outq)

      if (sed) then

        do n = 1, nps

          call new (id, attr(n), outc(n))

        end do

      end if

      do i = 1, limit - 1

        qrain = 0.
        qu_sum = 0.

        if (sed) then
          do n = 1, nps
            qs_sum(n) = 0.
          end do
        end if

        do j = 1, nu

          if (upr(j) .gt. 0) then

            call get (upr(j), i, qru)

            qrain = qrain + qru

          end if

          call get (upq(j), i, qu)

          qu_sum = qu_sum + qu

          if (sed) then

            do n = 1, nps

              call get (upc(j,n), i, cu)

              qs_sum(n) = qs_sum(n) + cu * qu ! sed. discharge

            end do

          end if

        end do

        call store (outr, i, qrain)
        call store (outq, i, qu_sum)

        if (sed) then
          do n = 1, nps
            call store (outc(n), i, qs_sum(n)/qu_sum)
          end do
        end if

      end do

      return

      end
