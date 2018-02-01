C ---------------file  of miscell. subprograms, including:
C  Name      type         function
C  -------  ----------    ----------
C  CLERK    subroutine    handles storage of data with internal file and
C                         multiple entries for reading and writing
C  new      entry
C  old      entry
C  Look     entry         like "old", but does not free the block
C  Store    entry
C  get      entry
C  stora    entry
C  geta     entry
C  orphan   entry         checks for unused element runoff
C  fmt10    subroutine    chooses size 10 float format for variable
C  streal8  subroutine    returns a real*8 variable from a chaacter string
C  ITER     subroutine    does generalized Newton-Raphson iteration
C  
C ----------------------------------------------------------------     
C   CLERK:  Code type: FORTRAN subroutine with additional entry points

C   Programmed by: C. Unk-Rich, revised by R.E. Smith 20.02.2002

C   Date: 1/4/95

C   Description:

C     Provides storage for element outputs such as discharge, sediment
C     concentration, etc. These values are then available as input to
C     downstream elements, output, etc. The storage entity is a a binary
C     scratch file conceptually partitioned into blocks equal in length to
C     the number of time steps. Calling programs use "new" to allocate a free
C     block. The maximum number of blocks is currently set at 1000. "Old" is
C     used to obtain the number of an existing block associated with a
C     particular identifier and attribute. When "old" is called, the block is
C     freed to be overwritten, allowing output values to be written one step
C     behind input values being read. "Store" and "get" are used to store and
C     retrieve the values for a given time step or index in the specified
C     block.

C     A pair of entry points, "stora" and "geta", are used to store and
C     retrieve contributing areas based on the element identifiers.

C     Another entry point, "orphan," was added to check for unprocessed blocks
C     at the end of the run (excluding the last element). The presence of an
C     unprocessed implies the element was not connected to the system -- an
C     "orphan."
C------------------------------------------------------------------------------

C  Arguments:

C  -- new, old --

C     id            int           element id number,
C     attr          char*2        attribute code denoting what is being stored,
C                                 i.e., runoff, sediment, etc.,
C     positn         int           storage block number,

C  -- new --

C     ierr          int           = 0 if block is found,
C                                 = 1 if not found,

C  -- store, get --

C     indx         int           storage location within a block,
C     value         real          value to be stored or retrieved,
C     positn         int           storage block number,

C  -- stora, geta --

C     id            int           element id number,
C     area          real          contributing area (sq.m or sq.ft).
C------------------------------------------------------------------------------


      subroutine clerk (lim)


      integer :: i, k, ierr, indx, positn, limit, id, norphs, icap, lim

      integer :: idir(1000), unconn(1000), list(1000)

      character(len=2) :: attr, adir(1000), free

      real, save :: areas(1000),  storall(1000000) !, adir(100)

      save  limit, icap, idir, list, adir

      parameter (free = '**')

      open (44, form = 'unformatted', access = 'direct', recl = 4,
     &      status = 'scratch')


      limit = lim
C                                    initialize attribute directory as all free
      do i = 1, 1000

        adir(i) = free
        idir(i) = 0

      end do

      do i = 1, 1000

        list(i) = 1.E8
        areas(i) = -1.1

      end do
C                                                 capacity of storall in blocks
      icap = 1000000 / limit

      return

C------------------------------------------------------------------------------


      entry new (id, attr, positn)


      do i = 1, 1000

        if(adir(i) .eq. free) go to 10

      end do

      call errxit ('CLERK', 'directory full')

10    idir(i) = id
      adir(i) = attr
      positn = i

      return

C------------------------------------------------------------------------------

      entry old (id, attr, positn, ierr)


      ierr = 0

      do i = 1, 1000

        if(idir(i) .eq. id .and. adir(i) .eq. attr) go to 20

      end do

      ierr = 1

      return

20    positn = i
C                                           block is now free to be overwritten
      adir(i) = free

      return

C------------------------------------------------------------------------------

      entry look (id, attr, positn, ierr)
C this entry is like OLD except the directory location is not freed:

      ierr = 0
      i = 1
C
      do while(i .le. 1000)
C
        if(idir(i) .eq. id .and. adir(i) .eq. attr) then
          positn = i
          return
        end if
C
        i = i + 1
      end do

      ierr = 1

      return


C------------------------------------------------------------------------------


      entry store (positn, indx, value)

C
      if (positn .le. icap) then
C                                                          put value in storall
        i = (positn - 1) * limit + indx
        storall(i) = value
C        write(77,'(" storing:",3i3,g13.3)')positn,indx,i,value

      else
C                                                                 write to file
        i = (positn - icap - 1) * limit + indx
        write (44, rec = i) value

      end if

      return

C------------------------------------------------------------------------------


      entry get (positn, indx, value)


      if (positn .le. icap) then
C                                                        get value from storall
        i = (positn - 1) * limit + indx
        value = storall(i)

      else
C                                                                read from file
        i = (positn - icap - 1) * limit + indx
        read (44, rec = i) value

      end if

      return

C------------------------------------------------------------------------------


      entry stora (id, areain)


      do i = 1, 1000

        if (areas(i) .le. -1.) go to 40

      end do

      call errxit ('CLERK','stora: array full')

40    list(i) = id
      areas(i) = areain

      return

C------------------------------------------------------------------------------


      entry geta (id, area)


      do i = 1, 1000

        if (list(i) .eq. id) go to 50

      end do

      call errxit ('CLERK','geta: id not found')

50    area = areas(i)
      areas(i) = -1.1

      return

C------------------------------------------------------------------------------


      entry orphan (idl, norphs, unconn)
C                          loop to look for unconnected, "orphaned" elements
      k = 0
C 
      norphs = 0

      do i = 1,1000
C** exclude urban interunit flows:   4/2003        
        if(idir(i) .ne. idl .and. adir(i) .eq. 'RO' .and. idir(i) .ne.
     & -888) then
C  outflow has not been read
          k = k+1
          unconn(k) = idir(i)
          norphs = norphs + 1
        end if
      end do


      return

      end
C-----------------------------------------------------------------------
      subroutine fmt10 (var, ref, string, n)


      integer fmt, n, j

      real :: var, valu, mult, ref

      character string*10


      valu = var

C                                          compute # of decimal places based on
C                                            7 sig. fig.'s from reference value
      if (ref .ge. 1.0) then

        j = 6 - int (alog10 (ref))

      else

        j = 7

      end if

      if (valu .ge. 1.E9 .or. -valu .ge. 1.E8) then
C                                                      use exponential notation
        assign 185 to fmt

      else if (j .ge. 7) then

        assign 105 to fmt

      else if (j .eq. 6) then

        assign 115 to fmt

      else if (j .eq. 5) then

        assign 125 to fmt

      else if (j .eq. 4) then

        assign 135 to fmt

      else if (j .eq. 3) then

        assign 145 to fmt

      else if (j .eq. 2) then

        assign 155 to fmt

      else if (j .eq. 1) then

        assign 165 to fmt

      else if (j .le. 0) then

        assign 175 to fmt

        mult = 1.

        if (var .ge. 1.E8) then
C                                                  remove nonsignificant digits
          mult = 100.

        else if (var .ge. 1.E7) then

          mult = 10.

        end if

        valu = mult * float (int (var / mult))

      end if

      write (string, fmt) valu
C                                     determine position of first nonblank char
      do n = 1, 10

        if (string(n:n) .ne. ' ') go to 20

      end do

20    return

105   format (f10.7)
115   format (f10.6)
125   format (f10.5)
135   format (f10.4)
145   format (f10.3)
155   format (f10.2)
165   format (f10.1)
175   format (f10.0)
185   format (e10.4)

      end
C-----------------------------------------------------------------------
C   Code type: Fortran subroutine STreal8

C   Compiler: Microsoft FORTRAN v5.1

C   Programmed by: C. Unk-Rich

C   Date: 3/15/95

C   Description: Decodes a character string into an 8-byte (double precision)

C                real number. A decimal point does not have to appear in the

C                string.

C   Arguments:

C     string    char*(*)     character string representing the numeric value,

C     var8      real*8       8-byte real number to be returned,

C     ierr      integer      0 = no error,

C                            1 = error (most likely an invalid character).
C-----------------------------------------------------------------------------


      subroutine streal8 (string, var8, ierr)


      integer ierr, i, j, k, m, ival

      real(kind=8) :: var8, base, signv

      character string*(*), c*1


      ierr = 1

      i = len (string)


      do j = i, 1, -1
C                                                                      find end
        if (ichar (string(j:j)) .gt. 32) go to 10

      end do

10    continue


      do i = 1, j
C                                                                find beginning
        if (ichar (string(i:i)) .gt. 32) go to 20

      end do

20    signv = 1.D0

      c = string(i:i)

C                                                               check for - / +
      if (c .eq. '-') then

        signv = -1.D0

        i = i + 1

      else if (c .eq. '+') then

        i = i + 1

      end if

      m = j + 1

C                                                  m is location of decimal pt.
      do k = i, j


        if (string(k:k) .eq. '.') then

          m = k

          go to 30

        end if

      end do

30    var8 = 0.D0

      base = 1.D0

C                                      decode digits to the left of decimal pt.
      do k = m - 1, i, -1

        ival = ichar (string(k:k)) - 48
C                                                              non-numeric char
        if (ival .lt. 0 .or. ival .gt. 9) return

        var8 = var8 + base * real(ival,4)

        base = base * 10.D0

      end do

      base = 0.1D0

C                                     decode digits to the right of decimal pt.
      do k = m + 1, j

        ival = ichar (string(k:k)) - 48
C                                                              non-numeric char
        if (ival .lt. 0 .or. ival .gt. 9) return

        var8 = var8 + base * real(ival,4)

        base = base * 0.1D0

      end do

      ierr = 0

      return

      end
C ----------------------------------------------------------------------
C   Code type: FORTRAN subroutine ITER

C   Compiler: Fortran77

C   Programmed by: C. Unkrich

C   Date: 10/93

C   Description:

C     Newton-Raphson iterative procedure to locate the root of
C     "errfct", where "errfct" is an equation written in residual
C     form, i.e. f(x) = 0. Successive estimates are required to
C     bracket the root; an estimate which falls outside the previous
C     interval will be recomputed based on a bisection of the
C     interval. The routine will return when the residual or
C     difference between successive estimates is smaller than the
C     tolerance. These test criteria are scaled by the magnitude of
C     the estimate when it is greater than one, so that the
C     relative accuracy is more or less independent of scale.
C------------------------------------------------------------------------------

C  Arguments:

C    errfct     --           external subroutine which computes the residual
C                            and derivative of the error function at an
C                            assumed root,

C         x     real         when called, contains an initial estimate of the
C                            root, on exit, the final estimate

C      xmin     real         lower bound on x,

C      xmax     real         upper bound for the initial interval,

C      ierr     integer      0 = no error,
C                            1 = did not converge.
C
C      trace    character    location of source of call
C------------------------------------------------------------------------------
C
      subroutine iter (errfct, x, xmin, xmax, ierr, trace)
C
C
      integer i, ierr
      character(LEN=6) trace

      data tol, rtol /1.E-7,1.e-4/

      xlow = xmin
      xhigh = xmax
      ierr = 1
      call errfct (xlow, fxlow, deriv)

      call errfct (xhigh, fxhigh, deriv)

C      if(xlow. gt. xhigh .and. fxlow .gt. fxhigh) then
      if(xlow. gt. xhigh) then
        write(77, 977) trace, xlow, xhigh
        stop ' Limits contradiction in ITER.  see outfile for details'
      end if
  977  format(' Lower bound higher than upper bound in call to ITER'/
     &       ' from ',a6,/'  xlow = ',g13.4,',  xhigh = ',g13.4)
      if(xlow. eq. xhigh) then
        x = xlow
        ierr = 0
        return
      end if

      x1 = x

      call errfct (x1, fx, deriv)
C                                                                iteration loop
C------------------------------------------------------------------------------

      do i = 1, 50

        if (deriv .ne. 0. .and. i .le. 10) then
C                                                       compute newton estimate
          x2 = x1 - fx / deriv

          if (x2 .le. xlow .or. x2 .ge. xhigh) then
C                                                                        bisect
            x2 = (xlow + xhigh) / 2.

          end if

        else
C                                                                        bisect
          x2 = (xlow + xhigh) / 2.

        end if
C                                            get function value of new estimate
        ofx = fx
  10    call errfct (x2, fx, deriv)
C
  900 format(' itr:', i2,3g12.4)

        test1 = abs(fx)
C                                             scale the test criteria if x2 < 1
        if (abs(x2) .gt. 1000.*tol)   test1 = abs (fx / x2)
C!!  modify denominator in case x2 is zero:  12/02
        xdenom = max(abs(x1),abs(x2))
        if(xdenom .le. 1.e-9) then
C          write(77,'(a6,i2,6g12.3)') trace,i,x,x1,x2,xmin,xmax
C          stop ' error stop in iter '
          test2 = 2.*rtol
        else
          test2 = abs ((x2 - x1) / xdenom)
        end if
        dfx = abs(ofx - fx) 
C                                                              convergence test
        if (test1 .lt. tol .or. test2 .lt. rtol .or. (dfx .lt. tol
     &  .and. abs(x2) .lt. tol)) then

          x = x2
C          ierr = 0
          ierr = -i
          return

        end if
C                                                              bracket the root
        if (fx * fxhigh .gt. 0.) then

          xhigh = x2

        else

          xlow = x2

        end if

        x1 = x2

      end do
C                                                            end iteration loop
C------------------------------------------------------------------------------

      return

      end
