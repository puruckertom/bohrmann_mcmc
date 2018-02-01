C   Code type: Fortran subroutine with additional entry points

C   Compiler: Fortran77

C   Programmed by: C. Unkrich

C   Date: 3/15/95

C   Description:

C     * THIS VERSION IS NOT CASE SENSITIVE

C     Reads a block of labeled data from a text file, passing the
C     data to calling programs when requested. Each input block
C     begins with a 'begin blockname', appearing on its own line,
C     and terminates with an 'end blockname', also on its own line:
C
C              begin blockname
C                .
C                .
C                .
C
C              end blockname
C
C     Two types of constructs are recognized, 'list' and 'column'.
C     These types can be mixed within a block but not within a line.

C     LIST input is:
C
C                label = datum
C
C     or:        label = datum1, datum2, ...
C
C     where the data list may continue over successive lines. There
C     may be more than one list on a line.
C
C     The COLUMN type of construct has labels on one line followed
C     by at least one line of corresponding data:
C
C              label1      label2      ...
C
C              datum       datum       ...

C              datum       datum       ...
C                .           .
C                .           .
C
C     For both types of constructs, having more than one datum
C     assigned to a given label implies an indexed variable, and
C     the datum returned will be determined by its position in the
C     list or column, as indicated by the 'indx' argument passed
C     by the calling program.
C
C     Acceptable delimiters include blanks, commas and tabs, and
C     additional blanks/tabs and blank lines can be used to make the
C     file more readable. Labels must be strictly alphanumeric,
C     i.e., 0 - 9 and a - z  or A - Z, with at least one alphabetic
C     character. Numeric values can have a leading + or -, without
C     intervening blanks. A comment is preceded by an ' ' and can
C     follow the data on the same line or be on a separate line.

C     The data block is read into a buffer when 'reader' is called,
C     which returns the name of the block. Subsequent access to the
C     data is provided by six additional entry points:

C     getr8 (returns a double precision value)
C     getr4 (returns a real value)
C     geti4 (returns an integer*4 value)
C     getstr  (returns a character string and its length)
C------------------------------------------------------------------------------

C   Arguments:

C   -- reader --

C     filen           integer      file unit connected to the input file,

C     blockname       char*(*)     block name read from the file; returned to
C                                  calling program,

C     ierr            integer      0 = no error,
C                                  1 = end of file,
C                                  2 = invalid block structure,
C                                  3 = memory limit exceeded.

C  -- getr4, getr8, geti4, getstr --

C     lname           char*(*)     string used to find the desired label in the
C                                  data block; it can be abreviated but should
C                                  contain enough characters to avoid confusing
C                                  it with other labels in the same block,
C                                  otherwise it will return the 1st ocurrence.
C                                  Also, the label search is CASE SENSITIVE.

C     indx            integer      index indicating the desired list element or
C                                  column row,

C     var8 (rget8)    double       returned value,

C     var4 (rget4)    real         returned value,

C     ivar4 (iget4)   integer      returned value,

C     string (sget)   char*(*)     returned string,

C     leng (sget)     integer      length of the string,

C     ierr            integer      0 = no error,
C                                  1 = 'lname' not found,
C                                  2 = datum not found,
C                                  3 = invalid chars in numeric datum.
C------------------------------------------------------------------------------

C   Subroutines/functions:

C   streal8         (miscel.for)
C------------------------------------------------------------------------------


      subroutine reader (filen, blockname, ierr)



      logical alpha, good

      integer filen, i, j, k, m, n, l, indx, nx, ivar4, leng, flag,
     &        ierr, col

      real var4

      double precision var8, val8

      character blockname*(*), buffer*210, block*100000, c*1,
     &          comma*1, equals*1, bang*1, lname*(*), astrsk*1,
     &          delim*1, leader*1, trailer*1, string*(*)

      save n, block

      parameter (comma = ',', equals = '=', bang = '!', astrsk = '*')


      ierr = 1
C                                                                 read new line
C------------------------------------------------------------------------------

10    read (filen, 20, end = 50) buffer

20    format (a)
C                                                          convert to uppercase
      call upper (buffer)

      do i = 1, 200
C                                                      look for 'BEGIN ' string
        c = buffer(i:i)

        if (c .eq. bang) go to 10
C                                                                  comment line
        good = alpha (c)

        if (good) then
C                                                             check for 'BEGIN'
          j = i + 4

          if (buffer(i:j) .eq. 'begin' .or.
     &        buffer(i:j) .eq. 'BEGIN') then

            do k = j + 1, 200
C                                                          look for 'blockname'
              c = buffer(k:k)

              if (c .ne. bang) then

                good = alpha (c)

                if (good) then
C                                                             'blockname' found
                  do l = k + 1, 200

                    good = alpha (buffer(l:l))

                    if (.not. good) then
C                                                            return 'blockname'
                      m = l - 1
                      blockname = buffer(k:m)

                      go to 25
C                                                        continue reading block
                    end if
                  end do
                end if
              end if
            end do
          end if
C                                                          invalid string found
          blockname = ' '
          ierr = 2

          go to 50

        end if
      end do
C                                                                    blank line
      go to 10

25    m = 1

      ierr = 2
C                                                                 read new line
C------------------------------------------------------------------------------

30    read (filen, 20, end = 50) buffer
C                                                          convert to uppercase
      call upper (buffer)

      flag = 0

      delim = astrsk

      j = 1

40    do i = j, 200
C                                             find first char of next substring
        c = buffer(i:i)

        if (c .eq. bang) go to 30
C                                                                       comment
        good = alpha (c)

        if (good) go to 45
C                                                                    first char
        if (c .eq. equals) delim = equals
C                                                     substring preceded by '='
      end do
C                                                                   end of line
      go to 30

45    if (flag .eq. 0) then
C                                     first substring - check for 'END ' string
        k = i + 3

        if (buffer(i:k) .eq. 'end ' .or. buffer(i:k) .eq. 'END ') then

          ierr = 0

          go to 50
C                                                                        return
        end if

        flag = 1

      end if
C                                             delimiter (comma, equals, astrsk)
      block(m:m) = delim

      do j = i + 1, 200
C                                                   find last char of substring
        c = buffer(j:j)
        good = alpha (c)
        if (.not. good) go to 47

      end do

47    m = m + 1
      n = j - i + m - 1

      if (n .ge. 99999) then
C                                                         memory limit exceeded
        ierr = 3

        go to 50
C                                                                        return
      end if

      block(m:n) = buffer(i:j)
C                                                 put substring into data block
      m = n + 1
      delim = comma

      go to 40
C                                                       look for next substring

50    return
C------------------------------------------------------------------------------


      entry getr8 (lname, indx, var8, ierr)
C                                                       return double precision


      flag = 1

      go to 55

C------------------------------------------------------------------------------


      entry getr4 (lname, indx, var4, ierr)
C                                                                   return real


      flag = 2

      go to 55

C------------------------------------------------------------------------------


      entry geti4 (lname, indx, ivar4, ierr)
C                                                              return integer*4


      flag = 3

      go to 55

C------------------------------------------------------------------------------


      entry getstr (lname, indx, string, leng, ierr)
C                                                            return char string


      flag = 5

C                                      look for substring corresponding to lname
C------------------------------------------------------------------------------

55    if (indx .eq. 0) then
C                                  non-indexed variable equivalent to index = 1
        nx = 1

      else

        nx = indx

      end if
C                                             convert keyword lname to uppercase
      call upper (lname)

      ierr = 2
      col = 1
      l = len (lname)
      k = l - 1

      do i = 2, n - k

        j = i + k

        if (block(i:j) .eq. lname) then
C                                                               substring found
          im1 = i - 1

          leader = block(im1:im1)
C                                                              must be preceded
C                                                                by a delimiter
          if (leader .eq. comma .or. leader .eq. astrsk) then

            j = j + 1

            do m = j, n
C                                                       find trailing delimiter
              trailer = block(m:m)

              if (trailer .eq. equals) then
C                                                                     list mode
                j = m

                go to 60

              else if (trailer .eq. comma .or.
     &                       trailer .eq. astrsk) then

                j = m - 1

                go to 70
C                                                                   column mode
              end if
            end do
          end if
        end if

        if (block(i:i) .eq. astrsk) then

          col = 1
C                                                         keep track of columns
        else if (block(i:i) .eq. comma) then

          col = col + 1

        end if

      end do

      ierr = 1
C                                                end of block - label not found
      return

60    do m = 1, nx

C                                    move to list element indicated by index --
C                                                    count substring delimiters
C------------------------------------------------------------------------------

        i = j + 1

        do j = i, n

          c = block(j:j)

          if (c .eq. comma) then

            go to 66

          else if (c .eq. astrsk) then
C                                                                end of line --
            if (m .eq. nx) then
C                                                         last value is the one
              go to 80

            else

              return
C                                                               value not found
            end if

          else if (c .eq. equals) then
C                                                     beginning of next list --
            return
C                                                               value not found
          end if

        end do
C                                                end of block - value not found
        if (m .lt. nx) return

66      continue

      end do
C                                                              decode substring
      go to 80

C                                            move to line indicated by index --
C                                                           count newline chars
C------------------------------------------------------------------------------

70    do m = 1, nx

        i = j + 1

        do j = i, n

          if (block(j:j) .eq. astrsk) go to 72

        end do

        if (m .le. nx) return
C                                                  end of block: line not found
72      continue

      end do
C                                            move to correct column by counting
C                                                          substring delimiters
      do m = 1, col

        i = j + 1

        do j = i, n

          c = block(j:j)

          if (c .eq. comma) then

            go to 74

          else if (c .eq. astrsk) then
C                                                                end of line --
C                                                       last column is the one?
            if (m .eq. col) then

              go to 80

            else

              return
C                                                               value not found
            end if
          end if
        end do
C                                                 end of block: value not found
        if (m .lt. col) return

74      continue

      end do

80    ierr = 0

      j = j - 1

      if (flag .eq. 5) then
C                                                      return string and length
        string = block(i:j)
        leng = j - i + 1

      else

        call streal8 (block(i:j), val8, ierr)
C                                                  decode into double precision
        if (ierr .gt. 0) then

          ierr = 3

          return

        end if

        if (flag .eq. 1) then
C                                                              double precision
          var8 = val8

        else if (flag .eq. 2) then
C                                                                          real
          var4 = sngl (val8)

        else
C                                                                     integer*4
          ivar4 = idint (val8)

        end if

      end if

      ierr = 0

      return

      end

C----------------------------------------------------------------------------

      logical function alpha (c)

C     .true.  c = {0-9, A-Z, ".", "-", ":", "\", "/", "_"}

C     .false.  c does not belong to the above set

      character*1 c
      integer i

      i = ichar (c)

      if ((i .le. 57 .and. i .ge. 45) .or. (i .le. 90 .and. i .ge. 65)
     &    .or. i .eq. 95 .or. i .eq. 58 .or. i .eq. 92)
     &then

        alpha = .true.

      else

        alpha = .false.

      end if

      return

      end

C----------------------------------------------------------------------

      subroutine upper (string)
C                                              converts string to all uppercase
      character*(*) string

      do j = 1, len(string)

        i = ichar (string(j:j))

        if (i .ge. 97 .and. i .le. 122) then

          i = i - 32
          string(j:j) = char(i)

        end if

      end do

      return

      end

