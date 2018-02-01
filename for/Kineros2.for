C   Code type: FORTRAN main program

C   Compiler: Standard Fortran 77/90

C   Programmed by: C. Unkrich

C   Date: 1/9/95.  modified 2/02

C   Description:

C     The main program unit for KINEROS2. Its job is to open the input files,
C     open the output file , get the elapsed time to the end of the simulation
C     and the time step from the user, and compute the number of time steps.
C     It then gets the characteristic length and the units specifier from the
C     parameter file. It calls a subroutine to initiallize storage for element
C     output and calls a subroutine to read the rainfall data. It also asks
C     the user if sediment is to be simulated. If so, it calls a subroutine to
C     get global sediment parameters and do any preliminary calculations that
C     are necessary. It then loops to read successive parameter blocks,
C     calling the routing subroutine indicated by the block name, passing the
C     file unit numbers for the parameter and output files. The loop continues
C     until the end of the parameter file is reached, at which time a summary
C     of the event is written to the output file.
C------------------------------------------------------------------------------

C   GLOBAL input Block (input file labels in parenthesis):

C     clen (C)       real       characteristic length, m or ft,

C     c (U)          char*1     units, 'M' = metric [units=1], 'E' = english,

C     ChRain (CHR)   char*1     rain on channels, Y or N

C   -- sediment --

C     temp (T)       real       water temperature, C or F,

C     dpd(5) (DI)    real       particle diameters, mm or in,

C     rho(5) (DE)    real       particle densities.
C------------------------------------------------------------------------------

C  Subroutines/functions (file):

CXC     banner         (banner.for)
C     clerk          (miscel.for)
C     reader         (reader.for)
C     getr4          (reader.for)
C     sed00          (kinsed.for)
C     rain           (rain.for)
C     wrt00          (writer.for)
C     event               "
C     errxit         (K2RunW.for)
C     plane          (plane.for)
C     channl         (channel.for)
CXC     pipe           (pipe.for)
C     pond           (pond.for)
C     adder          (adder.for)
C     inject         (inject.for)
C     urban          (urban.for)
C------------------------------------------------------------------------------
C
      subroutine K2main
C 
      use runpars  ! this communicates with graph, k2run, chann, etc.
C                    includes files = file unit numbers
      use elpars

      integer  i, ierr,  sfile, nord !nprint, limit, units, nps

      logical ::  op , rarea

!*8/2005 need to compute kinematic viscosity (using temperature) for laminar flow option in plane

	real :: temp, xnu

      character(LEN=1)   :: c
      character(LEN=10)  :: typl, dummy
      character(LEN=150) :: filnam
C 
      dimension sumbal(10), flag(5), vsl(5), dtm(4), nord(5) 
     &           !  ,qbal(10), dpd(5),rho(5)  into module


      equivalence (typl, dummy)

      data sfile /119/

C                                                     initialize global volumes
      do j = 1, 10
        sumbal(j) = 0.
      end do
C                                                  compute number of time steps
      limit = int (tfin / delt) + 1
C                                             convert time increment to seconds
      delt = delt * 60.
C!!    12/02
      depmax = 0.                             ! initialize
C                                                   read global parameter block
C------------------------------------------------------------------------------

      call reader (files(1), dummy, ierr)
      if (ierr .gt. 0) stop ' error - invalid global block'

      call getr4 ('C', 0, clen, ierr)
C                                                                   char length
      if (ierr .ne. 0) stop ' error - char. length not found'

      call geti4 ('NE', 0, nele, ierr)
C!!  set zero, not 1, when unspecified  18/12/02                       no of elements
      if (ierr .ne. 0) nele = 0
C                                                               system of units
C------------------------------------------------------------------------------

      call getstr ('U', 0, c, l, ierr)

      if (ierr .ne. 0) stop ' error - units not specified'

      if (c .eq. 'm' .or. c .eq. 'M') then
C                                                                        metric
        units = 1
        dlab = 'm.'
        dilb = 'mm'
        conv = 1000.
        wt = 1000.           !  kg/m3 water density
        arconv = 10000.      !  m^2 per ha.
        wconv = 1000.        !  kg/Tonne
C                          these factors go into runpars
      else if (c .eq. 'e' .or. c .eq. 'E') then
C                                                                       english
        units = 2
        dlab = 'ft'
        dilb = 'in'
        conv = 12.
        wt = 62.4    
        arconv = 43560.
        wconv = 2000.

      else

        stop ' error - units not recognized'

      end if
C
      call progress(1)
C
      call getstr ( 'CHR', 0, c, l, ierr)    ! rain on channel option
      
      chrain = .false.             ! this is default if CHRain not found
      
      if(ierr .eq. 0) then
        if(c .eq. 'y' .or. c .eq. 'Y') then
          chrain = .true.
        end if
      end if

!*8/2005 need temperature to compute kinematic viscosity for laminar flow option in plane

      call getr4 ('T', 0, temp, ierr)

      if (ierr .gt. 0) temp = 0.

      if (temp .gt. 0.) then

        if (units .eq. 2) temp = (temp - 32.) * 5. / 9.

        xnu = .0000017756 / (1.+ 0.03368 * temp + 0.000221 * temp*temp) ! m**2/s
C                                                            
        if (units .eq. 2) xnu = xnu * 10.76496 ! convert to ft**2/s

	else

        xnu = 0.

	end if
C
C optional initial water depth at surface:
C
      if(vAPI) then              ! value of api found.  Look for thbi
        write(files(3),95) api, dilb   ! using input units
        api = api/conv         ! converted to meters or feet 
  95  format(//' Run uses API option with value of ',f6.2,1x,a2)
C
      Else  
C no api.  Use SAT for entire profile
C
        api = -0.001
      End If

      if (sed) call sed00 ( vsl, nord) !units, delt, limit, nps, dpd, rho,

C                                                                    initialize
C                                                               output routines
C
      call wrt00 (files(3)) !, units, limit, delt, nps, dpd, rho, cour)

C                                                              allocate storage
      call clerk (limit)

C                                                            read rainfall data
      call rain (files(2))  !, units, limit, delt

      do j = 1, 4
        dtm(j) = delt
      end do
      nelt = 0
      ltab = 0
C                                              begin element processing loop...
C------------------------------------------------------------------------------

40    call reader (files(1), typl, ierr)


      if (ierr .gt. 0) then
        if (ierr .eq. 1) go to 50
C                                                                   end of file
        stop ' error - invalid input block'
C                                                                   input error
      end if

      diag = .false.
      nprint = 0
      nchan = 1     ! default
      nelt = nelt + 1
      ltab = ltab + 1         ! these will not agree only when there is a
C                               compound channel with overbank infil
C                                                                  print option
C------------------------------------------------------------------------------

      call geti4 ('PR', 0, nprint, ierr)

      if (ierr .eq. 0) then
        nlpr = nprint         ! pass this on thru module runpar to run prog

        if (nprint .eq. 4) then
C                                                      print to diagnostic file
          inquire (99, opened = op)
          if (.not. op) open (99, file = 'diagno.out',
     &                        status = 'unknown')
          write(99,1009) fname(1),fname(2)  ! how do we get these without module?
 1009 format(2x,"Diagnostic output for run with:"/5x,"Parameter File: ",
     %a60,/5x,"     Rain File: ",a60/) 
          diag = .true.
          nprint = 2   ! this allows printout of hydrograph for elements of interest

        else if (nprint .eq. 3) then
C                                                         open spreadsheet file
          call getstr ('FI', 0, filnam, j, ierr)

          if (ierr .gt. 0) then

C            write (*, 5)
C5           format (' error: PRINT = 3 & no file specified')
            stop ' error: PRINT=3, but no file specified in PAR file '

          else
C!!  add error directive
            open (sfile, file = filnam, status = 'unknown', err=90)

          end if

        else if (nprint .eq. 5) then
C                                                         open spreadsheet file
          call getstr ('FI', 0, filnam, j, ierr)

          if (ierr .gt. 0) then

            stop ' error: PRINT=5, but no file specified in PAR file '

          else

            open (xfile, file = filnam, status = 'unknown', err=90)

          end if
        end if
      end if
C
C------------------------------------------------------------------------------
      Rarea = .true.

      if (typl(1:5) .eq. 'PLANE') then
C                                                                        plane
        call plane (dtm, xnu)
!     &              qbal, delt, limit, clen, units, sed, diag, cour,

      else if (typl(1:7) .eq. 'CHANNEL') then
C                                                                        channel
        call channl (dtm)
        Rarea = chrain
C
      else if (typl(1:4) .eq. 'PIPE') then
C                                                                        pipe
        call pipe (dtm)                      

      else if (typl(1:6) .eq. 'INJECT') then
C                                                                        inject
        call inject (nord) !, delt,qbal, limit, sed, nps, rho

      else if (typl(1:4) .eq. 'POND') then
C                                                                        pond
        call pond (vsl) !qbal, delt, limit, units, sed, diag, nps, rho,

      else if (typl(1:5) .eq. 'URBAN') then                       
C** 4/03                                                                 urban
        call urban (dtm)      !qbal, diag,

      else if (typl(1:5) .eq. 'ADDER') then                       

        call adder ()

      else
C                                                             unrecognized type
        call errxit (typl, 'invalid element')

      end if

      if (typl(1:5) .ne. 'ADDER') then                       

        If(Rarea) sumbal(1) = sumbal(1) + qbal(1)
C                                                contribution to global volumes
        do i = 2, 9
          sumbal(i) = sumbal(i) + qbal(i)
        end do
C                                                            write element info
        call writer (sfile)

      end if

      go to 40
C                                             ...end of element processing loop
C------------------------------------------------------------------------------

C                                                                 final outflow
  50  sumbal(10) = qbal(10)
C                                                          write  event summary
      call event (sumbal, dtm)
C
      return
C!! new error printout when file unopenable
  90  continue
      write(files(3),'(" For print option 3, unable to open file
     & specified as: ",a60 )') filnam
      stop ' unable to open file given for output option 3 '

      end
