! Front end program unit to run K2 from a command shell (MS-DOS, Unix, etc)

	program K2shell

	use runpars
	use messagew
	use multip

  integer, parameter :: PROMPT = 1
  integer, parameter :: READ_FROM_FILE = 2
  integer, parameter :: NONE = 3

	integer ::            file_error, option

  logical ::            file_exists, char2log

  character(len=1)   :: answer, cour_char, sed_char, tabl_char
  character(len=20)  :: version, tfin_char, delt_char, rmks_char, rmn_char, api_char
  character(len=20)  :: rmcv_char, rmg_char, rmin_char, rmcoh_char, rmspl_char
  character(len=150) :: parfile, rainfile, outfile, multfile ! long file paths for AGWA

10 format(a) ! for general character i/o

  version = '3.2 (Dec 2003)'

! Write banner to screen

  write(*, 10)
  write(*, 10) '  KINEROS2'
  write(*, 10) '  Kinematic Runoff and Erosion Model'
  write(*, 10) '  Version '//version
  write(*, 10) '  U. S. Department of Agriculture'
  write(*, 10) '  Agricultural Research Service'
  write(*, 10)

! Read event-specific info, with two options:

! 1) Use info from the previous run saved in kin.fil
! 2) Prompt

  option = PROMPT

  inquire(file = 'kin.fil', exist = file_exists)

  if(file_exists) then
    write(*, "('$ Repeat previous run? ')")
		read(*, 10) answer

		if(char2log(answer)) then ! read event info from file kin.fil

		  open(4, file = 'kin.fil', status = 'old', iostat = file_error)
			if(file_error .ne. 0) call errxit('K2shell', "Can't open kin.fil")

      option = READ_FROM_FILE

    else

      option = PROMPT

    end if

  else

    option = PROMPT

	end if

  select case(option)

  case(READ_FROM_FILE)

    read(4, 10) parfile
    read(4, 10) rainfile
    read(4, 10) outfile
    read(4, 10) title
    read(4, 10) tfin_char
    read(4, 10) delt_char
    read(4, 10) cour_char
    read(4, 10) sed_char
    read(4, 10) multfile

!   For backwards compatibility, don't require or expect to read
!   the tabular summary option or API in kin.fil

    read(4, 10, iostat = file_error) tabl_char
    if(file_error .ne. 0) tabl_char = 'n'

    read(4, 10, iostat = file_error) api_char
    if(file_error .ne. 0) api_char = 'n'

    close(4)

  case(PROMPT)

    write(*, "('$            Parameter file: ')")
    read(*, 10) parfile
    write(*, "('$             Rainfall file: ')")
    read(*, 10) rainfile
    write(*, "('$               Output file: ')")
    read(*, 10) outfile
    write(*, "('$               Description: ')")
    read(*, 10) title
    write(*, "('$            Duration (min): ')")
    read(*, 10) tfin_char
    write(*, "('$           Time step (min): ')")
    read(*, 10) delt_char
    write(*, "('$ Courant Adjustment? (y/n): ')")
		read(*, 10) cour_char
    write(*, "('$           Sediment? (y/n): ')")
		read(*, 10) sed_char
    write(*, "('$   Multipliers? (y/n/file): ')")
    read(*, 10) multfile
    write(*, "('$    Tabular Summary? (y/n): ')")
    read(*, 10) tabl_char
    write(*, "('$ API Initializing? (val/n): ')")
    read(*, 10) api_char

  end select

! Open parameter, rainfall and output files - rainfall file is treated as
! optional (file unit set to -1 tells subroutine rain that no file exists
! and to create a single raingage entry with zero depths)

  open(files(1), file = parfile, status = 'old', iostat = file_error)
	if(file_error .ne. 0) call errxit('K2shell', "Can't open "//parfile)

  open(files(2), file = rainfile, status = 'old', iostat = file_error)
	if(file_error .ne. 0) then
    if(len_trim(rainfile) .eq. 0) then ! no rainfall file specified     
      files(2) = -1
    else
      call errxit('K2shell', "Can't open "//rainfile)
    end if
  end if

  open(files(3), file = outfile, status = 'unknown', iostat = file_error)
	if(file_error .ne. 0) call errxit('K2shell', "Can't open "//outfile)

! Decode character input into floating and logical values (module runpars)

  read(tfin_char , *) tfin
  read(delt_char , *) delt
  cour = char2log(cour_char)
  sed = char2log(sed_char)
  tabl = char2log(tabl_char)
  vAPI = char2log(api_char(1:1))
  if(.not. vAPI) then
    api = 0.
  else
    read(api_char, *) api
  end if

! Three options for parameter multipliers:

! 1) Prompt
! 2) Use file specified by user or info from the previous run saved in mult.fil
! 3) No multipliers 

  if(len_trim(multfile) .gt. 1) then ! user specified a file

    open(4, file = multfile, status = 'old', iostat = file_error)
		if(file_error .ne. 0) call errxit('K2shell', "Can't open "//multfile)

    option = READ_FROM_FILE

  else
    if(char2log(multfile)) then ! user wants multipliers

      inquire(file = 'mult.fil', exist = file_exists)

      if(file_exists) then
        write(*, "('$ Use previous multipliers? ')")
	  	  read(*, 10) answer

		    if(char2log(answer)) then ! read multipliers from mult.fil

		      open(4, file = 'mult.fil', status = 'old', iostat = file_error)
			    if(file_error .ne. 0) call errxit('K2shell', "Can't open mult.fil")

          option = READ_FROM_FILE

        else

          option = PROMPT

        end if

      else

        option = PROMPT

      end if

    else ! user doesn't want multipliers

      option = NONE

    end if
  end if

  select case(option)

  case(READ_FROM_FILE)

    read(4, 10) rmks_char
    read(4, 10) rmn_char
    read(4, 10) rmcv_char
    read(4, 10) rmg_char
    read(4, 10) rmin_char
    read(4, 10) rmcoh_char
    read(4, 10) rmspl_char
    close(4)

  case(PROMPT)

    write(*, "('$            Ks multiplier (0-1): ')")
    read(*, 10) rmks_char
    write(*, "('$ Manning/Chezy multiplier (0-1): ')")
    read(*, 10) rmn_char
    write(*, "('$            CV multiplier (0-1): ')")
    read(*, 10) rmcv_char
    write(*, "('$             G multiplier (0-1): ')")
    read(*, 10) rmg_char
    write(*, "('$  Interception multiplier (0-1): ')")
    read(*, 10) rmin_char
    write(*, "('$ Soil cohesion multiplier (0-1): ')")
    read(*, 10) rmcoh_char
    write(*, "('$   Rain splash multiplier (0-1): ')")
    read(*, 10) rmspl_char

  case(NONE)

    rmks_char  = '1.0'
    rmn_char   = '1.0'
    rmcv_char  = '1.0'
    rmg_char   = '1.0'
    rmin_char  = '1.0'
    rmcoh_char = '1.0'
    rmspl_char = '1.0'

  end select

! Decode multipliers into floating values (module multip)

  read(rmks_char, *) rmks
  read(rmn_char, *) rmn
  read(rmcv_char, *) rmcv
  read(rmg_char, *) rmg
  read(rmin_char, *) rmin
  read(rmcoh_char, *) rmcoh
  read(rmspl_char, *) rmspl

! Save event info to kin.fil

  lpar = len_trim(parfile)
  lrain = len_trim(rainfile)
  lout = len_trim(outfile)
  lmult = len_trim(multfile)

  open(4, file = 'kin.fil', status = 'unknown')

  write(4, 10) parfile(1:lpar)
  write(4, 10) rainfile(1:lrain)
  write(4, 10) outfile(1:lout)
  write(4, 10) title
  write(4, 10) tfin_char
  write(4, 10) delt_char
  write(4, 10) cour_char
	write(4, 10) sed_char
  write(4, 10) multfile(1:lmult)
  write(4, 10) tabl_char
  write(4, 10) api_char
  close(4)

! Save multipliers to mult.fil

  open(4, file = 'mult.fil', status = 'unknown')

  write(4, 10) rmks_char
  write(4, 10) rmn_char
  write(4, 10) rmcv_char
  write(4, 10) rmg_char
  write(4, 10) rmin_char
  write(4, 10) rmcoh_char
  write(4, 10) rmspl_char
  close(4)

! Also write info to output file

  if(files(2) .eq. -1) rainfile = 'None'

   write(files(3), 20) version, title, parfile(1:lpar), rainfile(1:lrain), tfin_char, &
                       delt_char, cour_char, sed_char, multfile(1:lmult), tabl_char, api_char

20 format(' KINEROS2 Version ',a//                                                   &
          ' Title: ',a//                                                             &
          ' Parameter File Used....... ',a/                                          &
          ' Rainfall File Used........ ',a//                                         &
          ' Length of Run, minutes.... ',a/                                          &
          ' Time Step, minutes........ ',a/                                          &
          ' Use Courant criteria?..... ',a/                                          &
          ' Simulate Sed. Transport?.. ',a/                                          &
          ' Multiplier file (if any).. ',a/                                          &
          ' Tabular Summary?.......... ',a/                                          &
          ' API Initializing?......... ',a)

   write(files(3), 30) rmks_char, rmn_char, rmcv_char, rmg_char, &
                       rmin_char, rmcoh_char, rmspl_char

30 format(/' Multipliers Used:'//                                &
           ' Saturated Conductivity... ',a/                      &
           ' Manning n................ ',a/                      &
           ' CV of Ksat............... ',a/                      &
           ' Capillary Drive Coeff.... ',a/                      &
           ' Intercepted Depth........ ',a/                      &
           ' Sediment Cohesion Coeff.. ',a/                      &
           ' Sediment Splash Coeff.... ',a/)

! Finally, run kineros2

  call K2main

end


logical function char2log(char)

! Converts character values Y or N into false or true logical

	character(LEN=1) :: char

	if(char .eq. 'Y' .or. char .eq. 'y') then
		char2log = .true.
	else
		char2log = .false.
	end if
	return
end

subroutine errxit (id, message)

! Write error message to stdout & output file

   use runpars

	 character(len=*) id, message

   write(files(3), 10) id, message
	 write(*, 10)        id, message
10 format(//, ' error - ', a, /, 1x, a)

	 stop ' '
end


subroutine progss (id)

! Report overall computational progress - note that the compiler
! should be configured to interpret the first character of an output
! record as a carriage control indicator (fortran carriage control).
! For Compaq Visual Fortran, use the "/ccDefault: Fortran" or "/vms"
! option.

  integer :: nv

  character(len=*) :: id

  write (*, "('+ Processing ', a, '             ')") id

  return

entry progress(nv)

! This is a dummy version of a routine in Roger's Lahey Win32 version.
! Since it is called from K2main before the element loop, I'll use it
! to insert a line before the progress report line.

  write (*, '(/)')

  return
end
