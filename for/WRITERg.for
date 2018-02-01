C ---------------file of output subprograms, including:
C  Name       type        function
C  -------  ----------    ----------
C  writer   subroutine    handles writing
C   qwrt      entry
C   swrt      entry
C   event     entry         see below
C  blank1   subroutine    writes a blank line to the output file
C  blank2   subroutine    Write a blank line to the spreadsheet file.
C ------------------------------------------------------------------------
C   Code type: FORTRAN subroutine + 4 entry points

C   Programmed by: C. Unkrich
C   Revised: R. E. Smith  2/2002

C   Date: 7/6/95

C   Description:

C     -- writer --

C     For NPRINT > 0, writes an element volume and sediment yield summary to
C     the output file. The summary includes volumes in cu.m or cu.ft and mm or
C     inches over area.  For planes and urban elements area is the element
C     area; for channels, ponds, and pipes it is the upstream contributing
C     area.  Peak flow in mm/hr or in/hr for all elements refers to
C     contributing area.  In the summary, rain volume is less interception.

C     For NPRINT = 1, only the summary is written.  
C
C     For NPRINT = 2 or 4, an 80-column listing with total sediment discharge
C     is also written. For nprint = 3, a separate file is written which
C     includes sediment discharge by particle class, with labels enclosed in
C     quotes and columns delimited by commas, for convenient import into
C     spreadsheet programs. Called from the main program after each element
C     is processed.

C     -- qwrt --

C     Transfers summary volume balance information and discharge storage
C     locations to writer's local variables. Called from plane, channl, pipe,
C     pond, inject, infil, etc.

C     -- swrt --

C     Transfers summary sediment balance information and sediment storage
C     locations to writer's local variables. Called from sedsum.

C     -- event --

C     Called by the main program to write the end-of-event volume and sediment
C     yield summary to the output file and screen. Volumes are written as both
C     depth (mm or in) over basin area and cubic meters/feet. Sediment yield
C     is given as weight per hectare or acre. The weight measure is in
C     kilograms or pounds if the total yield is less than one ton (1000 kg or
C     2000 lb), otherwise it is in tons. A volume balance error is computed
C     and rounded to the nearest whole percent, unless it is less than one
C     percent, in which case it is noted as such, to imply that it is
C     insignificant. When the error is between 1 - 5%, it is simply reported,
C     and when it is greater than 5% the statement is concluded with "!" to
C     imply that there is a potential problem.

C     In the volume summary, only nonzero volumes are reported. The numerical
C     format for volumes and yields is adjusted for magnitude, such that there
C     are 7 significant digits (roughly the precision) with a maximum of 7
C     decimal places.

C     If there is more than one particle class, sediment yield is given by
C     class, along with the equivalent percentage by weight of each class
C     compared to the total yield.

C     -- qwrt00 --

C     Transfers global values to and initiallizes local variables used by
C     writer or event. This routine is called only once by the main  program.
C------------------------------------------------------------------------------

C   Entry points:

C     qwrt000     reads & initializes global variables related to flow routing,

C     qwrt        retains volume balance components and storage locations for
C                 element discharge,

C     swrt        retains sediment balance components and storage locations for
C                 element sediment discharge,

C     event       writes end-of-event summary,
C
C     wrt00  
C------------------------------------------------------------------------------

C   Arguments:

C     idstr       char*(*)      element identifying string,

C     k           int           = 1 for all elements except:
C                               = 2 for compound channels,

C     ltyp       int           index for type of element:
C                               = 0 for a plane,
C                               = 1 for a channel,
C                               = 2 for a pipe,
C                               = 3 for a pond,
C                               = 4 for a compound channel,
C                               = 5 for an injection site,
C                               = 6 for an urban element,
C                               = 7 for overbank part of channel

C     qpk         real          peak discharge (cu.m/s or cu.ft/s),

C     tpk         real          time to peak (sec),

C     vi          real          total inflow volume (cu.m or cu.ft),

C     sumwsi      real          total sediment in (equivalent volume of water),

C     sumwse      real          sediment deposited (+) or eroded (-),

C     sumwss      real          suspended sediment,

C     sumwso      real          sediment out (by particle class),

C     area        real          contributing area (sq. m or sq. ft),

C     storq       int           location of element discharge values,

C     storc       int           location of sediment concentration,

C     qbal        real          array with global volume balance components;
C                               qbal(1)  = area                  (sq.m/sq.ft),
C                               qbal(2)  = rainfall              (cu.m/cu.ft),
C                               qbal(3)  = baseflow              (cu.m/cu.ft),
C                               qbal(4)  = injection             (cu.m/cu.ft),
C                               qbal(5)  = initial pond storage  (cu.m/cu.ft),
C                               qbal(6)  = plane infiltration    (cu.m/cu.ft),
C                               qbal(7)  = channel infiltration  (cu.m/cu.ft),
C                               qbal(8)  = interception          (cu.m/cu.ft),
C                               qbal(9)  = storage               (cu.m/cu.ft),
C                               qbal(10) = outflow               (cu.m/cu.ft),
C                               qbal(11) = overbank area               "
C                               qbal(12) = overbank rainfall           "
C                               qbal(13) = channel bottom infil        "

C     file        int           output file unit,

C     sfile       int           spreadsheet file unit

C     nprint      int           print option,

C     dtm(4)      real          smallest dt to satisfy the Courant condition
C                               during each quarter of the simulation,

C     units       int           1 = metric, 2 = english,

C     lim         int           number of time steps,

C     del         real          user-specified time step (sec),

C     np          int           number of particle classes,

C     dpd         real          particle diameters (m or ft),

C     rho         real          particles densities,

C     storr       int           location of rain rate,

C     qrmax       real          maximum rain rate (cu.m/s or cu.ft/s),

C     courc       log           .true. = time step is being adjusted.
C------------------------------------------------------------------------------

C   Subroutines:

C     fmt10       (miscel.for)
C     orphan      (miscel.for\clerk)
C------------------------------------------------------------------------------

      subroutine writer ( sfile )

      use elpars
      use runpars
C                                                                     arguments
C------------------------------------------------------------------------------

C      logical courc

      integer :: fileo, storq, storc, typ, sfile, storr !, id
      integer, dimension(1000) :: unconn

      dimension storq(nchan), storc(5,nchan), sumwso(5), dtm(4) 
     &            ,sumbal(10)

      character idstr*(*)


C                                                               local variables
C------------------------------------------------------------------------------

      logical errout, fill2, fillov, lowfin(2), twol, chtyp

      integer i, j, l, m, n, k, file1, file2, ierr, 
     &         flag, norphs, ii, outr, idl, fintgs !, ltyp, lfh 
      integer, save :: outc(5,2), outq(2)
C                                                   8
      integer,dimension(10) :: minl, nx, je = (/7,9,8,8,8,6,7,7,7,8/)
      dimension  qso(5,2), pso(5), label(9), q(2),
     &           floss(2)
C** dp() moved to here
      real,dimension(5),save :: wso, wsj, dp
C
      character(LEN=9),dimension(10) :: fmdat = ',1x,f10.2'
      character(LEN=6),dimension(10) :: fmhed = ',1x,A9'
      Character(LEN=9),dimension(3,10) :: hedstr = '         '
      Character(LEN=14) :: charid, head(3)
      Character(LEN=33) :: fminj = '(a14,t50,f10.2                   '
      Character(LEN=120) :: fmthl, fmtdat
      Character(LEN=9) :: blank9 = '         '
      real, dimension(10) :: dval
      logical :: winj
C
      character qlab1*2, qlab2*5, qlab3*7, qlab4*5, slab1*5, slab2*7,
     &          slab3*4, alab*2, label*21, wlab*2, buffer*400,
     &          string*10, header1*67, header2*74, lines*71, c1*1, c2*2,
     &          c3*3, idbuff*20, arlab*4
      character(len=139), save :: colunits
      character(len=139), save :: head1, head2
      character(LEN=26), save :: str2a = '  Upper Layer:    Subsoil ', 
     &                     str2b = 'IniAvail  Over FC  Infil. '
      character(len=10), save :: fmtstr, fmttab
      character(LEN=11), save :: empty
C      character(len=11),dimension(15) :: tabstrng
      character(len=1) :: aster = ' '
      character(LEN=2) :: sedlab
      character(LEN=3) :: adjnote = '   '
      character(LEN=8) :: adds1 = '        ', adds2 = '        '
      character(len=8),dimension(0:7) :: eltype = (/'   Plane',
     &  ' Channel','   Pipe ','   Pond ','Compound','  Inject',
     &  '  Urban ','Overbank'/)
      common /writ/ file1, file2

      data label(1) /'             Rainfall'/
      data label(2) /'             Baseflow'/
      data label(3) /'            Injection'/
      data label(4) /' Initial pond storage'/
      data label(5) /'   Plane infiltration'/
      data label(6) /' Channel infiltration'/
      data label(7) /'         Interception'/
      data label(8) /'              Storage'/
      data label(9) /'              Outflow'/

      data header1 /'             Water balance                         
     &Sediment balance'/

      data header2 /'  Elapsed Time      Rainfall      Outflow      Outf
     &low      Total Sediment'/

      data lines /' --------------------------------------        ------
     &------------------'/

      data  errout, fill2, fillov, lowfin(1), lowfin(2) /5*.false./
C
      save outr, idbuff, qp, tp, coarea, conar, vin, rfmax, alab, qlab1,
     &     qlab2, qlab3, qlab4, slab1, slab2, wlab, fill2, filltv, dt,
     &     filtm, lowfin, floss, fillov, filtv, filtvo, filtmo, slab3,
     &     qsp, tsp, wsi, wse, wss, arlab, dref
C
      head1(1:96)= '   ID   Element        Areas           Inflow   Rain
     &fall    Outflow    Peak    Total    Initial '
      head1(123:139) = '                '
      Head2(123:139) = '                '
      head2(1:96)= '         Type     Element Cumulated                
     &                   Flow    Infil     Water  '
C
      empty = '           '

C------------------------------------------------------------------------------

C  area over which rain is accumulated
      areaR = qbal(1) 
      ltov = ltab
      if(ltyp .eq. 2 .or. ltyp .eq. 5) then
        areaR = 0.  ! no rain for pipes or injects
        sumtab(ltab)%twola = .false.
C
      else if(ltyp .eq. 1 .or. ltyp .eq. 4 ) then  ! channels
        if(ltyp .eq. 4) then    !           compound
          if(qbal(6) .gt. 1.e-6) then        ! overbank has infil
            ltov = ltab + 1
            AreaR = AreaR + qbal(11)
            sumtab(ltov)%idel = id
            sumtab(ltov)%are = qbal(11)
            sumtab(ltov)%volrn = qbal(12)
            sumtab(ltov)%ftot = qbal(6)
            sumtab(ltov)%itype = 7
          else
            qbal(1) = qbal(1) + qbal(11)
            qbal(2) = qbal(2) + qbal(12)
          end if
        end if
        sumtab(ltab)%abot = barea
        sumtab(ltab)%vbot = qbal(13)
        if(.not. chrain) areaR = 0.
C
      end if
C
      vr = qbal(2) - qbal(8)
      vf = qbal(6) + qbal(7)
      area1 = qbal(1)           ! this is element area. 
C
      sumtab(ltab)%cumare = coarea            ! these two are used in final graf
      sumtab(ltab)%itype = ltyp
      if(tabl) then
        sumtab(ltab)%idel = id
        sumtab(ltab)%are = area1
        sumtab(ltab)%volro = qbal(10)
        sumtab(ltab)%ropeak = qp * conar
        sumtab(ltab)%volrn = vr
        sumtab(ltab)%volin = vin
        if(ltyp .ne. 4) then
          sumtab(ltab)%ftot = vf
        else 
          sumtab(ltab)%ftot = qbal(7)
        end if

      end if
C
      if (nprint .eq. 0) return

      file2 = sfile
C - - - - - - - - - - - - - - - - - - -  nprint > 0: write volume & sediment summary
C                                   m^2 or ft^2
      chtyp = .false.
      if(ltyp .eq. 1 .or. ltyp .eq. 2 .or. ltyp .eq. 4) chtyp = .true.
C
      if(area1 .lt. coarea .and. .not. chtyp) aster = '*'   ! comes via entry qwrt
      if (area1 .eq. 0.) area1 = coarea

      call blank1
      call write1 (' '//idbuff)

      call blank1
C                                                             contributing area
      area2 = coarea / arconv                              !   in hectares;

      call fmt10 (area2, area2, string, j)

      m = 35 - j
      buffer(1:m) = ' Contributing area = '//string(j:10)//' '//alab

      call write1 (buffer(1:m))

      call blank1
C
      if(rsoil) then
        write(buffer, "(t44,'theta   rel.sat.')") 
        call write1 (trim(buffer))
        write(buffer, "('  surface initial water content =         ',
     &        f6.4,3x,f6.4)") thinu, satin
        call write1 (trim(buffer))
        write(buffer, "('  estimated wilting point water content = ',
     &        f6.4,3x,f6.4)") thwilt, satwl
        call write1 (trim(buffer))
        write(buffer, "('  estimated field capacity =              ',
     &  f6.4,3x,f6.4)") thfc, satfc
        call write1 (trim(buffer))
	end if

      call blank1 
C     
      call fmt10 (qp, qp, string, j)

      m = 33 - j
      buffer(1:m) = ' Peak flow = '//string(j:10)//' '//qlab3

      if (coarea .gt. 0.) then
C                                                        convert cu ?/s -> ?/hr
        qp2 = qp * conar

        call fmt10 (qp2, qp2, string, j)

        n = m + 19 - j
        buffer(m+1:n) = '('//string(j:10)//' '//qlab4//')'
        m = n

      end if

      tp2 = tp / 60.

      write (string, fmtstr) tp2

      do j = 1, 10
        if (string(j:j) .ne. ' ') go to 10
      end do

10    n = m + 20 - j
      buffer(m+1:n) = ' at '//string(j:10)//' min'

      call write1 (buffer(1:n))

      call blank1
C
      if(fill2) then
        buffer(1:37) = '           Upper soil saturated with '
        write(string,'(f6.1,4x)') filtv
        buffer(38:51) = string(1:7)//qlab1//'. at '
        write(string,'(f6.1,4x)') filtm
        buffer(52:63) = string(1:6)//' min.'
        call write1 (buffer(1:63))
        call blank1
        fill2 = .false.
      end if      

      if(fillov) then
        buffer(1:37) = ' Overbank: Upper soil saturated with '
        write(string,'(f6.1,4x)') filtvo
        buffer(38:51) = string(1:7)//qlab1//'. at '
        write(string,'(f6.1,4x)') filtmo
        buffer(52:63) = string(1:6)//' min.'
        call write1 (buffer(1:63))
        call blank1
        fillov = .false.
      end if
C
      ilst = 1
      ib = 0
      if(ltyp .eq. 4) then
        ilst = 2
        ib = 3
      end if
C
      do ie = 1,ilst
        if(lowfin(ie)) then
          ity = ltyp + (ie-1)*ib
          write(string,'(f6.2,1x)') floss(ie)
          buffer(1:9) = string(1:7)//qlab1
          buffer(10:43) = ' Infiltration into lower layer of '
          write(string,'(a8)') eltype(ity)
          buffer(44:51) = string(1:8)
          call write1 (buffer(1:51))
          call blank1               ! writes blank line = skip line
          lowfin(ie) = .false.
        end if
      end do
C
      if (sed) then

        qsp2 = qsp * wt

        call fmt10 (qsp2, qsp2, string, j)

        m = 48 - j
        buffer(1:m) = ' Peak sediment discharge = '//string(j:10)
     &                //' '//slab3//' at '

        tsp2 = tsp / 60.
        write (string, fmtstr) tsp2

        do j = 1, 10

          if (string(j:j) .ne. ' ') go to 20

        end do

20      n = m + 15 - j
        buffer(m+1:n) = string(j:10)//' min'

        call write1 (buffer(1:n))

        call blank1

      end if
C                                                            write header lines
      if(chtyp) then
        me = 26
        buffer(1:26) = header1(6:31)
        if(sed) then
          me = 52
          buffer(27:52) = header1(42:67)
        end if
        call write1(buffer(1:me))
C
        buffer(1:26) = lines(1:26)
        buffer(27:33) = '       '
        me = 26
        if(sed) then
          me = 57
          buffer(34:57) = lines(48:71)
        end if
        call write1(buffer(1:me))
C
      else
        j = 26
        if (sed) j = 67
        call write1 (header1(1:j))
C
        j = 39
        if (sed) j = 71
        call write1 (lines(1:j))
      end if

      if (sed) then
C                                               find maximum sediment component
        pwsi = wsi * wt
        pwse = wse * wt
        pwss = wss * wt
        twso = 0.

        do n = 1, nps
          twso = twso + wso(n)
        end do

        pwso = twso * wt
        wmax = amax1 (pwso, pwsi, pwss, abs (pwse))
C for graph scaling:
        gpwso = pwso
        sedlab = wlab
        if(pwso .lt. 0.02 .and. wlab .eq. 'kg') then ! adj. units for small sed
          gpwso = pwso*1000.
          sedlab = 'gm'
        else if(pwso .gt. 10000. .and. wlab .eq. 'kg') then
          gpwso = pwso/1000.
          sedlab = 'TN'
        end if
Cgr:
        write(sumsed,'(f8.2,1x,a2)') gpwso, sedlab
        sumsed(13:20) = 'sediment'
        sumsed = adjustl (sumsed)
C
      end if
C                                                     find max volume component
      vmax = amax1 (vin, vr, qbal(6), qbal(7), qbal(9))
C!!
      if(ltyp .eq. 5) vmax = qbal(10)

      call fmt10 (vr, vmax, string, j)
C                                                         write rainfall volume
      buffer(1:25) = '   Rain: '//string//' '//qlab2

      if(.not. chtyp) then
        if (areaR .gt. 0.) then

          v2max = vmax / areaR * conv
          vr2 = vr / areaR * conv

          call fmt10 (vr2, v2max, string, j)

        else

          string = '      ****'

        end if

        buffer(26:39) = ' '//string//' '//qlab1
        m = 39
      else
        v2max = vmax/area1*conv             ! scaling value for channels
      
        m = 25
      end if

      if (sed) then
C                                                         write sediment inflow
        call fmt10 (pwsi, wmax, string, j)
        mp = m+1
        m = mp + 31
        buffer(mp:m) = '               In: '//string//' '//wlab

      end if

      call write1 (buffer(1:m))

C ----------------------------------

      call fmt10 (vin, vmax, string, j)
C                                                           write inflow volume
      buffer(1:25) = ' Inflow: '//string//' '//qlab2

      if(.not. chtyp) then
        if (coarea .gt. 0.) then

          vin2 = vin / coarea * conv

          call fmt10 (vin2, v2max, string, j)

        else

          string = '      ****'

        end if

        buffer(26:39) = ' '//string//' '//qlab1
        buffer(40:40) = aster
        m = 40
      else
        buffer(26:26) = ' '
        m = 26
      end if

      if (sed) then
C                                                              write deposition
        call fmt10 (pwse, wmax, string, j)
        mp = m+1
        m = mp+30
        buffer(mp:m) = '       Deposited: '//string//' '//wlab

      end if

      call write1 (buffer(1:m))

C ------------------------------------ next line:
C                                                      write infiltrated volume

      call fmt10 (vf, vmax, string, j)

      buffer(1:25) = ' Infilt: '//string//' '//qlab2

      if(.not. chtyp) then
        if (area1 .gt. 0.) then
C
          areai = area1
          adjnote = '   '
C        if(ltyp .eq. 1 .or. ltyp .eq. 4) then
C for channels with small botton widths, use effective average wetted
C perimeter to calculate infiltration depth:
C          adjnote = '(a)'
C        end if

          vf2 = vf / areai * conv

          call fmt10 (vf2, v2max, string, j)

        else

          string = '      ****'

        end if
        m = 42
        buffer(26:m) = ' '//string//' '//qlab1//adjnote

      else
        m = 28
        buffer(26:m) = '   '
      end if

      if (sed) then
C                                                      write suspended sediment
        call fmt10 (pwss, wmax, string, j)
        mp = m+1
        m = mp+28
        buffer(mp:m) = '     Suspended: '//string//' '//wlab

      end if

      call write1 (buffer(1:m))

C ------------------------------------------- next line:       write storage
C                                                          
      call fmt10 (qbal(9), vmax, string, j)

      buffer(1:25) = ' Stored: '//string//' '//qlab2

      if(.not. chtyp) then

        if (area1 .gt. 0.) then

          vs2 = qbal(9) / area1 * conv

          call fmt10 (vs2, v2max, string, j)

        else

          string = '      ****'

        end if

        m = 39
        buffer(26:m) = ' '//string//' '//qlab1
      else
        m = 25

      end if

      if (sed) then
C
        sumtab(ltab)%sedout = pwso
        wms = wmax
C                                                        write sediment outflow
        call fmt10 (pwso, wms, string, j)
        mp = m+1
        m = mp+31
        buffer(mp:m) = '              Out: '//string//' '//wlab

      end if

      call write1 (buffer(1:m))

C -----------------------------------next line:         write outflow volume
C                                                 
      call fmt10 (qbal(10), vmax, string, j)

      buffer(1:25) = '    Out: '//string//' '//qlab2

      if (coarea .gt. 0.)  then
        vo2 = qbal(10) / coarea * conv
        call fmt10 (vo2, v2max, string, j)
      else 
        vo2 = 0.
        string = '      ****'
      end if
      buffer(26:39) = ' '//string//' '//qlab1
      buffer(40:40) = aster
        m = 40

      if(chtyp) then

        m = 25

      end if
Cgr:
      write(sumro,'(f8.2,1x,a2)') vo2, qlab1
      sumro(13:18) = 'runoff'
      sumro = adjustl (sumro)

      if (sed) then
C                                                              compute sediment
C                                                             balance error (%)
        supply = wsi

        if (wse .lt. 0.) supply = supply - wse

        if (supply .gt. 1.E-6) then

          errs = (wsi - wse - wss - twso) * 100. / supply

        else

          errs = 0.

        end if

        write (string, 35) errs
        mp = m + 1            ! 41
        m = mp + 30           !71
        buffer (mp:m) = '           Error: '//string//' %'

      end if

      call write1 (buffer(1:m))

C -----------------------------  compute element water balance error (%)

      if (vr + vin .gt. 1.E-4) then
C                            infil               stor      out               
        errq = (vin + vr - qbal(6) - qbal(7) - qbal(9) - qbal(10))
     &         * 100. / (vin + vr)
      else

        errq = 0.0

      end if

      write (string, 35) errq
35    format (f10.2)
      m = 21
      buffer (1:m) = '  Error: '//string//' %'
      if(aster .eq. '*') then
        buffer (22:69) = 
     & '                 [*: based on contributing area]'
        m = 69
      end if

      call write1 (buffer(1:m))

      if(adjnote .eq. '(a)') write(77,'(18x,"[",a3," = adjusted for mean 
     & total infiltrating width ]")')adjnote

      call blank1

      if (cour .and. .not. errout .and. qbal(10) .gt. 0.) then
C                                compute outflow volume based on user time step
C------------------------------------------------------------------------------

        qvdelt = 0.

        do i = 2, limit -1
C
          do j = 1, nchan

            call get (outq(j), i, qi)
            qvdelt = qvdelt + qi

          end do
C
        end do

        qvdelt = 2. * qvdelt

        do j = 1, nchan

          call get (outq(j), 1, qi)

          qvdelt = qvdelt + qi

        end do

        do j = 1, nchan

          call get (outq(j), limit, qi)

          qvdelt = qvdelt + qi

        end do

        qvdelt = delt * qvdelt / 2.
C                                                           check if difference
C                                                            is greater than 1%
        if (abs ((qvdelt - qbal(10)) / qbal(10)) .gt. .01)
     &  errout = .true.

      end if

      ltab = ltov
      if (nprint .eq. 1) then

        return

      end if
C=======================================================================

      if (nprint .eq. 2) then
C                                                 printable (< 80 columns) file
        j = 54
        if (sed) j = 74

        call write1 (header2(1:j))

        m = 54
        buffer(1:m) = '           min         '//qlab4//
     &                 '        '// qlab4//'      '//qlab3

        if (sed) then

          m = 74
          buffer(55:m) = '                '//slab3

        end if

        call write1 (buffer(1:m))

        call blank1

      else if (nprint .eq. 3) then
C                                          spreadsheet (comma/quote delimiters)
        m = 55
        buffer(1:m) =
     &  '    "Time",    "Rainfall",     "Outflow",     "Outflow"'

        if (sed) then

          m = 70
          buffer(56:m) = ',   "Total Sed"'
C                                                                particle sizes
          do n = 1, nps

            l = m + 1
            m = m + 15

            call fmt10 (dp(n), dref, string, i)

            buffer(l:m) = ',"'//string//qlab1//'"'

          end do
        end if

        if (ltyp .eq. 4) then
C                                         labels for main and overbank sections
          l = m + 1
          m = m + 30

          buffer(l:m) = ',        "Main",        "Main"'

          if (sed) then

            l = m + 1
            m = m + 15

            buffer(l:m) = ',   "Total Sed"'
C                                                                particle sizes
            do n = 1, nps

              l = m + 1
              m = m + 15

              call fmt10 (dp(n), dref, string, i)

              buffer(l:m) = ',"'//string//qlab1//'"'

            end do
          end if

          l = m + 1
          m = m + 30

          buffer(l:m) = ',    "Overbank",    "Overbank"'

          if (sed) then

            l = m + 1
            m = m + 15

            buffer(l:m) = ',   "Total Sed"'
C                                                                particle sizes
            do n = 1, nps

              l = m + 1
              m = m + 15

              call fmt10 (dp(n), dref, string, i)

              buffer(l:m) = ',"'//string//qlab1//'"'

            end do
          end if
        end if
C                                                          write labels to file
        call write2 (buffer(1:m))

        m = 55
        buffer(1:m) = '   "(min)",     "('//qlab4//')",     "('//qlab4//
     &                ')",   "('//qlab3//')"'

        if (sed) then
C                                                    write labels for sed units
          m = 70
          buffer(56:m) = ',      "('//slab3//')"'

          do n = 1, nps

            l = m + 1
            m = m + 15
            buffer(l:m) = ',      "('//slab3//')"'

          end do
        end if

        if (ltyp .eq. 4) then
C                                          units for main and overbank sections
          do j = 1, 2

            l = m + 1
            m = m + 30
            buffer(l:m) = ',     "('//qlab4//')",   "('//qlab3//')"'

            if (sed) then

              l = m + 1
              m = m + 15
              buffer(l:m) = ',      "('//slab3//')"'

              do n = 1, nps

                l = m + 1
                m = m + 15
                buffer(l:m) = ',      "('//slab3//')"'

              end do
            end if
          end do
        end if

        call write2 (buffer(1:m))

      end if
C                                                             ...loop over time
C------------------------------------------------------------------------------

      do i = 1, limit

        qi = 0.
C                                                           obtain discharge(s)
        do j = 1, nchan

          call get (outq(j), i, q(j))

          qi = qi + q(j)

        end do

        if (sed) then

          tqso = 0.

          do n = 1, nps

            do j = 1, nchan

              call get (outc(n,j), i, cs)

              qso(n,j) = cs * q(j) * rho(n) * wt   ! this would be mass/time
              tqso = tqso + qso(n,j)

            end do
          end do
C          sd(i) = tqso  ! for graphing
        end if
C                                                                       time...
C------------------------------------------------------------------------------

        time = float (i - 1) * dt
C        tgr(i) = time
C
        write (string, fmtstr) time

        if (nprint .eq. 2) then
C                                                              80-column format
          buffer(1:14) = '    '//string

        else if (nprint .eq. 3) then
C                                                            spreadsheet format
          buffer(1:10) = string

        end if
C                                                 get rainfalls (mm/hr, in/hr)...
C------------------------------------------------------------------------------

        qr = 0.
        rf = 0.

        if (outr .gt. 0) then

          call get (outr, i, qr)

          if (conar .gt. 0.) then
            rf = qr * conar
          end if
C          rgr(i) = rf

        end if

        call fmt10 (rf, rfmax, string, j)

        if (nprint .eq. 2) then
C                                                              80-column format
          buffer(15:28) = '    '//string

        else if (nprint .eq. 3) then
C                                                            spreadsheet format
          buffer(11:25) = ',    '//string

        end if
C                                                   discharge (mm/hr, in/hr)...
C------------------------------------------------------------------------------
        if (conar .gt. 0.) then

          qa = qi * conar
C          qd(i) = qa                                    ! for graphing
          call fmt10 (qa, qp2, string, j)

        else

          string = '**********'

        end if

        if (nprint .eq. 2) then
C                                                              80-column format
          buffer(29:41) = '   '//string

        else if (nprint .eq. 3) then
C                                                            spreadsheet format
          buffer(26:40) = ',    '//string

        end if
C                                                discharge (cu.m/s, cu.ft/s)...
C------------------------------------------------------------------------------

        call fmt10 (qi, qp, string, j)

        if (nprint .eq. 2) then
C                                                              80-column format
          buffer(42:54) = '   '//string

        else if (nprint .eq. 3) then
C                                                            spreadsheet format
          buffer(41:55) = ',    '//string

        end if

        if (sed) then
C                                                         sediment discharge...
C------------------------------------------------------------------------------

          call fmt10 (tqso, qsp2, string, j)
C                                                                total sediment
          if (nprint .eq. 2) then
C                                                              80-column format
            buffer(55:74) = '          '//string

          else if (nprint .eq. 3) then
C                                                         spreadsheet format --
C                                                   discharge by particle class
            m = 70
            buffer(56:m) = ',    '//string

            do n = 1, nps

              l = m + 1
              m = m + 15
              dqso = qso(n,1)
              if (ltyp .eq. 4) dqso = dqso + qso(n,2)

              call fmt10 (dqso, qsp2, string, j)

              buffer(l:m) = ',    '//string

            end do
          end if
        end if

        if (nprint .eq. 2) then
C                                                       write to 80-column file
          call write1 (buffer(1:m))

        else if (nprint .eq. 3) then

          if (ltyp .eq. 4) then
C                                     columns for main and overbank sections...
C------------------------------------------------------------------------------

            do ii = 1, 2

              if (conar .gt. 0.) then

                qa = q(ii) * conar
C                                                      discharge (mm/hr, in/hr)
                call fmt10 (qa, qp2, string, j)

              else

                string = '**********'

              end if

              l = m + 1
              m = m + 15
              buffer(l:m) = ',    '//string
C                                                   discharge (cu.m/s, cu.ft/s)
              call fmt10 (q(ii), qp, string, j)

              l = m + 1
              m = m + 15
              buffer(l:m) = ',    '//string

              if (sed) then
C                                                      total sediment discharge
                tqso = 0.

                do n = 1, nps

                  tqso = tqso + qso(n,ii)

                end do

                call fmt10 (tqso, qsp2, string, j)

                l = m + 1
                m = m + 15
C                                                             by particle class
                buffer(l:m) = ',    '//string

                do n = 1, nps

                  l = m + 1
                  m = m + 15

                  call fmt10 (qso(n,ii), qsp2, string, j)

                  buffer(l:m) = ',    '//string

                end do
              end if
            end do
          end if
C                                                     write to spreadsheet file
          call write2 (buffer(1:m))

        end if
      end do

      close (file2)

      call blank1

      return

C------------------------------------------------------------------------------

      entry qwrt (idl, k, idstr, qpk, tpk, vi, area,
     &            storq, typ, storr, qrmax)

       
      if(k .eq. 0) then     ! report filling of upper layer
C
        if(storr .eq. 1) then
          fill2 = .true.
          filtv = qpk*conv
          filtm = qrmax/60.

        else if(storr .ge. 2) then 
          fillov = .true.  ! overbank fill
          filtvo = qpk*conv
          filtmo = qrmax/60.
C
        end if
C
      else if(k .lt. 0) then  ! report loss to lower layer
C
        lowfin(storr) = .true.
        floss(storr) = qpk*conv
C
      else
        id = idl
        idbuff = idstr     ! string(16) for id no and type
        qp = qpk
        tp = tpk
        vin = vi
        coarea = area

        if (coarea .gt. 0.) then

          conar = conv * 3600. / coarea
          rfmax = qrmax * conar

        else

          conar = 0.
          rfmax = 0.

        end if

        do i = 1, k

          outq(i) = storq(i)

        end do

        if(typ .eq. 6) ltyp = typ
        outr = storr
      end if

      return

C------------------------------------------------------------------------------


      entry swrt ( sumwso, sumwsi, sumwse, sumwss, qpk, tpk, storc)


      wsi = sumwsi
      wse = sumwse
      wss = sumwss
      qsp = qpk
      tsp = tpk

      do n = 1, nps

        wso(n) = sumwso(n)
C                                                        keep track of injected
C                                                                      sediment
        if (ltyp .eq. 5) wsj(n) = wsj(n) + wso(n)    ! wso already multiplied by rho

        do i = 1, nchan

          outc(n,i) = storc(n,i)

        end do

      end do

      sed = .true.

      return

C------------------------------------------------------------------------------
C------------------------------------------------------------------------------

      entry event (sumbal, dtm)   !

      call blanks
C                                                check for unconnected elements
C------------------------------------------------------------------------------

      call orphan (id, norphs, unconn)

      if (norphs .gt. 0) then
C                                                                        report
        buffer(1:41) = ' Input error - unconnected element: ID = '
C**  shouldn't go to (norphs+1), since unconn will then be undefined  4/2003
        do i = 1, norphs     !+ 1
C                                                     don't report last element
          if (unconn(i) .ne. id) then

            write (buffer(41:46), '(i6)') unconn(i)

            call writes (buffer(1:46))

            call blanks

          end if

        end do

        return

      end if
C                                       write volumes to output file and screen
C------------------------------------------------------------------------------
      if(nele .gt. 0 .and. nele. ne. nelt) then
        write (file1,'(" Number of elements processed [",i4,"] does not 
     &equal NELE [",i4,"] in input file. ")') nelt, nele
      else
        write(file1,'(i4," elements Processed"/)') nelt
      end if

C                                                                watershed area
      coarea = sumbal(1)

      call writes (' Event Volume Summary:')

      call blanks
C                                                       compute total volume in
      qin = 0.

      do i = 2, 5

        qin = qin + sumbal(i)

      end do
C                                                      compute total volume out
      qout = 0.

      do i = 6, 10

        qout = qout + sumbal(i)

      end do

      flag = 0
      qmax = amax1 (qin, qout)

      if (coarea .gt. 0.) qa = qin / coarea * conv
C                                                          write nonzero volume
C                                                                    components
      do l = 2, 10

        gv = sumbal(l)

        if (gv .gt. 0.) then
C                                                   compute volume / basin area
          if (coarea .gt. 0.) then

            gv2 = gv / coarea * conv

            call fmt10 (gv2, qa, string, len)

          else

            string = '**********'

          end if

          j = l - 1
          buffer(1:33) = label(j)//'  '//string

          if (flag .eq. 0) then
C                                                        units (first row only)
            buffer(34:41) = ' '//qlab1//'     '

          else

            buffer(34:41) = '        '

          end if
C                                                         volume in cubic units
          call fmt10 (gv, qmax, string, len)

          buffer(42:51) = string
          m = 51

          if (flag .eq. 0) then
C                                                        units (first row only)
            buffer(52:57) = ' '//qlab2
            m = 57
            flag = 1

          end if

          call writes (buffer(1:m))

        end if

      end do
C                                                     compute vol balance error
      if (qin .gt. 0.) then

        err = 100. * abs (qin - qout) / qin

        if (err .lt. 1.) then

          ierr = 0

        else

          ierr = nint (err)

        end if

      else

        ierr = -1

      end if

      if (ierr .ge. 0) then
C                                                           write balance error
C------------------------------------------------------------------------------

        call blanks

        buffer(1:41) = ' Error (Volume in - Volume out - Storage)'
        m = 53

        if (ierr .lt. 1) then

          buffer(42:53) = ' < 1 percent'

        else if (ierr .lt. 5) then

          write (c1, '(i1)') ierr
          buffer(42:53) = ' = '//c1//' percent'

        else if (ierr .lt. 10) then

          write (c1, '(i1)') ierr
          buffer(42:54) = ' = '//c1//' percent!'
          m = 54

        else if (ierr .lt. 100) then

          write (c2, '(i2)') ierr
          buffer(42:56) = ' = '//c2//' percent!!'
          m = 56

        else

          if (ierr .ge. 1000) ierr = 999
          write (c3, '(i3)') ierr
          buffer(42:58) = ' = '//c3//' percent!!!'
          m = 58

        end if

        call writes (buffer(1:m))

      end if

      call blanks
C                                                          write time step info
C------------------------------------------------------------------------------

      if (cour) then

        call writes
     &  (' Time step was adjusted to meet Courant condition')

        call blanks

        if (errout) then
          call writes (' Time step is too large: element hydrographs are
     & not well-represented!')

          call blanks

        end if

      else
C                                                    rank minimum dt's for each
C                                              quarter from smallest to largest
        do i = 1, 4
          do j = 1, 3

            if (dtm(j+1) .lt. dtm(j)) then
              temp = dtm(j)
              dtm(j) = dtm(j+1)
              dtm(j+1) = temp
            end if

          end do
        end do

        if (dtm(1) .lt. delt) then

          buffer(1:39) = ' Time step distribution (100,75,50%) = '
          m = 40

          do i = 1, 3

            dtm(i) = dtm(i) / 60.

            call fmt10 (dtm(i), delt, string, j)

            n = m + 12 - j
            buffer(m:n) = string(j:10)//', '
            m = n + 1

          end do

          n = m + 1
          buffer(m-2:n) = ' min'

          call writes (buffer(1:n))

          call blanks

        end if

      end if
C                                             compute area in hectares or acres
C------------------------------------------------------------------------------

      coarea = coarea / arconv

      call fmt10 (coarea, coarea, string, j)

      m = 38 - j
      buffer(1:m) = ' Total watershed area = '//
     &               string(j:10)//' '//alab

      call writes (buffer(1:m))

      call blanks

      if (sed) then
C                                               compute total injected sediment
C------------------------------------------------------------------------------

        twsj = 0.

        do i = 1, nps

          twsj = twsj + wsj(i)

        end do

        if (twsj .gt. 1.E-7) then

          twsj = twsj * wt

          if (twsj .le. wconv) then

            wconv = 1.
            flag = 0

          else
C                                                               convert to tons
            twsj = twsj / wconv
            flag = 1

          end if

          call fmt10 (twsj, twsj, string, j)

          m = 33 - j

          buffer(1:m) = ' Injected sediment = '//string(j:10)//' '

          if (flag .eq. 0) then
C                                                                 primary units
            n = m + 2
            buffer(m+1:n) = wlab

          else
C                                                                          tons
            n = m + 4
            buffer(m+1:n) = 'tons'

          end if

          call writes (buffer(1:n))

          call blanks

        end if

        if (coarea .gt. 0.) then
C                                                        compute total sediment
C------------------------------------------------------------------------------

          twso = 0.

          do i = 1, nps

            twso = twso + wso(i)

          end do

          if (twso .lt. 1.E-7) return
C                                                     convert to kg/ha or lb/ac
          twso = twso * wt / coarea

          if (twso .le. wconv) then

            wconv = 1.
            flag = 0

          else
C                                                 convert to tons/ha or tons/ac
            twso = twso / wconv
            flag = 1

          end if

          call fmt10 (twso, twso, string, j)

          m = 30 - j
          buffer(1:m) = ' Sediment yield = '//string(j:10)//' '

          if (flag .eq. 0) then
C                                                                 primary units
            n = m + 6
            buffer(m+1:n) = slab1

          else
C                                                                          tons
            n = m + 8
            buffer(m+1:n) = slab2

          end if

          call writes (buffer(1:n))

          if (nps .gt. 1) then
C                                      compute sediment yield by particle class
C------------------------------------------------------------------------------

            do j = 1, nps
C                                                          convert to wt / area
              qso(j,1) = wso(j) * wt / coarea / wconv
C                                                            compute percentage
C                                                                      of total
              pso(j) = qso(j,1) / twso * 100.

            end do

            call blanks

            call writes (' Sediment yield by particle class:')

            call blanks

            buffer(1:19) = ' Particle size ('//qlab1//')'
            m = 19
C                                                                particle sizes
            do i = 1, nps

              l = m + 1
              m = m + 12

              call fmt10 (dp(i), dref, string, j)

              buffer(l:m) = '  '//string

            end do

            call writes (buffer(1:m))

            if (flag .eq. 0) then

              buffer(1:19) = '      Yield ('//slab1//')'

            else

              buffer(1:19) = '    Yield ('//slab2//')'

            end if

            m = 19
C                                              sediment yield by particle class
            do i = 1, nps

              l = m + 1
              m = m + 12

              call fmt10 (qso(i,1), twso, string, j)

              buffer(l:m) = '  '//string

            end do

            call writes (buffer(1:m))
C                                                     percent by particle class
C            write (*, 75) (pso(i), i = 1, nps)
            write (file1, 75) (pso(i), i = 1, nps)
75          format ('   % of total yield', 5(6x, f6.2))

          end if
        end if
      end if
C - - - -- - - - - - - - - - - - - -  optional table of element data
      if(tabl) then
C
        call blanks
        call blanks
        write(file1,'(20x," Tabular Summary of Element Hydrologic Compon
     &ents"/)')
C! this complete table writing section has been redesigned  12/02  RES
        head(1)= '   ID  Element'
        head(2)= '        Type  '
        head(3) = '              '
        aremx = 0.
        tarmx = 0.
        rfamx = 0.
        oqmax = 0.
        tsdmx = 0.
        finmx = 0.
        winj = .false.
        twol = .false.
C                              preliminary scan of tab stored values
        do jt = 1,ltab
          if(.not. twol .and. sumtab(jt)%twola) then
            twol = .true.
          end if
          if(sumtab(jt)%are .gt. aremx) aremx = sumtab(jt)%are
          if(sumtab(jt)%cumare .gt. tarmx) tarmx = sumtab(jt)%cumare
          if(sumtab(jt)%volrn .gt. rfamx) rfamx = sumtab(jt)%volrn
          if(sumtab(jt)%volro .gt. oqmax) oqmax = sumtab(jt)%volro
          if(sumtab(jt)%ftot .gt. finmx) finmx = sumtab(jt)%ftot
          if(sed) then
            if(sumtab(jt)%sedout .gt. tsdmx) tsdmx = sumtab(jt)%sedout
          end if
          if(sumtab(jt)%itype .eq. 5) winj = .true.
        end do
C get number of digits necessary
        nma = fintgs(tarmx)
        nme = fintgs(aremx)
        nmr = fintgs(rfamx)
        nmq = fintgs(oqmax)
        nmf = fintgs(finmx)
        if(sed) nms = fintgs(tsdmx)
C    
        if(units .eq. 1) then  ! metric
          ndp = 2
          ndq = 3
        else
          ndp = 1
          ndq = 2
        end if
C    
        nco = 8
        lfd = 14
C    
        do jc=1,8
          do Ln=1,3
            write(fmhed(jc)(6:6),'(i1)') je(jc)
          end do
          minl(jc) = je(jc) 
C    
          select case (jc)
          case(1)                    ! element area
            lneed = nme + 1 + ndp   !
            hedstr(1,1) = 'Element  '
            hedstr(2,1) = ' Area    '
            Hedstr(3,1) = '  '//arlab//'   '
            write(fmdat(1)(9:9),'(I1)') ndp
          case(2)                    !  cum area
            lneed = nma + 1 + ndp
            write(fmdat(2)(9:9),'(I1)') ndp
            hedstr(1,2) = 'Cumulated'
            hedstr(2,2) = '  Area   '
            hedstr(3,2) = '  '//arlab//'   '
          case(3,5,7)             ! inflow, outflow, total infil
            lneed = nmq + 1 + ndq
            write(fmdat(jc)(9:9),'(I1)')  ndq
            If(jc==3) hedstr(1,3) = ' Inflow  '
            hedstr(2,jc) = blank9
            If(jc==5) then
              hedstr(1,5) = 'Outflow  '
              ntbinj = lfd+2
              write(fminj(7:8),'(I2)') ntbinj
            end if
            If(jc==7) then
              lneed = nmf + 1 + ndq
              hedstr(1,7) = ' Total  '
              hedstr(2,7) = ' Infil  '
            end if
            hedstr(3,jc) = ' '//qlab2//'   '
          case(4)                              ! Rainfall vol
            lneed = nmr + 1 + ndq
            write(fmdat(jc)(9:9),'(I1)')  ndq
            hedstr(1,4) = 'Rainfall '
            hedstr(2,4) = blank9
            hedstr(3,4) = ' '//qlab2//'   '
          case(6)                              ! peak in mm(in)/hr
            lneed = 4 + ndq
            write(fmdat(jc)(9:9),'(I1)')  ndq
            hedstr(1,6) = ' Peak    '
            hedstr(2,6) = ' Flow    '
            hedstr(3,6) = ' '//qlab4//'     '
          case(8)                              ! initial theta
            lneed = 6
            hedstr(1,8) = 'Initial  '
            hedstr(2,8) = ' Water   '
            hedstr(3,8) = 'Content  '
            fmdat(8)(9:9) = '4'
          end select
C    
          if(lneed .gt. minl(jc)) then
            nx(jc) = lneed - minl(jc) + 1
            minl(jc) = lneed
          else
            nx(jc) = 1
          end if
          lfd = lfd + minl(jc) + 1
C    
          write(fmhed(jc)(2:2),'(I1)') nx(jc)
          write(fmhed(jc)(6:6),'(I1)') je(jc)
          write(fmdat(jc)(6:7),'(I2)') minl(jc)
C    
C          write(77,'(" fm",3i3,2a10)')jc,minl(jc),nx(jc),fmdat(jc)
        end do
C    
        fmthl(1:4) = '(A14'
        fmtdat(1:4) = '(A14'
        nhe = 4
        nde = 4
        do i=1,8
          nhs = nhe + 1
          nhe = nhe + 6
          fmthl(nhs:nhe) = fmhed(i)
          nds = nde + 1
          nde = nde + 9
          fmtdat(nds:nde) = fmdat(i)
        end do
C    
        if(winj) then                 ! prepare format for writing inject
          write(fminj(11:12),'(I2)')minl(5) 
          fminj(15:23) = fmdat(6)
          fminj(24:32) = fmdat(7)
          fminj(33:33) = ')'
        end if
C    
        if(twol) then                           ! add subsoil infil:
          nco = 9
          fmhed(nco)(6:6) = '7' 
          hedstr(1,9) = 'Subsoil  '
          hedstr(2,9) = ' Infil.  '
          hedstr(3,9) = '   '//qlab1//'    '
          nhs = nhe + 1
          nhe = nhe + 6
          fmthl(nhs:nhe) = fmhed(nco)
          fmdat(nco)(6:7) = '07'
          nds = nde + 1
          nde = nde + 9        
          fmtdat(nds:nde) = fmdat(nco)
        end if
C    
        if(sed) then                            ! add sed yield
          nco = nco + 1
          je(nco) = 8
          jc = nco
          minl(jc) = je(jc) + 1
          lneed = nms + 2 + ndq
          if(lneed .gt. minl(jc)) then
            nx(jc) = lneed - minl(jc) + 1
            minl(jc) = lneed
          else
            nx(jc) = 1
          end if
C    
          write(fmhed(jc)(2:2),'(I1)') nx(jc)
          write(fmhed(jc)(6:6),'(I1)') je(jc)
          write(fmdat(jc)(2:2),'(I1)') nx(jc)
          write(fmdat(jc)(6:7),'(I2)') minl(jc)
C    
          hedstr(1,jc) = 'Sediment '
          hedstr(2,jc) = ' Yield   '
          hedstr(3,jc) = '   '//wlab//'    '
          nhs = nhe + 1
          nhe = nhe + 6
          fmthl(nhs:nhe) = fmhed(nco)
          nds = nde + 1
          nde = nde + 9        
          fmtdat(nds:nde) = fmdat(nco)
        end if
        nhs = nhe + 1
        fmthl(nhs:nhs) = ')'
C
        nds = nde + 1
        fmtdat(nds:nds) = ')'
C                                        put headings into outfile
        do i=1,3
          write(file1,trim(fmthl)) head(i),
     &        (hedstr(i,jc)(1:je(jc)),jc=1,nco)
        end do
C        write(77,'(a90)') fmtdat
C                     first go thru list and report planes, injects, and urban els.
        do jp = 1,ltab
          itu = sumtab(jp)%itype
          if(itu .eq. 0 .or. itu .eq. 6 .or. itu .eq. 5) then
            write(charid(1:6),'(I5,1x)') sumtab(jp)%idel
            charid(7:14) = eltype(itu)
            if(itu .eq. 5) then
              write(file1,fminj) charid,sumtab(jp)%volro,
     &           sumtab(jp)%ropeak
            else
              dval(1) = sumtab(jp)%are
              dval(2) = sumtab(jp)%cumare
              dval(3) = sumtab(jp)%volin
              dval(4) = sumtab(jp)%volrn
              dval(5) = sumtab(jp)%volro
              dval(6) = sumtab(jp)%ropeak
              if(itu .eq. 3) then
                dval(7) = 0.
                dval(8) = 0.
              else
                dval(7) = sumtab(jp)%ftot
                dval(8) = sumtab(jp)%thst
              end if
              if(sumtab(jp)%twola) then
                dval(9) =  sumtab(jp)%flowr * conv
              else if (twol) then   ! this element has only one layer
                dval(9) = -0.0
              end if
              if(sed) then
                dval(nco) = sumtab(jp)%sedout
              end if
              write(file1,trim(fmtdat)) charid,(dval(i),i=1,nco)
            end if
          end if
        end do
C     now go thru list for channels, pipes,
        do jp = 1, ltab
          itu = sumtab(jp)%itype
          if(itu .eq. 1 .or. itu .eq. 2 .or. itu .eq. 4 .or. itu .eq. 7)
     &      then
            write(charid(1:6),'(I5,1x)') sumtab(jp)%idel
            charid(7:14) = eltype(itu)
            if(itu .ne. 2) then
              dval(1) = sumtab(jp)%are
            else                ! pipe
              dval(1) = 0.
            end if
            if(itu .eq. 7) then   ! overbank of channel
              do i=2,6
                dval(i) = 0.
              end do
            else
              dval(2) = sumtab(jp)%cumare
              dval(3) = sumtab(jp)%volin
              dval(4) = sumtab(jp)%volrn
              dval(5) = sumtab(jp)%volro
              dval(6) = sumtab(jp)%ropeak
            end if
            if(itu .eq. 2) then   ! pipe case
               dval(7) = 0.
               dval(8) = 0.
            else
              dval(7) = sumtab(jp)%ftot
              dval(8) = sumtab(jp)%thst
            end if
            if(sumtab(jp)%twola) then  !  pipe would not pass this test
              dval(9) = sumtab(jp)%flowr * conv
            else if( twol) then
              dval(9) = -0.
            end if
            if (sed) then
              dval(nco) = sumtab(jp)%sedout
            end if
            write(file1,trim(fmtdat)) charid,(dval(m),m=1,nco)
          end if
        end do
C     lastly look for ponds
        do jp = 1, ltab
          itu = sumtab(jp)%itype
          if(itu .eq. 3 ) then
            write(charid(1:6),'(I5,1x)') sumtab(jp)%idel
            charid(7:14) = eltype(itu)
            dval(1) = sumtab(jp)%are
            dval(2) = sumtab(jp)%cumare
            dval(3) = sumtab(jp)%volin
            dval(4) = sumtab(jp)%volrn
            dval(5) = sumtab(jp)%volro
            dval(6) = sumtab(jp)%ropeak
            dval(7) = sumtab(jp)%ftot
            dval(8) = -0.0
            if(twol) then
              dval(9) = -0.0
            end if
            if(sed) then
              dval(nco) = sumtab(jp)%sedout
            end if
            write(file1,trim(fmtdat)) charid,(dval(m),m=1,nco)
          end if
        end do
C
      end if
C
C
      return

C------------------------------------------------------------------------------


      entry wrt00 (fileo) !, lim, del, np, dpd, rh)


C                                             save arguments as local variables
      file1 = fileo
C 
      dt = delt / 60.    ! time step in minutes

C      nps = np

      if (units .eq. 1) then
C                                                                  metric units
!        conv = 1000.
!        wt = 1000.           ! kg/m3 water density
!        arconv = 10000.
!        wconv = 1000.
        wlab = 'kg'
        alab = 'ha'
C!!
        arlab = ' m^2'
        qlab1 = 'mm'
        qlab2 = 'cu m '
        qlab3 = 'cu m /s'
        qlab4 = 'mm/hr'
        slab1 = 'kg/ha'
        slab2 = 'tons/ha'
        slab3 = 'kg/s'
C!!
        fmttab = '(1x,f10.2)'
C  table output option units column
        colunits = '                     m^2       m^2        m^3      m 
     &^3        m^3      mm/h     m^3     Content '!   mm      mm      '
C
      else
C                                                                 english units
!        conv = 12.
!        wt = 62.4    
!        arconv = 43560.
!        wconv = 2000.
        wlab = 'lb'
        alab = 'ac'
C!!  add new label
        arlab = 'ft^2'
        qlab1 = 'in'
        qlab2 = 'cu ft'
        qlab3 = 'cu ft/s'
        qlab4 = 'in/hr'
        slab1 = 'lb/ac'
        slab2 = 'tons/ac'
        slab3 = 'lb/s'
C!!
        fmttab = '(1x,f10.1)'
        colunits = '                    ft^2      ft^2       ft^3      f 
     &t^3       ft^3     in/h     ft^3    Content '!  in      in  '
C
      end if
C
      if(sed) then
C                                               convert particle diams to mm or
       dmin = 1.E6
C                                                   inches & note smallest size
       do j = 1, nps

         wsj(j) = 0.
         dp(j) = dpd(j) * conv
         if (dp(j) .lt. dmin) dmin = dp(j)

       end do

       dref = dmin * 1.E6
C
      end if
C                                                 determine time format from dt
      if (dt .ge. 10.) then
        fmtstr = '(f10.0)'
C        assign 105 to ifmt

      else if (dt .ge. 1.) then
        fmtstr = '(f10.1)'
C        assign 115 to ifmt

      else
        fmtstr = '(f10.2)'
C        assign 125 to ifmt

      end if

C105   format (f10.0)
C115   format (f10.1)
C125   format (f10.2)

      return

      end

C
C======================================================================
      integer function fintgs(vmax)
C finds digits necessary for printing based on max value
      integer :: m
      real :: vmax
      if(vmax .gt. 0.) then
        m = int(log10(vmax)) + 1
      else
        m = 1
      end if
      fintgs = m
      return
      end
C------------------------------------------------------------------------------


      subroutine blank1

C   Write a blank line to the output file.

      call write1(' ')

      return

      end

C------------------------------------------------------------------------------


      subroutine blank2

C   Write a blank line to the spreadsheet file.


      integer file1, file2

      common /writ/ file1, file2


      write (file2, 5)
5     format (' ')

      return

      end

C------------------------------------------------------------------------------


      subroutine blanks

C   Write a blank line to the output file & screen.


      integer file1, file2

      common /writ/ file1, file2


      write (*, 5)
      write (file1, 5)
5     format (' ')

      return

      end

C------------------------------------------------------------------------------


      subroutine write1 (string)

C   Write "string" to the output file.


      integer file1, file2
      logical opn
      character string*(*)

      common /writ/ file1, file2

      inquire(file1, opened = opn)
      if(opn) write (file1, 5) string
5     format (a)

      return

      end

C------------------------------------------------------------------------------


      subroutine write2 (string)

C   Write "string" to the spreadsheet file.


      integer file1, file2

      character string*(*)

      common /writ/ file1, file2


      write (file2, 5) string
5     format (a)

      return

      end

C------------------------------------------------------------------------------


      subroutine writes (string)

C   Write "string" to the output file & screen.


      integer file1, file2

      character string*(*)

      common /writ/ file1, file2


      write (*, 5) string
      write (file1, 5) string
5     format (a)

      return

      end

