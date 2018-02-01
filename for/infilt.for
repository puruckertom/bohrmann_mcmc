C     
C This file contains infiltration related subroutines and functions
C   version uses API if specified
C
C   Code type: Fortran subroutine with additional entry points

C   Compiler: Lahey Fortran 90  This version with modules.

C   Programmed by: R.E. Smith/C. Unkrich

C   Origin Date: 10/93

C   Last modification: 5/05 by RES.  make good for tiny or zero th in Rectr

C   Description:

C     Computes an average infiltration rate for each spatial node
C     using a three-parameter model developed by Parlange et al
C     (1982). Modifed for treating general rain patterns, surface
C     head effect, and wet initial conditions, 6/92. Modified for a
C     two-layer soil profile, 10/93. A value of zero thickness for
C     the upper layer indicates a single layer. A value of zero for
C     the lower capillary parameter signals a constant rate through
C     the lower layer when water reaches the interface. Revised again
C     to handle parallel processing for compound channel sections.
C     Revised 10/01 to correct interactions between effects of CV, 
C     surface relief, and two layer systems under complex rainfalls.
C     NOTE: 99 is assumed the number of the diagnostic output file, if
C     logical DIAG is .T.
C------------------------------------------------------------------------------

C   Entry points:

C     infil0       obtains parameter values from input block and initializes
C                  variables,
C------------------------------------------------------------------------------

C   Arguments:
C
C     afac(20)     real        in: relative area covered by flowing water
C                              out: relative area to which negative vav is to
C                              be applied.  the rest of the area infiltrates at
C                              rate rf
C
C     rfi          real        rainfall rate
C
C     dtt          real        time step over which computations apply
C
C     h1(20,k)     real        flow depths (m or ft),

C     qu(20,k)     real        potential surface supply rate (m/s, ft/s), 
C                              upstream or lateral

C     dtt          real        time step (s),

C     vav(20,k)    real        dt-average infiltration rate (m/s or ft/s),
C
C     ichan       integer      1 except for overbank
C

C------------------------------------------------------------------------------

C   Input Block (input file labels in parenthesis):

C     sk1 (KS)           real      Ks, upper layer, mm/hr or in/hr,

C     sk2                          Ks, lower layer,

C     g1 (G)             real      capillary drive, upper layer, mm or in,

C     g2                           capillary drive, lower layer,

C     al1 (DI)           real      pore size distribution index, upper layer,

C     al2                          pore size distribution index, lower layer,

C     por1 (PO)          real      porosity, upper layer,

C     por2                         porosity, lower layer,

C     rock1 (RO)         real      volumetric rock fraction, upper layer,

C     rock2                        volumetric rock fraction, lower layer,

C     sint(SA):> thin    real      fraction of surface saturated pore space, initially

C     cv (CV)            real      coefficient of variation of Ks, {0,1},
C                                  [less than 0.1 treated as 0]
C     zfl (THI)          real      thickness of upper soil layer, m or ft,
C
C GLOBAL:
C     units              int       1 = metric,
C                                  2 = english,

C     nchan            int         1 for planes and main channels, 2 for call
C                                  from overbank area
C
C     nk                 int       number of spatial nodes,

C     msg              char*(*)    identifies the current element in case an error
C                                  is encountered,

C     api                real     anticedent precip depth [Optional]: used in an initial 
C                                  surface block of water content=thin and depth zb to 
C                                  make a water content=api, present at beginning of run.  

C     diagn              log       .true. indicates diagnostic output is desired,
C------------------------------------------------------------------------------

C   Subroutines/Functions:

C     errvf        (infilt.for)
C     errly             "
C     rectr             "
C     setg              "
C     th2oth1           "
C     th1oth2           "
C     rkohy             "
C     fpond             "
C     funk              "
C     flams             "
C     fgex              "
C     favdI             "    change in storage in dt at rate fc
C     fIsots            "    I*(t*)
C     vcurf             "
C     vfm               "
C     fcapv             "
C     iter         (miscel.for)
C     getr4        (reader.for)
C     qprw         (writer.   )  send information on soil sat. to writer
C     errxit       (K2RunW.for)
C - - - - - - - - - - - - - - - - - - - - - - - -- -- - - - - - - - --- - - -
C
C     Working Variables
C
C     thiu    effective current initial water content at pulse front
C     fi()    infiltrated water above fmin counted as prewet depth
C     fs()    overall infiltrated depth
C     fv()    infiltrated depth of pulse during r > K
C     fr()    infiltrated water subject to redistribution during hiatus
C     f2n     infiltrated depth at beginning of time step
C     f2b     infiltrated depth of main pulse, end of dtt
C     f2p     infiltrated depth of rewet pulse
C     ffl     infiltrated depth above soil interface
C     cumcm   max. capacity of upper layer:  cumcm >= ffl
C     gu      effective G depending on location of wetting front.
C     fkp     net infiltration due to flux of initial water
C     qki     relative flux of initial pulse during rewet
C     ref     effective surface supply rate including surface depth
C     erfi    rainrate equivalent to maintain initial soil water
C     rhc     coverage factor, based on % area covered by flowing water
C               -related to afac
C     tho()    water content at surface
C     zp()     depth of rewet pulse
C     zb()     depth of main wetting pulse
C     vmin     aerial effective final f for infilt at surface
C     vmu      value of vmin used for base or rewet calculation
C     twolay   logical .T. when there are two soil layers
C     subcon   logical .T. when lower layer has Ks2 .le. 2.*Ks1
C     lowet    logical .T. when wetting reaches layer interface
C     topfil   logical .T. when upper layer is saturated.

C------------------------------------------------------------------------------

      subroutine infilt (afac, rfi, h1, qu, dtt, vav, ichan )

      use multip          ! this module for LaheyF version only
      use itpars
      use elpars
      use runpars
C                                                                     arguments
C------------------------------------------------------------------------------

      logical :: diagn, op

      integer ::  nx, ichan  !, k
      real, dimension(20) :: afac
      real, dimension(20,nchan) :: h1, vav, qu

C                                                               local variables
C------------------------------------------------------------------------------
C
      character msg*(*)
      character*6 trace, source, fmode, set*2
      character(LEN=11) :: notice
      Character(LEN=7) :: calldat

      logical :: twolay, ifl, topfil(20,2), subcon(2), lowet(20,2)
     &, notify(2), dbg  !, diag

      integer :: ierr, nk, i, j, niter, idu(2)  

      dimension fs(20,2), fv(20,2), thiu(20,2), tho(20,2), fr(20,2),
     &          fi(20,2), fkp(20,2), qki(20,2), fcr(20,2), zp(20,2),
     &          zb(20,2), gu(2), cumcm(2), thin(2), ziu(20), 
     &          def1(2), def2(2), por1(2), por2(2), t(2)   !vlim(2), 

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j
     
C-----------------------------------------------------------------------

      external errvf, errly
C
C-----------------------------------------------------------------------------
      
      save  fs, fv, thiu, tho, fr, fi, fkp, qki, fcr, zp, zb, 
     &      gu, cumcm, thin, ziu, def1, def2, por1, por2, t, 
     &      nk, jd, notice, sumrf, sumif, topfil, subcon, lowet, 
     &      notify, niter, po, idu, defl, rxv
C------------------------------------------------------------------------------
      i = ichan
C   
      t(i) = t(i) + dtt
      trp = t(i)/60.
C      
      trace = 'infilt'
C
        do j = 2, nk
C                                                loop thru 2 -> nk locations...
C------------------------------------------------------------------------------
          rfj = rfi
          rfdt = rfj*dtt
          ht = h1(j,i)
          if(j .eq. jd) sumrf = sumrf + rfdt
          dtu = dtt

          if (sk1(i) .le. 0.) then
C                                           impermeable soil option
            vav(j,i) = 0.
            afac(j) = 1.0
            fmode = 'imperv'
            rxv = 0.           ! in case of diag PO

          else if (g1(i) .eq. 0.) then    ! fixed loss option
C   note this option does not now recognise two-layer systems (could be added)
C                                                   
            refs = rfj + qu(j,i)            ! constant rate, but
            rst = rfj / sk1(i)              ! varies with rst when CV > 0.1
            vfmf = 1.0
            if(cv(i) .gt. 0.05) vfmf = vfm(rst,cv(i))
            vmin = sk1(i) * vfmf
            if(vmin .lt. rfj) afac(j) = 1.0
            rhc = afac(j)
C
            vav(j,i) = 0.
            if(refs .gt. 1.e-8 .or. ht .gt. 1.e-6)
     &         vav(j,i) = min(refs/afac(j), vmin)   ! adjust for microtopog during recess
            fmode = 'fix_fc'
            fs(j,i) = fs(j,i) + vav(j,i)*dtu !**
C
C  normal option
          else
            sku = sk1(i)
            defs = def1(i)
C
            if(lowet(j,i)) then
              fup = ffl(j,i)
C 
            else
              fup = 0.
            end if
C
            if(subcon(i)) then  ! subsurface control condition:
              if(lowet(j,i)) then
                sku(j) = sk2(i)
                if((ffl(j,i) .ge. cumcm(i)) 
     &              .and. .not. topfil(j,i)) then
                  topfil(j,i) = .true.
                  tho(j,i) = ths1(i)
C  
                  if(j .eq. nk .and. .not. notify(i)) then
                    dum = 1.0
                    izr = 0
                    vfill = cumcm(i) - apa
                    call qwrt (idl, izr, trace, vfill, dum, dum,
     &       dum, idu, 0, i, t(i))  ! pass cumcm and t to writer
C                    if(diag .and. j .eq. jd) write(99,*) t(i)
                    notify(i) = .true.
                  end if  
                end if
C!! new option for sealed lower bound
                if(sk2(i) .le. 0.) then !**
                  if(.not. topfil(j,i)) then !**
                    ffl(j,i) = ffl(j,i) + rfdt + ht !**
                    ffl(j,i) = min(ffl(j,i),cumcm(i)+.001) !**
                  end if !**
                  vav(j,i) = 0. !**
                  afac(j) = 1.0 !**
                  fmode = 'imperv' !**
                  rxv = 0.           ! in case of diag PO
                  go to 100 !**
                end if !**
              else 
                sku(j) = sk1(i)
              end if
            else
C                     surface layer control conditions
              if (fs(j,i) .gt. cumcm(i)) then   ! was fv and fcr  12/99
C  this assumes twolay, since otherwise cumcm is too large          
                sku(j) = vmin2(i)
                defs = def2(i)
              end if
            end if
            rst = rfj / sku(j)
C
            ref = rfj
C                                         capillary drive includes surface head
            apu = defs * (ht + gu(i))
            if( apu .le. 0.) then
              iunit = 77              ! output file should always be opened
              inquire (99, opened = op)
              if(op) iunit = 99
              write(iunit,"(' APU negative in INFILT: '/,i3,(7g12.4))") 
     &           j, g1(i), zfl(i), trp, ht, gu(i), thiu(j,i), defs
              stop    ' negative value of cap. drive in INFILT for this
     &element '
            end if
C  
            rhc = 1.0
C                        factor applicable to spatial distrib. of Ksat:
CCC            if(topfil(j,i)) afac(j) = 1.0   ! new 6/02
            rhv = 0.1 + 0.9*afac(j)   !assumes afac always represents relative coverage
            vfmf = rhv + (1.-rhv) * vfm(rst,cv(i))
C 
            vmin = sku(j) * vfmf      ! new method. sku is found from setg
C                        when surface is covered with water, Keff = Kmean
            if((qu(j,i)+ht) .gt. 0.) then
C!!  the following is new criterea [10/01]
              famx = 10.
              if(fs(j,i) .le. 1.e-5) then
                famx = vmin*2.*(1.+sqrt(0.5*apu/dtu/vmin))
              end if
              viu = fs(j,i) - fup
              if(zp(j,i) .gt. 0.) viu = fv(j,i) - fup
              fap = favdI(apu, vmin, viu, dtu)/dtu                     ! favdI
              if(lowet(j,i)) fap = fap + (cumcm(i) - ffl(j,i))/dtu    ! 10/02
              if(rfj .le. fap ) then !  .and. rfj .le. vav(j,i)
   ! treat surface depth as part of potential:
                ref = rfj + ht/dtu
                if(qu(j,i) .le. 0. .and. ht .gt. 0.) ref = ref + ht/dtu
                ref = min(ref,1.01*fap,famx)
                rhc = afac(j) ! reduction when infil from runon
C                                                flow area or during recession
              else
               afac(j) = 1.
              end if
C              
            else
C 
              fap = -0.
              afac(j) = 1.0  ! this is returned for use with ql
            end if
C
            if(diag .and. j .eq. jd) write(99,877) t(i)/60., ref*po,
     &  fap*po, vmin*po, ziu(j), zb(j,i), fr(j,i), ht !, zb(j,i)
 877  format(/'   time        ref         fpot       vmin        ziu    
     &      zb          fr           ht'/9g12.4)
C
C   + + + + + + + + + +  +  +  +  +  +  +   +   +   +   +    +    +    +
            if (ref .lt. vmin ) then
C                                                  ponding NOT possible
              vav(j,i) = ref
              refdt = ref*dtu
              if (fv(j,i) .gt. 0.) then
C                                                 begining of redistribution --
C                                                  low rainfall and positive fv
                fr(j,i) = fv(j,i) + fi(j,i) + fr(j,i) 
                fup = 0.
                zp(j,i) = 0.
                fkp(j,i) = 0.
                fv(j,i) = 0.
                fi(j,i) = 0.
                topfil(j,i) = .false.
              end if

              fadd = rhc*refdt + (1.-rhc)*rfdt
              reffa = fadd/dtu
              if (fr(j,i) .gt. 0.) then
C                                                           redistribute during
C                                                                this time step
                fmode = ' redis'
                thor = tho(j,i)
                if(diag .and. thor .lt. 0.) then
                  write(99,'(" rerr:",a6,i2,4g12.4)') fmode, j, fr(j,i) 
                end if
C!!  new trace variable
                write(calldat,'(" r",i2,i3)')j,itm
C!!     new trace added to parameter list
                call rectr (fr(j,i), tho(j,i), thin(i), reffa, tho2, 
     &             calldat)
C
                fs(j,i) = fs(j,i) + fadd
                thiu(j,i) = tho(j,i)
C!!  add check  12/02
                deltho = tho(j,i) - thin(i)
                if(deltho .gt. 1.e-6) zb(j,i) = fs(j,i)/deltho
C!! end add
                ziu(j) = zb(j,i)
C  depth to which rewet pulse must go to find initial conditions:
                if (twolay(i)) then
                  fcr(j,i) = zfl(i)*(ths1(i) - thiu(j,i))
                  if(lowet(j,i)) then
C
                    if(.not. topfil(j,i)) 
     &                   ffl(j,i) = zfl(i)*(tho(j,i) - thin(i))
                    zb(j,i) = zfl(i) + (fs(j,i) - 
     &                        ffl(j,i))/(tho2 - thin2(i))
                  else  
                    ffl(j,i) = ffl(j,i) + fadd
                    if(zb(j,i) .ge. zfl(i)) lowet(j,i) = .true.
                  end if
                end if   ! end twolayer case
                srat = (thiu(j,i) - thr1(i)) / (ths1(i) - thr1(i))

                if (srat .le. 0.) then    !     upper layer relative redistr. flux

                  qki(j,i) = 0.

                else

                  qki(j,i) = srat**eps1(i)

                end if
C                  if(diag .and. j .eq. jd) then
C                    write(99,'(" redist ",8g12.4)') 
C     &    rfj*po, for, fr(j,i), thor, tho(j,i), tho2, qki(j,i)
C                  end if
C
              else    !  no redistributible volume:
C                        prewet -- low rf, no fr or fv:  accumulate fi()
C
                fmode = ' prewt'
                if(fi(j,i) .gt. 0. .and. tho(j,i) .gt. thin(i)) then
                  zz = fi(j,i)/(tho(j,i) - thin(i))
                else
                  tho(j,i) = thin(i) + .001  ! need a kickoff value for RECTR
                  zz = 0.001
                end if
                set = 'pr'
                thiv = thiu(j,i)
C!! twolay test new 12/02
                if(twolay(i) .and. zz .ge. zfl(i)) then
C                  ffl(j,i) = (tho(j,i) - thin(i))*zfl(i)
                  thiv = th2oth1(thiu(j,i))
                  if(subcon(i)) then
                    lowet(j,i) = .true. ! control in lower layer due to init. wetting
                    zb(j,i) = zz
                  end if
                end if
C here the value of sku and apu can change when zz reaches zfl
C!!  test for sk > 0.
                if(sku(j) .gt. 0.) !**
     &           call setg(thiv, zz, fi(j,i), 0., gu(i), i, sku(j),set)
C 
                srat = (thiu(j,i) - thr1(i)) / (ths1(i) - thr1(i))
C
                if (srat .le. 0.) then

                  qki(j,i) = 0.

                else

                  qki(j,i) = srat**eps1(i)   ! this is dimensionless

                end if

                thur = tho(j,i)
C                thiur = thur  ! temporary
                f2t = fi(j,i)   ! + ref*dtu
C
                if(diag .and. thur .lt. 0.) then
                  write(99,'(" rerr:",a6,i2,4g12.4)') fmode, j, f2t 
                end if
C!!
                write(calldat,'(" p",i2,i3)')j,itm
C!! new trace parameter in call list
                call rectr(f2t, thur, thin(i), reffa, tho2, calldat)  ! growing initial plse
C
                tho(j,i)= amin1 (thur, 0.99*ths1(i))
C
                if( .not. twolay(i)) then
                  thiu(j,i) = 0.6 * tho(j,i) + 0.4 * thin(i)
                else
                  thiu(j,i) = tho(j,i)
                end if
                ziu(j) = zz
                dth = thiu(j,i) - thin(i)
                if(f2t .gt. 1.e-8 .and. dth .gt. 0.) ziu(j) = 
     &                  f2t/(thiu(j,i) - thin(i))
C
                if (twolay(i)) fcr(j,i) = zfl(i)*(ths1(i) - thiu(j,i))
                fi(j,i) = fi(j,i) + fadd
                fs(j,i) = fs(j,i) + fadd
                if (subcon(i)) then
                  if( topfil(j,i) ) then
                    ffl(j,i) = zfl(i)*(tho(j,i) - thin(i))
                  else
                    ffl(j,i) = fs(j,i)
                  end if
                end if
C
C        if(diag .and. j .eq. jd) write(99,'(" prew ",i2,7g13.4)')
C     &  jd, f2t, fi(j,i), thiur, thur, zz, ziu(j), reffa*po

              end if
C   low rf case redistributed.  Now cycle:
              go to 100

            end if
C--------------------------------------------------------------------------
C   higher rf case ************************ potential for ponding
            if (ref .gt. vmin .and. vav(j,i) .le. vmin) vav(j,i) = ref
C    
            fbl = fs(j,i) - fup
            f2l = max(fbl,0.)   !fbl
C!!  change: compare with small value, not zero
C From version dated 3/30/2005: if (fi(j,i) .gt. 0.001) then
            if (fi(j,i) .gt. 0.0005) then  ! api .lt. 0.5mm ignored
C!! new 12/02: value of storage until prewet pulse is passed:
C From version dated 3/30/2005: capr = (ths1(i) - tho(J,i))*zb(j,i)
              capr = (ths1(i) - tho(J,i))*ziu(j)       ! zb replaced w. ziu
              if(subcon(i) .and. lowet(j,i)) then
               fi(j,i) = 0.
               fv(j,i) = fs(j,i)
               qki(j,i) = 0.
C!! new test at start for negligible fi
              else if(ref*dtu .lt. capr) then
               zb(j,i) = ziu(j)
               ziu(j) = 0. !**
               fr(j,i) = fr(j,i) + fi(j,i)  ! put initial input water into redist.
               fi(j,i) = 0.
               bav = ref
               f2b = fbl + ref*dtu
               fmode = 'rewet '
               go to 75                     ! go to redistr. section
C!!  Start new section added since version dated 3/30/2005:
              else  ! input bypasses initial wetting depth !**
                f2l = f2l + fi(j,i) !**
                fi(j,i) = 0. !**
                thiu(j,i) = thin(i) !**
              end if !**
C!!  new - include small fi from low api in f2l !**
            else !**
C!!              f2l = f2l + fi(j,i) !**
              if (fi(j,i) .gt. 0.) f2l = f2l + fi(j,i) ! mod by CLU 7/09
              ziu(j) = 0. !**
              fi(j,i) = 0. !**
C              thiu(j,i) = thin(i)
C!!  End new section added since version dated 3/30/2005:
            end if
C                                                        end preponding section
            thiv = thin(i)
            hv = ht                                 ! new 6/02
C!! add twolay test 12/02
            if(twolay(i) .and. zb(j,i) .ge. zfl(i)) then
              thiv = thin2(i)
              thi2 = thiv
CCC              if(topfil(j,i)) hv = ht + 0.5*zfl(i)  ! proposed 6/02   
            end if
C
            set = 'cd'
            fmode = 'base  '
            call setg(thiv, zb(j,i), fs(j,i), hv, gu(i), i,sku(j), set)
C
            vmb = sku(j)
C                                     set appropriate capillary drive parameter
C!! add twolay test 12/02
            if(twolay(i) .and. zb(j,i) .gt. zfl(i)) then
C                                                      preponding for base case
              
              if (g2(i) .le. 0.0001) then   ! fixed Ks lower layer
C                                                         fixed lower loss rate
                f2b = f2l + vmin2(i) * dtu

                go to 75

              end if
C
            end if
C
C
            If(.not. lowet(j,i)) then  ! single soil or else unfilled:
              if (ref .gt. vmin .and. tho(j,i) .lt. ths1(i)) then
                f2t = fs(j,i) !- fkp(j,i) + (rfj - vmin*qki(j,i))*dtu
                if (f2t .lt. 0.) f2t = 0.
C
                thur = max(tho(j,i),(thin(i)+.01))
                f2o = f2t
C!! new trace parameter for rectr call  18/12/02
                write(calldat,'(" b",i2,i3)')j,itm
C                                                               estimate rising
C                                                            surface saturation
                call rectr(f2t, thur, thin(i), ref, tho2, calldat)
 !
                tho(j,i) = amin1 (thur, ths1(i))
C
              end if
            else
              f2r = ffl(j,i)/cumcm(i)
              thotry = f2r*(ths1(i) - thin(i)) + thin(i)
              if(thotry .ge. tho(j,i)) tho(j,i) = thotry
            End If
C

            ifl = .false.

            Vol = f2l
            vst = vol/apu  ! dimensionless depth
            vmu = vmin
C      
            call fcapv (Vol, fo, ifl, dum)  ! uses vmu for Ksat
C                                 now calculate base infiltration, original
C                                    initial conditions
            fil = fo   !  needed in errvf
C
            If(vst .gt. 0.05 .or. cv(i) .ge. 0.1) then
C!! new calculator for estimating delta F:
              Vx = favdI(apu,vmin,vol,dtu)
              f2b = f2l + Vx
              vl = f2l + 1.01*Vx
C!!              vl = f2l+ 1.01*fil*dtu
            else
C
              tsl = vst + exp(-vst) - 1.  ! parlange/ dimensionless
              tsn = tsl + dtu*vmu/apu
              vsn = tsn + sqrt(2.*tsn) ! philip
              f2b = vsn*apu 
              vl = 1.5*f2b
            end if
C                                                              vl = upper limit
            f2m = f2l + vmin * dtu
C                                                             f2m = lower limit
C
      if((f2m .lt. 0. .or. f2l .lt. 0. .or. vl .lt. 0.) .and. diag) then
         write(99,'(" b",i3,4g13.4)') j, fs(j,i), fup, ffl(j,i), fv(j,i)
        end if
C From version dated 3/30/2005: if (f2b .lt. 1.e-5) f2b = 1.e-5
C!!            if (f2b .lt. 1.e-6) f2b = 1.e-6 !**
C    errvf uses vmu for effective Ksat:
C   iterate: try f2b with limits vl(max) and f2m(min): 
            trial = f2b
            source = 'infiln'

            call iter (errvf, f2b, f2m, vl, ierr, source)

            if (ierr .gt. 0) then
              ftav = (f2b - f2l)/dtu
              write(77,788)ierr, j, t(i)/60., f2l, trial, f2m, vl, f2b,  
     & vmu*po, ftav*po, vst
 788  format(" infilt err:",2I3,f7.1,8g12.4)
              stop ' error - no convergence (errvf) base'
            else if(j .eq. jd) then
              niter = niter - ierr
            end if
C
            if(diag .and. j .eq. jd) then
              ftav = (f2b - f2l)/dtu
             write(99,789) j, trial, f2l, f2m, vl, f2b, ftav*po, fil*po,
     &          apu
 789  format(" t, etc.:",I5,8g12.4) !,f7.1
            end if
C  compare capacity with available:
            flm = f2l + ref * dtu
C
            if (f2b .gt. flm) f2b = flm
            bav = (f2b - f2l) / dtu
C
C this is the depth assuming wetting to saturation: (not true after lowet)
            divz = tho(j,i) - thin(i) !**
C From version dated 3/30/2005: if(divz .le. 0.) then
C From version dated 3/30/2005:   stop ' infilt: divz zero '
C From version dated 3/30/2005: end if
C From version dated 3/30/2005: zb(j,i) = f2b  / divz ! - fkp(j,i)
            if(divz .gt. 0.) zb(j,i) = f2b  / divz ! - fkp(j,i) !**
C
            If(twolay(i)) then
C
              if(subcon(i)) then
C
                if(lowet(j,i) ) then
                  zb(j,i) = zfl(i) + (f2b )/(thm2(i)-thin2(i)) !- fkp(j,i)
                else if(zb(j,i) .gt. zfl(i)) then !**
                  lowet(j,i) = .true.
                end if
C
              else if (f2b .gt. cumcm(i)) then
C                          ! pulse beyond crust into lowerlayer
                if (g2(i) .le. 0.) then

                  trial = (f2b - cumcm(i)) / (f2b - f2l)*sk2(i)*dtu
                  if (trial .le. cumcm(i)) f2b = cumcm(i)

                else

                  zb(j,i) = zfl(i) + (f2b - cumcm(i)) /
     &                    (thm2(i) - thin2(i))
                end if
C
              end if
C     
            End If
C
  75        continue
C
CCC         if(lowet(i) .and. ref .gt. sk1(i))  ...to be developed
            if (fr(j,i) .gt. 0.) then
C               ************* rewet or secondary pulse ***********
              zs = zp(j,i) / zb(j,i)
              thiv = thiu(j,i) 
              fmode(5:6) = 'rw'
C
              gb = gu(i)
              vqu = vmin  ! vqu is used for rewet pulse
              if(zs .lt. 1. .and. zp(j,i) .lt. zfl(i)) then
                fupr = 0. 
                vqu = sk1(i)*vfmf
              else if( subcon(i)) then  
                if(zp(j,i) .gt. zfl(i)) then
                  fv(j,i) = f2l - fkp(j,i)
                  fr(j,i) = 0.  ! to bypass rewet infil calcs.
C  
                end if
              end if
              vqki = vqu * qki(j,i)
              f2p = f2l + ref*dtu
C              
C              if(diag .and. j .eq. jd) then
C                write(99,'(" rw ",7g12.4)') fr(j,i), vqu, vqki, bav, ref
C     & , zb(j,i), thiv
C              end if
C
              if(bav .lt. vqki .or. vqu .gt. ref ) then !.or. fr(j,i) .le. 0. 
                zp(j,i) = zb(j,i) + 1.e-8
C
              else  if(vqu .lt. ref) then
C!! add twlolay test 12/02                         create transitional drive term
                if (twolay(i) .and. zp(j,i) .ge. zfl(i)) then    
C                                       ! wetting into second layer
                  if (g2(i) .eq. 0.) then
C                                                       fixed lower loss rate
                    f2n = f2l + vmin2(i) * dtu
                    if(rfj .gt. vmin2(i)) afac(j) = 1.0
                    go to 100

                  end if

                  thi2 = th2oth1 (thiu(j,i))
                  thiv = thi2
C           
                end if
C
                set = 'rw'
                write(set,'(I2)') j
                call setg (thiv, zp(j,i), fv(j,i), ht, gu(i), i, vmb,
     &                set)
C 
                vqu = vmb * vfmf      ! ? new
                if (abs (vmin - vqu) .le. 0.001 * vqu) then ! vmb is sku found from setg above

                  fth = (ths1(i) - thiu(j,i)) / (ths1(i) - thin(i))

                  if (fth .lt. 0.15) then
                    gama = 1. - 0.15 * zs
                    if (gama .le. 0.5) gama = 0.5
                    zs = zs * gama
                  end if

                  exfun = 0.

                  if (zs .lt. 0.999 .and. zs .gt. 0.) exfun = zs**1.5
C    effective net drive value
                  ge = gu(i) + (gb - gu(i)) * exfun

                  if( ge .lt. 0.) ge = 0. !** added since version dated 3/30/2005

                  apu = (ge / gu(i)) * apu
!**                  if( apu .le. 0.) then
!**                   iunit = 77   ! output file should always be opened
!**                   inquire (98, opened = op)
!**                   if(op) iunit = 98
!**            write(iunit,"(' APU negative in infilt-rw: ',2i3,(6g13.5))") 
!**     &           j, i, t(i),zb(j,i),ht, gb, gu(i), thiv 
!**                 stop ' negative value of cap. drive in INFILT/rw for th
!**     &is element '
!**                  end if

                end if

                f2l = fv(j,i) - fkp(j,i)  ! - fupr
                vmu = vqu * (1. - qki(j,i))
C                                                             
                Vol = f2l

C From version dated 3/30/2005: vls = Vol/apu
                if( apu .gt. 0.) then !**
                  vls = Vol/apu
                else !**
                  vls = 1.E10 !**
                end if !**

                ifl = .false.

                call fcapv (Vol, fo, ifl, dum)  ! uses vmu set just above.
C                               Now make explicit estimate of f2p (vmin > vmu)
                fil = fo
C                 
                if(vls .ge. 0.05 .or. cv(i) .gt. 0.1) then
C!!  replacement estimator for delta(F)                
                  Vx = favdI(apu, vmu, vol, dtu) !**
                  f2p = f2l + Vx !**
C!!                  vl = f2l + 1.01*fil*dtu   !!  3/04 - better estimate of upper lim   
                  vl = f2l + 1.02*vx !**
C
                else
                  tsl = vls + exp(-vls) - 1.  ! parlange/ dimensionless

C From version dated 3/30/2005: tsn = tsl + dtu*vmu/apu
                  if( apu .gt. 0.) then !**
                    tsn = tsl + dtu*vmu/apu
                  else !**
                    tsn = tsl !**
                  end if !**

                  vns = tsn + sqrt(2.*tsn) ! philip
                  f2p = vns*apu 
                  vl = 1.5*f2p
                end if
C
                f2m = f2l + vmu * dtu
C iterate: try f2p within limits vl(hi) and f2m(lo): 
C                                                                   rewet pulse
      if((f2m .lt. 0. .or. f2l .lt. 0. .or. vl .lt. 0.) .and. diag) then
        write(99,'(" p",i3,6g13.4)') j, fs(j,i), fupr, ffl(j,i), fv(j,i)
     & , f2l
        end if
C                                                          infiltration depth
                source = 'infilr'
                trial = f2p
!                if( f2m .ge. f2p) f2m = 0.5*f2p !**
!                if( vl .le. f2p) vl = 2.*f2p !**
C
                call iter (errvf, f2p, f2m, vl, ierr, source)
C                                        ervf uses vmu for the effective Ksat
                if (ierr .gt. 0) then
         write(77,969) j,ref*po, f2l,f2p,f2m,vl,trial,vns,tsn
                  stop ' error - no convergence (errvf) rewet'
                end if
                if (f2p .lt. f2l) then
         write(77,969) j,ref*po, f2l,f2p,f2m,vmu,vl,bav*po,vqki*po
  969  format(i3,8g13.4)
                  stop ' error - infilt/iter: Fnew less than Flast'
                end if

                dft = ref * dtu
                flm = f2l + dft - vqki * dtu
                if (flm .lt. f2p) f2p = flm
                rwav = (f2p - f2l) / dtu + vqki
C 
C            if(diag .and. j .eq. jd) then    
C        write(99,790) j, f2m,vl,f2p,fv(j,i), vmu*po,
C     &  rwav*po,fkp(j,i),zp(j,i)
 790  format(" t, rew:",I5,8g12.4)
C            end if
C                                                       rewet flux should never
C                                                        be less than base flux

                if (rwav .lt. bav) f2p = f2l + (bav - vqki) * dtu

                fkp(j,i) = fkp(j,i) + vqki * dtu  !sum of vol moving thru init. pulse
                zp(j,i) = f2p / (thm1(i) - thiu(j,i)) !assumes rewet pulse is saturated
C!!
                if (twolay(i) .and. zp(j,i) .gt. zfl(i)) then

                  if (g2(i) .eq. 0.) then

                  trial = (f2p - fcr(j,i)) / (f2p - f2l) * sk2(i) * dtu
                    if (trial .le. fcr(j,i)) f2p = fcr(j,i)
                    zp(j,i) = zfl(i) + .01

                  else

                    zp(j,i) = zfl(i) + (f2p - fcr(j,i)) /
     &                      (thm2(i) - thi2)
                  end if
                end if
              end if    

              if(zp(j,i) .ge. zb(j,i)) then
C                                                           rewet pulse catches
C                                                                original pulse
                fr(j,i) = 0.
                fkp(j,i) = 0.
                f2l = fs(j,i) - fup
                thiu(j,i) = thin(i)
                f2n = f2b
                qki(j,i) = 0.
                zp(j,i) = 0.
                fcr(j,i) = cumcm(i)

              else

                f2n = f2p

              end if
C end rewet pulse section
            else 
C                            no rewet pulse
              zp(j,i) = 0.
              f2n = f2b

            end if
C                                  For two-layer system, use "excess" to
C                                  fill upper layer until saturated:
            rxv = 0.
            vavu = (f2n - f2l) / dtu + vmin * qki(j,i)
C   
            rmb = max((rfj-vavu),0.)*dtu
            vav(j,i) = vavu
            if(subcon(i)) then
              if( lowet(j,i) .and. .not. topfil(j,i)) then
                rxv = dtu*(ref- vavu)*rhc + rmb*(1.-rhc)
                if(rxv .gt. 0.) then
C      fill upper layer storage from excess at interface
                  defl = cumcm(i) - ffl(j,i)
                  if(rxv .le. defl) then
                    ffl(j,i) = ffl(j,i) + rxv
                    vav(j,i) = ref   ! this is infil rate at the surface
                  else
                    rxv = defl
                    ffl(j,i) = cumcm(i) + 1.e-8
                    vav(j,i) = vavu + defl/dtu
                  end if
                  tho(j,i) = thin(i) + ffl(j,i)/zfl(i)
                end if
              else if(.not. lowet(j,i)) then
                rfin = min(rfdt,(vavu*dtu))
                ffl(j,i) = ffl(j,i) + rhc*(f2n-f2l) + (1.-rhc)*rfin
              end if
C                     ! single layer system or else top already full
            end if
C                     weight increase of storage based on rhc
      fv(j,i) = rxv + (fup +f2l) + rhc*(f2n -f2l) + fkp(j,i) + 
     & (1.-rhc)*(rfdt - rmb)
      fs(j,i) = rxv + (fup+fbl) + rhc*(f2n -f2l) + (1.-rhc)*(rfdt-rmb)  

          end if
C          
100       continue
C
          if(sk1(i) .gt. 0.) then
            if (diag .and. j .eq. jd) then   ! put any diagnostics here:
              sumif = sumif + vav(j,1)*dtu
C              rmu = ref*po
              if(twolay(i)) then
                write(99,198) 
     & lowet(j,i), topfil(j,i), qki(j,i), tho(j,i), thiu(j,i), defl, rxv
              else
                write(99,'(" tho/thiu: ",2g12.4)') tho(j,i), thiu(j,i)
              end if
              write(99,199)
              write(99,'(f6.1,I3,i4,1x,a6,8g12.4)') trp,jd,niter,fmode,
     &  rfj*po,vav(j,1)*po,vmin*po,fs(j,i),fv(j,i),ffl(j,i),zp(j,i) ! , sumrf, sumif
C
            end if
          end if

        end do
C                                                           ...end 1 -> nk loop
C------------------------------------------------------------------------------
C
        if(diag) then
          sumfs = 0.
          nkm = nk - 1
          do j = 2,nk
            sumfs = sumfs +  fs(j,i) 
          end do
          sumfs = qbal(1)*sumfs/real(nkm,4)
          write(99,'(" net inf:", g12.4)') sumfs
        end if
  198 format(" lowet topfil    qki      tho         thiu         defl      
     &    rxv"/3x,L2,4x,L2,5g12.4)
  199  format(" infl:t  j iter  mode    rf         f          Ks(eff)      
     &   Fs         Fv          F(1)        zp")
C
      if(itm .ge. itlim ) then
        sumtab(ltab)%flowr = 0.
        if(diag.and. sk1(i) .gt. 0.) write(99,'(" Tot fs:"/,(8g12.4))') 
     &   (fs(j,1),j=2,nk)
        nkm = nk-1
        if(subcon(i)) then
          finlot = 0
C       
          do j = 2, nk 
            finlo = 0.
            if(lowet(j,i)) finlo = max((fs(j,i) - ffl(j,i)),0.)
            finlot = finlot + finlo
          end do

          finlom = finlot/real(nkm,4)
C
          if(finlom .gt. 1.e-6) then
            jmns = -1
C           write(99,'(" finlo ",3g13.4)') finlot, fs(2,1), ffl(2,1)
            call qwrt (idl, jmns, trace, finlom, dum, dum,
     &       dum, idu, 0, i, t(i))  ! pass cumcm and t to writer
          end if
C
          ni = ltab + i - 1
          sumtab(ni)%flowr = finlom
          if(diag)   write(99,'(" fl:",8g12.4)') (ffl(j,1),j=2,nk)
        end if
C here one could recalc and get thend = (fs - finlo)/zfl
      end if
C
C------------------------------------------------------------------------------
      return

C------------------------------------------------------------------------------


      entry infil0 (iv, ki, nx, msg, diagn, sat)
C                                                    get / initialize variables
      trace = 'infl0 '
      t = (/0.,0./)
      idl = iv
      idu(ki) = iv
      niter = 0
      notice = 'element    '
      rsoil = .true.
      write(notice(9:11),'(I3)') idl
      notify = .false.
      sumrf = 0.
      sumif = 0.
C
      diag = diagn
      dbg = diagn

      if (ki .eq. 1) then              ! plane or main channel

        nk = nx

        if (units .eq. 1) then

    !          conv = 1000.
          con2 = 1.

        else

    !          conv = 12.
          con2 = 3.28

        end if
C                                                                values assumed
!        smax = 0.95
        po = conv*3600.  ! printout coefficient to get mm/h(in./hr) 
        alf = 0.8        ! from m/.sec (ft/sec)

      end if
C 
      jd = 5    ! node for which diagnostics will be written
C
      if(jd .gt. nk) jd = nk
      i = ki    ! this is nchan for initialization
      zfl(i) = 0.
      nvi = ltab + i - 1
C                                                             Ks of upper layer
C                                                              (mm/hr or in/hr)
      call getr4 ('KS', 1, sk1(i), ierr)

      if (ierr .gt. 0) sk1(i) = 0.
C                                                           impervious surface?
      if (sk1(i) .eq. 0.) then
        rsoil = .false.
        sumtab(nvi)%twola = .false.
        sumtab(nvi)%thst = 0.
        return
      end if

      sk1(i) = rmks * sk1(i) / (3600. * conv)    ! units are m(f)/sec
C                                                          capillary suction of
C                                                        upper layer (mm or in)
      call getr4 ('G', 1, g1(i), ierr)

      if (ierr .gt. 0) g1(i) = 0.
      if(rmg .le. 0.) stop ' zero multiplier'

      g1(i) = rmg * g1(i) / conv                       ! convert to m or ft

      if (g1(i) .gt. 0.) then
C                                                         pore size dist. index
C                                                                of upper layer
        call getr4 ('DI', 1, al1(i), ierr)

        if (ierr .gt. 0) call errxit
     &  (notice, 'pore size dist. index (DI, upper layer) not found')
        if(al1(i) .gt. 1.5) call errxit(notice,' pore size dist. index c
     &annot be > 1.5')

        call getr4 ('PO', 1, por1(i), ierr)
C                                                       porosity of upper layer
        if (ierr .gt. 0) call errxit
     &  (notice, 'error - porosity (PO, upper layer) not found')

        call getr4 ('RO', 1, rock1(i), ierr)
C                                                             vol. rock content
C                                                                of upper layer
        if (ierr .gt. 0) rock1(1) = 0.
C                                                           residual saturation
        call getr4('SR', 1, ser1, ierr)
C
        if(ierr .gt. 0) then
          ser1 = -1.
        end if
C                                                            maximum saturation
        call getr4('SM', 1, smax, ierr)
C
        if(ierr .gt. 0) then
          smax = 0.95
        end if
C
        ths1(i) = por1(i) * smax
C
        if(SAT .lt. 0.) then

          call getr4 ('SA', 0, sint, ierr)
C                                                            initial saturation
          if (ierr .gt. 0) then
C          
            call errxit
     &    (notice, 'initial soil saturation (SA) not specified')
          end if
C        
          else
          
            sint = sat
          
        end if

        if (sint .eq. smax) then

          g1(i) = 0.
          
        else if (sint .gt. smax) then

          call errxit (notice, 'initial soil saturation (SA)
     & cannot be greater than the maximum (SMAX)')

        end if

      else
        rsoil = .false.
        ths1(i) = .4                               ! for g1 = 0.
        sint = 0.2
        thin(i) = .16

      end if

      call getr4 ('CV', 0, cv(i), ierr)
C                                                     coeff. of variation of Ks
      if (ierr .gt. 0) cv(i) = 0.

      cv(i) = cv(i) * rmcv   ! CV multiplier, default 1.0; from module 'multip'
        

      call getr4 ('THI', 0, zfl(i), ierr)
C                                                            thickness of upper
C                                                              layer (mm or in)
      if (ierr .gt. 0) zfl(i) = 0.
C!!  2 lines moved for default:  3/05
      g2(i) = 0.
      g2u(i) = 0.
       
      if (zfl(i) .eq. 0.) then
C                                                                  single layer
        thin2(i) = 0.
        ths2(i) = 0.
        twolay(i) = .false.
        subcon(i) = .false.
        sumtab(nvi)%twola = .false.
      else
C                                                                    two layers
        zfl(i) = zfl(i) / conv ! now in meters
        if((units .eq. 1 .and. zfl(i) .lt. 1.e-2) .or.  ! 1 cm limit
     &       units .ne. 1 .and. zfl(i) .lt. .05) 
     &     call errxit(notice,' upper soil layer too thin for use: check 
     & units ')
C
        twolay(i) = .true.
        depmax = max(zfl(i),depmax)
        sumtab(nvi)%twola = .true.

        call getr4 ('KS', 2, sk2(i), ierr)

        if (ierr .gt. 0) sk2(i) = 0.
C!!                                                             Ks of lower layer
C!!        if (sk2(i) .le. 0.) then
C                                                        impervious lower layer
C!!   corrected 11/03 was 0., but that causes problems
C!!          g2(i) = 10.                                       

C!!        else
        if ( sk2(i) .gt. 0.000001) then
C                                                                      pervious
          sk2(i) = rmks * sk2(i) / (3600. * conv)

          call getr4 ('G', 2, g2(i), ierr)
C                                                capillary drive of lower layer
          if (ierr .gt. 0) g2(i) = 0.
          g2(i) = rmg * g2(i) / conv


          if (g2(i) .gt. 0.) then
C                                                         pore size dist. index
C                                                                of lower layer
            call getr4 ('DI', 2, al2(i), ierr)

            if (ierr .gt. 0) call errxit
     &     (notice, 'pore size dist. index (DI, lower layer) not found')
            if(al2(i) .gt. 1.5) call errxit
     &      (notice,' pore size dist. index cannot be >1.5')

            call getr4 ('PO', 2, por2(i), ierr)
C                                                       porosity of lower layer
            if (ierr .gt. 0) call errxit
     &      (notice, 'porosity (PO), lower layer) not found')
C
            call getr4 ('SR', 2, ser2, ierr)
            if(ierr .gt. 0) ser2 = -1.

            call getr4 ('RO', 2, rock2(i), ierr)
C                                                      vol. rock of lower layer
              if (ierr .gt. 0) rock2(i) = 0.
              ths2(i) = smax * por2(i)

          else if (sk2(i) .ge. sk1(i)) then

            call errxit
     &      (notice, 'fixed lower seepage rate > Ksat of upper soil')

          end if
        end if
      end if
C                                                       initial calculations...
C------------------------------------------------------------------------------

C                                                                values assumed
      if (g1(i) .gt. 0.) then
C!                  estimating function, g in m or ft
        if(ser1 .lt. 0.) ser1 = 0.1 + 0.3/(1.+(con2/g1(i))**4.)**.25  
        thr1(i) = ser1 * por1(i)             ! could be *ths1(i)
        dth1 = ths1(i) - thr1(i)
C!        thin(i) = sint * dth1 + thr1(i) ! this needs revision
        thin(i) = sint * por1(i) ! revised 9-3-2009 by CLU
        if(thin(i) .lt. thr1(i)) thin(i) = thr1(i) ! revised 9-3-2009 by CLU
C
        sumtab(nvi)%thst = thin(i)
C
        def1(i) = (ths1(i) - thin(i)) * (1. - rock1(i))
        if(al1(i) .le. 0.) al1(i) = 0.2  ! insert default value for lambda if necessary
        if(al1(i) .gt. 1.5) al1(i) = 1.5
C            estimate a value for pb:
        pb1(i) = g1(i) * (2. + 5. * al1(i)) / (3.4 + 3. * al1(i))  ! rev 4/02
        eps1(i) = (2. + 3. * al1(i)) / al1(i)
        erfi(i) = sk1(i)*sint**eps1(i)
C        write(77,'(" rfi ",g12.4)') erfi(i)
        thwilt = thr1(i) + dth1*(pb1(i)/153.)**al1(i)
        thfc = thr1(i) + dth1*(pb1(i)/3.33)**al1(i)
C
        if(twolay(i)) then
C                    this is max water held in upper layer at field capacy:
          sumtab(nvi)%flcap = zfl(i)*(thfc - thin(i)) * (1.-rock1(i))
C
        end if
        thm1(i) = ths1(i)

      end if

      if(twolay(i) .and. g2(i) .gt. 0.) then

        pb2(i) = g2(i) * (2. + 5. * al2(i)) / (3.4 + 3. * al2(i))  ! rev 4/02
C             ! estimating function, g in m
        if(ser2 .lt. 0.) ser2 = 0.1 + 0.3/(1.+(con2/g2(i))**4.)**.25  
        thr2(i) = por2(i)*ser2  
        thin2(i) = th2oth1 (thin(i))
        def2(i) = (ths2(i) - thin2(i)) * (1. - rock2(i))
        if(al2(i) .le. 0.) al2(i) = 0.2  ! insert default value for lambda if necessary
        if(al2(i) .gt. 1.5) al2(i) = 1.5
        eps2(i) = (2. + 3. * al2(i)) / al2(i)
C!!  new lines 3/05
      else if (g2(i) .le. 0.) then !**
        thr2(i) = 0. !**
        def2(i) = 0. !**

      end if
C
      gu(i) = g1(i)
      thm1(i) = ths1(i)
      cumcm(i) = 1000.

      if (twolay(i)) then
C                                                        analyze flow potential
        cumcm(i) = zfl(i) * (ths1(i) - thin(i))
C
        if(api .gt. 0.) then
          sumtab(nvi)%pored = (zfl(i)*(ths1(i)-thwilt) - api)*conv
        else
          sumtab(nvi)%pored = zfl(i)*def1(i)*conv   !
!* comment out following line as j is not defined
!*          ziu(j) = 0.   !! new
        end if
C
        vmin2(i) = sk2(i)
C!! test also for impervious lower layer 
        if(g2(i) .gt. 0. .and. sk2(i) .gt. 0.) then

          hc = zfl(i) * (1. - sk2(i) / sk1(i))   ! can be large
          thm2(i) = ths2(i)

          if (sk1(i) .lt. sk2(i)) then
C!!                                                            get vmin and thm's
C                                                              for crusted layered case
            vmu = sk1(i)
            defl = 0.
            rxv = 0.
C
            rc(i) = sk1(i) * (1. + alf / (exp (alf *
     &              cumcm(i) / g1(i)) - 1.))
C!!  abs() added
            hc = abs(hc)
            if (hc .gt. pb2(i))  then

              source = 'errlay'
              call iter (errly, hc, pb2(i), 100. * pb2(i), ierr, source)
C    this iteration finds effective K for crusted case:  vmu
              if (ierr .gt. 0) stop ' error - no convergence (errly)'

              vmin2(i) = vmu   ! effective composite K for crust
              thm2(i) = thr2(i) + (ths2(i) - thr2(i)) *
     &                  (pb2(i) / hc)**al2(i)

              if (hc .gt. pb1(i)) thm1(i) = thr1(i) + (ths1(i) -
     &                            thr1(i)) * (pb1(i) / hc)**al1(i)

              g2u(i)= pb2(i)*flams(al2(i),thm2(i),ths2(i),thr2(i),trace)

            else

              g2u(i) = g2(i) - hc

            end if

            rat = sk2(i) / sk1(i)
            bf = (alog (vmin2(i) / sk1(i)) + 1.) / (alog (rat) + 1.)

            if (bf. lt. 1.0) then

              wcf(i) = bf / (1. - bf)

              if (wcf(i) .gt. 2. * rat) wcf(i) = 2. * rat

            else

              wcf(i) = 2. * rat

            end if

          else  ! sk1 > sk2
C                                                       restrictive lower layer
            g2u(i) = g2(i) !  - hc
            vmu = sk2(i)
            if(sk1(i) .gt. 2.*sk2(i)) subcon(i) = .true.
C
          end if
C!! fixed ks lower layer
        else !**
           g2u(i) = g2(i) !**
           if(sk1(i) .gt. 2.*sk2(i))  subcon(i) = .true. !**
        end if
C!!    11/03
      else
C                                                               no second layer
        vmu = sk1(i)
        thm1(i) = ths1(i)
        zfl(i) = 1000.
        cumcm(i) = 1000.
        sk2(i) = sk1(i)

      end if
C
      satin = sint
      fsz = 0.
      thinu = thin(i)
C!! corrected position, 9/03      
      zd = 0.
      apa = api
      satfc = (thfc-thr1(i))/dth1 !**
      satwl = (thwilt-thr1(i))/dth1 !**
      if(api .gt. 0.) then
        thdel = max(0.,(thin(i)-thwilt))
        zd = zfl(i)
        if(thdel .gt. 0.) then  ! if thdel .le. 0. there is no pulse
          zt =  api/thdel
          if(zt .ge. zd) then  ! set apa limited by interface
            apa = zd*thdel
            zd = zfl(i) + .001
            thin2(i) = th2oth1(thwilt)
            def2(i) = (ths2(i) - thin2(i)) * (1. - rock2(i))
            if(twolay(i)) then  ! redefine with thin, no api 
              sumtab(nvi)%pored = zfl(i)*(ths1(i)-thin(i)) *conv
            end if
          else  !  set bottom of pulse < interface
            zd = zt
            thin(i) = min(thin(i),thwilt)
C  this operates below api level:
            if(twolay(i)) then   ! these have been set, but need be changed
C this is defined as field capacity in mm LESS initial storage:
             sumtab(nvi)%flcap = zfl(i)*(thfc-thwilt)*(1.-rock1(i)) -apa
              cumcm(i) = zfl(i) * (ths1(i) - thin(i))
              thin2(i) = th2oth1(thin(i))
              def2(i) = (ths2(i) - thin2(i)) * (1. - rock2(i))
            end if
            fsz = api
          end if
        end if
      end if
C
      do j = 1, nk
        sku(j) = sk1(i)
        tho(j,i) = thinu
        thiu(j,i) = thinu
        fs(j,i) = fsz
        fr(j,i) = 0.
        fi(j,i) = 0.
!* the next line was commented, out, but without it ziu may not be initialized
        ziu(j) = 0.
        if(zd .lt. zfl(i) .and. thdel .gt. 0.) then
C          zb(j,i) = zd
          ziu(j) = zd
C          fr(j,i) = apa
          fi(j,i) = apa
        end if
        ffl(j,i) = fs(j,i)
        fv(j,i) = 0.
    ! new
        qki(j,i) = 0.
        fkp(j,i) = 0.
        fcr(j,i) = cumcm(i)
        zb(j,i) = 0.
        zp(j,i) = 0.
        topfil(j,i) = .false.
        lowet(j,i) = .false.

      end do
C
      if(diagn) then
        write(99,999) vmu*po, cumcm(i), g1(i), g2(i), g2u(i),gu(i),cv(i)
        write(99,998) thin(i),thin2(i),ths1(i),ths2(i),zfl(i),zd,thr1(i)
      end if
 999  format(2x,' vmu          cumc         G1          G2         G2U 
     &        GU         CV(K)'/8g12.4)
 998  format(2x,' thi          thi(2)      ths1        ths2       zfl
     &    zd         thr1 '/8g12.4)
      return

      end
C------------------------------------------------------------------------------


      function th2oth1 (th)

C   Computes saturation in SUBso7il to match head in SURFACE soil

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      logical :: twolay, dbg

      integer i

      thrto = (th - thr1(i)) / (ths1(i) - thr1(i))
      if (thrto .lt. 0.) then
        thrto = 0.
        h1 = 1.e5
      else
        h1 = pb1(i)/thrto**(1./al1(i))
      end if
      if(pb2(i) .le. 0.) then !**
        th2oth1 = thr2(i) !**
      else if(h1 .gt. pb2(i)) then !**
        th2oth1 = thr2(i) + (ths2(i)-thr2(i))*(pb2(i)/h1)**al2(i)
      else
       th2oth1 = 0.99*ths2(i)
      end if
      return

      end
C------------------------------------------------------------------------------


      function th1oth2 (th)

C   Computes saturation in SURFACE soil to match head in SUBsoil


      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      logical :: twolay, dbg

      integer i

      thrto = (th - thr2(i)) / (ths2(i) - thr2(i))
      if(thrto .le. 0.) then
        thrto = 0.
        h2 = 1.e5
      else
        h2 = pb2(i)/thrto**(1./al2(i))
      end if
      if (h2 .gt. pb1(i)) then
        th1oth2 = thr1(i) + (ths1(i)-thr1(i))*(pb1(i)/h2)**al1(i)
      else
        th1oth2 = 0.99*ths1(i)
      end if
      return

      end
C------------------------------------------------------------------------------

      subroutine errvf (F2n, ferr, dferr)

C   Computes the residual and derivative at F2 of the error function form of
C     the 3-parameter infil equation for infiltrated volume as a f(I). 
C     This is found by an external iterating function ITER()
C   Cum. Infil. and infiltrability fil at beginning of step, F2l, are
C      assumed known.
C!!  changed ?/04 using vmu rather than vmin in invariant Ks ferr expressions
C
      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

      logical :: twolay, dbg, ifn

      integer i, idl, j

      data ifn/.true./
C
      if(cv(i) .ge. 0.1) then
C                        ! gets fv2 with spatial capacity function:
        call fcapv (F2n, fv2, ifn, dfdI)  ! finds fv2
        ferr = F2n - F2l - 0.5 * dtu * (fv2 + fil)
        dferr = 1. - 0.5 * dtu * dfdI
C
      else  ! simpler and more precise for no spatial variance case:
C
        xtermL = alf*F2L/apu
        xtermN = alf*F2N/apu
        if(xtermL .gt. 80. .or. xtermN .gt. 80.) then
          ferr = F2N - F2L - vmu*dtu
          dferr = 1.0
          if(dbg) write(77,99) idl, j, trp, F2L, apu, F2n
        else
          FG2 = exp(xtermN)
          ratl = (FG2 - 1. + alf)/(exp(alf*F2l/apu) - 1. + alf)
          if(ratl .le. 0. .or. F2l .lt. 0.) then
            write(77,99) idl, j, ratl, F2l, apu, F2n
  99  format(' neg. F2L or logf invalid in errvf: ',2i4,4g14.4)
            stop ' invalid argument for logf '
          end if
          ferr = F2n - F2l - apu*log(ratl) - vmu*dtu*(1.- alf)
          dferr = 1.0 - alf*FG2/(FG2 - 1. + alf)
        end if
C
      end if
      return
C
      end
C------------------------------------------------------------------------------


      subroutine rectr (fv, th, thi, rf, th2, origin)

C   Computes the redistribution (reduction in saturation) of wetted soil pulse,
C   which is assumed to elongate with time, according to Smith, et al,
C   WRR (1992), modified for a two-layer profile.

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      logical twolay, op, dbg

      integer   :: i, md, nloop, ncuts
      character(LEN=4) :: note, look
      character(LEN=6) :: trace
      character(LEN=7) :: origin

      trace = 'rect  '
      note = 'norm'
      look = 'norm'
      il = 1
      fl = fv
      thu =  th
      if (th .le. thi .and. rf .lt. erfi(i)) then

        this = thi
        go to 90

      end if

      zn = fv / (th - thi + .001)

      if (twolay(i) .and. zn .ge. zfl(i)) then                ! two-layer case

        if (g2(i) .eq. 0.) then
C                                               fixed flux - lower soil
          dfl = sk2(i) * dtu
          fl = fv - dfl
          this = th - dfl / zfl(i)

          if (this .lt. thi) then

            dfl = 0.
            this = thi

          end if

          go to 90

        end if

        il = 2
        th1 = th
        alm = al2(i)
        ths = ths2(i)
        thr = thr2(i)
        pb = pb2(i)
        thu = th2oth1 (th)
        if(thu .lt. thr) stop ' thu<thr in 2lay '
        thit = thin2(i)

      else
C                                            redistribute in upper soil
        il = 1
        alm = al1(i)
        pb = pb1(i)
        ths = ths1(i)
        thr = thr1(i)
        thu = th
        if( thu .lt. thr) stop 'thu<thr in 1lay'
        thit = thi

      end if

      fac = 2.
      if (rf .gt. 1.E-4) fac = 1.5
      if(rf .lt. -1.E-4) fac = 3.
C
      thwt = 0.99 * thu    ! corrected 11/99
C
      trace(5:6) = 'i '
      write(trace(6:6),'(i1)') il

      fog = flams (alm, thwt, ths, thr, trace)    

      gup = pb * fog
      if(thwt .gt. ths .or. thwt .lt. thr) then
        iunit = 77   ! output file should always be opened
        inquire (99, opened = op)
        if(op) iunit = 99
        write(iunit,'(" rect bounds ",i3,5g13.4)') il, thr, thu,ths, gup 
      end if

 
 47   continue
      this = thu  
      t = 0.
      odtt = dtu
      nloop = 0
C                  Runga-Kutta move thru adjustment of saturation
      do while (t .lt. dtu)
        trace(5:6) = '01'
        if(this .le. thr) then
          st = vmin*((thit - thr)/(ths - thr))**eps1(i)
          write(77,277)origin,nloop,fl,thr,thit,thu,this,rf,st,
     &        fk1, fk2, fk3, fk4, th2, th3, th4  !,  avg0, avg1, avg3
  277     format(1x,a7,i3,7g12.4/7g12.4/7g12.4)
          note = 'diag'
          exit
        end if
C
        avg0 = flams (alm, this, ths, thr, trace) * pb          !flams

        fk1 = funk (fl, rf, this, avg0, fac, thit, il, note)           ! funk
C set time step
        dtt = 0.002/abs(fk1)
        if(dtt .gt. 1.5*odtt) dtt = 1.5*odtt  ! don't grow too fast
        if((t+dtt) .gt. dtu) dtt = dtu - t + 1.e-4

        thls  = this
        ncuts = 0
C!!
 50     continue
C!!
        th2 = this - 0.5 * dtt * fk1
        th2 = max(th2,thit)
        odtt = dtt          ! added 10/02

        trace(5:6) = '02'
        avg1 = flams (alm, th2, ths, thr, trace) * pb           !flams

        hft = fl + 0.5 * rf * dtt

        fk2 = funk (hft, rf, th2, avg1, fac, thit, il, note)           ! funk

        th3 = this - 0.5 * dtt * fk2
        th3 = max(th3,thit)

        fk3 = funk ( hft, rf, th3, avg1, fac, thit, il, note)          ! funk

        th4 = this  - dtt * fk3
        th4 = max(th4, thit)

        trace(5:6) = '04'
        avg3 = flams (alm, th4, ths, thr, trace) * pb           ! flams

        ft = fl + rf * dtt

        fk4 = funk (ft, rf, th4, avg3, fac, thit, il, note)            ! funk

        dfac = (fk1 + 2.*(fk2 + fk3) + fk4)/6.
        this = this - dtt * dfac
!! new 4/05
        if(this .lt. thr) then
          if(ncuts .lt. 10) then
            dtt = dtt*0.5
            ncuts = ncuts + 1
            this = thls
            go to 50
          else
            stop ' unable to resolv small th in RECTR with 10 time step 
     &cuts '
          end if
        end if
!!
        fl = ft
        t = t + dtt
        nloop = nloop + 1
        if(dbg .and. nloop .gt. 100) then
          if(rf .gt. vmin .and. this .ge. 0.99*ths) go to 90
          write(99,'(A7,i3,3f9.4/5g13.4/10x,4g12.4)') origin, nloop,
     &     dtt,t,ft,   this, fk1, fk2, fk3, fk4,  th2, th3, th4, dfac 
          if(nloop .gt. 200) stop ' too many loops in rectr'
        end if


      end do
      if(this .lt. thit .or. this .gt. ths) then
        if(this .lt. thr) stop ' debug rect stop'
      end if

90    continue

      th = this

      if (il .ge. 2) then !th = th1oth2 (this)                        !th1oth2
        th2 = th
        th = th1 + this - thu   ! needed for vol. balance
      end if
      fv = fl

      return

      end
C------------------------------------------------------------------------------

      function flams (a, thw, ths, thr, trace)

C   Integral of kr(psi) up to a limit thw, scaled on the value of PB
C   Revised version good up to h = 0.  value of thw first interpreted
C   into equivalent h, and parameterized values of integral used.

C            a   is lambda
C            ths and thr are sat and resid thetas, respectively
C
      character(len=6) :: trace
      logical :: op

      if(thw .gt. ths) thw = ths
      s = (thw - thr) / (ths - thr)
      dp = 1.0 !**
      if(a .lt. 0.2) dp = 0.96 + 0.2*a   ! empirical second order improvement
C
      GP = GofP(a) !**                                                 ! GofP

      if (s .gt. 1.e-10) then
C
        hs = hsofs( s, a )                                           ! hsofs
        be = BParm( a ) 
C!!                                             ! BParm
        if(abs(hs -.05) .le. 1.e-10) then
          flams = GP
        else
          flams = Gp*(1.- dp*(1. +(Aparm(a)/(hs-.05))**be)**(-1./be)) !  AParm()
        end if
      else
C
C        iunit = 77   ! output file should always be opened
        inquire (99, opened = op)
        if(op) then
          iunit = 99
          Write(iunit,"(' Se too small in flams from ',a6,/(5g13.4))") 
     &        trace,  thw, ths, s, a
C        stop ' neg flams '
        end if
  !
        flams = GP*(1.- dp)

      end if

      return

      end

C ----------------------------------------------------------------------

      function hsofs( Se, alam )

C  Transitional Brooks-Corey function to obtain abs(scaled h),
C  assuming C = 5 and Sh/Pe = .05:
C
      If(Se .le. 0.  .or. Se .gt. 1.0) then
        stop ' Call to hsofs with Se out of range'
      Else
        C = 5.
        pw = C/alam
        pwi = 1./C
        term = 1./Se**pw
        hsofs = .05 + (term - 1.0)**pwi
      End If
      Return
      end
C
C ---------------------------------------------------
C
      function GofP( alam )
C
      if(alam .gt. 0.) then
        Gofp = (3.4 + 3.*alam)/(2. + 5.*alam)
      else
        stop ' call to Gofp with alam out of range'
      end if
      return
      end
C ---------------------------------------------------
C
      function Aparm( alam)
C
      if(alam .gt. 0.) then
        Aparm = 0.75*(2. + 3.*alam)/(1. + 3.*alam)
      else
        stop ' call to Aparm with alam out of range'
      end if
      return
      end
C ---------------------------------------------------
C
      function Bparm (alam)
C
      if(alam .gt. 0.) then
        Bparm = 1.8 + 2.1*alam
      else
        stop ' call to Bparm with alam out of range'
      end if
      return
      end

C------------------------------------------------------------------------------

      function funk (f, r, thw, gu, fac, thi, il, str)

C   Decay of water content as block of infiltrated water redistributes:
C   d(thet)/dt

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      character(len=4) :: str
      logical :: twolay, dbg

      integer i, il

      thwt = thw

      if (il .eq. 1) then

        eps = eps1(i)
        ths = ths1(i)
        thr = thr1(i)

      else

        eps = eps2(i)
        ths = ths2(i)
        thr = thr2(i)

      end if

      if (thwt .lt. thr) thwt = thr
      thwu = thwt - thi
C
      thwr = thwt - thr
      if (thwr .lt. 0.) thwr = 0.
C      
      phi = ths - thr
      bf = 0.85
C
      fu = f 
      tm1 = (thwr/phi)**eps
C
      if(fu .gt. 0. .and. thwu .gt. 1.e-4) then
       zuinv = (thwu + .0001)/(fu + .0001)
      else
       zuinv = 5.  ! this is d(theta) divided by infil depth
      endif
C
      tm2 = bf * fac * gu * zuinv
      funk = zuinv * (vmin * (tm1 + tm2) - r)
C
      if(str .eq. 'diag') then
        write(77,'(" funk",5g12.4/5g12.4)')funk, thwu, thwt, fu, zuinv, 
     &   tm1, tm2, gu, r
      end if
C
      return

      end
C------------------------------------------------------------------------------


      function fpond (rr)

C   Current ponding F as a function of current rainrate RR, excluding constant
C   total due to KI, if any.  Use 3=parameter infiltrability function.

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      logical :: twolay, ifn

      integer i


      if (rr .le. vmin) then

        fpond = 1.E6

      else

        fpond = apu / alf * alog ((rr - vmin
     &          * (1. - alf)) / (rr - vmin))

      end if

      return

      end
C------------------------------------------------------------------------------

      function rkohy (p, al, h, ifd, dkdh)

C   Finds simple relative K as f(H) and dK/dH, using Brooks and Corey.
C   Assumes head and entry head P are both positive.


      if (h .lt. p) then

        rkohy = 1.
        if (ifd .gt. 0) dkdh = 0.

      else

        eta = 2. + 3. * al
        rh = h / p
        rkohy = rh**(-eta)
        if (ifd .gt. 0) dkdh = -eta / p * rh**(-eta - 1.)

      end if

      return

      end
C------------------------------------------------------------------------------

      subroutine setg (thi, zb, fv, ht, gef, k, sk, set)

C   Sets appropriate values of G and AP depending on location
C   of wetting front, for 1 or 2 soil layers

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

      logical :: twolay, dbg, op
      character*6 trace, set*2

      integer i, k, idl

      trace = 'sg   0'
      trace(3:4) = set
      if (.not. twolay(i) .or. zb .lt. zfl(i)) then
C                                                              upper layer only
        sk = sk1(i)

        gef = g1(i) - pb1(i) * flams (al1(i),thi,ths1(i),thr1(i),trace)

        apu = (gef + ht) * (ths1(i) - thi) * (1. - rock1(i))
C  
        if(gef .le. 0.) then
           iunit = 77   ! output file should always be opened
           inquire (98, opened = op)
           if(op) iunit = 98
          write(iunit,"(' Gu negative in setg: '/a6,i3,(5g13.4))") 
     &           trace, i, gef, pb1(i),thi, g1(i), ht 
          stop ' negative value of cap. drive in setg for this element '
        end if
      else

        if (sk1(i) .lt. sk2(i)) then
C                                                        upper layer control --
C                                                          matching value for g
          cumi = zfl(i) * (ths1(i) - thi) * (1. - rock1(i))
          defth2 = (thm2(i) - thi2) * (1. - rock2(i))
          trace(5:6) = ' a'
          g2e = g2u(i) - pb2(i)*flams(al2(i),thi2,ths2(i),thr2(i),trace)

          gfm = alf * cumi / (1. + alog (alf *
     &          vmin2(i) / (rc(i) - vmin2(i)))) / defth2
          sk = vmin2(i)
          rx = fv / cumi - 1.
          if (rx .lt. 0.) rx = 0.
          omwt = exp (-wcf(k) * rx / 16.7)
          gef = omwt * gfm + (1.- omwt) * g2e
C          thwf = thm2(i)
          apu = (gef + ht) * defth2
          if(gef .le. 0.) then
           iunit = 77   ! output file should always be opened
           inquire (98, opened = op)
           if(op) iunit = 98
           write(iunit,"(' Gu negative in setg: '/a6,i3,(5g13.4))") 
     &           trace, i, gef, thi, gfm, zb, zfl(i)
          stop ' negative value of cap. drive in setg for this element '
          end if

        else
C                                             possible lower layer control
          sk = sk2(i)
C!!           setg can be called when zz exceeds zc above inpervious rock, 
C!!             so a check is needed here:
          if(sk2(i) .le. 0.) then !**
            apu = 1.       ! dummy !**
            return !**
          end if !**
C!!             end add  11/03
          if (g2(i) .ge. 0.0001) then
C 
            trace(5:6) = ' b'
            gef = g2u(i) -pb2(i)*flams(al2(i),thi,ths2(i),thr2(i),trace)

            apu = (gef + ht) * (thm2(i) - thi2) * (1. - rock2(i))  ! ht new, 6/02

            if( apu .le. 0.) then
              iunit = 77   ! output file should always be opened
              inquire (98, opened = op)
              if(op) iunit = 98
          write(iunit,"(' APU negative in INFILT/s: '/a6,i3,(5g13.4))") 
     &           trace, i, gef, thi, thi2, thm2(i), gef 
              stop ' negative value of cap. drive in setg for this
     &element '
            end if
          end if
        end if
      end if

      return

      end
C------------------------------------------------------------------------------


      subroutine errly (hc, ferr, dferr)

C   External function for iteratively finding subcrust head and final f
C   for flow thru Brooks-Corey crust.  Adapted for use by KINEROS2.

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      logical :: twolay, dbg

      integer i
C  values of H and PB are presumed both positive transformed

      hct = hc
      itrr = 0
      dh = abs(hct / 50.)
      vsum = 0.
      ssum = 0.
      hn = hct

      rk2 = rkohy (pb2(i), al2(i), hct, 1, drkdh)

      dkdh = drkdh * sk2(i)
      vkh2 = rk2 * sk2(i)
      vmu = vkh2

      rk1c = rkohy (pb1(i), al1(i), hct, 0, dum)

      vkh1c = rk1c * sk1(i)

10    continue

      ho = hn
      hn = hn - dh

      if (hn .lt. 0.) hn = 0.

      hb = 0.5 * (hn + ho)

      rk1 = rkohy (pb1(i), al1(i), hb, 0, dum)

      vkh1 = sk1(i) * rk1
      term = dh * vkh1 / (vkh2 - vkh1)
      sterm = term / (vkh2 - vkh1)

      if (term .le. 0.) then

        ferr = 10. * hc
        dferr = 50.

        return

      end if

      vsum = vsum + term
      ssum = ssum + sterm
      itrr = itrr + 1

      if (hn .gt. 0.01) go to 10

      ferr = vsum - zfl(i)
      dferr =  vkh1c / (vkh2 - vkh1c) - dkdh * ssum
      vmu = vkh2

      return

      end
C------------------------------------------------------------------------------
C
      function fgex (I, g, a)
C  dimensionless infiltrability as fn of cum I, 3-par model, no CV.
C  This is the upper envelope of all CV-affected infiltrability curves.
      real I,Is
C
      Is = I/g !**
      if( Is .le. 1.e-8) then
        vpart = 1.e8
      else if(Is .le. 1.e-5) then
        vpart = 1./Is
      else !**
        vpart = a / (exp (a * Is) - 1.)
      end if
      
      fgex = 1. + vpart

      return

      end
C!!-------------------------------------------------------------new 12/02
      function favdI(Apu, AK, Iv, dt)
C finds increment of potential infil depth in interval dt
      real :: Iv, Is1, Is2, Ino, Isz, gam = 0.8
C Iv could be intent IN or returned as new value of cumI      
C
C  get scaling parameters
      Ino = Apu
      Tn = Ino/AK
C  scale Iv
      If( Ino .le. 1.1e-6) then
        write(77,'(3g13.3)') Iv,apu
C        stop ' invalid arg. in favDI '
        FavdI = AK*dt
        return
      end if
      Is1 = Iv/Ino
      if(Is1 .gt. 100.) then
        FavdI = AK*dt
        return
      end if
      ts1 = (Is1 -log((exp(gam*Is1) - 1. + gam) / gam)) / (1.- gam) ! 3-par funct
      Isz = fIsots(ts1,gam)
C  step scaled time
      ts2 = ts1 + dt/Tn
C  approximate decay function
      Is2 = fIsots(ts2, gam)
C  rescale deltaI:
      favdI = (Is2 - Isz)*Ino
C
      return
      end
C * --------------------------------------------------------------------
      function fIsots(ts, gam)
C
C   explicit Cum function of Smith, with omega interpreted from gam. 
C   ts is scaled time, infiltrability from ponding 
C   gam is three-par. weighting function: 1 = S-P, 0 = G-A
C
      om = 1./(8.- 6.*gam*gam)
C
      if(ts .le. 0.) then
        fIsots = 0.
      else
        t8 = ts/8.
        tb = 2.*ts
        term1 = ts*(0.75 - om)
        term2 = sqrt(ts*ts + 8.*ts)/4.
        term3 = sqrt(t8) + sqrt(t8+1.)
        term4 = sqrt(om*om*ts*ts + ts/2.)
        term5 = om*sqrt(tb) + sqrt(om*om*tb + 1.)
        fIsots = term1 + term2 - 2.*log(term3) + term4 + 
     &           0.5*log(term5)/om
      end if
      return
      end
C------------------------------------------------------------------------------


      function vcurf (rs, cov)

C   Finds infilt curve coeff, cc, as a function of Cov(Ksat) and R*


      if (cov .lt. 0.1) then

        vcurf = 21.

        return

      else if (rs .ge. 0.1) then

        rsf = 1. - exp (-0.6 * (rs - .1))
        cvf = 0.75 / cov**1.3

      else

        rsf = 0.
        cvf = 1.

      end if

      vcurf = 1.+ rsf * cvf

      if (vcurf .lt. 1.) stop ' error - infilt: curve coeff. < 1'

      return

      end
C------------------------------------------------------------------------------


      function vfm (rs, cov)

C   Finds dimensionless reduction factor to apply to ks to get areal fmin


      if (rs .lt. 0.2) then

         vfm = 1.

      else if (cov .lt. 0.1) then
C
          vfm = 1.0

      else

        p = 1.8 / cov**.85
        vfm = (1. + (1. / rs)**p)**(-1. / p)

      end if

      return

      end
C------------------------------------------------------------------------------


      subroutine fcapv (Ic, fu, ifl, dfdi)

C   Finds infiltration capacity fu for randomly heterogeneous area, and for
C   flag IFL true, finds derivative w/r cumulative infilt Ic

C   Parameters:
C  in common /layr/:
C     ref     rainrate [L/T]
C     vmu     effective final f 
C     cc      curvature coefficient
C
C     Ic     cumul. infil depth, [L]
C     fu     infil capacity at depth ic [L/T]
C     ifl    flag signaling the requirement for a derivative value
C     dfdi   derivative of fv w/r Ic

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &       pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

      logical :: twolay, dbg, ifl

      integer i, idl

      real ic, fu

C                                                               define function
      fexi (al, fu, ag) = exp (al * fu / ag)

      if( vmu .le. 0.) then !**
        fu = 0. !**
        dfdi = 0. !**
        return !**
      end if !**

      res = ref / vmu

      if (Ic .le. 1.E-6 .and. cv(i) .lt. 0.1) then

        fu = vmu * fgex(Ic, apu, alf)
        if(ifl) then
          fe = fexi(alf, Ic, apu)
          dfdi = -vmu*alf*alf*fe/apu/(fe-1.)**2
        end if

      else

        if (res .le. 1.) then !**     ! use .le. rather than .lt.

          fu = ref
          dfdi = 0.

        else

          rat = Ic / apu
          cc = vcurf(res, cv(i))

          if (rat * alf .ge. 38.) then

            fs = 1.
            fe = -10.
            d1 = 0.

          else

            fe = fexi (alf, Ic, apu)

            rsm = res - 1.

            if (ifl)  d1 = vmu * rsm**2 / apu

            core = fe - 1.

            if (core .le. 0.) then      ! new test

              fs = 1. + rsm

            else

              vu = rsm * core / alf
              test = cc * alog (vu)

              if (test .lt. 39. .and. cc .lt. 20.) then

                if (test .gt. -40) then

                  fs = 1. + rsm * (1.+ vu**cc)**(-1. / cc)

                  if (ifl) d1 = d1 * vu**(cc - 1.) /
     &                          (1. + vu**cc)**(1. + 1. / cc)

                else

                  fs = 1. + rsm
                  d1 = 0.

                end if

              else  !  insignificant CV:

                fs = 1. + alf / (fe - 1.)
                d1 = vmu*alf*alf/apu/(fe - 1.)**2

              end if
            end if
          end if

          fu = fs * vmu
          if(ifl) dfdi = -d1 * fe

C         write(99,'(" test ",7g13.4)') Ic, fu, dfdi, core, apu, vu
        end if
      end if

      return

      end
