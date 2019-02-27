program deltaG
implicit none
integer  ,  parameter                   ::      dp = selected_real_kind(15, 307)
integer                                 ::      error_flag                                !for qvib subroutine checking
integer                                 ::      alloc_err                                 !Deallocation status for allocatable arrays
integer                                 ::      s                                         !The number of reactants, products, or transition states        
integer                                 ::      i                                         !variable for reactants, products and TS loops     
integer                                 ::      k                                         !variable for Temperature loop
integer                                 ::      nTemp                                     !number of different Temperatures
character(20)                           ::      filename                                  !The variable for the name of input file
real(dp) ,  allocatable, dimension (:)  ::      freq, qvibr,lnqvibr,qvibp,lnqvibp         !frequencies, vibrational partition function of reactants, products
real(dp) ,  allocatable, dimension (:)  ::      qvibts,lnqvibts                           !for transition state
real(dp) ,  allocatable, dimension (:)  ::      scfr,scfp,scfts                           !SCF energy of reactants, products & TS
real(dp) ,  allocatable, dimension (:)  ::      ZPCr,ZPEr,ZPCp,ZPEp,ZPCts,ZPEts           !zero point corrected energy/ zp energy of reactants.product & TS
real(dp)                                ::      h, kb, c,qvib,ZPE,scf,LT,HT,IT,Temp
real(dp) ,  allocatable, dimension (:)  ::      T                                         !Temperature array
real(dp) ,  allocatable, dimension (:)  ::      delE, delEF, delER                       !reaction energy, forward and reverse activation barrier
real(dp) ,  allocatable, dimension (:)  ::      Ertot,Eptot,Etstot                        !total energy of reactants,products, and transition state
real(dp) ,  allocatable, dimension (:)  ::      qvibrtot,qvibptot,qvibtstot               !total vibrational frequency of reactants,product and TS
real(dp) ,  allocatable, dimension (:)  ::      qelec, qelecF, qelecR                     !Reaction electronci partition function, forward and reverse
real(dp) ,  allocatable, dimension (:)  ::      Keq,kf,kr                                 !equilibrium constant, forward and reverse rate constant
real(dp) ,  allocatable, dimension (:)  ::      delG,delGF,delGR                          !free energy of reaction and forward and reverse free energy barrier
real(dp) ,  allocatable, dimension (:)  ::      Gr1,Gr2,Gp1,Gp2,Gts1,Gts2                 !Gibbs free energy for all materials
real(dp) ,  allocatable, dimension (:)  ::      qvibr1,qvibr2,qvibp1,qvibp2,qvibts1,qvibts2 ! vibrational partition functions
real(dp)                                ::      cutoff                                    !cutoff value for frequencies

h  = 4.13566845e-15                                                                       !in eV/K
kb = 8.61733502e-5                                                                        !in eV/K
c  = 299792458                                                                            !light speed m/s                                                                              
!print the information about this code
call help()

!get the number of Reactants, Products, or Transition States
write(*,*) ' Please type the number of Species (each of Reactants, Products or, &
            TS you have:'
read (*,*) s

!get Temperatures
write(*,*) 'Please type the temperature range you wish in the following order:'
write(*,*) 'Lower value of T. Higher value of T. Temperature interval'
read (*,*) LT, HT, IT

if (LT > HT .or. LT < 0 .or. HT < 0 .or. IT < 0)  then                    !FLAG!
    write(*,*) " ERROR: Your lower T is higher than your higher T or one of , &
                 LT, HT, or IT is negative!!!" 
    write(*,*) " The program is going to exit"
    call EXIT(0)
end if


nTemp=((HT-LT)/IT)+1

if (nTemp == 0)  then                    !FLAG!
    write(*,*) " ERROR: Your have defined no Temperrature, nothing to do!!!"
    write(*,*) " The program is going to exit"
    call EXIT(0)
end if



allocate (T(nTemp))
do i=1,nTemp
   T(i)=LT+(i-1)*IT
end do

write(*,*) " Please enter you cutoff value for frequencies, we usually use 50 or 100 "
read (*,*) cutoff

!allocatable arrays for calculation in different T
allocate(delE(nTemp))
allocate(delEF(nTemp))
allocate(delER(nTemp))
allocate(Ertot(nTemp))
allocate(Eptot(nTemp))
allocate(Etstot(nTemp))
allocate(qvibrtot(nTemp))
allocate(qvibptot(nTemp))
allocate(qvibtstot(nTemp))
allocate(qelec(nTemp))
allocate(qelecF(nTemp))
allocate(qelecR(nTemp))
allocate(Keq(nTemp))
allocate(Kf(nTemp))
allocate(Kr(nTemp))
allocate(delG(nTemp))
allocate(delGF(nTemp))
allocate(delGR(nTemp))
allocate(Gr1(nTemp))
allocate(Gr2(nTemp))
allocate(Gp1(nTemp))
allocate(Gp2(nTemp))
allocate(Gts1(nTemp))
allocate(Gts2(nTemp))
allocate(qvibr1(nTemp))
allocate(qvibr2(nTemp))
allocate(qvibp1(nTemp))
allocate(qvibp2(nTemp))
allocate(qvibts1(nTemp))
allocate(qvibts2(nTemp))

!loop over all calculations for different temperature
do k = 1,nTemp

!Allocatable arrays for reactants
allocate(qvibr(s))
allocate(lnqvibr(s))
allocate(scfr(s))
allocate(ZPEr(s))
allocate(ZPCr(s))   

!assign the temperature to variable Temp
        Temp = T(k)
!loop over reactants
        do i = 1,s
                write(filename,100) i
                100 format('Reactant-', I1)
                call q(filename,Temp,qvib,ZPE,scf)
                if ( error_flag == 1  ) exit           !FLAG
                qvibr(i) = qvib
                scfr(i)  = scf
                ZPEr(i)  = ZPE
                ZPCr(i)  = scfr(i)+ZPEr(i)
                lnqvibr(i)= log(qvibr(i))
        end do
Ertot(k)           =  SUM(ZPCr)
qvibrtot(k)        =  PRODUCT(qvibr)
Gr1(k)             =  ZPCr(1)-kb*T(k)*lnqvibr(1)
Gr2(k)             =  ZPCr(2)-kb*T(k)*lnqvibr(2)
qvibr1(k)          =  qvibr(1)
qvibr2(k)          =  qvibr(2)

deallocate(qvibr,lnqvibr,scfr,ZPEr,ZPCr, stat = alloc_err)

!allocatable arrays for products
allocate(qvibp(s))
allocate(lnqvibp(s))
allocate(scfp(s))
allocate(ZPEp(s))
allocate(ZPCp(s))

!loop over products
        do i = 1,s
                write(filename,120) i
                120 format('Product-', I1)
                call q(filename,Temp,qvib,ZPE,scf)
                if ( error_flag == 1  ) exit           !FLAG
                qvibp(i) = qvib
                scfp(i)  = scf
                ZPEp(i)  = ZPE
                ZPCp(i)  = scfp(i)+ZPEp(i)
                lnqvibp(i)= log(qvibp(i))
        end do
Eptot(k)            =  SUM(ZPCp)
qvibptot(k)         =  PRODUCT(qvibp)
Gp1(k)             =  ZPCp(1)-kb*T(k)*lnqvibp(1)
Gp2(k)             =  ZPCp(2)-kb*T(k)*lnqvibp(2)
qvibp1(k)          =  qvibp(1)
qvibp2(k)          =  qvibp(2)

deallocate(qvibp,lnqvibp,scfp,ZPEp,ZPCp, stat = alloc_err)

!allcatable arrays for transition states
allocate(qvibts(s))
allocate(lnqvibts(s))
allocate(scfts(s))
allocate(ZPEts(s))
allocate(ZPCts(s))

!loop over transition states
       do i = 1,s
                write(filename,140) i
                140 format('TS-', I1)
                call q(filename,Temp,qvib,ZPE,scf)
                if ( error_flag == 1  ) exit           !FLAG
                qvibts(i) = qvib
                scfts(i)  = scf
                ZPEts(i)  = ZPE
                ZPCts(i)  = scfts(i)+ZPEts(i)
                lnqvibts(i)= log(qvibts(i))
        end do
Etstot(k)           =  SUM(ZPCts)
qvibtstot(k)        =  PRODUCT(qvibts)
Gts1(k)             =  ZPCts(1)-kb*T(k)*lnqvibts(1)
Gts2(k)             =  ZPCts(2)-kb*T(k)*lnqvibts(2)
qvibts1(k)          =  qvibts(1)
qvibts2(k)          =  qvibts(2)

deallocate(qvibts,lnqvibts,scfts,ZPEts,ZPCts, stat = alloc_err)


!calculate all reaction properties

delE(k)    =       Eptot(k)   -  Ertot(k)
delEF(k)   =       Etstot(k)  -  Ertot(k)
delER(k)   =       Etstot(k)  -  Eptot(k)  
qelec(k)   =       exp ( - delE(k)   / kb / Temp )
qelecF(k)  =       exp ( - delEF(k)  / kb / Temp )
qelecR(k)  =       exp ( - delER(k)  / kb / Temp )
Keq(k)     =       qelec(k)   *  qvibptot(k)     /  qvibrtot(k)
kf(k)      =       kb      *  Temp    *  qelecF(k)       *  qvibtstot(k)    /  qvibrtot(k)  / h
kr(k)      =       kb      *  Temp    *  qelecR(k)       *  qvibtstot(k)    /  qvibptot(k)  / h
delG(k)    =      -kb      *  Temp    *  log(Keq(k))
delGF(k)   =      -kb      *  Temp    *  log( kf(k) * h / (kb * Temp ))
delGR(k)   =      -kb      *  Temp    *  log( kr(k) * h / (kb * Temp ))

end do
!end of the loop for Temparature


!Wrtie the data in files named 'Rxn-Energy & Kinetics & Gibbs'
open    ( unit = 500, file = 'Rxn-Energy', status = 'new')
write    (500, 400)
400 format ( 10x, "T          :     Temperature in Kelvin",/, 10x "deltaE     :     Reaction energy in eV", &
             /,10x,"deltaEF    :     Reaction Forward Energy, Zero Point Corrected, in eV", &
             /,10x,"deltaER    :     Reaction Reverse Energy, Zero Point Corrected, in eV", &
             /,10x,"deltaG     :     Reaction Free Energy in eV", &
             /,10x,"deltaGF    :     Reaction Forwaed Free Energy in eV", &
             /,10x,"deltaGR    :     Reaction Reverse Free Energy in eV",///)

write   ( 500, 501)
501 format ( 3x, "T(K)", 11x, "deltaE", 16x,"deltaEF",&
             16x, "deltaER" , 16x, "deltaG"     ,&
             16x, "deltaGF",16x, "deltaGR")
write   (500,502)
502 format (2x, "=============================================",&
                "================================================", &
            "===================================================")
do k=1,nTemp
write   (500, 503) T(k), delE(k), delEF(k), delER(k), delG(k), delGF(k), delGR(k)
503 format (2x, F6.2, 3x , 6(ES20.13,3x))
end do
close(500)

open    ( unit = 5000, file = 'Kinetics', status = 'new')
write    (5000, 4000)
4000 format ( 10x, "T          :     Temperature in Kelvin", &
             /,10x,"Keq        :     Equilibrium Constant", &
             /,10x,"kf         :     Forward Reaction Constant, 1/s", &
             /,10x,"kr         :     Reverse Reaction Constant, 1/s",///)

write   ( 5000, 5010)
5010 format ( 3x, "T(K)", 11x, "Keq", 16x,"kf",&
             16x, "kr")
write   (5000,5020)
5020 format (2x, "=============================================",&
                "================================================")
do k=1,nTemp
write   (5000, 5030) T(k), Keq(k), kf(k), kr(k)
5030 format (2x, F6.2, 3x , 3(ES20.13,3x))
end do
close(5000)

open    ( unit = 50000, file = 'Gibbs', status = 'new')
write    (50000, 40000)
40000 format ( 10x,"T              :     Temperature in Kelvin", &
             /,10x,"Gr1,Gr2        :     Gibbs Free Energy of Reactan-1 and 2, eV ", &
             /,10x,"Gts1,Gts2      :     Gibbs Free Energy of Product-1 and 2, eV ", &
             /,10x,"Gp1,Gp2        :     Gibbs Free Energy of TS-1 and 2, eV ",///)

write   ( 50000, 50100)
50100 format ( 3x, "T(K)", 11x, "Gr1", 20x, "Gr2", 20x,"Gts1",&
             20x, "Gts2", 20x, "Gp1", 20x, "Gp2")
write   (50000,50200)
50200 format (2x, "=====================================================================",&
                "===========================================================================")
do k=1,nTemp
write   (50000, 50300) T(k), Gr1(k), Gr2(k), Gts1(k), Gts2(k), Gp1(k), Gp2(k)
50300 format (2x, F6.2, 3x, 6(ES20.13,3x))
end do
close(50000)

open    ( unit = 50001, file = 'Qvib', status = 'new')
write   (50001, 40001)
40001 format ( 10x," T                 :     Temperature in Kelvin", &
             /,10x," qvibr1, qvibr2    :     Vibrational partition functions of Reactan-1 and 2, eV ", &
             /,10x," qvibp1, qvibp2    :     vibrational partition functions of Product-1 and 2, eV ", &
             /,10x," qvibts1, qvibts2  :     vibrational partition functions of TS-1 and 2, eV ",///)

write   ( 50001, 50101)
50101 format ( 3x, "T(K)", 11x, "qvibr1",18x, "qvibr2", 15x,"qvibts1",&
             17x, "qvibts2", 17x, "qvibp1", 17x, "qvibp2")
write   (50001,50201)
50201 format (2x, "=====================================================================",&
                "===========================================================================")
do k=1,nTemp
write   (50001, 50301) T(k), qvibr1(k), qvibr2(k), qvibts1(k), qvibts2(k), qvibp1(k), qvibp2(k)
50301 format (2x, F6.2, 3x, 6(ES20.13,3x))
end do
close(50001)


!deallocation of all arrays
deallocate(delE, delEF,delER,Ertot,Eptot,Etstot,qvibrtot,qvibptot,qvibtstot, stat = alloc_err)
deallocate(qelec, qelecF, qelecR,Keq,kf,kr,delG,delGF,delGR, stat = alloc_err)
deallocate(Gr1, Gr2, Gts1, Gts2, Gp1, Gp2, stat = alloc_err)
deallocate(qvibr1,qvibr2,qvibp1,qvibp2,qvibts1,qvibts2, stat = alloc_err)
deallocate(T, stat = alloc_err)

end program deltaG



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Subroutine q!!!!!!!!!!!!!!!!!!!!!!
!This subroutine calculate vibrational partition function and read scf frome file
subroutine q(filename,Temp,qvib,ZPE,scf)

implicit none

integer  ,  parameter                   ::      dp = selected_real_kind(15, 307)
real     ,  parameter                   ::      h  = 4.13566845e-15                   !in eV/K
real     ,  parameter                   ::      kb = 8.61733502e-5                    !in eV/K
real     ,  parameter                   ::      c  = 299792458                        !light speed m/s




!local variable
real(dp) ,  allocatable, dimension (:)  ::      freq
real(dp)                                ::      freqtot, cutoff
integer                                 ::      i,j,ierror,alloc_err,nfreq
integer                                 ::      error_flag

!parameter types and definition
character(20)    ,  intent(in)                  ::      filename 
real(dp)         ,  intent(in)                  ::      Temp 
real(dp)         ,  intent(out)                 ::      scf
real(dp)         ,  intent(out)                 ::      qvib
real(dp)         ,  intent(out)                 ::      ZPE

error_flag = 0

!read the # of freq in freq.dat file
nfreq=-1       ! The first line is Escf
open(1, file=filename)
do
 read (1, *, end=10)
 nfreq=nfreq+1
end do
10 close(1)



allocate (freq(nfreq))
   qvib         =       1.0
   freqtot      =       0.0
   open (unit=99, file = filename, status='old', action = 'read', iostat = ierror)
   read(99,*,iostat = ierror) scf 
       
       if (ierror /=0 .or. scf == 0) then                  !!!!FALG
          write(*,*) " Problem with reading file:  ", filename
          write(*,*) " You probably misspelled its name or it is empty"
          write(*,*) " The program is going to exit"
          error_flag = 1
          call EXIT(0)
        end if
  
         do i=1,nfreq
                read(99,*, iostat = ierror) freq(i)
                if (ierror /=0 .or. freq(i) == 0) exit
                if (freq(i) < cutoff) then
                        freq(i) = cutoff
                end if       
                freqtot=freqtot+freq(i)
                qvib=qvib*(1/(1-exp((-h*freq(i)*100*c)/(kb*Temp))))
        end do
        ZPE=0.5*h*100*c*freqtot
 close(99)
deallocate(freq, stat= alloc_err)
end subroutine q


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Subroutine Help!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine help()
integer         ::      status = 0
character(10)   ::      answer
write (*,950)
    950 format(//," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" ,&
                /," ! This code is written  by Mehdi Zare in 5/20/2018 for calculating      !" ,&
                /," ! reaction properties of Surface Reactions                              !" ,&
                /," ! This code works only for the same number of species of Reactanat,     !" ,&
                /," ! Product, and Transition State                                         !" ,&
                /," ! For Example, if you insert 2 as the number of species, then           !", &
                /," ! You need to have six input files for this code:                       !", &
                /," ! Reactant-1, Reactant-2, Product-1, Product-2, TS-1, TS-2:             !", &
                /," ! The first line of these six files is SCF energy and the rest are      !", &
                /," !                          frequencies                                  !", &
                /," ! HINT: DO NOT FORGET to remove IMAGINARY frequency of transition state !", &
                /," !                           !!!AGAIN!!!                                 !", &
                /," ! This code is written for the same number of species as reactants,     !", &
                /," ! products and transition states.                                       !", &
                /," ! If you have different number of species for them                      !", &
                /," ! you should make necessary changes to this code                        !", &
                /," !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//)
   
250  write (*,960)
        960 format (1X,"Do you still want to use this code? type yes or no",//)
        read (*,*) answer
        if (answer == 'no') then
          write(*,1010)
          1010 format(//,1X,"See You Later!",//)
          call exit (status)
          elseif (answer=='yes') then
          write (*,970)
          970 format (//,1X,"Great! Glad that you can use this code!",//)
          else
          write(*,980)
          980 format (//,1X,"Your answer is not recognized, please try again",//)
          go to 250
        end if                          
   End Subroutine help
