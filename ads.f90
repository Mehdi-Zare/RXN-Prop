program deltaG
implicit none
integer  ,  parameter                   ::      dp = selected_real_kind(15, 307)
integer                                 ::      error_flag                                !for qvib subroutine checking
integer                                 ::      alloc_err                                 !Deallocation status for allocatable arrays
integer                                 ::      s                                         !The number of reactants, products, or transition states        
integer                                 ::      i                                         !variable for reactants, products and TS loops     
real(dp)                                ::      freqtot, zpeGas
integer                                 ::      ierror,nfreq
integer                                 ::      k                                         !variable for Temperature loop
integer                                 ::      nTemp                                     !number of different Temperatures
character(20)                           ::      kind                                      !Reaction Kind: adsorption/desorption or surface reaction
character(20)                           ::      filename,lnq                              !The variable for the name of input file
real(dp) ,  allocatable, dimension (:)  ::      freq, qvibr,lnqvibr                       !frequencies, vibrational partition function of reactants, products
real(dp)                                ::      part,qvibp,lnqvibp,ZPEp,ZPCp,scfp
real(dp) ,  allocatable, dimension (:)  ::      scfr                                      !SCF energy of reactants, products & TS
real(dp) ,  allocatable, dimension (:)  ::      ZPCr,ZPEr                                 !zero point corrected energy/ zp energy of reactants.product & TS
real(dp)                                ::      Na,kbSI,M,N0,h, kb, c,qvib,ZPE,scf,LT,HT,IT,Temp
real(dp) ,  allocatable, dimension (:)  ::      T                                         !Temperature array
real(dp) ,  allocatable, dimension (:)  ::      qtotr1                                    !gas phase partition function
real(dp) ,  allocatable, dimension (:)  ::      delE,delG
real(dp) ,  allocatable, dimension (:)  ::      Ertot,Eptot,Etstot                        !total energy of reactants,products, and transition state
real(dp) ,  allocatable, dimension (:)  ::      qvibrtot,qvibptot,qvibtstot               !total vibrational frequency of reactants,product and TS
real(dp) ,  allocatable, dimension (:)  ::      qelec
real(dp) ,  allocatable, dimension (:)  ::      Keq,kf,kr
real(dp)                                ::      cutoff                                    !cutoff value for frequencies

Na      = 6.02214076e23
kbSI    = 1.380649e-23
N0      = 1.866010378e19
h       = 4.13566845e-15                                                                       !in eV/K
kb      = 8.61733502e-5                                                                        !in eV/K
c       = 299792458                                                                            !light speed m/s                                                                              

!print the information about this code
call help()


!get Temperatures
write(*,*) 'please insert the molecular weight of your adsorbate in g/mol'
read(*,*) M

if ( M < 0)  then                    !FLAG!
    write(*,*) " ERROR: Molecular weight cannot be negative!!!"
    write(*,*) " The program is going to exit"
    call EXIT(0)
end if


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



nTemp=((HT-LT)/IT)+1
allocate (T(nTemp))
do i=1,nTemp
   T(i)=LT+(i-1)*IT
end do

write(*,*) " Please enter you cutoff value for frequencies, we usually use 50 or 100 "
read (*,*) cutoff

!allocatable arrays for calculation in different T
allocate(delE(nTemp))
allocate(Ertot(nTemp))
allocate(Eptot(nTemp))
allocate(qtotr1(nTemp))
allocate(qvibrtot(nTemp))
allocate(qvibptot(nTemp))
allocate(qelec(nTemp))
allocate(Keq(nTemp))
allocate(kf(nTemp))
allocate(kr(nTemp))
allocate(delG(nTemp))

do k = 1, nTemp
        open (unit=990, file = 'lnq', status='old', action = 'read', iostat = ierror)
        read(990,*,iostat = ierror) part
         if (ierror /=0 .or. part == 0) then                  !!!!FALG
           write(*,*) " Problem with reading file: lnq "
           write(*,*) " probable causes: misspell its name, #of Temperature,&
                        greater number of partition functions you have in your lnq file, it is empty, you do not,&
                        have that file at all"
           write(*,*) " The program is going to exit"
           call EXIT(0)
        end if
          qtotr1(k) = exp(part)
end do
  close(990)

!loop over all calculations for different temperature
do k = 1,nTemp

!Allocatable arrays for reactants
s=2
allocate(qvibr(s))
allocate(lnqvibr(s))
allocate(scfr(s))
allocate(ZPEr(s))
allocate(ZPCr(s))   

!for Gas phase Reactant

!read the # of freq in freq.dat file
nfreq=-1        ! The first line is Escf 
open(1, file='Reactant-1')
do
 read (1, *, end=10)
 nfreq=nfreq+1
end do
10 close(1)

allocate (freq(nfreq))
   freqtot      =       0.0
   open (unit=99, file = 'Reactant-1', status='old', action = 'read', iostat = ierror)
   read(99,*,iostat = ierror) scfr(1)
        if (ierror /=0 .or. scfr(1) == 0) then                  !!!!FALG
           write(*,*) " Problem with reading file: Reactant-1 "
           write(*,*) " You probably misspelled its name or it is empty"
           write(*,*) " The program is going to exit"
           call EXIT(0)
        end if
          
        do i=1,nfreq
                read(99,*, iostat = ierror) freq(i)
                if (ierror /=0) exit
                if (freq(i) == 0) exit
                if (freq(i) < cutoff) then
                        freq(i) = cutoff
                end if
                freqtot=freqtot+freq(i)
        end do
        ZPEr(1)=0.5*h*100*c*freqtot
        zpeGas=ZPEr(1)    !this is just for printing the ZPE of Gas phase Reactant
        ZPCr(1)=scfr(1)+ZPEr(1)
        qvibr(1)=qtotr1(k)   ! this is the total partion function for gas phase
 
 close(99)


deallocate(freq, stat= alloc_err)


! for the second Reactant
Temp = T(k)
!for reactant #2 (free site)
        i=2
                write(filename,100) i
                100 format('Reactant-', I1)
                call q(filename,Temp,qvib,ZPE,scf)
                if ( error_flag == 1  ) exit           !FLAG
                qvibr(i) = qvib
                scfr(i)  = scf
                ZPEr(i)  = ZPE
                ZPCr(i)  = scfr(i)+ZPEr(i)
                lnqvibr(i)= log(qvibr(i))
        
Ertot(k)           =  SUM(ZPCr)
qvibrtot(k)        =  PRODUCT(qvibr)
deallocate(qvibr,lnqvibr,scfr,ZPEr,ZPCr, stat = alloc_err)
   
!for product
      i = 1
                write(filename,120) i
                120 format('Product-', I1)
                call q(filename,Temp,qvib,ZPE,scf)
                if ( error_flag == 1  ) exit           !FLAG
                qvibp = qvib
                scfp  = scf
                ZPEp  = ZPE
                ZPCp  = scfp+ZPEp
                lnqvibp= log(qvibp)
        
Eptot(k)            =  ZPCp
qvibptot(k)         =  qvibp

!calculate Rxn Properties
delE(k)    =       Eptot(k)   -  Ertot(k)
qelec(k)   =       exp ( - delE(k)   / kb / Temp )
Keq(k)     =       qelec(k)   *  qvibptot(k)     /  qvibrtot(k)
kf (k)     =       100000*sqrt(Na)/N0/sqrt(2*3.1415927*M*kbSI*Temp)
kr (k)     =       kf(k)/Keq(k)
delG(k)    =      -kb      *  Temp    *  log(Keq(k))

end do
!end of the loop for Temparature


!Wrtie the data in file named 'Rxn-Energy & Kinetics'
open    ( unit = 500, file = 'Rxn-Energy', status = 'new')
write    (500, 400)
400 format ( 10x, "T           :     Temperature in Kelvin",/, 10x "deltaE     :     Reaction energy in eV", &
             /,10x,"deltaG     :     Reaction Free Energy in eV", &
             /,10x,"Keq        :     Equilibrium constant", &
             /,10x,"Kf         :     Forward Rxn Constant, 1/s", &
             /,10x,"kr         :     Reverse Rxn constatn, 1/s",///)

write   ( 500, 501)
501 format ( 3x, "T(K)", 11x, "deltaE", 16x,"deltaG",&
             20x, "Keq",20x, "kf", 20x, "kr")
write   (500,502)
502 format (2x, "=============================================",&
                "================================================", &
            "===================================================")
do k=1,nTemp
write   (500, 503) T(k), delE(k), delG(k), Keq(k), kf(k), kr(k)
503 format (2x, F6.2, 3x , 6(ES20.13,3x))
end do
write (500,504)
504 format (2x, "=============================================",&
                "================================================", &
            "===================================================")
write (500,*) "ZPE_GasReactant=  " , zpeGas
close(500)

!deallocation of all arrays
deallocate(delE,delG,Ertot,Eptot,qvibrtot,qvibptot,stat = alloc_err)
deallocate(qelec,Keq,kf,kr, stat = alloc_err)
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


!read the # of freq in freq.dat file
nfreq=-1        ! The first line is Escf
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
                if (ierror /=0) exit
                if (freq(i) == 0) exit
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
   950 format(//,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" ,&
               /,"! This code is written  by Mehdi Zare in 6/4/2018 for calculating reaction properties of!" ,&
               /,"!                             adsorption / desorption Reaction                          !" ,&
               /,"! It is written for reaction:     A(g)   +   *   <==>   A*                              !" ,&
               /,"! You need to have four input files for this code:                                      !" ,&
               /,"!             Reactant-1, Reactant-2, Product-1: The first line of these                !" ,&
               /,"! three files is SCF energy and the rest are frequencies.                               !" ,&
               /,"!      HINT: DO NOT FORGET to remove required frequencies of gas phase (Reactant-1)     !" ,&
               /,"! And you also need lnq file which contains partition function for Reactant-1 at        !" ,&
               /,"! different temperatures.                                                               !" ,&
               /,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//)
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
