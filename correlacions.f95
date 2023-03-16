program correlacions
!Programa que calcula fdr e as correlaci�ns dos programas anteriores
!Carlos L�pez Bueno

use define_precision
use variables_comunes

implicit none

!Definimos os par�metros para as correlaci�ns
integer(kind=enteiro),parameter::korixe=10,max_n=5000,n_tau=300,nmax_fdr=1000
real(kind=doblep),parameter::dr=0.01d00
!korixe:cada cantos pasos toma a orixe
!max_n: n�mero m�ximo de execuci�ns para as que est� programado
!n_tau: n�mero m�ximo de pasos en tempo
!nmax_fdr: n�mero m�ximo de valores da funci�n de distribuci�n radial

!defino as demais variables
integer(kind=enteiro)::np,ktotal,kpaso,k,n,i,j,l,m,kfdr,k_tau
!np: n�mero de part�culas na rede
!k::variable empregada para o conteo do bucle
!ktotal:n�mero de pasos totais
!kpaso:valor que determina cada cantos pasos gardamos os valores obtidos
!n: n�mero de veces que se executa (ktotal/kpaso)
!i:variable para as definici�ns dos vectores rfdr e t_tau (impl�cito)
!j,l,m:variables para os bucles

real(kind=doblep),dimension(max_n,npmax)::rx,ry,rz,vx,vy,vz
!Posici�ns e velocidades
real(kind=doblep),dimension(n_tau+1):: t_tau,desplcm,corr_v
!
real(kind=doblep),dimension(nmax_fdr+1)::rfdr,fdr
!
real(kind=doblep):: dis_max,factor,gas_ideal
!
real(kind=doblep)::epot,ecin,et,deltat
!epot,ecin,et: enerx�as do sistema
!
real(kind=doblep)::rrx,rry,rrz,rijx,rijy,rijz,dis,riix,riiy,riiz
!posici�ns das part�culas (escalares), escalares do vector rij e distancia entre part�culas

character(len=25):: name1,name2,name3,name4,name5
!name1: Ficheiro cos datos da simulaci�n
!name2: Ficheiro coa configuraci�n do sistema
!name3: Ficheiro cos datos de rv
!name4: Ficheiro cos datos de fdr
!name5: Ficheiro cos datos das correlaci�ns

integer (kind=enteiro)::rv
!Valor empregado para ler os datos (access direct)
rv=6*8*npmax+8



write(*,*) 'Escriba o nome do ficheiro onde se atopan os datos da simulacion:'
read(*,9000) name1

9000	format(a25)

!Lemos os datos que precisamos da simulaci�n
open(90,file=name1)
read(90,*) np,pl,pli,rc,rc2
read(90,*) vol, dens,ktotal,kpaso,dt
read(90,9000) name2
read(90,9000) name2
close(90)

!Calculamos os diferentes valores necesarios para a simulaci�n
n=ktotal/kpaso
deltat=dble(kpaso)*dt
dis_max=pl/2.d00-dr/2.d00
factor=4.d00*pi*dens/3.d00
rfdr=(/(dble(i)*dr,i=0,nmax_fdr)/)
t_tau=(/(dble(i)*deltat,i=0,n_tau)/)

!Pido o nome do ficheiro onde se atopan os datos de rv
write(*,*) 'Indique o nome do ficheiro onde se almacenan os datos de rv'
read(*,9000) name3

!E leo os datos deste
open(91,file=name3,form='unformatted',access='direct',recl=rv)
do i=1,n
  read(91,rec=i) rx(i,:),ry(i,:),rz(i,:),vx(i,:),vy(i,:),vz(i,:)
end do 
close(91)

!po�emos a cero os acumuladores
desplcm=0.d00
corr_v=0.d00
fdr=0.d00

!E comezamos cos c�lculos. En primeiro lugar obte�amos fdr

do l=1,n !Eleximos cada paso gravado
 do i=1,np-1
  !Por optimizaci�n temporal extraemos os valores precisos dos vectores (todo � similar � subrutina potencial)
  rrx=rx(l,i)
  rry=ry(l,i)
  rrz=rz(l,i)
  
  do j=i+1,np
    !calculamos as distancias
    rijx=rrx-rx(l,j)
    rijy=rry-ry(l,j)
    rijz=rrz-rz(l,j)
    
    !E reducimos a distancia para as imaxes m�is achegadas
    rijx=rijx-pl*dnint(rijx*pli)
    rijy=rijy-pl*dnint(rijy*pli)
    rijz=rijz-pl*dnint(rijz*pli)

    !finalmente calculamos a distancia
    dis=dsqrt(rijx*rijx+rijy*rijy+rijz*rijz)

    !Tan s� calculamos ata L/2 da caixa
    if (dis<dis_max) then
      kfdr=nint(dis/dr)
      fdr(kfdr)=fdr(kfdr)+2.d00
     end if
   end do 
 end do 
end do 

!Calculamos o valor medio
fdr=fdr/(dble(n)*dble(np))

!Dividimos polo sistema de gas ideal
do i=1,nmax_fdr
gas_ideal=factor*((rfdr(i)+dr/2.d00)**3-(rfdr(i)-dr/2.d00)**3)
fdr(i)=fdr(i)/gas_ideal
end do 

!Pido o nome do ficheiro onde quero gardar as fdr
write(*,*) 'Indique o nome do ficheiro onde desexa gardar os datos de fdr'
read(*,9000) name4

!E finalmente gravo os valores obidos
open(92,file=name4)
do i=0,nmax_fdr
  if(fdr(i)>=1.d-12) then
    write(10,9001) rfdr(i),fdr(i)
  end if
end do 
close(92)

9001	format (1pe13.6,2x,e13.6)

!Pasamos ao c�lculo das correlaci�ns, empregando as f�rmulas obtidas en teor�a

k=0
do l=1,n-n_tau,korixe
  do m=l,l+n_tau
    k_tau=m-l
    do i=1,np
      riix=rx(m,i)-rx(l,i)
      riiy=ry(m,i)-ry(l,i)
      riiz=rz(m,i)-rz(l,i)
      desplcm(k_tau)=desplcm(k_tau)+riix*riix+riiy*riiy+riiz*riiz
      corr_v(k_tau)=corr_v(k_tau)+vx(m,i)*vx(l,i)+vy(m,i)*vy(l,i)+vz(m,i)*vz(l,i)
    end do 
  end do
  k=k+1 
end do 

!Calculamos valores medios
desplcm=desplcm/(dble(np)*dble(k))
corr_v=corr_v/(dble(np)*dble(k))

!Pedimos o nome do ficheiro onde gravar os datos das correlaci�ns
write(*,*) 'Indique o nome do ficheiro onde desexa gardar os datos das correlacions:'
read(*,9000) name5

!Xa para rematar, gravamos os datos
open(93,file=name5)
write(93,9002) (t_tau(i),desplcm(i),corr_v(i),i=0,n_tau)
close(93)

9002	format(1pe13.6,2x,e13.6,2x,e13.6)


!E damos o programa por finalizado

stop
end correlacions