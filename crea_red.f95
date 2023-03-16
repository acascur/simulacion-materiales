program crea_red

!Programa que determina a configuración inicial dun sistema de partículas con estructura fcc no conxunto microcanónico.
!Autor: Carlos López Bueno

use variables_comunes !Usa o módulo variables comunes
use define_precision !Usa o módulo define precisión
implicit none

integer(kind=enteiro):: i,j,k,np,idum
!define os contadores dos bucles (i,j,k) e o número de partículas da caixa.
!idum e o valor que debemos introducir para iniciar a función random

real(kind=doblep), dimension(npmax):: rx,ry,rz,vx,vy,vz,ax,ay,az
!vectores posición, velocidade e aceleración das partículas

real(kind=doblep):: dnpmax, et, ecq,epot,dfiv,d2fiv,random,ecin,px,py,pz,fixE
!dnpmax: npmax convertido a real en dobre precisión
!et: enerxía total
!ecq: enerxía cinética correxida
!epot: enerxía potencial
!dfiv: derivada da enerxía potencial
!d2fiv: derivada sergunda da enerxía potencial
!random: número aleatorio entre 0 e 1
!ecin: enerxía cinética calculada a partir das velocidades
!px,py,pz: momento total do sistema de partículas
!fixE:factor empregado para a corrección das velocidades fixando a enerxía total do sistema

character(len=25):: name1, name2
!fichero con los datos de la simulación (name1) y con los datos obtenidos de posiciones, velocidades y aceleraciones (name2)

write(*,*) "Indique a densidade:"
read(*,*)dens
write(*,*) "Indique o radio de corte do potencial:"
read(*,*)rc
write(*,*)"Indique o nome do ficheiro no que desexa gardar os datos da simulacion (.dat):"
read(*,9000)name1
write(*,*)"Indique o nome do ficheiro no que desexa gardar a configuracion do sistema (.dat):"
read(*,9000)name2

9000	format(a25) !Especifica o formato para as cadenas de caractéres

!calculamos algunhas variables que empregaremos posteriormente
dnpmax=dble(npmax) 	!Convirte npmax en real de dobre precisión
vol=dnpmax/dens   	
pl=(vol)**(1.d00/3.d00) 
pli=1.d00/pl   	
rc2=rc*rc    		
dl=pl/dble(numk)	
dl12=dl/2.d00	

np=0 !Iniciamos o contador de partículas a cero para así levar conta do número de partículas que temos colocadas
 do i=0,numk-1
   do j=0,numk-1
     do k=0,numk-1					!Cada un deste bucles coloca as partículas nunha red cúbica simple. Combinando 4 con orixes diferentes, obtemos a rede fcc desexada.
       np=np+1						
       rx(np)=0.d00+dble(i)*dl		
       ry(np)=0.d00+dble(j)*dl		
       rz(np)=0.d00+dble(k)*dl
     end do
   end do
 end do
 do i=0,numk-1
   do j=0,numk-1
     do k=0,numk-1
       np=np+1
       rx(np)=0.d00+dble(i)*dl
       ry(np)=dl12+dble(j)*dl
       rz(np)=dl12+dble(k)*dl
     end do
   end do
 end do
  do i=0,numk-1
   do j=0,numk-1
     do k=0,numk-1
       np=np+1
       rx(np)=dl12+dble(i)*dl
       ry(np)=0.d00+dble(j)*dl
       rz(np)=dl12+dble(k)*dl
     end do
   end do
 end do
   do i=0,numk-1
   do j=0,numk-1
     do k=0,numk-1
       np=np+1
       rx(np)=dl12+dble(i)*dl
       ry(np)=dl12+dble(j)*dl
       rz(np)=0.d00+dble(k)*dl
     end do
   end do
 end do
 write(*,*) "Numero de particulas colocadas na rede fcc:",np
 

call potencial(np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)
!Chamamos a subrutina potencial que, aportandolle a posición das partículas que acabamos de obter e o número de partículas, danos a aceleración e o potencial total, así comoa as súas derivadas.

!Agora imos asignar as velocidades de cada unha das partículas de xeito aleatorio empregando a función random aportada polo profesor
!Esta función precisa que lle asignemos un valor inicial, o que nos queiramos sempre que sexa enteiro (idum)

write(*,*) "Teclee un numero enteiro calquera"
read(*,*)idum

!Como queremos velocidades (reducidas) comprendidas entre -1 e 1 e a función random danos valores entre 0 e 1, multiplicamos por dous e restamos un a este valor, asegurando así que esté comprendido neste intervalo
do i=1,np
    	vx(i)=2.d00*random(idum)-1.d00
        vy(i)=2.d00*random(idum)-1.d00  !Con este bucle asignamos as tres compoñentes da velocidade a cada partícula
        vz(i)=2.d00*random(idum)-1.d00
end do

!Calculamos agora os momentos totales (masa unidade).

px=sum(vx)
py=sum(vy)
pz=sum(vz)

write(*,*)"Momentos antes da correcion:",px,py,pz !Escribimos os momentos en pantalla, observando que estes non son cero

!Dado que sabemos que o momento total debe ser nulo, pero non temos a certeza de que isto pase ao xenerar os números aleatorios xa que temos un número finito destes, empregamos a seguinte correción nas velocidades

px=px/dnpmax
py=py/dnpmax
pz=pz/dnpmax
    	
vx=vx-px
vy=vy-py
vz=vz-pz

!Trivialmente vese que agora o momento total si deberá ser exactamente cero

px=sum(vx)
py=sum(vy)
pz=sum(vz)

write(*,*) "Momentos despois da correccion:",px,py,pz !Escribimos de novo os momentos para comprobar que agora sí son nulos.

!Calculamos a enerxía cinética do xeito habitual
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!A enerxía total neste momento será
write(*,*) "Enerxia total do sistema:",ecin+epot

!Agora, como estamos no microcanónico debemos fixar a enerxía total do sistema. Evidentemente de partida esta non será a que temos neste momento, senon que debemos variar a velocidade das partículas para que asi a suma da cinética e a potencial tome o valor desexado
write(*,*) "Indique o valor desexado da enerxia do sistema"
read(*,*) et

!Fixamos o valor da enerxía cinética
ecq=et-epot

!Para asegurarnos, poñemos un aviso en caso de que esta enerxía sexa negativa, cousa que nunca debería acontecer
if (ecq<=0.d00) then
      write (*,*)"Error: enerxia cinetica negativa"
      stop
end if

!Agora calculamos o factor polo que deberemos multiplicar as velocidades para que a enerxía cinética sexa a desexada
fixE=dsqrt(ecq/ecin)
!Comprobase de xeito trivial que isto é así
vx=fixE*vx
vy=fixE*vy
vz=fixE*vz
!Posto que simplemente aparece o mesmo factor en todas as compoñentes da velocidade, o momento total seguirá sendo nulo.

!Calculamos agora a enerxía cinética, que loxicamente deberá coincidir con ecq
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!Finalmente, volcamos os datos obtidos aos diferentes ficheiros
open(90,file=name1)
write(90,9001) np,pl,pli,rc,rc2
write(90,9002) vol,dens
write(90,9003) ecin+epot,ecin,epot
write(90,9000) name1
write(90,9000) name2
close(90)
!Este ficheiro estará en ASCII

!Debido á gran cantidade de datos de posicións, velocidades e aceleracións dos que dispoñemos, grabamos un ficheiro binario (sen formato) con estes)

open(91,file=name2,form='unformatted')
write(91) rx,ry,rz,vx,vy,vz,ax,ay,az
close(91)
   
!Por último especificamos os formatos empregados (o 9000 xa foi definido con anterioridade)
9001	format(i4,2x,1pe19.12,3(2x,e19.12))
9002	format(1pe19.12,2x,e19.12)
9003	format(1pe19.12,2x,e19.12,2x,e19.12)



end program crea_red