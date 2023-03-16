program crea_red

!Programa que determina a configuraci�n inicial dun sistema de part�culas con estructura fcc no conxunto microcan�nico.
!Autor: Carlos L�pez Bueno

use variables_comunes !Usa o m�dulo variables comunes
use define_precision !Usa o m�dulo define precisi�n
implicit none

integer(kind=enteiro):: i,j,k,np,idum
!define os contadores dos bucles (i,j,k) e o n�mero de part�culas da caixa.
!idum e o valor que debemos introducir para iniciar a funci�n random

real(kind=doblep), dimension(npmax):: rx,ry,rz,vx,vy,vz,ax,ay,az
!vectores posici�n, velocidade e aceleraci�n das part�culas

real(kind=doblep):: dnpmax, et, ecq,epot,dfiv,d2fiv,random,ecin,px,py,pz,fixE
!dnpmax: npmax convertido a real en dobre precisi�n
!et: enerx�a total
!ecq: enerx�a cin�tica correxida
!epot: enerx�a potencial
!dfiv: derivada da enerx�a potencial
!d2fiv: derivada sergunda da enerx�a potencial
!random: n�mero aleatorio entre 0 e 1
!ecin: enerx�a cin�tica calculada a partir das velocidades
!px,py,pz: momento total do sistema de part�culas
!fixE:factor empregado para a correcci�n das velocidades fixando a enerx�a total do sistema

character(len=25):: name1, name2
!fichero con los datos de la simulaci�n (name1) y con los datos obtenidos de posiciones, velocidades y aceleraciones (name2)

write(*,*) "Indique a densidade:"
read(*,*)dens
write(*,*) "Indique o radio de corte do potencial:"
read(*,*)rc
write(*,*)"Indique o nome do ficheiro no que desexa gardar os datos da simulacion (.dat):"
read(*,9000)name1
write(*,*)"Indique o nome do ficheiro no que desexa gardar a configuracion do sistema (.dat):"
read(*,9000)name2

9000	format(a25) !Especifica o formato para as cadenas de caract�res

!calculamos algunhas variables que empregaremos posteriormente
dnpmax=dble(npmax) 	!Convirte npmax en real de dobre precisi�n
vol=dnpmax/dens   	
pl=(vol)**(1.d00/3.d00) 
pli=1.d00/pl   	
rc2=rc*rc    		
dl=pl/dble(numk)	
dl12=dl/2.d00	

np=0 !Iniciamos o contador de part�culas a cero para as� levar conta do n�mero de part�culas que temos colocadas
 do i=0,numk-1
   do j=0,numk-1
     do k=0,numk-1					!Cada un deste bucles coloca as part�culas nunha red c�bica simple. Combinando 4 con orixes diferentes, obtemos a rede fcc desexada.
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
!Chamamos a subrutina potencial que, aportandolle a posici�n das part�culas que acabamos de obter e o n�mero de part�culas, danos a aceleraci�n e o potencial total, as� comoa as s�as derivadas.

!Agora imos asignar as velocidades de cada unha das part�culas de xeito aleatorio empregando a funci�n random aportada polo profesor
!Esta funci�n precisa que lle asignemos un valor inicial, o que nos queiramos sempre que sexa enteiro (idum)

write(*,*) "Teclee un numero enteiro calquera"
read(*,*)idum

!Como queremos velocidades (reducidas) comprendidas entre -1 e 1 e a funci�n random danos valores entre 0 e 1, multiplicamos por dous e restamos un a este valor, asegurando as� que est� comprendido neste intervalo
do i=1,np
    	vx(i)=2.d00*random(idum)-1.d00
        vy(i)=2.d00*random(idum)-1.d00  !Con este bucle asignamos as tres compo�entes da velocidade a cada part�cula
        vz(i)=2.d00*random(idum)-1.d00
end do

!Calculamos agora os momentos totales (masa unidade).

px=sum(vx)
py=sum(vy)
pz=sum(vz)

write(*,*)"Momentos antes da correcion:",px,py,pz !Escribimos os momentos en pantalla, observando que estes non son cero

!Dado que sabemos que o momento total debe ser nulo, pero non temos a certeza de que isto pase ao xenerar os n�meros aleatorios xa que temos un n�mero finito destes, empregamos a seguinte correci�n nas velocidades

px=px/dnpmax
py=py/dnpmax
pz=pz/dnpmax
    	
vx=vx-px
vy=vy-py
vz=vz-pz

!Trivialmente vese que agora o momento total si deber� ser exactamente cero

px=sum(vx)
py=sum(vy)
pz=sum(vz)

write(*,*) "Momentos despois da correccion:",px,py,pz !Escribimos de novo os momentos para comprobar que agora s� son nulos.

!Calculamos a enerx�a cin�tica do xeito habitual
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!A enerx�a total neste momento ser�
write(*,*) "Enerxia total do sistema:",ecin+epot

!Agora, como estamos no microcan�nico debemos fixar a enerx�a total do sistema. Evidentemente de partida esta non ser� a que temos neste momento, senon que debemos variar a velocidade das part�culas para que asi a suma da cin�tica e a potencial tome o valor desexado
write(*,*) "Indique o valor desexado da enerxia do sistema"
read(*,*) et

!Fixamos o valor da enerx�a cin�tica
ecq=et-epot

!Para asegurarnos, po�emos un aviso en caso de que esta enerx�a sexa negativa, cousa que nunca deber�a acontecer
if (ecq<=0.d00) then
      write (*,*)"Error: enerxia cinetica negativa"
      stop
end if

!Agora calculamos o factor polo que deberemos multiplicar as velocidades para que a enerx�a cin�tica sexa a desexada
fixE=dsqrt(ecq/ecin)
!Comprobase de xeito trivial que isto � as�
vx=fixE*vx
vy=fixE*vy
vz=fixE*vz
!Posto que simplemente aparece o mesmo factor en todas as compo�entes da velocidade, o momento total seguir� sendo nulo.

!Calculamos agora a enerx�a cin�tica, que loxicamente deber� coincidir con ecq
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!Finalmente, volcamos os datos obtidos aos diferentes ficheiros
open(90,file=name1)
write(90,9001) np,pl,pli,rc,rc2
write(90,9002) vol,dens
write(90,9003) ecin+epot,ecin,epot
write(90,9000) name1
write(90,9000) name2
close(90)
!Este ficheiro estar� en ASCII

!Debido � gran cantidade de datos de posici�ns, velocidades e aceleraci�ns dos que dispo�emos, grabamos un ficheiro binario (sen formato) con estes)

open(91,file=name2,form='unformatted')
write(91) rx,ry,rz,vx,vy,vz,ax,ay,az
close(91)
   
!Por �ltimo especificamos os formatos empregados (o 9000 xa foi definido con anterioridade)
9001	format(i4,2x,1pe19.12,3(2x,e19.12))
9002	format(1pe19.12,2x,e19.12)
9003	format(1pe19.12,2x,e19.12,2x,e19.12)



end program crea_red