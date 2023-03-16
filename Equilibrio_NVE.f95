program Equilibrio_NVE

use define_precision
use variables_comunes

implicit none

integer(kind=enteiro):: np,ndatos,k,ktotal,kpaso
!np: número de partículas na rede
!ndatos:: valor que empregaremos no bucle para saber o número de datos da enerxía que foron gardados
!k::variable empregada para o conteo do bucle
!ktotal:número de pasos totais
!kpaso:valor que determina cada cantos pasos gardamos os valores obtidos

real(kind=doblep),dimension(npmax)::rx,ry,rz,vx,vy,vz,ax,ay,az
!posicións velocidades e aceleracións das partículas

real(kind=doblep):: et,epot,dfiv,d2fiv,ecin,t,etd

!et: enerxía total
!epot: enerxía potencial
!dfiv: derivada da enerxía potencial
!d2fiv: derivada sergunda da enerxía potencial
!ecin: enerxía cinética
!t: tempo, empregaráse para gardalo xunto cos datos e posteriormente realizar unha gráfica
!etd:Enerxía total desexada, é dicir, a do sistema (empregarase para a corrección debida ao salto producido inicialmente.

character (len=25):: name1,name2,name3,name4
!name1:ficheiro cos datos da simulacion
!name2:ficheiro cos datos de r,v,a
!name3:ficheiro para gardar as enerxías no intervalo de tempo que definamos
!name4:datos da nova simulación

write(*,*) 'Inserte o nome do ficheiro de datos da simulacion:'
read(*,9000) name1

9000	format(a25)
!Lemos os datos da simulación realizada para a obtención da configuración inicial (crea_red)
open(90,file=name1)
read(90,*) np,pl,pli,rc,rc2
read(90,*) vol,dens
read(90,*) et,ecin,epot
read(90,9000) name2
read(90,9000) name2
close(90)

!Poñemos que a enerxía total desexada é a que lemos neste instante e así almacenamos a variable
etd=et

!Lemos os datos da configuración inicial obtidos
open(91, file=name2,form='unformatted')
read(91) rx,ry,rz,vx,vy,vz,ax,ay,az
close(91)

!Facemos unha petición dos datos restantes
write(*,*) 'Indique o numero de pasos a realizar'
read(*,*) ktotal
write(*,*) 'Indique cada cantos pasos desexa gardar os datos das enerxias:'
read(*,*) kpaso
write(*,*) 'Indique o paso de tempo (dt) que desexa para a simulacion:'
read(*,*) dt
write(*,*) 'Indique o ficheiro no que desexa gardar os datos das enerxias:'
read(*,9000) name3
write(*,*) 'Indique o nome do novo ficheiro no que desexa que se garden os novos datos da simulacion:'
read(*,9000) name4

!calculamos algunhas variables que empregaremos posteriormente   	 	    		
dt12=0.5d00*dt
dt2=0.5d00*dt*dt


!En primeiro lugar, posto que a perfecta colocación das partículas fai que, cando comeza a súa evolución temporal, se produza un salto na enerxía
!Isto implica que, a enerxía á que se vai equilibrar o sistema non vai ser exactamente a que nos marcamos, e posto que estamos no microcanónico, debemos correxir o erro
!Polo tanto, evolucionamos o sistema nun pequeno intervalo temporal, o equivalente a 1000 pasos (0.1 s), para posteriormente correxir a enerxía cinética e deixar o sistema evolucionar ata o equilibrio

do k=1,1000
  call verlet(np,rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,dfiv,d2fiv) 
end do
et=ecin+epot
write(*,*) 'Enerxia total antes da correccion:', et

!Chamamos a subrutina que correxirá esta enerxía
call corr_enerx(np,vx,vy,vz,etd,ecin,epot)
et=ecin+epot
write(*,*) 'Enerxia total despois da correccion:', et

!A enerxía xa é a desexada. Agora pasamos a equilibrar o sistema

!Creamos un bucle para que repita o algoritmo de verlet o número desexado de veces, chegando ao equilibrio
ndatos=0
do k=1,ktotal
call verlet(np,rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,dfiv,d2fiv) 
et=ecin+epot


!Imos agora impoñer a condición de que, cada kpaso, se almanece o valor da enerxía
	if (mod(k,kpaso)==0) then
	t=dble((ndatos)*kpaso)*dt  
	ndatos=ndatos+1
	open(93,file=name3,access='append') !Poññemos access append para que así escriba ao final do archivo e non elimine os datos anteriores
	write(93,9003)  t,epot,ecin,et
	close(93)
	end if

!Xa temos r,v,a para as partículas un tempo dt posterior.
end do

!Almacenamos os datos do sistema no equilibrio
open(92,file=name4)
write(92,9001) np,pl,pli,rc,rc2
write(92,9002) vol,dens,ktotal,kpaso,dt
write(92,9000) name1
write(92,9000) name2
write(92,9000) name3
write(92,9000) name4
close(92)
9001	format(i4,2x,1pe19.12,3(2x,e19.12))
9002	format(1pe19.12,2x,e19.12,2x,i7,2x,i4,2x,e19.12)

!Ao finalizar o bucle teremos o sistema en equilibrio (sempre e cando elixamos intervalos temporais axeitados)

9003 	format (1pe13.6,2x,e13.6,2x,e13.6,2x,e13.6)

!Para rematar, almacenamos os datos de r,v,a no equilibrio

open(91, file=name2,form='unformatted')
write(91) rx,ry,rz,vx,vy,vz,ax,ay,az
close(91)

stop
end program
