subroutine verlet (np,rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,dfiv,d2fiv)

use variables_comunes
use define_precision

implicit none

real (kind=doblep),dimension(npmax),intent(inout):: rx,ry,rz,vx,vy,vz,ax,ay,az
integer (kind=enteiro),intent(in):: np
real (kind=doblep), intent(out):: ecin,epot,dfiv,d2fiv
!Definimos as variables que emprega a subrutina verlet, así como se son de entrada ou saida

!Empregando o algoritmo de verlet, calculamos as posicións e as velocidades. Estas últimas deben ser calculadas en dúas etapas, xa que necesitamos a aceleración en t e en t+dt
!A aceleración calcúlase coa subrutina potencial unha vez que temos r(t+dt)

!posicións en t+dt:
rx=rx+vx*dt+ax*dt2
ry=ry+vy*dt+ay*dt2
rz=rz+vz*dt+az*dt2

!primeira parte do cálculo das velocidades
vx=vx+ax*dt12
vy=vy+ay*dt12
vz=vz+az*dt12

!Agora calculamos as aceleracións en t+dt para as novas posicións empregando a subrutina potencial
call potencial(np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)

!E finalmente calculamos a segunda parte da velocida, co que nos queda a expresión completa para v en t+dt
vx=vx+ax*dt12
vy=vy+ay*dt12
vz=vz+az*dt12

!Para finalizar, calculamos a nova enerxía cinética do sistema
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!regresamos ao programa principal
return
end subroutine verlet 
