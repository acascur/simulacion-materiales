subroutine verlet (np,rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,dfiv,d2fiv)

use variables_comunes
use define_precision

implicit none

real (kind=doblep),dimension(npmax),intent(inout):: rx,ry,rz,vx,vy,vz,ax,ay,az
integer (kind=enteiro),intent(in):: np
real (kind=doblep), intent(out):: ecin,epot,dfiv,d2fiv
!Definimos as variables que emprega a subrutina verlet, as� como se son de entrada ou saida

!Empregando o algoritmo de verlet, calculamos as posici�ns e as velocidades. Estas �ltimas deben ser calculadas en d�as etapas, xa que necesitamos a aceleraci�n en t e en t+dt
!A aceleraci�n calc�lase coa subrutina potencial unha vez que temos r(t+dt)

!posici�ns en t+dt:
rx=rx+vx*dt+ax*dt2
ry=ry+vy*dt+ay*dt2
rz=rz+vz*dt+az*dt2

!primeira parte do c�lculo das velocidades
vx=vx+ax*dt12
vy=vy+ay*dt12
vz=vz+az*dt12

!Agora calculamos as aceleraci�ns en t+dt para as novas posici�ns empregando a subrutina potencial
call potencial(np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)

!E finalmente calculamos a segunda parte da velocida, co que nos queda a expresi�n completa para v en t+dt
vx=vx+ax*dt12
vy=vy+ay*dt12
vz=vz+az*dt12

!Para finalizar, calculamos a nova enerx�a cin�tica do sistema
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!regresamos ao programa principal
return
end subroutine verlet 
