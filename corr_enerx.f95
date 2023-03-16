subroutine corr_enerx(np,vx,vy,vz,etd,ecin,epot)

use variables_comunes
use define_precision

implicit none

integer (kind=enteiro),intent(in)::np
real(kind=doblep), dimension(np), intent(inout)::vx,vy,vz
real(kind=doblep),intent(inout):: ecin
real(kind=doblep),intent (in):: etd,epot
real(kind=doblep):: ecind,factor,px,py,pz,dnp
!ecind: Enerxía cinética desexada
!facotr: Factor de corrección das velocidades
!px,py,pz: Momentos
!dpn:np en dobre precisión

dnp=dble(np)
!En primeiro lugar aseguramos que o momento sexa nulo
px=sum(vx)
py=sum(vy)
pz=sum(vz)
!Correximos as velocidades
px=px/dnp
py=py/dnp
pz=pz/dnp
vx=vx-px
vy=vy-py
vz=vz-pz

ecind=etd-epot
 !Poñemos una aviso por se nos da unha enerxía cinética negativa
 
if (ecind<=0.d00) then
	   	write (*,*) "Error: enerxia cinetica negativa"
	   	stop
end if

!Calculamos o factor de corrección
factor=dsqrt(ecind/ecin)

!Multiplicamos as velocidades por este factor

vx=factor*vx
vy=factor*vy
vz=factor*vz

!Comprobamos que o momento sexa nulo
px=sum(vx)
py=sum(vy)
pz=sum(vz)

write(*,*) 'Momento tras a correccion:',px,py,pz

!E finalmente calculamos a enerxía cinética
ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

return
end subroutine corr_enerx

