module variables_comunes

use define_precision !Usamos o módulo define precisión
implicit none

integer(kind=enteiro), parameter:: npmax=500, numk=5 
!definimos o número máximo de partículas na caixa (npmax) e o número de celdas fcc de cada lado da caixa
real(kind=doblep):: pi=3.14159265359 !definimos el número pi para tenerlo en memoria
real(kind=doblep):: pl, pli, vol, dens, rc, rc2, dl, dl12, fac, dt, dt2,dt12
!pl: lado da caixa
!pli: inverso de pl
!vol: volume da caixa
!dens: densidade de partículas na caixa
!rc: radio de corte do potencial.
!rc2: rc elevado ao cadrado
!dl:lonxitude dunha celda
!dl12: metade de dl
!fac: factor empregado para a corrección da enerxía potencial
!dt: paso empregado para o algoritmo verlet
!dt2: 1/2*dt*dt
!dt12: 1/2*dt 
!Estes dous últimos empreganse no algoritmo verlet

real(kind=doblep):: corr_ener=0, corr_sum_rup=0, corr_sum_r2upp=0
!definimos os parámetros para a corrección de enerxías potenciais e as súas derivadas e iniciámolas a cero

end module variables_comunes