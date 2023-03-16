module variables_comunes

use define_precision !Usamos o m�dulo define precisi�n
implicit none

integer(kind=enteiro), parameter:: npmax=500, numk=5 
!definimos o n�mero m�ximo de part�culas na caixa (npmax) e o n�mero de celdas fcc de cada lado da caixa
real(kind=doblep):: pi=3.14159265359 !definimos el n�mero pi para tenerlo en memoria
real(kind=doblep):: pl, pli, vol, dens, rc, rc2, dl, dl12, fac, dt, dt2,dt12
!pl: lado da caixa
!pli: inverso de pl
!vol: volume da caixa
!dens: densidade de part�culas na caixa
!rc: radio de corte do potencial.
!rc2: rc elevado ao cadrado
!dl:lonxitude dunha celda
!dl12: metade de dl
!fac: factor empregado para a correcci�n da enerx�a potencial
!dt: paso empregado para o algoritmo verlet
!dt2: 1/2*dt*dt
!dt12: 1/2*dt 
!Estes dous �ltimos empreganse no algoritmo verlet

real(kind=doblep):: corr_ener=0, corr_sum_rup=0, corr_sum_r2upp=0
!definimos os par�metros para a correcci�n de enerx�as potenciais e as s�as derivadas e inici�molas a cero

end module variables_comunes