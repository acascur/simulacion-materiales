subroutine potencial (np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)

use define_precision
use variables_comunes

implicit none

real(kind=doblep):: rrx,rry,rrz,rijx,rijy,rijz
!rrx,rry,rrz: posición de cada partícula, empregase para o bucle para que non teña que acceder ao vector cada vez que o precise
!rijx,rijy,rijz: distancias ao átomo correspondente máis próximo
real(kind=doblep), dimension(npmax),intent(in):: rx,ry,rz
integer(kind=enteiro), intent(in):: np
real(kind=doblep), dimension(npmax),intent(out):: ax,ay,az
real(kind=doblep),intent(out)::epot,dfiv,d2fiv
!definimos as variables que se empregan na subrutina, indicando se estas son de entrada ou de saida.
integer(kind=enteiro):: i,j
!Variables empregadas para os bucles
real(kind=doblep):: a2,a6,a12,dnpmax,fmod,aux1,dist2
!fmod: módulo da forza entre partículas
!a2,a6,a12,aux1: factores do potencial de Lennard-Jones
!dist2: distancia ao cadrado entre partículas

dnpmax=dble(npmax) 	!Convirte npmax en real de dobre precisión

 epot=0.d00
 dfiv=0.d00
 d2fiv=0.d00
 ax=0.d00
 ay=0.d00
 az=0.d00
 !Para asegurarnos de que os acumuladores non teñen un valor distinto de cero na memoria, reiniciámolos a este valor

  do i=1,np-1
        rrx=rx(i)
        rry=ry(i) !Sacamos fóra do bucle de j estes pasos para reducir o tempo do programa e facelo máis eficiente
        rrz=rz(i)
        do j=i+1,np
          rijx=rrx-rx(j)
          rijy=rry-ry(j)
          rijz=rrz-rz(j)

          rijx=rijx-pl*dnint(rijx*pli)		
          rijy=rijy-pl*dnint(rijy*pli)
          rijz=rijz-pl*dnint(rijz*pli)
!Algunha partículas non están dentro do radio de corte, que serán as únicas que consideremos para o cálculo da enerxía potencial, pero si teñen unha partícula equivalente dentro deste, que se situará a unha distancia igual a un número enteiro de veces o lado da caixa.
!Isto corríxese coas tres ecuacións anteriores.
          
          dist2=rijx*rijx+rijy*rijy+rijz*rijz !Calculamos a distancia entre as particulas

          

          if (dist2<=rc2)then        !Só calculamos a interacción entre as partículas se estas se atopan dentro da esfera con radio o radio de corte.
            a2=1.d00/dist2			!Con estas tres liñas calculamos o potencial para cada par de partículas. 
            a6=a2*a2*a2
            a12=a6*a6
            epot=epot+a12-a6		!Calcula a interacción entre a partícula i e a j se j se atopa dentro do radio de corte. Isto vaise engadindo en cada iteracción ao acumulador do potencial total
            aux1=-2.d00*a12+a6		!Derivada primeira do potencial. Eliminamos os factores comúns para evitar operacións innecesarias no bucle, engadíndose ao final nunha das correccións
            dfiv=dfiv+aux1						
            d2fiv=d2fiv+26.d00*a12-7.d00*a6		!Derivada segunda do potencial
            fmod=-aux1*a2						
            ax(i)=ax(i)+fmod*rijx				 
            ay(i)=ay(i)+fmod*rijy				
            az(i)=az(i)+fmod*rijz   !Cálculamos as aceleracións
            ax(j)=ax(j)-fmod*rijx
            ay(j)=ay(j)-fmod*rijy
            az(j)=az(j)-fmod*rijz
          end if
        end do
      end do
      !Xa temos os potenciais e as aceleracións calculadas, salvo por correccións e factores que engadimos a continuación
fac=pi*dnpmax*dnpmax/(vol*rc**3)					
corr_ener=8.d00*fac*(1.d00/(3.d00*rc**6)-1.d00)/3.d00 !Cálculo da corrección na enerxía
corr_sum_rup=16.d00*fac*(-2.d00/(3.d00*rc**6)+1.d00)   !Cálculo da corrección na primeira derivada da enerxía potencial
corr_sum_r2upp=16.d00*fac*(26.d00/(3.d00*rc**6)-7.d00)  !Cálculo da correción na segunda derivada da enerxía potencial
epot=4.d00*epot+corr_ener   !Enerxía potencial correxida
dfiv=24.d00*dfiv+corr_sum_rup !Derivada da enerxía potencial correxida
d2fiv=24.d00*d2fiv+corr_sum_r2upp !Derivada segunda da enerxía potencial correxida

ax=24.d00*ax 
ay=24.d00*ay
az=24.d00*az
!Aceleracións correxidas cos factores correspondentes    
dfiv=dfiv/(3.d00*vol)
d2fiv=(d2fiv-2.d00*dfiv)/(9.d00*vol*vol)
!Derivadas do potencial correxidas cos factores correspondentes

return  

end subroutine potencial