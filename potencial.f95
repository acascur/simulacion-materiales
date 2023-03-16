subroutine potencial (np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)

use define_precision
use variables_comunes

implicit none

real(kind=doblep):: rrx,rry,rrz,rijx,rijy,rijz
!rrx,rry,rrz: posici�n de cada part�cula, empregase para o bucle para que non te�a que acceder ao vector cada vez que o precise
!rijx,rijy,rijz: distancias ao �tomo correspondente m�is pr�ximo
real(kind=doblep), dimension(npmax),intent(in):: rx,ry,rz
integer(kind=enteiro), intent(in):: np
real(kind=doblep), dimension(npmax),intent(out):: ax,ay,az
real(kind=doblep),intent(out)::epot,dfiv,d2fiv
!definimos as variables que se empregan na subrutina, indicando se estas son de entrada ou de saida.
integer(kind=enteiro):: i,j
!Variables empregadas para os bucles
real(kind=doblep):: a2,a6,a12,dnpmax,fmod,aux1,dist2
!fmod: m�dulo da forza entre part�culas
!a2,a6,a12,aux1: factores do potencial de Lennard-Jones
!dist2: distancia ao cadrado entre part�culas

dnpmax=dble(npmax) 	!Convirte npmax en real de dobre precisi�n

 epot=0.d00
 dfiv=0.d00
 d2fiv=0.d00
 ax=0.d00
 ay=0.d00
 az=0.d00
 !Para asegurarnos de que os acumuladores non te�en un valor distinto de cero na memoria, reinici�molos a este valor

  do i=1,np-1
        rrx=rx(i)
        rry=ry(i) !Sacamos f�ra do bucle de j estes pasos para reducir o tempo do programa e facelo m�is eficiente
        rrz=rz(i)
        do j=i+1,np
          rijx=rrx-rx(j)
          rijy=rry-ry(j)
          rijz=rrz-rz(j)

          rijx=rijx-pl*dnint(rijx*pli)		
          rijy=rijy-pl*dnint(rijy*pli)
          rijz=rijz-pl*dnint(rijz*pli)
!Algunha part�culas non est�n dentro do radio de corte, que ser�n as �nicas que consideremos para o c�lculo da enerx�a potencial, pero si te�en unha part�cula equivalente dentro deste, que se situar� a unha distancia igual a un n�mero enteiro de veces o lado da caixa.
!Isto corr�xese coas tres ecuaci�ns anteriores.
          
          dist2=rijx*rijx+rijy*rijy+rijz*rijz !Calculamos a distancia entre as particulas

          

          if (dist2<=rc2)then        !S� calculamos a interacci�n entre as part�culas se estas se atopan dentro da esfera con radio o radio de corte.
            a2=1.d00/dist2			!Con estas tres li�as calculamos o potencial para cada par de part�culas. 
            a6=a2*a2*a2
            a12=a6*a6
            epot=epot+a12-a6		!Calcula a interacci�n entre a part�cula i e a j se j se atopa dentro do radio de corte. Isto vaise engadindo en cada iteracci�n ao acumulador do potencial total
            aux1=-2.d00*a12+a6		!Derivada primeira do potencial. Eliminamos os factores com�ns para evitar operaci�ns innecesarias no bucle, engad�ndose ao final nunha das correcci�ns
            dfiv=dfiv+aux1						
            d2fiv=d2fiv+26.d00*a12-7.d00*a6		!Derivada segunda do potencial
            fmod=-aux1*a2						
            ax(i)=ax(i)+fmod*rijx				 
            ay(i)=ay(i)+fmod*rijy				
            az(i)=az(i)+fmod*rijz   !C�lculamos as aceleraci�ns
            ax(j)=ax(j)-fmod*rijx
            ay(j)=ay(j)-fmod*rijy
            az(j)=az(j)-fmod*rijz
          end if
        end do
      end do
      !Xa temos os potenciais e as aceleraci�ns calculadas, salvo por correcci�ns e factores que engadimos a continuaci�n
fac=pi*dnpmax*dnpmax/(vol*rc**3)					
corr_ener=8.d00*fac*(1.d00/(3.d00*rc**6)-1.d00)/3.d00 !C�lculo da correcci�n na enerx�a
corr_sum_rup=16.d00*fac*(-2.d00/(3.d00*rc**6)+1.d00)   !C�lculo da correcci�n na primeira derivada da enerx�a potencial
corr_sum_r2upp=16.d00*fac*(26.d00/(3.d00*rc**6)-7.d00)  !C�lculo da correci�n na segunda derivada da enerx�a potencial
epot=4.d00*epot+corr_ener   !Enerx�a potencial correxida
dfiv=24.d00*dfiv+corr_sum_rup !Derivada da enerx�a potencial correxida
d2fiv=24.d00*d2fiv+corr_sum_r2upp !Derivada segunda da enerx�a potencial correxida

ax=24.d00*ax 
ay=24.d00*ay
az=24.d00*az
!Aceleraci�ns correxidas cos factores correspondentes    
dfiv=dfiv/(3.d00*vol)
d2fiv=(d2fiv-2.d00*dfiv)/(9.d00*vol*vol)
!Derivadas do potencial correxidas cos factores correspondentes

return  

end subroutine potencial