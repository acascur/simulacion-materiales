!Programa de dinámica molecular. Carlos López Bueno

program dinamica_molecular

use define_precision
use variables_comunes

implicit none

!Programa que realiza a simulación e obtén os valores medios das diferentes variables termodinámicas
!Executaremos este programa un total de 10 veces para obter 10 datos e, posteriormente, no programa resultados_finais, calcularemos simplemente a media destes valores.
!Emprega as subrutinas anteriormente creadas verlet e potencial, así como os módulos variables_comunes e define_precision

integer(kind=enteiro):: np,ndatos,k,ktotal,kpaso
!np: número de partículas na rede
!ndatos:: valor que empregaremos no bucle para saber o número de datos da enerxía que foron gardados
!k::variable empregada para o conteo do bucle
!ktotal:número de pasos totais
!kpaso:valor que determina cada cantos pasos gardamos os valores obtidos

real(kind=doblep),dimension(npmax)::rx,ry,rz,vx,vy,vz,ax,ay,az
!posicións velocidades e aceleracións das partículas

real(kind=doblep):: dnp,et,epot,dfiv,d2fiv,ecin,inv_ecin,t,dktotal
!dnpmax: np convertido a real en dobre precisión
!et: enerxía total
!epot: enerxía potencial
!dfiv: derivada da enerxía potencial
!d2fiv: derivada sergunda da enerxía potencial
!ecin: enerxía cinética
!dktotal: ktotal convertido a real en dobre precisión
!inv_ecin: Inversa da enerxía cinética

real(kind=doblep):: grad_lib,temperatura,aux1,aux2,capv,capp,cesp,cesv,presion,gamma,comp_s,comp_t,alfa_p,alfa_s,alfa_e !alfa_ef??
!grad_lib:grados de liberdade
!temperatura: temperatura do sistema
!aux1,aux2: factores que empregaremos para evitar cálculos innecesaios
!capv,capp: capacidade calorífica a volume/presión constante
!cesp,cesv: calor específico a presión/volume constante
!presion: presión do sistema
!gamma: gamma de Grunesein
!comp_s: factor de compresibilidade adiabática
!comp_t: factor de compresibilidade isotérmico
!alfa_p,alfa_s,alfa_e: coeficientes de compresión

character (len=25):: name1, name2, name3,name4,name5
!name1:Ficheiro cos datos da simulación
!name2:Ficheiro coa configuración do sistema (r,v,a)
!name3:Ficheiro de gravación das enerxías
!name4:Ficheiro cos novos datos da simulación
!name5:Ficheiro de datos de r e v para empregar posteriormente nas correlacións

real(kind=doblep):: ecin_med,inv_ecin_med,epot_med,et_med,dfiv_med,d2fiv_med,dfiv_ecin_med,dfiv2_ecin_med
!Acumuladores para os valores medios de todas as magnitudes a calcular


write(*,*) 'Indique o nome do ficheiro onde estan os datos da simulacion:'
read(*,9000) name1

9000 format(a25) !Formato para as cadeas de caracteres

!Lemos os datos da simulación
open(90,file=name1)
read(90,*) np,pl,pli,rc,rc2
read(90,*) vol, dens,ktotal,kpaso,dt
read(90,9000) name2
read(90,9000) name2
close(90)

!Agora lemos as r,v,a das partículas
open(91,file=name2,form='unformatted')
read(91) rx,ry,rz,vx,vy,vz,ax,ay,az
close(91)

!calculamos os pasos necesarios
dnp=dble(np)
dktotal=dble(ktotal)
dt12=0.5d00*dt
dt2=0.5d00*dt*dt
grad_lib=3.d00*dnp-3.d00
aux1=1.d00-2.d00/grad_lib
aux2=grad_lib/2.d00-1.d00

!Poñemos os acumuladores para as medias a cero

ecin_med=0.d00
inv_ecin_med=0.d00
epot_med=0.d00
et_med=0.d00
dfiv_med=0.d00
d2fiv_med=0.d00
dfiv_ecin_med=0.d00
dfiv2_ecin_med=0.d00

!Agora comezamos cos cálculos. Estes son moi similares ao programa de equilibración

!Pedimos o nome do arquivo no que queremos almacenar os datos r,v cada kpaso e o arquivo para os datos finales
write(*,*) 'Indique o nome do arquivo onde desexa almacenar os datos das enerxias:'
read(*,9000) name3
write(*,*) 'Indique o nome do arquivo onde desexa almacenar os datos finais da simulacion:'
read(*,9000) name4
write(*,*) 'Indique o nome do arquivo onde desexa gardar os datos de r e v (para empregar posteriormente nas correlacions):'
read(*,9000) name5

!En primeiro lugar inicializamos o número de datos gravados a 0
ndatos=0

!Abrimos o ficheiro onde se almacenan os datos e facemos o bucle que nos dará os acumuladores para os valores medios

do k=1,ktotal
  call verlet(np,rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,dfiv,d2fiv)
  et=ecin+epot
  inv_ecin=1.d00/ecin

  !Valores medios
  ecin_med=ecin_med+ecin
  inv_ecin_med=inv_ecin_med+inv_ecin
  epot_med=epot_med+epot
  et_med=et_med+et
  dfiv_med=dfiv_med+dfiv
  d2fiv_med=d2fiv_med+d2fiv
  dfiv_ecin_med=dfiv_ecin_med+dfiv*inv_ecin
  dfiv2_ecin_med=dfiv2_ecin_med+dfiv*dfiv*inv_ecin

  !Indicamos que grave os datos cando desexamos
  if (mod(k,kpaso)==0) then
    ndatos=ndatos+1
    t=dble((ndatos-1)*kpaso)*dt
    open(92,file=name3,access='append')
    write(92,9001) t,et,ecin,epot
    close(92)
    open(94,file=name5,form='unformatted',access='append')
    write(94) rx,ry,rz,vx,vy,vz
    close(94)
  end if

  end do !bucle finalizado

9001	 format(1pe13.6,2x,e13.6,2x,e13.6,2x,e13.6)

!Agora calculamos os valores medios
ecin_med=ecin_med/dktotal
inv_ecin_med=inv_ecin_med/dktotal
epot_med=epot_med/dktotal
et_med=et_med/dktotal
dfiv_med=dfiv_med/dktotal
d2fiv_med=d2fiv_med/dktotal
dfiv_ecin_med=dfiv_ecin_med/dktotal
dfiv2_ecin_med=dfiv2_ecin_med/dktotal

!Gravamos a última configuración

open(91,file=name2,form='unformatted')
write(91) rx,ry,rz,vx,vy,vz,ax,ay,az
close(91)

!Finalmente calculamos todas as magnitudes termodinámicas coas relacións coñecidas

temperatura=2.d00*ecin_med/grad_lib
presion=dnp*temperatura/vol-dfiv_med

capv=1.d00/(1.d00-aux1*ecin_med*inv_ecin_med)
cesv=capv/dnp
    
gamma=1.d00/cesv+vol*aux2*(dfiv_med*inv_ecin_med-dfiv_ecin_med)

comp_s=dnp*temperatura*(1.d00+2.d00*gamma-1.d00/cesv)/vol+vol*d2fiv_med
comp_s=comp_s-vol*aux2*(dfiv2_ecin_med-2.d00*dfiv_med*dfiv_ecin_med+dfiv_med*dfiv_med*inv_ecin_med)
comp_t=comp_s-temperatura*capv*gamma*gamma/vol
comp_s=1.d00/comp_s
comp_t=1.d00/comp_t

capp=capv*comp_t/comp_s
cesp=capp/dnp

alfa_p=capv*gamma*comp_t/vol
alfa_e=1.d00/(presion*vol/capv-gamma*temperatura)
alfa_s=-1.d00/(gamma*temperatura)	

!Xa para rematar gardamos os datos obtidos

open(93,file=name4,access='append')
write(93,9002) name1,np,vol,et_med
write(93,9003) pl,rc
write(93,9004) ktotal,kpaso,dt
write(93,9005) ecin_med,epot_med
write(93,9006) temperatura,presion
write(93,9007) capv,cesv
write(93,9008) capp,cesp
write(93,9009) gamma
write(93,9010) comp_s,comp_t
write(93,9011) alfa_p,alfa_e,alfa_s
write(93,9012) ecin_med/dnp,epot_med/dnp,et_med/dnp
close(93)

9002	format(a25,'np=',i4,2x,'vol=',1pe13.6,2x,'Et=',e13.6)
9003	format('lado caixa=',1pe13.6,2x,'radio de corte=',e13.6)
9004 	format('pasos=',i7,2x,'gardo tras=',i3,2x,'dt=',1pe13.6)
9005	format('Ec=',1pe13.6,2x,'Ep=',e13.6)
9006	format('Temperatura=',1pe13.6,2x,'Presion=',e13.6)
9007	format('Cv=',1pe13.6,2x,'cv=',e13.6)
9008	format('Cp=',1pe13.6,2x,'cp=',e13.6)
9009	format('gamma=',1pe13.6)
9010	format('comp_s=',1pe13.6,2x,'comp_t=',e13.6)
9011	format('alfa_p=',1pe13.6,2x,'alfa_e=',e13.6,2x,'alfa_s=',e13.6)
9012	format('Ec por particula=',1pe13.6,2x,'Ep por particula=',e13.6,2x,'Et por particula=',e13.6)

!O programa está finalizado. Executarase un total de 10 veces e despois realizaráse a media entre os datos obtidos
!Con isto conseguiremos uns valores finais para as variables termodinámicas, así como o seu erro.


end program dinamica_molecular