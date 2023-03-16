program medias

use define_precision

implicit none



integer (kind=enteiro):: nsim,np,ktotal,kpaso,i
character (len=25)::name1,name2,name3
real(kind=doblep),dimension(10)::et,epot,ecin,temperatura,presion,capv,comp_s,gamma,alfa_e
real(kind=doblep):: vol,pl,rc,dt,cesv,capp,cesp,comp_t,alfa_s,alfa_p,etp,ecp,epp
real(kind=doblep)::et_media,epot_media,ecin_media,temp_media,p_media,capv_media,comp_s_media,gamma_media,alfa_e_media
real(kind=doblep)::et2_media,epot2_media,ecin2_media,temp2_media,p2_media,capv2_media,comp_s2_media,gamma2_media,alfa_e2_media
real(kind=doblep)::err_et,err_epot,err_ecin,err_temp,err_p,err_capv,err_comp_s,err_gamma,err_alfa_e


write(*,*)'Indique o nome do ficheiro onde estan os datos finais das simulacions:'
read(*,9000) name1
write(*,*)'Indique o numero de simulacions realizadas:'
read(*,*) nsim
write(*,*) 'Indique o nome do ficheiro onde desexa gardar os datos finais:'
read(*,*) name2

9000	format(a25)

!Lemos os datos

open(90,file=name1)

do i=1,nsim
  read(90,9002) name3,np,vol,et(i)
  read(90,9003) pl,rc
  read(90,9004) ktotal,kpaso,dt
  read(90,9005) ecin(i),epot(i)
  read(90,9006) temperatura(i),presion(i)
  read(90,9007) capv(i),cesv
  read(90,9008) capp,cesp
  read(90,9009) gamma(i)
  read(90,9010) comp_s(i),comp_t
  read(90,9011) alfa_p,alfa_e(i),alfa_s
  read(90,9012) ecp,epp,etp
end do
 
  close(10)


9002	format(a25,3x,i4,6x,1pe13.6,5x,e13.6)
9003	format(11x,1pe13.6,17x,e13.6)
9004 	format(6x,i7,13x,i3,5x,1pe13.6)
9005	format(3x,1pe13.6,5x,e13.6)
9006	format(12x,1pe13.6,10x,e13.6)
9007	format(3x,1pe13.6,5x,e13.6)
9008	format(3x,1pe13.6,5x,e13.6)
9009	format(6x,1pe13.6)
9010	format(7x,1pe13.6,9x,e13.6)
9011	format(7x,1pe13.6,9x,e13.6,9x,e13.6)
9012	format(17x,1pe13.6,19x,e13.6,19x,e13.6)


!Calculamos as medias
et_media=sum(et(1:nsim))/dble(nsim)
epot_media=sum(epot(1:nsim))/dble(nsim)
ecin_media=sum(ecin(1:nsim))/dble(nsim)
temp_media=sum(temperatura(1:nsim))/dble(nsim)
p_media=sum(presion(1:nsim))/dble(nsim)
capv_media=sum(capv(1:nsim))/dble(nsim)
gamma_media=sum(gamma(1:nsim))/dble(nsim)
comp_s_media=sum(comp_s(1:nsim))/dble(nsim)
alfa_e_media=sum(alfa_e(1:nsim))/dble(nsim)

!Calculamos a media dos cadrados
et2_media=sum(et(1:nsim)*et(1:nsim))/dble(nsim)
epot2_media=sum(epot(1:nsim)*epot(1:nsim))/dble(nsim)
ecin2_media=sum(ecin(1:nsim)*ecin(1:nsim))/dble(nsim)
temp2_media=sum(temperatura(1:nsim)*temperatura(1:nsim))/dble(nsim)
p2_media=sum(presion(1:nsim)*presion(1:nsim))/dble(nsim)
capv2_media=sum(capv(1:nsim)*capv(1:nsim))/dble(nsim)
gamma2_media=sum(gamma(1:nsim)*gamma(1:nsim))/dble(nsim)
comp_s2_media=sum(comp_s(1:nsim)*comp_s(1:nsim))/dble(nsim)
alfa_e2_media=sum(alfa_e(1:nsim)*alfa_e(1:nsim))/dble(nsim)

!E finalmente calculamos os erros
err_et=2.d00*dsqrt(et2_media-et_media*et_media)/dble(nsim-1)
err_epot=2.d00*dsqrt(epot2_media-epot_media*epot_media)/dble(nsim-1)
err_ecin=2.d00*dsqrt(ecin2_media-ecin_media*ecin_media)/dble(nsim-1)
err_temp=2.d00*dsqrt(temp2_media-temp_media*temp_media)/dble(nsim-1)
err_p=2.d00*dsqrt(p2_media-p_media*p_media)/dble(nsim-1)
err_capv=2.d00*dsqrt(capv2_media-capv_media*capv_media)/dble(nsim-1)
err_gamma=2.d00*dsqrt(gamma2_media-gamma_media*gamma_media)/dble(nsim-1)
err_comp_s=2.d00*dsqrt(comp_s2_media-comp_s_media*comp_s_media)/dble(nsim-1)
err_alfa_e=2.d00*dsqrt(alfa_e2_media-alfa_e_media*alfa_e_media)/dble(nsim-1)

!Por último calculamos os datos restantes a partir dos valores obtidos
comp_t=1/comp_s_media-temp_media*capv_media*gamma_media*gamma_media/vol
comp_t=1/comp_t

capp=comp_t*capv_media/comp_s_media

alfa_p=capv_media*gamma_media*comp_t/vol
alfa_s=-1.d00/(gamma_media*temp_media)

!Para rematar gardamos os datos obtidos

open(91,file=name2,access='append')
write(91,9013) np,dble(np)/vol,pl,vol
write(91,9014) '         Et=',et_media,err_et,'Et por particula=',et_media/dble(np),err_et/dble(np)
write(91,9014) '       Ecin=',ecin_media,err_ecin,'Ecin por particula=',ecin_media/dble(np),err_ecin/dble(np)
write(91,9014) '       Epot=',epot_media,err_epot,'Epot por particula=',epot_media/dble(np),err_epot/dble(np)
write(91,9015) 'Temperatura=',temp_media,err_temp
write(91,9015) '    Presion=',p_media,err_p
write(91,9014) '         Cv=',capv_media,err_capv,'                cv=',capv_media/dble(np),err_capv/dble(np)
write(91,9016) '         Cp=',capp,'                cp=',capp/dble(np)
write(91,9015) '         ks=',comp_s_media,err_comp_s
write(91,9017) '         kt=',comp_t
write(91,9015) '      Gamma=',gamma_media,err_gamma
write(91,9015) '     Alfa_e=',alfa_e_media,err_alfa_e
write(91,9017) '     Alfa_s=',alfa_s
write(91,9017) '     Alfa_p=',alfa_p       

9013	format('SIMULACION MICROCANONICO: np=',i4,'    densidade=',1pe13.6,'    lado caixa=',e13.6,2x,'    volume=',e13.6)
9014	format(a12,1pe13.6,' err=',e13.6,5x,a19,e13.6,' err=',e13.6)
9015	format(a12,1pe13.6,' err=',e13.6)
9016	format(a12,1pe13.6,23x,a19,e13.6)
9017	format(a12,1pe13.6)




end program medias