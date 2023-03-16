module define_precision

implicit none

integer, parameter :: enteiro=SELECTED_INT_KIND(9)      
integer, parameter :: doblep=SELECTED_REAL_KIND(15,307)
!definimos os dous tipos de números que imos empregar (enteiros e números con doble precisión)


end module define_precision