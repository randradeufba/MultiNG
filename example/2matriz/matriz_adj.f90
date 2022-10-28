!Program matriz a adjacência 

! a partir de um de terminado sigma, gera a matriz de adjacência a partir de uma matriz de identidade

!------------------------------

parameter(npm=10000)
real*4 am(npm,npm)                  
integer*4 bm(npm,npm)
real lar(npm), y
character*100 entrada
character*100 saida1,saida2

!--------------------------------

 open (unit=3, file='matriz_adj.dat') !abrindo arquivo

!---------------------------------
!lendo arquivo de dados      
!---------------------------------
read (3,*)n                !numero de vertices
read (3,*)y                 !controle
read (3,500)entrada          !matriz de similaridade
read (3,500)saida1            !matriz de adj
read (3,500)saida2
!----------------------------------

open (unit=4,file=entrada)
open (unit=5,file=saida1)

!----------------------------------

open (unit=4,file=entrada)
     
do i = 1,n
read(4,*)(lar(j),j=1,n)
do j = 1,n
am(i,j) = lar(j)
enddo
enddo

close(unit=4)

!----------------------------------
!teste de leitura (excluir depois)
!----------------------------------
!do i = 1,n
!write(*,530)(am(i,j),j=1,n)
!enddo
!-----------------------------------
!transformando matriz em adjacencia
!-----------------------------------

open (unit=4,file=entrada)
do i = 1,n
 read(4,*)(lar(j),j=1,n)
 do j = 1,n
  am(i,j) = lar(j)
if (am(i,j) .GE. y) then 
  bm(i,j)=1
 else if(am(i,j) .LT. y) then 
  bm(i,j)=0
 else
  bm(i,j)=am(i,j)
end if
 enddo
enddo

!------------------------------------
!Saidas
!------------------------------------
write(*,*)
open(unit=5,file=saida1)
do i = 1,n
!write(*,510)(bm(i,j),j=1,n)
write(5,510)(bm(i,j),j=1,n)
enddo

close(unit=5)

write(*,*)
open(unit=6,file=saida2)
do i = 1,n
!write(*,510)(bm(i,j),j=1,n)
write(6,540)(bm(i,j),j=1,n)
enddo

close(unit=6)

write(*,*)'saida encontra-se no diretorio respectivo'

!------------------------------------
!formatações 
!------------------------------------

500   format(a30)
510   format(10000i1)
520   format(10000f5.2) !entrada com espaço
530   format(10000f4.2) !entrada sem espaço
540   format(10000i2)


!-----------------------------------

end program            

 
