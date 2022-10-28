c     programa multiG
c     constroi multiplex de n camadas a partir de redes com mesmo numero de nós
c     permite escolher se faz o acoplamento ou não

      parameter(npm=10000)

      integer am(npm,npm)
      integer lar(npm)
      integer nm, nm2, k, l, u1, u2
      character*100 entrada(npm)
      character*100 saida, saida1
      

c	entrada de dados
c==================================================================
      open (unit=3,file='multiplex.dat')
      read (3,*)nm,ncop, l
	nm2 = nm * l
      
      do k=1,l
       read (3,fmt='(a30)')entrada(k)
      enddo
      read (3,fmt='(a30)')saida, saida1

c==================================================================
c       Arquivo de entrada

c	nm: numero de nos em cada camada
c	ncop: valor do acoplamento (0 ou 1)
c       l: numero de camadas
c	entrada(k): guarda a matriz de adjacencia das redes
c       saida: multiplex sem espaçamento 
c       saida1: multiplex com espaçamento 
c==================================================================

      write(*,'(/,"nos = ",i4, 3x, "camadas = ",i2)') nm, l
      do k=1,l
      write(*,fmt='(/,"matriz",i2,2x,A)') k, entrada(k)
      enddo

c============================================================
c   monta a diagonal da supramatriz 

        do i = 1,nm2
	 do j = 1,nm2
	  am(i,j) = 0
	 enddo
	enddo

      do k = 1, l
        
      open (unit=3+k,file=entrada(k))
     
        do i = 1,nm
	 read(3+k,fmt='(1000I1)')(lar(j),j=1,nm)
	 do j =1,nm
	  if(lar(j).ne.1)lar(j)=0
	 enddo
	 do j = 1,nm
	  am(i+(k-1)*nm,j+(k-1)*nm) = (lar(j))
	 enddo
	enddo

        close(unit=3+k)
       
        enddo

c==================================================================
c     monta o coplamento entre as matrizes 

      do u1 = 1, l
      do u2= 1, l

       if(u1 .ne. u2)then
	do i = 1,nm
         am(i+(u1-1)*nm,i+(u2-1)*nm)=ncop
        enddo
       endif
      enddo
      enddo

c==================================================================
c    saidas (sem e com espaçamento)     
    
      open (unit=3+(l+1),file=saida)


      do i = 1,nm2
       write(unit=3+(l+1),fmt='(1000I1)')(am(i,j),j=1,nm2)
      enddo

      
      close(unit=3+(l+1))
 

c=====================================================================

      end program 

