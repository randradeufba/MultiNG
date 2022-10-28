c programa Caminhos minimos
c==============================================================================
c baseado em madchar13
c
c usa o algoritimo em madchar13 para obter os minimos caminhos entre dois nos
c de um multiplex. usa como entrada a matriz de adjacencia de um multiplex 
c determinada com aux�lio do programa multip02 
c
c calcula propriedades de redes a partir de sua matriz de adjacencia
c ou atraves de uma lista de pares de nodos conectados
c
c usa o algoritmo baseado no produto booleano de matrizes de adjacencia
c otimiza o produto de matrized mediante 
c (1) introdu��o de listas; (2) elimina��o de elementos j� calculados; 
c (3) identifica��o de n�s que atigiram o maximo valor de vizinhan�a  
c
c usa uma unica matriz (npm,npm) e um vetor de dimensao npm*(npm-1)/2
c usa o tri�ngulo superior para a matriz de adjacencia e 
c usa o tri�ngulo inferior para acumula��o das matrizes de vizinhan�a de 
c ordem superior
c
c n�o calcula espectro das matrizes de adjacencia
c
c calcula as seguintes medidas de rede:
c (a) coeficiente de aglomera��o de cada n�, coeficiente de aglomera��o 
c medio da rede para cada vizinhan�a de ordem superior
c (b) grau de cada no, grau medio da rede para cada vizinhan�a de ordem 
c superior
c (c) distancia minima media de cada n�, distancia minima media da rede
c (d) di�metro da rede
c (e) dist�ncia minima entre cada par de nos = elemento da matriz de 
c     vizinha�a 
c (f) coeficiente de assortatividade de cada no, coeficiente de assortati-
c     vidade medio da rede para cada vizinhan�a de ordem superior
c (g) dimens�o fractal da rede
c
c controla o calculo das medidas (a), (f) e (g) para vizinhan�as de 1 at� 
c limites superiores npa, npf e npg
c
c arquivos de saida: resultados para  valores de ipx=1,..,npx
c
c saida7: (a) coeficiente de clusteriza��o de cada n� e valor m�dio da rede
c saida8: (b) valor do grau de cada n�  e valor m�dio da rede
c saida9: (c,d) distancia minima m�dia por cada no, distancia minima media  
c         da rede e tambem di�metro da rede
c saida11:(e) matriz de vizinhan�a
c saida12:(f) coeficiente de assortatividade para as diversas matrizes
c saida13:(g) dados para calculo da dimens�o fractal
c saida14:(h) matriz de vizinhan�a da multiplex
c==============================================================================

c      use msimsl
c      use portlib
      parameter(nvm=300,npm=1000,nfaclis=50)
      integer*2 am(npm,npm)                 ! 
      integer*2 cm(npm*(npm+1)/2)
      integer lar(npm), pra(npm), pro(npm), kk(0:npm,0:nvm)
      real lis(nfaclis*npm,2)
      real ga(0:npm), kmedio(0:nvm),rk(0:nvm)
      real gamedio(0:nvm),Lmedio(0:npm),gam(npm,nvm),dfr(1000,0:2)
      character*100 entrada
      character*100 saida7,saida8,saida9,saida11,saida12,saida13,saida14

c==============================================================================
c identifica��o dos parametros
c==============================================================================
c nvm:                valor maximo do di�metro
c npm:                valor maximo do numero de nos
c nfaclis:            parametro interno usado na otimizacao do produto  
c                     booleano. se ultrapassado durante o calculo, mensagem 
c                     de erro � emitida  
c am(npm,npm):        matriz com informa��o sobre a rede: no tri�ngulo 
c                     superior guarda a matriz de adjacencia e no tri�ngulo 
c                     inferior acumula as matrizes  de vizinhan�a de ordem 
c                     superior
c cm(npm*(npm+1)/2):  vetor que guarda a matriz de vizinha�a
c lar(npm):           espa�o interno usado para leitura dos elementos de 
c                     matriz
c pra(npm):           espa�o interno usado para controle do produto booleano
c pro(npm):           espa�o interno usado para controle do produto booleano
c kk(0:npm,0:nvm:)    guarda o grau de cada no em cada vizinhan�a. na coluna 
c                     0 guarda a soma dos graus em cada vizinhan�a o que 
c                     equivale ao tamanho do cluster ao qual pertence
c lis(nfaclis*npm,2): espa�o interno usado para acelerar produto booleano 
c ga(0:npm):          valor do coeficiente de aglomera��o de cada no, obtido na 
c                     rotina rede1
c kmedio(0:nvm):      grau medio da rede em cada vizinhan�a 
c rk(0:nvm):                  grau de assortatividade m�dio da rede em cada vizinhanca
c gamedio(0:nvm):     valor medio do coeficiente de aglomera��o da rede em cada 
c                     vizinhan�a 
c Lmedio(0:npm):      distancia minima media de cada no
c gam(npm,nvm):       valor do coeficiente de aglomera��o de cada no em cada 
c                     vizinhan�a 
c dfr(1000,0:2):      dados para o calculo da dimensao fractal

c      tempo1 = cpsec()
c      write(*,*)'tinic=',tempo1
      
       open (unit=3, file='caminho.dat')

c==============================================================================
c      entrada de dados
c==============================================================================

10    read (3,*)irede
      if (irede.lt.0)stop
      read (3,*)nm,ns,np,npa,npf,npg
      read (3,*)ient,npula
      read (3,*)ncamada

	nms = nm/ncamada

c==============================================================================
c irede: identifica o tipo de rede com os seguintes valores
c        0-tridiagonal
c        1-small world
c        2-randomica
c        3-scale free
c        4-apolonia
c        5-le matriz em arquivo
c nm:    numero de nos na rede, no caso apoloniano, nm � inicialmente a gera��o
c ns:    numero de diferentes redes que ser�o analizadas
c np:    potencia booleana maxima
c npa:   potencia booleana maxima para o calculo do coeficiente de aglomera��o
c npf:   potencia booleana maxima para o calculo do grau de assortatividade
c npg:   controla calculo da dimensao fractal: 0, calcula; 1, n�o calcula
c ient:  identifica maneira de entrada de dados
c        0-matriz de adjacencia
c        1-lista
c npula: numero de linhas a ser pulada no arqui de entrada de dados, seja matriz
c        de adjacencia, seja lista

      read (3,500)entrada
      read (3,500)saida7,saida8,saida9,saida11,saida12,saida13,saida14

c==============================================================================
c arquivo de entrada
c==============================================================================
c entrada: se ient = 0 guarda a matriz de adjacencia durante todo o programa
c==============================================================================

      open (unit=4,file=entrada)
      open (unit=7,file=saida7)
      open (unit=8,file=saida8)
      open (unit=9,file=saida9)
      open (unit=11,file=saida11)
      open (unit=12,file=saida12)
      open (unit=13,file=saida13)

      open (unit=20,file='tempocpu.dat',access='append')

      intime = time()
      xtime = rand(intime)

c==============================================================================
c comeco do grand loop no numero de amostras
c==============================================================================
      do is = 1,ns

c==============================================================================
c coloca zero nas diversas vari�veis

       do i = 1,nm
        ga(i) = 0.
        Lmedio(i) = 0.
        do j = 1,nvm
         gam(i,j) = 0.
         kk(i,j) = 0
         kmedio(j) = 0.
         gamedio(j) = 0.
        enddo
       enddo

       npp = np

c==============================================================================
c gera a rede

       if (irede.eq.0) then
        call tridi(am,nm,npm)
       else if (irede.eq.1) then
        read (3,*)pc       
        call smalw(am,nm,npm,pc)
       else if (irede.eq.2) then 
        read (3,*)pc       
        call rando(am,nm,npm,pc)
       else if (irede.eq.3) then
        read (3,*)pc       
        call scafr(am,nm,npm)
       else if (irede.eq.4) then
        call apolo(am,nm,ng,npm)
       else if (irede.eq.5) then

        open (unit=4,file=entrada)
        do i = 1,npula
         read(4,510)
        enddo
        if(ient.eq.0)then
         do i = 1,nm
          read(4,510)(lar(j),j=1,nm)
          do j = 1,nm
           am(i,j) = lar(j)
          enddo
         enddo
        else if (ient.eq.1) then
        endif
        close(unit=4)

       endif

c==============================================================================
c escreve am(i,j) em cm(ll)

       do i = 1,nm-1
        do j = i+1,nm
         am(j,i) = am(i,j)
         ll = (2*nm-i)*(i-1)/2+j-i
         cm(ll) = am(i,j) 
        enddo
       enddo
       do i = 1,nm
        am(i,i) = 1
        pra(i) = 1
        pro(i) = 0
       enddo

c==============================================================================
c prepara a lista de elementos nao nulos de a(i,j)

       do i = 1,nm
        do j = 1,2
         lis(i,j) = 0.
        enddo
       enddo

       iq = 1

       do i = 1,nm
        lis(i,2) = iq
        do j = 1,nm
         if(am(i,j).eq.1)then
          lis(iq,1) = j
          iq = iq + 1
         endif
        enddo
       enddo

       lis(i,2) = iq

       if (iq.gt.nfaclis*npm)then
        write(*,*)'iq>nfaclis*npm, iq =',iq,'  nfaclis*npm =',nfaclis*npm
        stop
       endif

c==============================================================================
c calcula o numero de arestas e de nos conectados

       na = 0
       do i = 1,nm
        ic = 0
        do j = 1,nm
         if(am(i,j).gt.0)then 
          na = na + 1
          ic = 1
         endif
        enddo
        nc = nc + ic
       enddo
       na = na/2

c==============================================================================
c comeca o loop para calculo das diferentes matrizes mad no triangulo superior

       do i = 1,nm-1
        do ip = 1,np
         if(pra(i).eq.0.and.pro(i).lt.ip)goto 90 
c==============================================================================
c desvio para finalizar o programa

         if (ip.eq.np) go to 100
c==============================================================================
         pra(i) = 0

         do j = i+1,nm

          am(i,j) = 0

          if (am(j,i).eq.1) goto 126
            
          do k = int(lis(j,2)),int(lis(j+1,2)-1)
           k1 = lis(k,1)
           if(i.lt.k1) then
            if (cm((2*nm-i)*(i-1)/2+k1-i).gt.0)then
             if(cm((2*nm-i)*(i-1)/2+k1-i).le.ip)go to 124
            endif
           else if (i.gt.k1) then
            if (cm((2*nm-k1)*(k1-1)/2+i-k1).gt.0)then
             if(cm((2*nm-k1)*(k1-1)/2+i-k1).le.ip)go to 124
            endif
           endif
          enddo

          goto 125
124       am(i,j) = 1
          pra(i) = 1
125       continue
          ll = (2*nm-i)*(i-1)/2+j-i
          cm(ll) = am(i,j)*(ip+1)
          pro(j) = max(pro(j),cm(ll))

126      enddo

c==============================================================================
c transfere a matriz mad(ip) para o triangulo inferior 

         do j = i+1,nm
          am(j,i) = am(j,i) + am(i,j)
         enddo

c==============================================================================
c fim do loop em ip

        enddo

c==============================================================================
c fim do loop em i

90     enddo

c==============================================================================
c completa calculo dos resultados para a amostra is

100    continue

c==============================================================================
c calcula diametro npp

       npp = 0
       do i = 1,nm*(nm-1)/2
        npp = max(cm(i),npp)
       enddo

c==============================================================================
c copia o triangulo inferior no triangulo superior

       do i = 1,nm
        do j = i+1,nm
         am(i,j) = am(j,i)
         enddo
       enddo

c==============================================================================
c calcula o tamanho de cluster de cada no

       do i = 1,nm
        kk(0,i) = 0                                                                                
        do j = 1,nm
         kk(0,i) = kk(0,i) + am(j,i)
         enddo
       enddo

c==============================================================================
c calcula o caminho minimo medio de cada no e da rede

       xlmedio  = 0
	 do i = 1,nm
	  lmedio(i) = 0
        do j = 1,i-1
         ll = (2*nm-j)*(j-1)/2+i-j
         lmedio(i) = lmedio(i) + cm(ll)
        enddo
        lar(i) = 0
        do j = i+1,nm
         ll = (2*nm-i)*(i-1)/2+j-i
         lmedio(i) = lmedio(i) + cm(ll)
        enddo
	  xlmedio = xlmedio + lmedio(i)
       enddo

       xlmedio = xlmedio/(nm*(nm-1.))

c==============================================================================
c calcula o grau de cada n� de cada rede

       do ip = 1,npp
        do i = 1,nm
         kk(i,ip) = 0
        enddo
       enddo 

	 do i = 1,nm
        do j = 1,i-1
         ll = (2*nm-j)*(j-1)/2+i-j
         am(i,j) = cm(ll)
        enddo
        am(i,i) = 0
        do j = i+1,nm
         ll = (2*nm-i)*(i-1)/2+j-i
         am(i,j) = cm(ll)
        enddo
       enddo

       do i = 1,nm
        do j = 1,nm
         kk(i,am(i,j)) = kk(i,am(i,j)) + 1
        enddo
       enddo 

       do ip = 1,npp
        kmedio(ip) = 0
        do i = 1,nm
         kmedio(ip) = kmedio(ip) + kk(i,ip)
        enddo
        kmedio(ip) = kmedio(ip)/nm
       enddo
c==============================================================================
c calcula o numero de nos no maior cluster

       mcl = 0
       do i = 1,nm
        mcl = max(kk(0,i),mcl)
       enddo

c==============================================================================
c calcula a distancia minima m�dia da rede limitada ao(s) maior(es) cluster(s)

       xm2 = 0.
       if (mcl.gt.1)xm2 = 1./float((mcl-1)*mcl)
       imc = 0
       do i = 1,nm
        if (kk(0,i).eq.mcl)then
         imc = imc + 1
         xlmd = xlmd + lmedio(i) 
         ylmd = ylmd + lmedio(i)/kk(0,i)/kk(0,i)
         zlmd = zlmd + lmedio(i)/kk(0,i)/mcl
        endif
       enddo
       idege = imc/mcl
       xlmd = xlmd*xm2*mcl/imc

c==============================================================================
c calculo do coeficiente de clusteriza��o e de assortatividade das redes ip
       
       do ip = 1,npp
       
        if (ip.le.npa.or.ip.le.npf) then
        
         do i = 1,nm
          do j = 1,i-1
           am(i,j) = 0
      	 ll = (2*nm-j)*(j-1)/2+i-j
           if(cm(ll).eq.ip)am(i,j) = 1
          enddo
          am(i,i) = 0
          do j = i+1,nm
           am(i,j) = 0
           ll = (2*nm-i)*(i-1)/2+j-i
           if(cm(ll).eq.ip)am(i,j) = 1
          enddo
         enddo

c==============================================================================
c chama rotina clust: coeficiente de clusteriza��o gam da rede ip

         if (ip.le.npa) then

          call clust(am,nm,npm,ga)
          gamedio(ip) = ga(0)/nm
          do j = 1,nm
           gam(j,ip) = ga(j)
          enddo

         endif
c==============================================================================
c chama rotina assorta: grau de assortatividadeda rk matriz mad(ip)

         if (ip.le.npf) then

          call assorta(am,nm,npm,nvm,rrk,kk,ip)
          rk(ip) = rrk

         endif

        endif

       enddo
c==============================================================================
c calcula a dimens�o fractal

       if(npg.eq.0) then

        do i = 1,nm
         do j = 1,i-1
          ll = (2*nm-j)*(j-1)/2+i-j
          am(i,j) = cm(ll)
         enddo
         am(i,i) = 0
         do j = i+1,nm
          ll = (2*nm-i)*(i-1)/2+j-i
          am(i,j) = cm(ll)
         enddo
        enddo

        call fracnet(am,nm,npm,npp,dfr)

       endif

c==============================================================================
c saida dos resultados
c==============================================================================

c==============================================================================
c escreve coeficiente de clusteriza��o em saida7

       write(7,500)entrada
       write(7,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na
       write(7,*)'coeficiente de clusteriza��o C(L), L=1,...,',npp

       if(irede.lt.5)write(7,*)' pc = ',pc

       write(7,560)((k),k=1,npp+1)
       write(7,540)(gamedio(k),k=1,npp+1)
       write(7,*)'C para cada n� e valor de L'
       do j = 1,nm
        write(7,520)j,(gam(j,k),k=1,npp+1)
       enddo

c==============================================================================
c escreve distribui��o de nos em saida8

       write(8,500)entrada
       write(8,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na
       write(8,*)'numero de n�s k(L), L=1,...,',npp

       if(irede.lt.5)write(8,*)' pc = ',pc

       write(8,560)((k),k=1,npp+1)
       write(8,570)(kmedio(k),k=1,npp+1)
       write(8,*)'k para cada n� e valor de L'
       do j = 1,nm
        write(8,520)j,(float(kk(j,k)),k=0,npp+1)
       enddo
c==============================================================================
c escreve di�metro, distancia minima m�dia por cada no e a distancia minima 
c media da rede em saida9

       xlmedio = 0.
       nm2 = nm*nm
       
       do i = 1,nm
        xlmedio = xlmedio + Lmedio(i)/nm2
       enddo

       write(9,500)entrada
       write(9,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na

       if(irede.lt.5)write(9,*)' pc = ',pc
       write(9,*)'di�metro da rede'
       write(9,520)npp
       write(9,*)'dist�ncia minima m�dia da rede'
       write(9,600)xlmedio
       write(9,*)'dist�ncia minima m�dia de cada n�'
       do i = 1,nm
        write(9,520)i,Lmedio(i)/nm
       enddo

c==============================================================================
c escreve a matriz de vizinhan�a da rede em saida11 

       open(unit=11,file=saida11)

       write(11,500)entrada
       write(11,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na,
     . '  diametro = ',npp
       write(11,*)'matriz de vizinhan�a ^m '

       if(irede.lt.5)write(11,*)' pc = ',pc
        
       do i = 1,nm
        write(11,550)(am(i,j),j=1,nm)									 
       enddo
c	write(*,*)am(1,2),am(1,130),am(129,2),am(129,130)
c==============================================================================
c escreve coeficientes de assortatividade em saida12

       write(12,500)entrada
       write(12,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na
       write(12,*)'coeficientes de assortatividade A(L), L=1,...,',npp

       if(irede.lt.5)write(12,*)' pc = ',pc

       write(12,560)((k),k=1,npp+1)
       write(12,590)(rk(i),i=1,npp+1)

c==============================================================================
c      escreve dados para c�lculo da dimens�o fractal em saida13

       if(npg.eq.0)then

        write(13,500)entrada
        write(13,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na
        write(13,*)'dados para c�lculo da dimens�o fractal n(L), L=1,...
     .  .,',npp

        if(irede.lt.5)write(13,*)' pc = ',pc

        write(13,*)'      log(n(L))     log(L) ' 

        do i = 1,npp
         write(13,580)dfr(i,1),dfr(i,2)
        enddo

       endif

c==============================================================================
c determina a matriz de vizinhan�a projetada da multiplex a partir da matriz de vizinhanca 
c da rede onde os n�s em cada camada � identificado individualmente 

       do i = 1,nms
        do j = 1,nms
         do k = 1,ncamada
          do l = 1,ncamada
           am(i,j) = min(am(i+(k-1)*nms,j+(l-1)*nms),am(i,j))
c	write(*,*)i,j,k,l,am(i,j)
c	write(*,*)i+(k-1)*nms,j+(l-1)*nms,am(i+(k-1)*nms,j+(l-1)*nms)
c	pause
          enddo
         enddo
        enddo
       enddo

       do i = 1,nm
        am(i,i) = 0
       enddo

c==============================================================================
c escreve a matriz de vizinhan�a projetada da multiplex em saida14 

       open(unit=14,file=saida14)

       write(14,500)entrada
       write(14,*)'is =',is,' irede =',irede,' nm =',nm,', e na =',na,
     . '  diametro = ',npp
       write(14,*)'matriz de vizinhan�a da multiplex '

       if(irede.lt.5)write(14,*)' pc = ',pc
        
       do i = 1,nms
        write(14,550)(am(i,j),j=1,nms)
       enddo

c==============================================================================
      enddo

      close(unit=4)
      close(unit=7)
      close(unit=8)
      close(unit=9)
      close(unit=11)
      close(unit=12)      
      close(unit=13)      
      close(unit=14)      

c      tempo2 = cpsec()
c      write(20,*)'is=',is,' irede=',irede,' nm=',nm, 
c     .' tfinal=',tempo2-tempo1

      goto 10

500   format(a30)
505   format(a13)
510   format(10000i1)
520   format(i5,300(2x,e13.7))     
530   format(i5,300(1x,i4))     
540   format(5h Cmed,300(2x,e13.7))     
550   format(8768i3)
560      format(5x,300(6x,i2,7x))
570   format(5h kmed,300(2x,e13.7))     
580   format(300(2x,e13.7))     
590   format(5h Amed,300(2x,e13.7))     
600   format(6h<dmin>,1x,e13.7)     
      end
            
c==============================================================================
c subroutine clust(a,n,np,nvm,ga,ip,ns)
c calcula o grau e o coeficiente de aglomera��o de uma rede
c==============================================================================

      subroutine clust(a,n,np,ga)
      integer*2 a(np,np)
      real ga(0:np)
      integer vvz(np)

c     calcula o coeficiente de aglomera��o de cada n�

      ga(0) = 0.

      do i = 1,n
       
       do j = 1,n  
        vvz(j) = -1
       enddo
      
       ordem = 0
       do j = 1,n
        if (a(i,j).eq.1) then
         ordem = ordem+1
         vvz(int(ordem)) = j
        endif
       enddo

       ctotal = ((ordem * (ordem-1))/2)
       sumg = 0
       
       do j = 1,int(ordem)
        li = vvz(j)
        do k = 1,n
         if (a(li,k).eq.1.and.a(i,k).eq.1.and.li.ne.k.and.k.ne.i)then
          sumg = sumg+1
         endif
        enddo
       enddo
      
       clocal = (sumg)/2      
       if (ctotal.gt.0) then
        ga(i) = clocal/ctotal
       else
        ga(i) = 0
       endif
       ga(0) = ga(0) + ga(i)
      
      enddo 

      return
      end

c fim da subroutine clust
c==============================================================================

c==============================================================================
c subroutine assorta(a,n,np,nvm,kk,ip,rk)
c calcula o grau de assortatividade de uma rede
c==============================================================================

      subroutine assorta(a,n,np,nv,rk,kk,ip)
      integer kk(0:np,0:nv)
      integer*2 a(np,np)

      s1 = 0.
      s2 = 0.
      s3 = 0.
      na = 0

      do i = 1,n-1
       do j = i+1,n
        if (a(i,j).eq.1) then
         s1 = s1 + kk(i,ip)+kk(j,ip)
         s2 = s2 + kk(i,ip)*kk(j,ip)
         s3 = s3 + kk(i,ip)*kk(i,ip)+kk(j,ip)*kk(j,ip)
         na = na + 1
        endif
       enddo
      enddo

      s1 = s1/na
      s2 = s2/na
      s3 = s3/na

c      calcula coeficiente de assortatividade
      
      rkk = (0.5*s3 - (0.5*s1)**2)
      
      if (rkk.ne.0) rkk = (s2 - (0.5*s1)**2)/rkk

      rk = rkk
 
      return
      end

c fim da subroutine assorta
c==============================================================================

c==============================================================================
c subroutine fracnet(a,n,np,diam,dfr)
c calcula a dimensao fractal da rede
c==============================================================================
      subroutine fracnet(a,n,np,diam,dfr)
      
      integer*2 a(np,np)
      integer diam,fim
      dimension dfr(1000,0:2),labelb(0:n,0:n),numberb(0:n)

      aux2 = diam
      aux = n

      do i = 1,diam
       dfr(i,1) = 0
       dfr(i,2) = 0
      enddo

      do iref = 1,int(aux2)

       do i = 0,n
        numberb(i) = 0
        do j = 0,n
         labelb(i,j) = 0
        enddo
       enddo

       nbox = 0

       do i = 1,n

        if (nbox.eq.0) then
         nbox = 1
          labelb(0,0) = 1
          numberb(0) = 1

        else

         j = 0
         fim = 0

         do while (j.lt.nbox.and.fim.eq.0)
          fim = 0
          do k = 0,numberb(j)-1
           if (a(i,labelb(k,j)).le.iref)then
            fim = fim + 1
           else
            go to 100
           endif
          enddo

100          continue

          if (fim.eq.numberb(j)) then

           fim = 1
           labelb(numberb(j),j) = i
           numberb(j) = numberb(j)+1

          else

           fim = 0

          endif

          j = j+1
      
         enddo 

         if (fim.eq.0) then

          labelb(0,nbox) = i
          numberb(nbox) = numberb(nbox) + 1
          nbox = nbox + 1

         endif

        endif

       enddo

       dfr(iref,1) = log(float(iref+1));
       dfr(iref,2) = log(float(nbox));

      enddo

      return

      end

c fim da subroutine fracnet
c==============================================================================

c==============================================================================
c subroutine tridi(a,n,np)
c gera uma rede tridiagonal
c==============================================================================

      subroutine tridi(a,n,np)
      integer*2 a(np,np)

      do i = 1,n
       do j = 1,n
        a(i,j) = 0
       enddo
      enddo

      do i = 1,n-1
       a(i,i+1) = 1
       a(i+1,i) = 1
      enddo

c      se os valores abaixo s�o 1, cadeia fechada, se s�o 0, cadeia aberta

      a(1,n) = 0
      a(n,1) = 0

      return
      end

c fim da subroutine tridi
c==============================================================================

c==============================================================================
c subroutine smalw(a,n,np,pc)
c gera uma rede de mundo pequeno
c==============================================================================

      subroutine smalw(a,n,np,pc)
      integer*2 a(np,np)
      
      do i = 1,n
       do j = 1,n
        a(i,j) = 0
       enddo
      enddo

      do i = 1,n-2
       a(i,i+1) = 1
       a(i+1,i) = 1
       a(i+2,i) = 1
       a(i,i+2) = 1
      enddo

      a(1,n) = 1
      a(n,1) = 1
      a(2,n) = 1
      a(n,2) = 1

      do i = 1,n
       x1 = rand(0)
       if (x1.lt.pc) then
        j1 = int(rand(0)*n+1)
        a(i,j1) = 1
        a(j1,i) = 1
       endif
      enddo

      do i = 1,n
       a(i,i) = 0
      enddo

      return
      end

c fim da subroutine smalw
c==============================================================================

c==============================================================================
c subroutine rando(a,n,np,pc)
c gera uma rede aleatoria
c==============================================================================

      subroutine rando(a,n,np,pc)
      integer*2 a(np,np)

      xlb = 2*n
      xup = n*sqrt(float(n))

c     aqui o numero de conexoes pode variar a cada amostra

      xt = rand(0) 
      nt = int(xlb + xt*(xup-xlb))

c     aqui o numero de conexoes � fixo = pc*nm
      if (pc.gt.0) nt = pc * n * (n-1) / 2

      do i = 1,n
       do j = 1,n
        a(i,j) = 0
       enddo
      enddo

      do i = 1,nt
       i1 = int(rand(0)*n+1)
       j1 = int(rand(0)*n+1)
       a(i1,j1) = 1
       a(j1,i1) = 1
      enddo

      do i = 1,n
       a(i,i) = 0
      enddo

      return
      end

c fim da subroutine rando
c==============================================================================

c==============================================================================
c subroutine scafr(ma,nvert,npm)
c gera uma rede livre de escala
c==============================================================================

      subroutine scafr(ma,nvert,npm)
      integer tgraus
      integer*2 ma(npm,npm)

c      o parametro m deve ser no minimo 4      

      parameter (m=4)
      integer grau(nvert)

      i = 0
      j = 0
      k = 0

c      nosaleat: numero de nos aleatorios iniciais 
c      nnos: numero total de arestas do grafo

      nosaleat = m
      nnos = 0
      tgraus = 0

      do i = 1,nvert
       grau(i) = 0
       do j = 1,nvert
        ma(i,j) = 0
       enddo
      enddo

c      gera pequena rede aleatoria inicial
      
      nosaleat = m
      nnos = 1
      
      do i = 2,nosaleat
       do k = 1,nnos
        j = int(rand(0)*nnos+1)
        if (i .ne. j .and. ma(i,j) .eq. 0) then
         ma(i,j) = 1
         ma(j,i) = 1
         grau(i) = grau(i)+1
         grau(j) = grau(j)+1
         tgraus = tgraus+2
        endif
       enddo
       nnos = nnos+1
      enddo

c      varre o restante do Grafo conectando os vertices

      do i = nosaleat+1,nvert

c      m arestas para cada n� novo

       k = m
       do while(k .gt. 0)
        j = int(rand(0)*(nnos)+1)

c      calcula probabilidade de um n� ser conectado;
        
        pr = float(grau(j))/tgraus
        
        if (pr.gt.rand(0).and.i.ne.j.and.ma(i,j).eq. 0) then
         ma(i,j) = 1
         ma(j,i) = 1
         grau(i) = grau(i)+1
         grau(j) = grau(j)+1
         k = k-1
         tgraus = tgraus+2
        endif
       
       enddo
       
       nnos = nnos+1
      enddo

      return

      end

c fim da subroutine scafr
c==============================================================================

c==============================================================================
c subroutine apolo(a,n,ngm,np)
c gera uma rede apoloniana
c==============================================================================

      subroutine apolo(a,n,ngm,np)

      integer ng(n) 
c      double precision a(np,np)
      integer*2 a(np,np)

      ngm = n

      do i = 1,ngm
       ng(i) = (3**(i-1)+5)/2
      enddo
      
      n = ng(ngm)

      do i = 1,ng(ngm)
       do j = 1,ng(ngm)
        a(i,j) = 0
       enddo
      enddo  
      
      do i = 1,4
       do j = i+1,4
        a(i,j) = 1
        a(j,i) = 1
       enddo
      enddo  

      do ig = 3,ngm
      
       ng0 = ng(ig)
       ng1 = ng(ig-1)
       ng2 = ng(ig-2)
       ng3 = 2*ng1 - 3
       
       do i = 1,ng1 - 3
        do j = 1,ng1 - 3
       
         a(ng1+i,ng1+j) = a(2+i,2+j)
         a(ng3+i,ng3+j) = a(2+i,2+j)
        
        enddo
       enddo
       
       do j = 1,ng1 - 3
       
        a(1,ng1+j) = a(1,2+j)
        a(ng1,ng1+j) = a(ng1,2+j)
        a(ng0,ng1+j) = a(2,2+j)
        a(2,ng3+j) = a(2,2+j)
        a(ng1,ng3+j) = a(ng1,2+j)
        a(ng0,ng3+j) = a(1,2+j)

        
        a(ng1+j,1) = a(1,ng1+j)
        a(ng1+j,ng1) = a(ng1,ng1+j)
        a(ng1+j,ng0) = a(ng0,ng1+j)
        a(ng3+j,2) = a(2,ng3+j)
        a(ng3+j,ng1) = a(ng1,ng3+j)
        a(ng3+j,ng0) = a(ng0,ng3+j)

       enddo

       a(ng0,1) = a(ng1,1)
       a(ng0,2) = a(ng1,2)
       a(ng0,ng1) = a(ng1,ng2)
       
       a(1,ng0) = a(ng1,1)
       a(2,ng0) = a(ng1,2)
       a(ng1,ng0) = a(ng1,ng2)

      enddo
       
500   format (10000i1)
510   format (1x,10000i1)

      end

c fim da subroutine apolo
c==============================================================================

