!     programa multiNG ultima vers‹o
!     calcula o dendograma de um multiplex a partir da eliminação sucessiva
!     dos nós com maior grau de betweeness 
!     não faz a eliminação das conexões entre as camadas
!     faz renumeracao dinamica dos nos, de modo que todos os nos de um cluster
!     sejam numerados contiguamente 
!     a mudança da numeração é feita simultaneamente para todas as camadas do multiplex
!     constroi arquivo para tracar dendograma no origin
!     calcula a modularidade Q definida por Newman e Girvan e generalizada por Mucha
!     entrada dos dados via matriz de vizinhança
!     usa rotina rannyu para gerar numeros aleatorios
!     duas maneiras de gerar sementes independentes, para windows e linux
!     igual a dendo1, mas calcula a distância entre duas matrizes de vizinhancas
!     após a eliminação de uma conexão

!     criado em 26/11/2019, a partir de mdendo2u(de 28/12/2018) - comando 510 format em i4 -
!     faz a gravacao em arquivos do resultado parcial apos a eliminacao de um numero escolhido de ligacoes, 
!     permitindo reiniciar o processo caso o computador falhe por falta de luz durante o processo. 
!     inclui a correção para o problema da não conexão de um cluster que é quebrado pelo primeiro link eliminado 

!      use portlib 
      parameter(npm=5000,nml=2000000)
      integer icc(npm),iord(npm,4),ila(npm),inc(npm)
	integer am(npm,npm),bm(npm,npm),lar(npm),Em(npm,npm),kk(npm)
      real yc(0:npm) !cm(npm,npm),yc(0:npm)
	real*8 rannyu ,cm(npm,npm),xb,xb1
      integer ie(nml),ig(0:nml),iv(2,nml),iq(npm),nlik(20)
	character*8 saida
	character*26 saida1,saida2,saida3,saida4,saida5,saida6,saida7
	character*26 entrada2,entrada3,entrada4,Entrada1
	common /rnyucm/ mk1,mk2,mk3,mk4,km1,km2,km3,km4

!==================================================================
!     entrada2: numeração dos nós - no inicio usa a numeração identidade. caso o programa seja
!     reinicializado, os dados devem ser os contidos no arquivo saida5='d2'//saida//'5.dat'
!     entrada3: matriz de vizinhança da rede já com os links eliminados. no inicio, esta entrada
!     eh identica a entrada 4. caso o programa seja reinicializado, os dados devem ser os contidos
!     no arquivo saida7='d2'//saida//'7.dat'
!     entrada4: matriz de vizinhança da rede sem os links eliminados. no inicio, esta entrada
!     eh identica a entrada 3. caso o programa seja reinicializado, os dados devem ser os contidos
!     no arquivo saida3='d2'//saida//'3.dat'
!     Entrada1: matriz de vizinhança da rede original sempre, mesmo se o programa for reinicializado
!     nm: numero de nos na rede
!     npula: numero de linhas antes dos elementos da matriz
!     ngrava: numero de ligações eliminadas necessarias para se proceder a uma noma
!             gravacao de dados para reinicializar o programa
!     iran: se 0, escolhe sempre o primeiro link com maximo grau de intermediacao
!           se 1, usa numero aleatorio para sortear o link 

!     saida1: informa, no minimo, a sequencia das conexoes eliminadas 
!             pode-se habilitar outras informaçoes sobre a evolucao do processo
!             contudo, para redes grandes, o arquivo pode tornar-se muito grande
!     saida2: informa o cluster ao qual cada nó pertence, para cada valor de itera
!             na ultima linha informa a numeracao inicial de cada nó
!     saida3: guarda a matriz de vizinhança com nós renumerados de acordo com a 
!             sequencia de eliminaçoes de conexoes
!     saida4: dados para o dendograma. 
!     saida5: coluna unica, com a mesma informacao que na ultima linha de saida2
!             usada para ilustrar dendograma feito com o origin
	
!     saida6: distancia entre duas matrizes de vizinhancas sucessivas em funcao de itera
!     saida7: guarda a matriz de vizinhança com links eliminados e nós renumerados de acordo com a 
!             sequencia de eliminaçoes de conexoes
!==================================================================
!     reinicializacao

!     caso o programa seja interrompido por algum motivo, o processo pode ser reinicializado 
!     tomado os seguintes passos:

!     1) renomear o arquivo saida5='d2'//saida//'5.dat' para ser utilizado como novo entrada2 
!     2) apagar a primeira linha deste arquivo
!     3) renomear o arquivo saida7='d2'//saida//'7.dat' para ser utilizado como novo entrada3
!     4) alterar o arquivo de entrada edendo2t.dat, passando npula de 3 para 1  
!     5) renomear o arquivo saida3='d2'//saida//'3.dat' para ser utilizado como novo entrada4
!     6) alterar o arquivo de entrada edendo2t.dat, passando npula de 3 para 1  
!     7) opcional: mudar o nome dos arquivos de saida (variável "saida") em edendo2t.dat
!==================================================================

!     constantes para a rotina rnyucm e nao podem ser trocados
      mk1=0
      mk2=1
      mk3=3513
      mk4=821

!      sementes do aleatorio, devem ser trocados para cada run

      km1=3577 
      km2=5735
      km3=4689
      km4=9087

!     entrada de dados
!==================================================================
      open (unit=3,file='multiNG.dat')

10    read (3,*,end=1000)nnos,npula,iseed,ngrava,ncamada,iran
	if(nnos.lt.0)stop
	nms = nnos
      nm = nms*ncamada

!     altera ou nao semente do aleatorio
!==================================================================

!     comandos para gerar semente aleatoria em linux
	if (iseed.eq.-1) then
      intime = time()
      km = intime/10000
      km = intime - km*10000
	else
	 km = iseed
	endif

!     comandos para gerar semente aleatoria em windows
!      call seed(iseed)

!      km = 10000*rand(0)

      km1 = km + km1
      km2 = (km + km2)*(km + km2)/km1
      km3 = (km + km3)*(km + km3)/km2
      km4 = (km + km4)*(km + km4)/km3

!     termina alteracao da semente do aleatorio
!     continua entrada de dados
!==================================================================
	read (3,500)Entrada1
	read (3,500)entrada2
	read (3,500)entrada3 
	read (3,500)entrada4
	read (3,505)saida
	saida1 = 'mdQu'//saida//'1.dat'
	saida2 = 'mdQu'//saida//'2.dat'
	saida3 = 'mdQu'//saida//'3.dat'
	saida4 = 'mdQu'//saida//'4.dat'
	saida5 = 'mdQu'//saida//'5.dat'
	saida6 = 'mdQu'//saida//'6.dat'
	saida7 = 'mdQu'//saida//'7.dat'

!     entrada da matriz am
!==================================================================
! 	open (unit=1,file=Entrada1)
	open (unit=2,file=entrada2)
	open (unit=4,file=entrada3)
	open (unit=16,file='tempdendo1s.dat',form='unformatted')
	 
	do i = 1,npula
	 read(4,*)
	enddo
	   
	nlink = 0

	do i = 1,nm
	 read(4,*)(lar(j),j=1,nm)
	 do j = 1,nm
	  am(i,j) = lar(j)
	  if (am(i,j).eq.1)nlink = nlink + 1
	 enddo
	enddo

      do k = 1,ncamada
	 nlik(k) = 0
	 do i = (k-1)*nms+1,k*nms
	  do j = (k-1)*nms+1,k*nms
	   if (am(i,j).eq.1) then 
	    nlik(k) = nlik(k) + 1
	   endif
	  enddo
	 enddo
	 nlik(k) = nlik(k)/2
	enddo

	close(unit=4)

	nlink = nlink/2
	nlink1 = nlink
	write (*,*)'nlink = ',nlink,(nlik(k),k=1,4)
	if (nlink.gt.nml)stop
!==================================================================
!     entrada da matriz bm

	open (unit=4,file=entrada4)
	 
	do i = 1,npula
	 read(4,*)
	enddo
	   
	do i = 1,nm
	 read(4,*)(lar(j),j=1,nm)
	 do j = 1,nm
	  bm(i,j) = lar(j)
	 enddo
	enddo

	close(unit=4)
!==================================================================
!     entrada da matriz Em

	open (unit=4,file=Entrada1)
	 
	do i = 1,npula
	 read(4,*)
	enddo
	   
	do i = 1,nm
	 read(4,*)(lar(j),j=1,nm)
	 do j = 1,nm
	  Em(i,j) = lar(j)
	 enddo
	enddo

	do i = 1,nm
	 kk(i) = 0
	 do j = 1,nm
	  if (Em(i,j).eq.1) kk(i) = kk(i) + 1
	 enddo
	 kk(i) = kk(i) - ncamada + 1
	enddo

	close(unit=4)
!==================================================================
!     grava memoria do ordenamento inicial dos sitios

	do i = 1,nm
 	 read(2,*)iord(i,4)
	enddo

	open (unit=5,file=saida1)
	open (unit=6,file=saida2)
	open (unit=8,file=saida4)
	open (unit=9,file=saida5)
	open (unit=10,file=saida6)

	write(5,*)' resultados de mdend2Qu com arquivo ',entrada4
	write(6,*)' resultados de mdend2Qu com arquivo ',entrada4
	write(9,*)' resultados de mdend2Qu com arquivo ',entrada4
	write(10,*)' resultados de mdend2Qu com arquivo ',entrada4

	close(unit=5)
	close(unit=6)
	close(unit=8)
	close(unit=9)
	close(unit=10)
!     começa o grande loop para eliminação dos links com maior betweenness
!==================================================================

	id1 = 1
	itera = 0
	
      do while(id1.gt.0)

	 open (unit=5,file=saida1,access='append')
	 open (unit=6,file=saida2,access='append')
	 open (unit=9,file=saida5,access='append')
	 open (unit=10,file=saida6,access='append')

	 itera = itera + 1

	 if (itera.gt.nml)then
	  write (*,*) 'itera>nml'
	  stop
	 endif

!     determina o diametro
!==================================================================

       do i = 1,nm
	  do j = 1,nm
	   cm(i,j) = 0.
	  enddo
	 enddo

       do i = 1,nm
	  if (am(i,i).ne.0)write(*,*)itera,i,am(i,i)
	  if (am(i,i).ne.0)write(5,*)itera,i,am(i,i)
	 enddo

	 id1 = 0
	 iedge = 0
       do i = 1,nm
	  do j = 1,nm
	   if (am(i,j).gt.id1) id1 = am(i,j)
	   if (am(i,j).eq.1)cm(i,j) = am(i,j)
	   if (am(i,j).eq.1)iedge = iedge + 1
	  enddo
	 enddo

	if (iedge.eq.nms*ncamada*(ncamada-1))id1=0 

!     identificação dos clusters antes da primeira eliminacao
!	identificação inicial dos clusters
!==================================================================
	 if (itera.eq.1)then
	  do i = 1,nm
	   icc(i) = 0
	   inc(i) = 0
	  enddo

	  id = 1	 					   
	  nc = 0

	  do i = 1,nm
	   if(nc.lt.nm) then
	    ic = 0
	    do j = 1,nm
	     if(icc(j).eq.0) then
	      if(am(i,j).ne.0.or.i.eq.j) then
	       icc(j) = id
	       nc = nc + 1
	       ic = 1
	      endif
	     endif
	    enddo
	    if (ic.eq.1) id = id + 1
	   endif
	  enddo

	  incm = 0
        do i = 1,nm
	   inc(icc(i)) = inc(icc(i)) + 1
	   incm = max(incm,icc(i))
	  enddo

	  nden = 0
	  do i = 1,incm
	   nden = nden + inc(i)*(inc(i)-1)/2
	  enddo

!     renumeracao dos sitios antes da primeira eliminacao
!==================================================================

	  do i = 1,nms
	   iord(i,1) = i
	   iord(i,2) = icc(i)
	  enddo

	  call orderset(iord,nms,4,2)

	  iwc = 0

	  do i = 1,nms
	   iord(i,3) = i
	   ila(i) = i
	  enddo

15	  iw = 0
	  iwc = iwc + 1
	   do i = 1,nms
	    if (iord(i,1).ne.iord(i,3)) then
	     i1 = iord(i,1)
	     i3 = iord(i,3)
	     j1 = i
	     j2 = ila(i1)
	     do j = 1,nms
	      if(j.ne.j2.and.j.ne.j1)then
	      iw = 1
	      do k1 = 1,ncamada
	      do k2 = 1,ncamada
	     mtemp = am(j1+(k1-1)*nms,j+(k2-1)*nms)
	     am(j1+(k1-1)*nms,j+(k2-1)*nms)=am(j2+(k1-1)*nms,j+(k2-1)*nms)
 	     am(j2+(k1-1)*nms,j+(k2-1)*nms)=mtemp
	     am(j+(k2-1)*nms,j1+(k1-1)*nms)=am(j1+(k1-1)*nms,j+(k2-1)*nms)
	     am(j+(k2-1)*nms,j2+(k1-1)*nms)=am(j2+(k1-1)*nms,j+(k2-1)*nms)
	     mtemp = bm(j1+(k1-1)*nms,j+(k2-1)*nms)
	     bm(j1+(k1-1)*nms,j+(k2-1)*nms)=bm(j2+(k1-1)*nms,j+(k2-1)*nms)
 	     bm(j2+(k1-1)*nms,j+(k2-1)*nms)=mtemp
	     bm(j+(k2-1)*nms,j1+(k1-1)*nms)=bm(j1+(k1-1)*nms,j+(k2-1)*nms)
	     bm(j+(k2-1)*nms,j2+(k1-1)*nms)=bm(j2+(k1-1)*nms,j+(k2-1)*nms)
	      enddo
	      enddo
	      endif
	     enddo
	     iord(ila(i1),3) = i3
	     iord(i,3) = i1
	     mtemp = ila(i1)
	     ila(i1) = ila(i3)
	     ila(i3) = mtemp

	    endif
	   enddo
	  if (iw.eq.1) goto 15
 
	 endif

!     comeca a calcular o grau de betweeness
!==================================================================

	 do id = 1,id1-1
	  i1 = id1 - id + 1 
	 
	  do i = 1,nm
         do j = i+1,nm
	    if(am(i,j).eq.i1) then
	     il = 0
	     do l = 1,nm
	      if (am(i,l).eq.1.and.am(j,l).eq.i1-1.or.
     .	   am(j,l).eq.1.and.am(i,l).eq.i1-1) il=il + 1
	     enddo
	     if (il.gt.0) then
	      do l = 1,nm
	       if (am(i,l).eq.1.and.am(j,l).eq.i1-1.or.
     .	    am(j,l).eq.1.and.am(i,l).eq.i1-1)then
	         cm(i,l) = cm(i,l) + (cm(i,j) + 1.)/il 
	         cm(j,l) = cm(j,l) + (cm(i,j) + 1.)/il 
		     cm(l,i) = cm(i,l)
	         cm(l,j) = cm(j,l)
		   endif
	      enddo
	     endif
	    endif
	   enddo
	  enddo

       enddo

!     corrige matriz de betweenness
!==================================================================

	 do i = 1,nm
	  do j = 1,nm
	   if (am(i,j).ne.1)cm(i,j)=0
	  enddo
	 enddo
	
	 do i = 1,ncamada-1
	  do j = i,ncamada
         do k = 1,nms
	    cm(k+(i-1)*nms,k+(j-1)*nms) = 0
	    cm(k+(j-1)*nms,k+(i-1)*nms) = 0
	   enddo
	  enddo
	 enddo

!     procura o sitio	com o maior grau de betweeness
!=====================================================

	 xb = 0. ; xb1 = 0.
	 inb = 0
	 do i = 1,nm-1
	  do j = i+1,nm
	   if (am(i,j).eq.1.and.cm(i,j).gt.0) then
	    xb1 = cm(i,j)
!!	    if (xb1.gt.xb) then
	    if (xb1-xb.gt.0.0001) then
	     xb = xb1
	     inb = 1
	     iv(1,inb) = i
	     iv(2,inb) = j

	    else if (abs(xb1-xb).lt.1e-4) then
	     inb = inb + 1
	     iv(1,inb) = i
	     iv(2,inb) = j

	    endif

	   endif
	  enddo
	 enddo

!     transfere a matriz am para cm
!=====================================================
	 do i = 1,nm
	  do j = 1,nm
	   cm(i,j) = am(i,j)
	  enddo
	 enddo

!     elimina um link
!=====================================================

	 if (inb.eq.1.or.iran.eq.0) then

	  inbb = 1
	  am(iv(1,inbb),iv(2,inbb)) = 0
	  am(iv(2,inbb),iv(1,inbb)) = 0
	  i1 = iv(1,inbb)
	  j1 = iv(2,inbb)

	 else if (inb.gt.1.and.iran.eq.1) then

	  xt = rannyu()
	  il = min(int(xt*inb)+1,inb)
	  am(iv(1,il),iv(2,il)) = 0
	  am(iv(2,il),iv(1,il)) = 0
	  i1 = iv(1,il)
	  j1 = iv(2,il)

	 endif

	 ie(1) = i1
	 ie(2) = j1
	 iv(1,1) = i1
	 iv(2,1) = j1
	 nie = 2
	 nv = 1

!     calcula os efeitos da eliminacao do link (i1,j1)
!=====================================================

	 write(5,*)' itera',itera,' link elim', i1,j1,
     . iord(i1,4),iord(j1,4)

	 iel = 1

	 do while (iel.ne.0)

	  iel = iel + 1
	  id = iel
	  nnie = 0
	  imid = 0

	  do ili = 1,nie
	   ia = ie(ili)
	   do l = 1,nm
	    if (am(ia,l).eq.id) then
	     ip = 0
	     do it = 1,nm
	      if(am(ia,it).eq.1.and.am(l,it).eq.id-1.
     .	  or.am(l,it).eq.1.and.am(ia,it).eq.id-1) ip = 1
	     enddo
	     if (ip.eq.0) then
	      am(ia,l) = 0
	      am(l,ia) = 0
	      nnie = nnie + 1
	      nv = nv + 1

            if (nv.gt.nml.or.nie+nnie.gt.nml) then
             write(5,*)'limite nml ultrapassado. nv, nie, nnie'
     .	   ,nv, nie, nnie
             stop
            endif
	
	      iv(1,nv) = min(ia,l)
	      iv(2,nv) = max(ia,l)
	      ie(nie+nnie) = l
	      imid = 1
	     endif
	    endif
	   enddo
	  enddo
	  nie = nie + nnie
	 if (imid.eq.0.and.id.gt.id1) ied = iel
	 if (imid.eq.0.and.id.gt.id1) iel = 0
	 enddo
 
!     reconstroi a matriz de vizinhança
!=====================================================

	 iel = ied + 1
	 iel = ied*0+2
	 do while (iel.ne.0)
!!	if(iel.le.1)write(*,*)'iel',iel
	  id = 0
	  do iiv = 1,nv
	   i1 = iv(1,iiv)
	   j1 = iv(2,iiv)
	   if(am(i1,j1).eq.0) then
	    id = 1
	    ip = 0
	    il = 1
          do while (ip.eq.0)
	     do it = 1,nm
	      if(am(i1,it).eq.il.and.am(j1,it).eq.iel-il.
     .      or.am(i1,it).eq.iel-il.and.am(j1,it).eq.il) ip = 1
	     enddo
	     if (ip.eq.1) am(i1,j1) = iel
	     if (ip.eq.1) am(j1,i1) = iel
	     il = il + 1  
	     if (il.gt.iel-1) ip = npm+1
	    enddo
	   endif
	  enddo
	  iel = iel + 1
	  if (iel.gt.2*ied.or.id.eq.0) iel = 0
	 enddo

!     calcula a distância entre duas matrizes sucessivas 
!==================================================================

	 dist = 0

	 do i = 1,nm-1
	  do j = i,nm
	   dist = dist + (cm(i,j) - am(i,j))**2
	  enddo
	 enddo

	 dist = sqrt(dist)
 
!     identificação dos clusters
!==================================================================
	 do i = 1,nms
	  icc(i) = 0
	  inc(i) = 0
	 enddo

	 id = 1
	 nc = 0

	 do i = 1,nms
	  if(nc.lt.nms) then
	   ic = 0
	   do j = 1,nms
	    if(icc(j).eq.0) then
	     ak = 0
	     do k = 1,ncamada
	      ak = ak + am(i+(k-1)*nms,j+(k-1)*nms)
	     enddo
	     if(ak.ne.0.or.i.eq.j) then
	      icc(j) = id
	      nc = nc + 1
	      ic = 1
	     endif
	    endif
	   enddo
	   if (ic.eq.1) id = id + 1
	  endif
	 enddo

	 incm = 0
       do i = 1,nms
	  inc(icc(i)) = inc(icc(i)) + 1
	  incm = max(incm,icc(i))
	 enddo

	 nden = 0
	 do i = 1,incm
	  nden = nden + inc(i)*(inc(i)-1)/2
	 enddo

!     renumera os sitios
!==================================================================

	 do i = 1,nms
	  iord(i,1) = i
	  iord(i,2) = icc(i)
	 enddo

	 call orderset(iord,nms,4,2)

	 iwc = 0

	 do i = 1,nms
	  iord(i,3) = i
	  ila(i) = i
	 enddo

20	 iw = 0
	 iwc = iwc + 1
	  do i = 1,nms
	   if (iord(i,1).ne.iord(i,3)) then
	    i1 = iord(i,1)
	    i3 = iord(i,3)
	    j1 = i
	    j2 = ila(i1)
	    do j = 1,nms
	     if(j.ne.j2.and.j.ne.j1)then
	     iw = 1
	     do k1 = 1,ncamada
	     do k2 = 1,ncamada
	     mtemp = am(j1+(k1-1)*nms,j+(k2-1)*nms)
	     am(j1+(k1-1)*nms,j+(k2-1)*nms)=am(j2+(k1-1)*nms,j+(k2-1)*nms)
 	     am(j2+(k1-1)*nms,j+(k2-1)*nms)=mtemp
	     am(j+(k2-1)*nms,j1+(k1-1)*nms)=am(j1+(k1-1)*nms,j+(k2-1)*nms)
	     am(j+(k2-1)*nms,j2+(k1-1)*nms)=am(j2+(k1-1)*nms,j+(k2-1)*nms)
	     mtemp = bm(j1+(k1-1)*nms,j+(k2-1)*nms)
	     bm(j1+(k1-1)*nms,j+(k2-1)*nms)=bm(j2+(k1-1)*nms,j+(k2-1)*nms)
 	     bm(j2+(k1-1)*nms,j+(k2-1)*nms)=mtemp
	     bm(j+(k2-1)*nms,j1+(k1-1)*nms)=bm(j1+(k1-1)*nms,j+(k2-1)*nms)
	     bm(j+(k2-1)*nms,j2+(k1-1)*nms)=bm(j2+(k1-1)*nms,j+(k2-1)*nms)
	     enddo
	     enddo
	     endif
	    enddo
	    iord(ila(i1),3) = i3
	    iord(i,3) = i1
	    mtemp = ila(i1)
	    ila(i1) = ila(i3)
	    ila(i3) = mtemp

	   endif
	  enddo
	 if (iw.eq.1) goto 20
 
	 write(6,511)itera,(iord(k,2),k=1,nms) 	  
       write(9,511)itera+0,(iord(k,4),k=1,nm) 	  
	 write(10,512)itera,dist,dist/max(1,nden)
	 write(16)itera,(iord(k,2),k=1,nms) 

!     grava resultados parciais 
!==================================================================

	 igrava = igrava + 1

	 if(igrava.eq.ngrava)then

	  open (unit=7,file=saida3)
        open (unit=9,file=saida5,access='append')
	  open (unit=11,file=saida7)

	  write(7,*)' matriz de vizinhança renumerada',itera
	  write(11,*)' matriz de vizinhança com link eliminado renumerada'
     . ,itera
	  igrava = 0

	  do i = 1,nm
	   write(7,510)(bm(i,j),j=1,nm)
	   write(11,510)(am(i,j),j=1,nm)
	  enddo
	 endif

       close(unit=5)
       close(unit=6)
       close(unit=7)
       close(unit=9)
       close(unit=10)
       close(unit=11)

      enddo

	open (unit=7,file=saida3)
	open (unit=9,file=saida5,access='append')
	open (unit=11,file=saida7)

	write(7,550)
	write(7,*)' resultados de mdend2Qu com arquivo ',entrada4
	write(7,*)' matriz de vizinhança renumerada',itera
	write(11,*)' resultados de mdend2Qu com arquivo ',entrada4
	write(11,*)' matriz de vizinhança com link eliminado renumerada'
     .,itera
	igrava = 0

	do i = 1,nm
	 write(7,510)(bm(i,j),j=1,nm)
	 write(11,510)(am(i,j),j=1,nm)
	enddo

      close(unit=7)
      close(unit=9)

      close(unit=16)

      nlink = itera - 1
!!	write(*,*)nlink,nlink1
      open (unit=6,file=saida2,access='append')

      write(6,510)
      write(6,511)100,(iord(k,4),k=1,nms) 	  
 
      close(unit=6)

!     calcula modularidade Q 
!==================================================================

!     identifica valores das coordenadas de cada cluster 
!==================================================================
	open (unit=9,file=saida5)  
	open (unit=16,file='tempdendo1s.dat',form='unformatted')

	read (9,*)
	write(11,*)

	ncomu =1

	do i = 1,nlink

	 do j = 1,nm
	  ie(j) = 0
	  lar(j) = 0
	 enddo

	 read (9,*)itera,(ie(k1),k1=1,nm)
	 read (16)itera,(lar(k1),k1=1,nms) 
	 do j = 1,nms
	  k = ie(j)
	  iq(k) = lar(j)
	 enddo

	 do j = 2,ncamada
	  do k = 1,nms
	   iq((j-1)*nms+k) = iq(k)
	  enddo
	 enddo

	 if (i.eq.1) then

	  qq = functionqq(Em,iq,nms,ncamada,nlink1,nlik,kk)
	  write(11,*)itera,qq, ncomu

	 elseif	(lar(nms).gt.ncomu) then

	  ncomu = lar(nms)
	  qq = functionqq(Em,iq,nms,ncamada,nlink1,nlik,kk)
	  write(11,*)itera,qq,ncomu

	 endif

	enddo

      close(unit=11)
      close(unit=9)
      close(unit=16)

!     constroi dendograma
!==================================================================

!     identifica valores das coordenadas de cada cluster 
!==================================================================

	open (unit=16,file='tempdendo1s.dat',form='unformatted')
	open (unit=17,file='tempdendo2s.dat',form='unformatted')

	itera = 0
	ig(0) = 1
	xlm = 0.5*(1.+ float(nms))

	do k = ii,j-1
	 yc(k) = xlm
	enddo

	write (17)itera,(yc(k1),k1=1,nms)

	do i = 1,nlink

	 read (16)itera,(lar(k1),k1=1,nms)
	 xlm = float(lar(1))
	 ii = 1
	 xne = 1.
	 ig(i) = 1
	 do j = 2,nms
	  if(lar(j).eq.lar(j-1))then

	   xlm = xlm*xne/(xne+1) + float(j)/(xne+1)
	   xne = xne + 1

	   if(j.eq.nms) then
	    do k = ii,j
	     yc(k) = xlm
	    enddo
	   endif

	  else

	   ig(i) = ig(i) + 1
	   do k = ii,j-1
	    yc(k) = xlm
	   enddo
	   xlm = float(j)
	   ii = j
	   xne = 1.
	   if(j.eq.nms) yc(k) = xlm

	  endif

	 enddo

	 write (17)itera,(yc(k1),k1=1,nms)

	 do k = 1,nms
	  yc(k) = 0.
	 enddo

	enddo

      close(unit=16)
      close(unit=17)

!     insere linhas extras nos valores de bifurcação e escreve dendograma
!==================================================================

      open (unit=8,file=saida4)
	open (unit=17,file='tempdendo2s.dat',form='unformatted')

	write(8,*)' resultados de mdend2Qu com arquivo ',entrada4

	do i = 0,nlink-1
	 read (17)itera,(yc(k1),k1=1,nms)
	 write(8,530)itera,(yc(k1),k1=1,nms)
	 if(ig(i).lt.ig(i+1))write(8,530)itera+1,(yc(k1),k1=1,nms) 
	enddo

	read (17)itera,(yc(k1),k1=1,nms)
	write(8,530)itera,(yc(k1),k1=1,nms)

      open (unit=9,file=saida5,access='append')
	
      do i = 1,nm
	 write(9,511)iord(i,4)
	enddo

	write(9,511)			    

      close(unit=8)
      close(unit=9)

	goto 10

500   format(a26)
505   format(a8)
510   format(10000i4)
511   format(i6,10000i4)
512   format(i6,2(2x,e13.6))
520   format(10000(2x,e11.3))     
530   format(i5,1x,10000(2x,e11.3))
540   format(/,2(i5,2x))
550   format(' resultados de mdend2Qu com arquivo ')     
600	format(3(i4,1x),e15.8,1x,3(i4,1x))	
1000	stop
      end
            

!======================================================
!======================================================
      subroutine orderset(set,n,m,icol)
!     it orders, descending, a set of n numbers, conserving the respective
!     positions of other sets
!     m indicates the number of columns in the set
!     icol indicates the column which will be considered for the ordering  
      parameter(npm=5000)

      integer set(npm,m),x1 
      itera = 0
10    icont = 0

      itera = itera + 1

      do  20  j= 1,n-1

!     to order in ascending order use the following command
      if (set(j,icol).gt.set(j+1,icol)) then

!     to order in descending order use the following command
!      if (set(j,icol).lt.set(j+1,icol)) then

      icont = 1
      do k = 1,m
      x1 = set(j,k)
      set(j,k) = set(j+1,k)
      set(j+1,k) = x1
      enddo

      endif

20    continue

      if (icont.ne.0) go to 10

      return
      
      end      


!======================================================
!======================================================
!
      function rannyu()
      implicit double precision (a-h,o-z)
!        real*8 rannyu, twom12
      parameter (twom12 = 1/4096.d0)
      common /rnyucm/ mk1,mk2,mk3,mk4,km1,km2,km3,km4
!
!     this is rannyu as modified by a. sokal 9/26/85.
!     it is linear congruencial with modulus m = 2**48,incremen_tempoc=1,
!     and multiplier a = (2**36)*mk1 + (2**24)*mk2 + (2**12)*mk3 + mk4.
!     the multiplier is stored in common (see subroutine setrn)
!     and is set to a = 31167285 (recommended by knuth, vol. 2,
!     2nd ed., p. 102).
!
      i1 = km1*mk4 + km2*mk3 + km3*mk2 + km4*mk1
      i2 = km2*mk4 + km3*mk3 + km4*mk2
      i3 = km3*mk4 + km4*mk3
      i4 = km4*mk4  +  1
      km4 = mod(i4, 4096)
      i3 = i3 + i4/4096
      km3 = mod(i3, 4096)
      i2 = i2 + i3/4096
      km2 = mod(i2, 4096)
      km1 = mod(i1 + i2/4096, 4096)
      rannyu = twom12*(dble(km1)+twom12*(dble(km2)+
     .         twom12*(dble(km3)+twom12*(dble(km4)))))
      return
      end

!======================================================
!======================================================

!======================================================
	function functionqq(am,iq,nms,ncamada,nlink,nlik,kk) 
      
      parameter(npm=5000,nmaxcomu=100)
	integer am(npm,npm)
      integer iq(npm),kk(npm),nlik(20)

	nm = nms*ncamada
	nlink0 = (nlink - ncamada*(ncamada-1)*nms/2)/ncamada
	nlink2  = (nlink - ncamada*(ncamada-1)*nms/2)
      qq = 0.
	qq1 = 0.
	qq2 = 0.
	qq3 = 0.

	cjxy = 1

      do k = 1,ncamada
	 do i = (k-1)*nms+1,k*nms
	  do j = (k-1)*nms+1,k*nms
	   if (iq(i).eq.iq(j)) then
	    if (am(i,j).eq.1) then
	     qq1 = qq1 + am(i,j) 
	    endif
	    qq2 = qq2  + 0.5*float(kk(i)*kk(j))/nlik(k)
	    if (i.eq.j.and.k.eq.1) then
	     qq3 = qq3 + cjxy 
	    endif
	   endif
	  enddo
	 enddo
	enddo

	qq1 = 0.5*qq1/nlink
	qq2 = 0.5*qq2/nlink
	qq3 = 0.5*qq3*ncamada*(ncamada-1)/nlink

	qq = qq1 - qq2 + qq3

	functionqq = qq

	return
	end
!======================================================



