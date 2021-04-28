C  ====================================================================

C       Copyright:
C       ---------
C       This program was written by 

C       Ernesto Martins
C       Departamento de Matematica, Universidade de Coimbra
C       Apartado 3008
C       3000 Coimbra, PORTUGAL
C       e-mail: eqvm@mat.uc.pt

C       The program is available on "as is" basis. 
C       It is not guaranteed to be free of bugs, and the author
C       assume no responsibility for any potential problems. 
 
C       The code may be used and MODIFIED for NOT-FOR-PROFIT use.
C       (note the "function clock", for instance)

C       THIS NOTICE MUST BE REMAINED.   

C       January 1997 version

C  ====================================================================

	function clock ()
	real a(2), c
	c = etime (a)
	clock = a(1)
	end



	subroutine quick
C
C   Iterative quick sort algorithm code.
C
	implicit integer ( a - z )

	parameter ( Mnodes = 1010, Marcs = 509545, mMarcs = -Marcs )

	common/inputnet/tail(Marcs), head(Marcs), cost(Marcs), n, m
	common/outputnet/up(Mnodes), down(Mnodes)
	common/stack/stackL(mMarcs:Marcs), stackR(mMarcs:Marcs)


	s = 1
	stackR(1) = m
	stackL(1) = 1
5	L = stackL(s)
	R = stackR(s)
	s = s - 1
10	i = L
	j = R
	half = ( L + R ) / 2
	meio = cost(half)
20	if ( cost(i) .lt. meio ) then
		i = i + 1
		go to 20
	end if
30      if ( cost(j) .gt. meio ) then
		j = j - 1
		go to 30
	end if
	if ( i .le. j ) then
		aux = head(i)
		head(i) = head(j)
		head(j) = aux
		aux = tail(i)
		tail(i) = tail(j)
		tail(j) = aux
		aux = cost(i)
		cost(i) = cost(j)
		cost(j) = aux
		i = i + 1
		j = j - 1
	end if
	if ( i .le. j ) go to 20
	if ( i .lt. R ) then
		s = s + 1
		stackL(s) = i
		stackR(s) = R
	end if
	R = j
	if ( L .lt. R ) go to 10
	if ( s .ne. 0 ) go to 5

	end


	subroutine Kruskal
C
C   Input data is in the form:
C	tail node of kth arc, head node of kth arc, distance of kth arc. 
C
C	Marcs      max number of arcs
C	Mnodes     max number of nodes
C
C   Data is sorted by cost arcs. An iterative QUICK sort algorithm is
C   used.
C
	implicit integer ( a - z )

	parameter ( Mnodes = 1010, Marcs = 509545 )

	common/inputnet/tail(Marcs), head(Marcs), cost(Marcs), n, m
	common/result/tree(Mnodes), solucao
	common/work1/rot(Mnodes), list(Marcs), inicio(Mnodes), fim(Mnodes)

	estao = 0
	k = 1
	vai = 0
5	ant = tail(k)
	suc = head(k)
	if ( rot(ant) + rot(suc) .eq. 0 ) then
		estao = estao + 1
		fim(estao) = k
		inicio(estao) = k
		rot(ant) = estao
		rot(suc) = estao
		vai = vai + 1
		tree(vai) = k
		solucao = solucao + cost(k)
	else if ( rot(ant) .eq. 0 ) then
		pos = rot(suc)
		rot(ant) = pos
		list(k) = inicio(pos)
		inicio(pos) = k
		vai = vai + 1
		tree(vai) = k
		solucao = solucao + cost(k)
	else if ( rot(suc) .eq. 0 ) then
		pos = rot(ant)
		rot(suc) = pos
		list(k) = inicio(pos)
		inicio(pos) = k
		vai = vai + 1
		tree(vai) = k
		solucao = solucao + cost(k)
	else if ( rot(suc) .ne. rot(ant) ) then
		if (rot(suc) .lt. rot(ant)) then
			pos = rot(suc)
			rotulo = rot(ant)			
		else
			pos = rot(ant)
			rotulo = rot(suc)
		end if
		no = inicio(pos)
10              rot(head(no)) = rotulo
                rot(tail(no)) = rotulo
		if ( no .ne. fim(pos) ) then
			no = list(no)
			go to 10
		end if
		list(fim(rot(ant))) = inicio(pos)
		fim(rot(ant)) = fim(pos)
		vai = vai + 1
		tree(vai) = k
		solucao = solucao + cost(k)
	end if
	if ( vai .lt. n - 1 ) then
		k = k + 1
		go to 5
	end if
	
	end



	program main

	implicit integer ( a - z )
	real cc1, cc2, clock, tempo


	parameter ( Mnodes = 1010, Marcs = 509545 )

	common/inputnet/tail(Marcs), head(Marcs), cost(Marcs), n, m
	common/result/tree(Mnodes), solucao
	common/work1/rot(Mnodes), list(Marcs), inicio(Mnodes), fim(Mnodes)

C
C   Input data from 'rede.dat' file.
C   Output data to 'ARVORE.res' file (the minimal spanning tree).
C   In the monitor it is displaied some complementary information:
C      Code name, execution time, solution cost, time for constructing data
C      (the execution time includes the time to construct data structure)
C
	open(1, file='rede.dat',status='unknown')
	open(2, file='ARVORE.res',status='unknown')

	read (1, *) n
	read (1, *) m
	if ( n .gt. Mnodes .or. m .gt. Marcs ) then
		print *, 'Max number of Nodes = ', Mnodes
		print *, 'Max number of  Arcs = ', Marcs
		stop
	end if
	do 10 k = 1, m
		read(1, *) tail(k), head(k), cost(k)
10		continue
	cc1 = clock( )
	call quick
	cc2 = clock( ) - cc1
	call Kruskal
	tempo = clock( ) - cc1
	print '(a24,f9.3,a17,i10,5x,f9.3, a)', ' kruskal_quick   time = ',
     &	tempo, ' sec, solution = ', solucao, cc2, ' sec '
	write(2, *) 'ARVORE'
	do 30 i = 1, n - 1
		k = tree(i)
30		write(2, *) k, tail(k)  , head(k), cost(k)

	end
