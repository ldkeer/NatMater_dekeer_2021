
!	Calculation of reaction probabilities P(mu) (not normalised)
 
	do mu=1,size(phi)
	P(mu)=hh(mu)*kmic(mu)
	enddo

!     Calculation of total reaction rate

	TotRxnRate=0.0D0
  	do mu=1,size(phi)
	TotRxnRate=
     &TotRxnRate
     &+P(mu)
      enddo

	if (TotRxnRate.EQ.0) then
	write(*,*) 'error, check input'
      pause
	go to 100
      endif

!     Update time
      Time=Time+1.0d0/TotRxnRate          
!      Time=Time-log(PRN2)/TotRxnRate
      
!	Explicit construction of cumulative probability function phi(µ) 

	phi(1)=P(1)
      do mu=1,size(phi)-1
	phi(mu+1)=
     &P(mu+1)
     &+phi(mu)
      enddo
      	
!     Selection of reaction channel µ

	call RANDOM_NUMBER(random)
      PRN=random*TotRxnRate
      
	muselected=0     
      left=0
	right=size(phi)

	do while (muselected.eq.0)
	dum1=ceiling(0.5D0*(left+right))
	if (PRN.GT.phi(dum1)) then
	left=dum1
	else
	right=dum1
	endif
	if(left.eq.(right-1)) then
	muselected=right
	endif
      enddo
      
