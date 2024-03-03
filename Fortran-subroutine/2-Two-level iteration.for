      subroutine vumat(
c Read only -
     1 nblock, ndir,  nshr,  nstatev,  nfieldv,  nprops,  lanneal,
     2 stepTime,  totalTime,  dt,  cmname,  coordMp,  charLength,
     3 props,  density,  strainInc,  relSpinInc,
     4 tempOld,  stretchOld,  defgradOld,  fieldOld,
     5 stressOld,  stateOld,  enerInternOld,  enerInelasOld,
     6 tempNew,  stretchNew,  defgradNew,  fieldNew,
c Write only - 
     1 stressNew,  stateNew,  enerInternNew,  enerInelasNew)

      include  'vaba_param.inc'
c     Define variable data type       
      real(8)::y,Cb,b,f,Cnc,Cnw,yc,db,CMt,dwz1,dwz2,dwz3,dwz,dcz1,
     *dcz2,dcz3,dcz,CK,Pt,d,deqps0,deqpsn,fdeqps0,fdeqps,ddeqps0,deqps,
     *e,xnu,S1,da,Ckc,Ckw,Cm,Ca,G,f0,f00,Ck0,Ck00,dc0,dw0,ds,Pt0,twomu,
     *alamda,vk,thremu,s11,s12,s13,s21,s22,s23,s31,s32,s33,sm,vmises,
     *peeqOld,peeq_rate,sig_y,dc,dw,sigdif,wc,x0,y0,dy0,deqps1,deqps2,
     *sigama
      integer::max_iter,j,n,i
c     All arrays dimensioned by (*) are not used  this algorithm
      dimension coordMP(nblock,*), charLength(nblock), props(nprops),
     *   density(nblock),strainInc(nblock,ndir+nshr),
     *   relSpinInc(nblock,nshr),  tempOld(nblock),
     *   stretchOld(nblock,ndir+nshr),
     *   defgradOld(nblock,ndir+nshr+nshr),
     *   fieldOld(nblock,nfieldv),  stressOld(nblock,ndir+nshr),
     *   stateOld(nblock,nstatev),  enerInternOld(nblock),
     *   enerInelasOld(nblock),  tempNew(nblock),
     *   stretchNew(nblock,ndir+nshr),
     *   defgradNew(nblock,ndir+nshr+nshr),
     *   fieldNew(nblock,nfieldv),
     *   stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     *   enerInternNew(nblock), enerInelasNew(nblock)
       parameter ( zero = 0.d0, one = 1.d0, two =2.d0,
     *             third = 1.d0 / 3.d0, half =0.5d0, vp5 = 1.5 )
c     Dislocation density constitutive variable
      e=props(1)      ! Young's Modulus
      xnu=props(2)    ! Poisson's Ratio
      S1=props(3)     ! Strain-independent Initial Stress                 
      da=props(4)     ! Dislocation Generation   
      db=props(5)     ! Dislocation Motion       
      Ckc=props(6)    ! Parameters of Dislocation Annihilation:kc          
      Ckw=props(7)    ! Parameters of Dislocation Annihilation:kw            
      Cnc=props(8)    ! Parameters of nc    
      Cnw=props(9)    ! Parameters of nw 
      Cm=props(10)    ! Strain Rate Sensitive Parameters:m  
      b=props(11)     ! Burgers Vector      
      CMt=props(12)   ! Taylor Factor:M        
      Ca=props(13)    ! Constant                
      G=props(14)     ! Shear Modulus               
      f0=props(15)    ! Initial Dislocation Cell Wall Volume Fraction 
      f00=props(16)   ! Saturated Dislocation Cell Wall Volume Fraction  
      yc=props(17)    ! Reference Decomposition Shear Strain Rate         
      Ck0=props(18)   ! Parameters related to plastic strain:K0
      Ck00=props(19)  ! Parameters related to plastic strain:K00                
      y=props(20)     ! Coefficient Characterizing Decay Rate              
      Cb=props(21)    ! Constant            
      dc0=props(22)   ! Initial dislocation density in the interior of the cell                
      dw0=props(23)   ! Initial dislocation density in the  cell walls            
      ds=props(24)    ! Initial Dislocation Cell Size   
      Pt0=props(25)   ! Initial Total Dislocation Density
      
      twomu = e / (one + xnu)
      alamda = xnu * twomu /( one-two * xnu)
      thremu = vp5 * twomu
      vk=e/(1-2*xnu) 
c     Initial time step determination
      if (stepTime .eq. zero) then         
          do k = 1, nblock              
c             Trial stress calculation
              trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
              stressNew(k,1) =  stressOld(k,1)
     *                       + twomu * strainInc(k,1) + alamda * trace
              stressNew(k,2) = stressOld(k,2)
     *                       + twomu * strainInc(k,2) + alamda * trace
              stressNew(k,3) = stressOld(k,3)
     *                       + twomu * strainInc(k,3) + alamda * trace
              stressNew(k,4)=stressOld(k,4) + twomu * strainInc(k,4)
              stressNew(k,5)=stressOld(k,5) + twomu * strainInc(k,5)
              stressNew(k,6)=stressOld(k,6) + twomu * strainInc(k,6)    
          enddo     
      else
          do k = 1, nblock
              trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
              s11 = stressOld(k,1)
     *            + twomu * strainInc(k,1) + alamda * trace
              s22 = stressOld(k,2)
     *            + twomu * strainInc(k,2) + alamda * trace
              s33 = stressOld(k,3)
     *            + twomu * strainInc(k,3) + alamda * trace
              s12=stressOld(k,4) + twomu * strainInc(k,4)
              s23= stressOld(k,5) + twomu * strainInc(k,5)
              s13 = stressOld(k,6) + twomu * strainInc(k,6)
              sm = (s11+s22+s33)/3
              s11=s11-sm
              s22=s22-sm
              s33=s33-sm
c             Calculate von Mises stress              
              vmises = sqrt( vp5 * ( s11 * s11 + s22 * s22 + s33 * s33 
     *               + two * s12 * s12 + two * s13 * s13 + two * s23 
     *               * s23 ) )
c             Parameter variable at current time step 
              peeqOld   = stateOld(k,1)
              peeq_rate = stateOld(k,2) 
              sig_y     = stateOld(k,3)
              f         = stateOld(k,4)
              Pt        = stateOld(k,5)
              d         = stateOld(k,6)
              dc        = stateOld(k,7)
              dw        = stateOld(k,8)
              
              if (peeqOld.eq.zero) then
                  sig_y=S1
                  Pt=Pt0
                  d=ds
                  dc=dc0
                  dw=dw0
                  f=2.5d-1
                  CK=100
                  peeq_rate=zero
              endif  
c             Yield condition determination
              sigdif = vmises - sig_y
c             Plastic deformation occurs
              if (sigdif.GT.0) then    
c                 Calculate plastic strain increment-1
                  deqps0=1.d-5
                  x0=1.d-50
c                 Tolerance
                  wc = 5.d-3
c                 Maximum number of iterations
                  max_iter = 1000  
c                 Original function       
                  y0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *               *((CMt*x0/(dt*yc))**(1.d0/Cm))+thremu*x0-vmises
c                 Derivative of the original function       
                  dy0 = CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *                * ((CMt*x0/(dt*yc))**((1.d0/Cm)-1))*(CMt/(yc*dt))
     *                + thremu

                  do j=1, max_iter
                      if (deqps0.GE.x0) then
c                         Original function
                          fdeqps0=S1+ CMt*G*b*Ca*(f*sqrt(dw)+(1-f)
     *                              * sqrt(dc))*((CMt*deqps0/(dt*yc))
     *                              ** (1.d0/Cm))+thremu*deqps0-vmises
c                         Derivative of the original function                      
                          ddeqps0=CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *                           *(1.d0/Cm)*((CMt*deqps0/(dt*yc))
     *                           **((1.d0/Cm)-1))*(CMt/(yc*dt))+thremu
                      else
                          fdeqps0=y0+dy0*(deqps0-x0)
                          ddeqps0=dy0
                      endif
                      deqps=deqps0-fdeqps0/ddeqps0
                      if (abs((deqps-deqps0)/deqps).gt.wc) then
                          deqps0=deqps
                      else 
                          exit
                      end if
                  enddo 
                  if (deqps.lt.0) then
                      deqps1=0
                      i=0
                      n=0
                      factor=1
                      deqps2=0
                      peeq_rate=0
                      stateOld(k,1)=peeqOld
                      stateOld(k,2)=peeq_rate
                      stateOld(k,3)=sig_y
                      stateOld(k,4)=f
                      stateOld(k,5)=Pt
                      stateOld(k,6)=d
                      stateOld(k,7)=dc
                      stateOld(k,8)=dw
                      stateOld(k,9)=CK
                      stateOld(k,10)=j
                      stateOld(k,11)=i
                      stateOld(k,12)=deqps1
                      stateOld(k,13)=deqps2 
                      stateOld(k,14)=n
           
                      stateNew(k,1)= stateOld(k,1)
                      stateNew(k,2)= stateOld(k,2)
                      stateNew(k,3)= stateOld(k,3)
                      stateNew(k,4)= stateOld(k,4)
                      stateNew(k,5)= stateOld(k,5)
                      stateNew(k,6)= stateOld(k,6)
                      stateNew(k,7)= stateOld(k,7)
                      stateNew(k,8)= stateOld(k,8)
                      stateNew(k,9)= stateOld(k,9)
                      stateNew(k,10)= stateOld(k,10)
                      stateNew(k,11)= stateOld(k,11)
                      stateNew(k,12)= stateOld(k,12)
                      stateNew(k,13)= stateOld(k,13)
                      stateNew(k,14)= stateOld(k,14)
                          go to 2
                  endif
c                 Plastic strain increment:deqps1
                  deqps1=deqps     
c                 Total number of iterations before meeting the tolerance condition:n
                  n=0
1                 continue
c                 Calculate dislocation wall volume fraction:f
                  f=f00+(f0-f00)*exp((-CMt*(peeqOld+deqps1))/y)
c                 Calculate varible:CK
                  CK=Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld+deqps1))
c                 Calculate the increment of dislocation density in the interior of the cell
                  dcz1=(da*CMt*deqps1*sqrt(dw))/(b*sqrt(3.d0))
                  dcz2=(db*CMt*deqps1*6.d0)/(b*(1-f)**(third)*d)
                  dcz3=dc*CMt*deqps1*Ckc
     *                *((CMt*deqps1/(dt*yc))**(1.d0/(-Cnc)))
                  dcz=dcz1-dcz2-dcz3
c                 Calculate the increment of dislocation density in the cell walls
                  dwz1=(CMt*deqps1*6.d0*db*((1.d0-f)**(2.d0/3.d0)))
     *                /(d*b*f)
                  dwz2=(sqrt(3.d0)*sqrt(dw)*db*(1.d0-f)*CMt*deqps1)
     *                /(b*f)
                  dwz3=CMt*deqps1*dw*Ckw
     *                *((CMt*deqps1/(dt*yc))**(1.d0/(-Cnw)))
                  dwz=dwz1+dwz2-dwz3
c                 Update the dislocation density in the interior of the cell
                  dc=dc+dcz
c                 Update the dislocation density in the cell walls
                  dw=dw+dwz
c                 Update total dislocation density
                  Pt=f*(dw)+(1-f)*(dc)
c                 Update dislocation cell size
                  d=CK/sqrt(Pt)  
c                 Calculate plastic strain increment-2
                  deqps0=deqps1
                  x0=1.d-50
                  y0=S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *              *((CMt*x0/(dt*yc))**(1.d0/Cm))+thremu*x0-vmises
      
                  dy0=CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *               *((CMt*x0/(dt*yc))**((1.d0/Cm)-1))
     *               *(CMt/(yc*dt))+thremu

                  do i=1, max_iter
                      if (deqps0.GT.x0) then
c                         Original function                          
                          fdeqps0=S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)
     *                           *sqrt(dc))*((CMt*deqps0/(dt*yc))
     *                           **(1.d0/Cm))+thremu*deqps0-vmises
c                         Derivative of the original function        
                          ddeqps0=CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *                           *(1.d0/Cm)*((CMt*deqps0/(dt*yc))
     *                           **((1.d0/Cm)-1))*(CMt/(yc*dt))+thremu
                      else
                          fdeqps0=y0+dy0*(deqps0-x0)
                          ddeqps0=dy0
                      endif
                      deqpsn=deqps0-fdeqps0/ddeqps0
                      if (abs((deqpsn-deqps0)/deqpsn).gt.wc) then
                          deqps0=deqpsn
                      else
                              exit
                      endif
                  enddo 
c                 Plastic strain increment:deqps2
                  deqps2=deqpsn             
c                 Tolerance of deqps1 and deqps2
                  sigama=5.d-3
                  if (abs((deqps1-deqps2)/deqps1).GE.sigama) then
                      n=n+1
                      deqps1=deqps2
                      go to 1
                  else
c                     Updating Equivalent Plastic Strain and Equivalent Plastic Strain Rate
                      peeqOld=peeqOld+deqps2
                      peeq_rate=deqps2/dt 
                  endif
c                 Update yield stress
                  sig_y = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *                  * ((CMt*deqps2/(dt*yc))**(1.d0/Cm))
c                 Calculate stress compensation factor
                  factor=(vmises-thremu*deqps2)/vmises 
c                 Assigning the end state of the current increment step(StateOld)
                  stateOld(k,1)=peeqOld
                  stateOld(k,2)=peeq_rate
                  stateOld(k,3)=sig_y
                  stateOld(k,4)=f
                  stateOld(k,5)=Pt
                  stateOld(k,6)=d
                  stateOld(k,7)=dc
                  stateOld(k,8)=dw
                  stateOld(k,9)=CK
                  stateOld(k,10)=j
                  stateOld(k,11)=i
                  stateOld(k,12)=deqps1
                  stateOld(k,13)=deqps2 
                  stateOld(k,14)=n
c                 Update variable data (StateNew)
                  stateNew(k,1)= stateOld(k,1)
                  stateNew(k,2)= stateOld(k,2)
                  stateNew(k,3)= stateOld(k,3)
                  stateNew(k,4)= stateOld(k,4)
                  stateNew(k,5)= stateOld(k,5)
                  stateNew(k,6)= stateOld(k,6)
                  stateNew(k,7)= stateOld(k,7)
                  stateNew(k,8)= stateOld(k,8)
                  stateNew(k,9)= stateOld(k,9)
                  stateNew(k,10)= stateOld(k,10)
                  stateNew(k,11)= stateOld(k,11)
                  stateNew(k,12)= stateOld(k,12)
                  stateNew(k,13)= stateOld(k,13)
                  stateNew(k,14)= stateOld(k,14)
c             If no plastic flow occurs  
              else
                  factor=1
                  deqps1=0
                  deqps2=0
                  peeq_rate=0
                  i=0
                  j=0
                  n=0

                  stateOld(k,1)=peeqOld
                  stateOld(k,2)=peeq_rate
                  stateOld(k,3)=sig_y
                  stateOld(k,4)=f
                  stateOld(k,5)=Pt
                  stateOld(k,6)=d
                  stateOld(k,7)=dc
                  stateOld(k,8)=dw
                  stateOld(k,9)=CK
                  stateOld(k,10)=j
                  stateOld(k,11)=i
                  stateOld(k,12)=deqps1
                  stateOld(k,13)=deqps2 
                  stateOld(k,14)=n

                  stateNew(k,1)= stateOld(k,1)
                  stateNew(k,2)= stateOld(k,2)
                  stateNew(k,3)= stateOld(k,3)
                  stateNew(k,4)= stateOld(k,4)
                  stateNew(k,5)= stateOld(k,5)
                  stateNew(k,6)= stateOld(k,6)
                  stateNew(k,7)= stateOld(k,7)
                  stateNew(k,8)= stateOld(k,8)
                  stateNew(k,9)= stateOld(k,9)
                  stateNew(k,10)= stateOld(k,10)
                  stateNew(k,11)= stateOld(k,11)
                  stateNew(k,12)= stateOld(k,12)
                  stateNew(k,13)= stateOld(k,13)
                  stateNew(k,14)= stateOld(k,14)
              endif
c             Update stress components
2             stressNew(k,1) = s11*factor+ sm  
              stressNew(k,2) = s22*factor+ sm  
              stressNew(k,3) = s33*factor+ sm  
              stressNew(k,4) = s12*factor  
              stressNew(k,5) = s23*factor  
              stressNew(k,6) = s13*factor
c             Update specific internal energy
              stressPower = half * (( stressOld(k,1) + stressNew(k,1) ) 
     *                    * strainInc(k,1)+( stressOld(k,2) 
     *                    + stressNew(k,2) ) * strainInc(k,2) 
     *                    +( stressOld(k,3) + stressNew(k,3) ) 
     *                    * strainInc(k,3) ) +( stressOld(k,4) 
     *                    + stressNew(k,4) ) * strainInc(k,4) 
     *                    +( stressOld(k,5) + stressNew(k,5) ) 
     *                    * strainInc(k,5) +(stressold(k,6) 
     *                    + stressnew(k,6) ) * straininc(k,6)
              enerInelasNew(k) = enerInelasOld(k)+stressPower/density(k)
c             Update non-elastic dissipation energy
              plasticWorkInc=sig_y*deqps2
              enerInelasNew(k)=enerInelasOld(k)+stressWorkInc/density(k)
          enddo   
      endif    
      return
      end subroutine
      
      
     
            
