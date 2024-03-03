      subroutine vumat(
c     Read only
     *nblock, ndir,  nshr,  nstatev,  nfieldv,  nprops,  lanneal,
     *stepTime,  totalTime,  dt,  cmname,  coordMp,  charLength,
     *props,  density,  strainInc,  relSpinInc,
     *tempOld,  stretchOld,  defgradOld,  fieldOld,
     *stressOld,  stateOld,  enerInternOld,  enerInelasOld,
     *tempNew,  stretchNew,  defgradNew,  fieldNew,
c     Write only
     1stressNew,  stateNew,  enerInternNew,  enerInelasNew)
      include  'vaba_param.inc'
c     Define variable data type
      real(8)::y,Cb,b,f,Cnc,Cnw,yc,db,CMt,dwz1,dwz2,dwz3,dwz,dcz1,
     *dcz2,dcz3,dcz,CK,Pt,d,deqps0,deqpsn,fdeqps0,fdeqps,ddeqps0,deqps,
     *e,xnu,S1,da,Ckc,Ckw,Cm,Ca,G,f0,f00,Ck0,Ck00,dc0,dw0,ds,Pt0,twomu,
     *alamda,vk,thremu,s11,s12,s13,s21,s22,s23,s31,s32,s33,sm,vmises,
     *peeqOld,peeq_rate,sig_y,dc,dw,sigdif,wc,x0,y0,dy0
      integer::max_iter,j
      
c     All arrays dimensioned by (*) are not used  this algorithm
      dimension coordMP(nblock,*), charLength(nblock), props(nprops),
     *density(nblock),strainInc(nblock,ndir+nshr),
     *relSpinInc(nblock,nshr),  tempOld(nblock),
     *stretchOld(nblock,ndir+nshr),
     *defgradOld(nblock,ndir+nshr+nshr),
     *fieldOld(nblock,nfieldv),  stressOld(nblock,ndir+nshr),
     *stateOld(nblock,nstatev),  enerInternOld(nblock),
     *enerInelasOld(nblock),  tempNew(nblock),
     *stretchNew(nblock,ndir+nshr),
     *defgradNew(nblock,ndir+nshr+nshr),
     *fieldNew(nblock,nfieldv),
     *stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     *enerInternNew(nblock), enerInelasNew(nblock)
      parameter ( zero = 0.d0, one = 1.d0, two =2.d0,
     *third = 1.d0 / 3.d0, half =0.5d0, vp5 = 1.5 )
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
      if (stepTime .eq. zero ) then
          do k = 1, nblock              
c             Trial stress calculation
              trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
              stressNew(k,1) =  stressOld(k,1)
     *                       + twomu * strainInc(k,1) + alamda * trace
              stressNew(k,2) = stressOld(k,2)
     *                       + twomu * strainInc(k,2) + alamda * trace
              stressNew(k,3) = stressOld(k,3)
     *                       + twomu * strainInc(k,3) + alamda * trace
              stressNew(k,4) = stressOld(k,4) + twomu * strainInc(k,4)
              stressNew(k,5) = stressOld(k,5) + twomu * strainInc(k,5)
              stressNew(k,6) = stressOld(k,6) + twomu * strainInc(k,6)
          end do
      else 
          do k = 1, nblock
              trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
              s11 = stressOld(k,1)
     *            + twomu * strainInc(k,1) + alamda * trace
              s22 = stressOld(k,2)
     *            + twomu * strainInc(k,2) + alamda * trace
              s33 = stressOld(k,3)
     *            + twomu * strainInc(k,3) + alamda * trace
              s12 = stressOld(k,4) + twomu * strainInc(k,4)
              s23 = stressOld(k,5) + twomu * strainInc(k,5)
              s13 = stressOld(k,6) + twomu * strainInc(k,6)      
              sm  = (s11+s22+s33)/3     
              s11 = s11-sm
              s22 = s22-sm
              s33 = s33-sm
c             Calculate von Mises stress
              vmises = sqrt( vp5 * ( s11 * s11 + s22 * s22 + s33 * s33 
     *               + two *s12 * s12 + two * s13 * s13 + two * s23 
     *               * s23 ))
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
                  sig_y = S1
                     Pt = Pt0
                     CK = CK0
                     dc = dc0
                     dw = dw0
                      f = f0
                      d = ds     
              endif  
c             Yield condition determination
              sigdif = vmises - sig_y
c             Plastic deformation occurs
              if (sigdif.GT.0) then    
c                 Calculate plastic strain increment(deqps)
                  deqps0 = 1.d-6
                  x0 = 1.d-50
c                 Tolerance
                  wc = 5.d-3
c                 Maximum number of iterations
                  max_iter = 1000                             
                  y0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *              *((CMt*x0/(dt*yc))**(1.d0/Cm))+thremu*x0-vmises   
                  dy0=CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *               *((CMt*x0/(dt*yc))**((1.d0/Cm)-1))*(CMt/(yc*dt))
     *               +thremu

                  do j = 1, max_iter
                      if (deqps0.GE.x0) then
c                         Original function
                          fdeqps0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)
     *                            * sqrt(dc))*((CMt*deqps0/(dt*yc))
     *                            **(1.d0/Cm))+thremu*deqps0-vmises
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
                      endif
                  enddo 

                  if (deqps.le.0) then
                      deqps=0
                      factor=1
                      peeq_rate=0
                      stateOld(k,1)=peeqOld
                      stateOld(k,2)=peeq_rate
                      stateOld(k,3)=sig_y
                      stateOld(k,4)=f
                      stateOld(k,5)=Pt
                      stateOld(k,6)=d
                      stateOld(k,7)=dc
                      stateOld(k,8)=dw

                      stateNew(k,1)=peeqOld
                      stateNew(k,2)=peeq_rate
                      stateNew(k,3)=sig_y
                      stateNew(k,4)=f
                      stateNew(k,5)=Pt
                      stateNew(k,6)=d
                      stateNew(k,7)=dc
                      stateNew(k,8)=dw
                      stateNew(k,9)=CK
                      stateNew(k,10)=j
                      stateNew(k,11)=deqps
                      go to 2
                  endif
c                 Calculate dislocation wall volume fraction:f
                  f=f00+(f0-f00)*exp((-CMt*(peeqOld+deqps))/y)
c                 Calculate varible:CK
                  CK=Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld+deqps))
c                 Calculate the increment of dislocation density in the interior of the cell
                  dcz1=(da*CMt*deqps*sqrt(dw))/(b*sqrt(3.d0))
                  dcz2=(db*CMt*deqps*6.d0)/(b*(1-f)**(third)*d)
                  dcz3=dc*CMt*deqps*Ckc*((CMt*deqps/(dt*yc))
     *                **(1.d0/(-Cnc)))
                  dcz=dcz1-dcz2-dcz3
c                 Calculate the increment of dislocation density in the cell walls
                  dwz1=(CMt*deqps*6.d0*db*((1.d0-f)**(2.d0/3.d0)))
     *                /(d*b*f)
                  dwz2=(sqrt(3.d0)*sqrt(dw)*db*(1.d0-f)*CMt*deqps)/(b*f)
                  dwz3=CMt*deqps*dw*Ckw*((CMt*deqps/(dt*yc))**(1.d0
     *                /(-Cnw)))
                  dwz=dwz1+dwz2-dwz3
c                 Update the dislocation density in the interior of the cell
                  dc=dc+dcz
c                 Update the dislocation density in the cell walls
                  dw=dw+dwz
c                 Update total dislocation density
                  Pt=f*(dw)+(1-f)*(dc)
c                 Update dislocation cell size
                  d=CK/sqrt(Pt)
c                 Update yield stress
                  sig_y = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*
     *                  ((CMt*deqps/(dt*yc))**(1.d0/Cm))
c                 Update plastic strain and plastic strain rate
                  peeqOld=peeqOld+deqps
                  peeq_rate=deqps/dt 
c                 Calculate stress compensation factor
                  factor=(vmises-thremu*deqps)/vmises
c                 Assigning the end state of the current increment step(StateOld)
                  stateOld(k,1)=peeqOld
                  stateOld(k,2)=peeq_rate
                  stateOld(k,3)=sig_y
                  stateOld(k,4)=f
                  stateOld(k,5)=Pt
                  stateOld(k,6)=d
                  stateOld(k,7)=dc
                  stateOld(k,8)=dw
c                 Update variable data (StateNew)
                  stateNew(k,1)=peeqOld
                  stateNew(k,2)=peeq_rate
                  stateNew(k,3)=sig_y
                  stateNew(k,4)=f
                  stateNew(k,5)=Pt
                  stateNew(k,6)=d
                  stateNew(k,7)=dc
                  stateNew(k,8)=dw
                  stateNew(k,9)=CK
                  stateNew(k,10)=j
                  stateNew(k,11)=deqps
c             If no plastic flow occurs
              else 
                  factor=1
                  deqps=0
                  peeq_rate=0
                  j=0

                  stateOld(k,1)=peeqOld
                  stateOld(k,2)=peeq_rate
                  stateOld(k,3)=sig_y
                  stateOld(k,4)=f
                  stateOld(k,5)=Pt
                  stateOld(k,6)=d
                  stateOld(k,7)=dc
                  stateOld(k,8)=dw

                  stateNew(k,1)=peeqOld
                  stateNew(k,2)=peeq_rate
                  stateNew(k,3)=sig_y
                  stateNew(k,4)=f
                  stateNew(k,5)=Pt
                  stateNew(k,6)=d
                  stateNew(k,7)=dc
                  stateNew(k,8)=dw
                  stateNew(k,9)=CK
                  stateNew(k,10)=j
                  stateNew(k,11)=deqps
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
     *                    * strainInc(k,1) + ( stressOld(k,2) 
     *                    + stressNew(k,2) ) * strainInc(k,2) 
     *                    +( stressOld(k,3) + stressNew(k,3) ) 
     *                    * strainInc(k,3) ) + ( stressOld(k,4) 
     *                    + stressNew(k,4) ) * strainInc(k,4) 
     *                    + ( stressOld(k,5) + stressNew(k,5) ) 
     *                    * strainInc(k,5) + (stressold(k,6) 
     *                    + stressnew(k,6) ) * straininc(k,6)
              enerInelasNew(k) =enerInelasOld(k)+stressPower /density(k)
c             Update non-elastic dissipation energy
              plasticWorkInc=sig_y*deqps
              enerInelasNew(k)=enerInelasOld(k)+stressWorkInc/density(k)
          enddo
      endif 
      return
      endsubroutine