      subroutine vumat(
c     Read only -
     *nblock, ndir,  nshr,  nstatev,  nfieldv,  nprops,  lanneal,
     *stepTime,  totalTime,  dt,  cmname,  coordMp,  charLength,
     *props,  density,  strainInc,  relSpinInc,
     *tempOld,  stretchOld,  defgradOld,  fieldOld,
     *stressOld,  stateOld,  enerInternOld,  enerInelasOld,
     *tempNew,  stretchNew,  defgradNew,  fieldNew,
c      Write only - 
     *stressNew,  stateNew,  enerInternNew,  enerInelasNew)
      include  'vaba_param.inc'
c     Define variable data type
      real(8):: y,Cb,b,f,Cnc,Cnw,yc,db,CMt,dwz1,dwz2,dwz3,dwz,dcz1,
     *dcz2,dcz3,dcz,CK,Pt,d,deqps0,deqpsn,fdeqps0,fdeqps,ddeqps0,deqps,
     *e,xnu,S1,da,Ckc,Ckw,Cm,Ca,G,f0,f00,Ck0,Ck00,dc0,dw0,ds,Pt0,twomu,
     *alamda,vk,thremu,s11,s12,s13,s21,s22,s23,s31,s32,s33,sm,vmises,
     *peeqOld,peeq_rate,sig_y,dc,dw,sigdif,wc,x0,y0,dy0,deqps1,deqps2,
     *dwz0,dcz0,f1,f2,f3,J11,J12,J13,J21,J22,J23,J31,J32,J33,invJ11,
     *invJ12,invJ13,invJ21,invJ22,invJ23,invJ31,invJ32,invJ33,detJ,
     *idetJ,wcmax3
      
      integer::L,max_iter
c     Definition matrix
      real(8):: J(3,3),invJ(3,3),cz(1,3),czz(3,1),hsjz(1,3),ff(3,1),
     *mbz(3,1),wcmax(1,3),wcmax2(3,1)
c All arrays dimensioned by (*) are not used  this algorithm
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
      if ( stepTime .eq. zero ) then    
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
          enddo
      else 
          do k = 1, nblock
              trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
              s11   = stressOld(k,1)
     *              + twomu * strainInc(k,1) + alamda * trace
              s22   = stressOld(k,2)
     *              + twomu * strainInc(k,2) + alamda * trace
              s33   = stressOld(k,3)
     *              + twomu * strainInc(k,3) + alamda * trace
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
              peeqOld=stateOld(k,1)
              peeq_rate=stateOld(k,2) 
              sig_y= stateOld(k,3)
              f =stateOld(k,4)
              Pt=stateOld(k,5)
              d =stateOld(k,6)
              dc=stateOld(k,7)
              dw=stateOld(k,8)
              L=stateOld(k,9)
              deqps=stateOld(k,10)
              
              if (peeqOld.eq.zero) then
                  sig_y=S1
                  Pt=Pt0
                  d=ds
                  dc=dc0
                  dw=dw0
                  f=2.5d-1
                  CK=1.d2       
              endif  
c             Yield condition determination
              sigdif = vmises - sig_y
c             Plastic deformation occurs
              if (sigdif.Gt.0) then    
                  deqps0=1.d-6
                  dwz0=1.d6
                  dcz0=1.d6
                  x0=1.d-50
c                 Tolerance
                  wc = 5.d-3
c                 Maximum number of iterations
                  max_iter=1000
c                 Solve nonlinear equations:[f1 f2 f3]
                  do L=1, max_iter  
c                     Initial value matrix
                      cz=reshape([deqps0,dwz0,dcz0],shape(cz))
                      czz=transpose(cz)
c                     Nonlinear equations:f1
                      f1 = S1+CMt*G*b*ca*((f00+(f0-f00)*exp((-CMt
     *                   * (peeqOld+deqps0))/y))*sqrt(dw+dwz0)+(1-(f00
     *                   + (f0-f00)*exp((-CMt*(peeqOld+deqps0))/y)))
     *                   * sqrt(dc+dcz0))*((CMt*deqps0/(dt*yc))
     *                   ** (1.d0/Cm))+thremu*deqps0-vmises
c                     Nonlinear equations:f2
                      f2 = (CMt*deqps0*6*db*((1-(f00+(f0-f00)
     *                   * exp((-CMt*(peeqOld+deqps0))/y)))**(2.d0/3)))
     *                   / (((Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld
     *                   + deqps0)))/sqrt((f00+(f0-f00)*exp((-CMt
     *                   * (peeqOld+deqps0))/y))*(dw+dwz0)+(1-(f00
     *                   + (f0-f00)*exp((-CMt*(peeqOld
     *                   + deqps0))/y)))*(dc+dcz0)))*b*(f00+(f0-f00)
     *                   * exp((-CMt*(peeqOld+deqps0))/y)))+(sqrt(3.d0)
     *                   * sqrt(dw+dwz0)*db*(1-(f00+(f0-f00)*exp((-CMt
     *                   * (peeqOld+deqps0))/y)))*CMt*deqps0)/(b*(f00
     *                   + (f0-f00)*exp((-CMt*(peeqOld+deqps0))/y)))
     *                   - CMt*deqps0*(dw+dwz0)*Ckw*((CMt*deqps0
     *                   / (dt*yc))**(1.d0/(-Cnw)))-dwz0
c                     Nonlinear equations:f3
                      f3 = (da*CMt*deqps0*sqrt(dw+dwz0))/(b*sqrt(3.d0))
     *                   - (db*CMt*deqps0*6)/(b*(1-(f00+(f0-f00)
     *                   * exp((-CMt*(peeqOld+deqps0))/y)))**(1.d0/3)
     *                   * ((Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld
     *                   + deqps0)))/sqrt((f00+(f0-f00)*exp((-CMt
     *                   * (peeqOld+deqps0))/y))*(dw+dwz0)+(1-(f00
     *                   + (f0-f00)*exp((-CMt*(peeqOld+deqps0))/y)))*(dc
     *                   + dcz0))))-(dcz0+dc)*CMt*deqps0*Ckc*((CMt
     *                   * deqps0/(dt*yc))**(1.d0/(-Cnc)))-dcz0
c                     Nonlinear equations:ff
                      hsjz=reshape([f1,f2,f3],shape(hsjz))
                      ff=transpose(hsjz)
c                     Jacobian matrix
                      J11 = thremu+CMt*Ca*G*b*((CMt*exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(dc+dcz0)**(1.d0/2)*(f0-f00))/y
     *                    - (CMt*exp(-(CMt*(deqps0+peeqOld))/y)
     *                    * (dw+dwz0)**(1.d0/2)*(f0-f00))/y)*((CMt
     *                    * deqps0)/(dt*yc))**(1.d0/Cm)-(CMt**(2)*Ca*G*b
     *                    * ((dc+dcz0)**(1.d0/2)*(f00+exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0-f00)-1)-(dw+dwz0)**(1.d0/2)
     *                    * (f00+exp(-(CMt*(deqps0+peeqOld))/y)
     *                    * (f0-f00)))*((CMt*deqps0)/(dt*yc))
     *                    ** ((1.d0/Cm)-1))/(Cm*dt*yc)

                      J12 = (CMt*Ca*G*b*(f00+exp(-(CMt*(deqps0+peeqOld))
     *                   / y)*(f0-f00))*((CMt*deqps0)/(dt*yc))
     *                   ** (1.d0/Cm))/(2*(dw+dwz0)**(1.d0/2))

                      J13 = -(CMt*Ca*G*b*(f00+exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0-f00)-1)*((CMt*deqps0)
     *                    / (dt*yc))**(1.d0/Cm))/(2*sqrt(dc+dcz0))

                      J21 = (6*CMt*db*((dw+dwz0)*(f00+exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0 - f00)) - (dc + dcz0)*(f00 
     *                    + exp(-(CMt*(deqps0+ peeqOld))/y)*(f0 - f00)
     *                    - 1))**(1.d0/2)*(1-exp(-(CMt*(deqps0 
     *                    + peeqOld))/y)*(f0 - f00)-f00)**(2.d0/3))/(b
     *                    * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0 
     *                    - f00))*(Ck00+ exp(-CMt*Cb*(deqps0 + peeqOld))
     *                    * (Ck0-Ck00)))- (3**(1.d0/2)*CMt*db*(dw +dwz0)
     *                    ** (1.d0/2)*(f00+ exp(-(CMt*(deqps0+peeqOld))
     *                    / y)*(f0 - f00) - 1))/ (b*(f00+exp(-(CMt
     *                    * (deqps0 + peeqOld))/y)*(f0 - f00)))
     *                    - (CMt*Ckw*(dw + dwz0))/((CMt*deqps0)/(dt*yc))
     *                    ** (1.d0/Cnw)+ (CMt**(2)*Ckw*deqps0*(dw+dwz0))
     *                    / (Cnw*dt*yc*((CMt*deqps0)/ (dt*yc))**(1.d0
     *                    / Cnw + 1))+(3*CMt*db*deqps0*((CMt*exp(
     *                    - (CMt*(deqps0+peeqOld))/y)*(dc+dcz0)
     *                    * (f0 -f00))/y - (CMt*exp(-(CMt*(deqps0 
     *                    + peeqOld))/y)*(dw+ dwz0)*(f0 - f00))/y)*(1 
     *                    - exp(-(CMt*(deqps0 + peeqOld))/ y)*(f0-f00)
     *                    - f00)**(2.d0/3))/(b*(f00+exp(-(CMt*(deqps0 
     *                    + peeqOld))/y)*(f0 -f00))*(Ck00 + exp(-CMt*Cb
     *                    * (deqps0 + peeqOld))*(Ck0 - Ck00))*((dw 
     *                    + dwz0)*(f00 + exp(-(CMt* (deqps0 + peeqOld))
     *                    / y)*(f0 - f00))-(dc + dcz0)*(f00 
     *                    + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0-f00)
     *                    - 1))**(1.d0/2))+ (3**(1.d0/2)*CMt**(2)*db
     *                    * deqps0*exp(-(CMt*(deqps0 + peeqOld))/y)*(dw 
     *                    + dwz0)**(1.d0/2)*(f0-f00))
     *                    / (b*y*(f00 + exp(-(CMt*(deqps0 + peeqOld))/y)
     *                    * (f0-f00)))+ (4*CMt**(2)*db*deqps0*exp(-(CMt
     *                    * (deqps0+peeqOld))/y)*(f0- f00)*((dw + dwz0)
     *                    * (f00 + exp(-(CMt*(deqps0+peeqOld))/y)
     *                    * (f0 - f00))-(dc + dcz0)*(f00 + exp(-(CMt
     *                    * (deqps0+ peeqOld))/y)*(f0 - f00) - 1))
     *                    ** (1.d0/2))/(b*y*(f00 + exp(-(CMt*(deqps0 
     *                    + peeqOld))/y)*(f0 - f00))*(Ck00 
     *                    + exp(-CMt*Cb*(deqps0 + peeqOld))*(Ck0-Ck00))
     *                    * (1 - exp(-(CMt*(deqps0 + peeqOld))/y)*(f0
     *                    - f00) -f00)** (1.d0/3))+(6*CMt**(2)*db
     *                    * deqps0*exp(-(CMt*(deqps0+ peeqOld))/y)*(f0
     *                    - f00)*((dw + dwz0)*(f00 + exp(-(CMt
     *                    * (deqps0 + peeqOld))/y)*(f0 - f00)) 
     *                    - (dc + dcz0)*(f00 + exp(-(CMt*(deqps0 
     *                    + peeqOld))/y)* (f0 - f00) -1))**(1.d0/2)*(1
     *                    - exp(-(CMt*(deqps0+peeqOld))
     *                    / y)*(f0 - f00) - f00)**(2.d0/3))/(b*y*(f00  
     *                    + exp(-(CMt*(deqps0+peeqOld))/y)*(f0 - f00))
     *                    ** (2)*(Ck00 + exp(-CMt*Cb*(deqps0+peeqOld))
     *                    * (Ck0 - Ck00)))+ (6*CMt**(2)*Cb*db*deqps0
     *                    * exp(-CMt*Cb*(deqps0+peeqOld))
     *                    * (Ck0 - Ck00)*((dw + dwz0)*(f00 + exp(-(CMt
     *                    * (deqps0+ peeqOld))/y)*(f0 - f00))-(dc+dcz0)
     *                    * (f00 + exp(-(CMt*(deqps0 + peeqOld))/y)
     *                    * (f0 - f00) -1))**(1.d0/2)*(1-exp(-(CMt
     *                    * (deqps0+peeqOld))/ y)*(f0- f00) -f00)
     *                    ** (2.d0/3))/(b*(f00+exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0 - f00))*(Ck00 + exp(-CMt*Cb
     *                    * (deqps0+ peeqOld))*(Ck0-Ck00))**(2)) 
     *                    - (3**(1.d0/2)*CMt**(2)*db
     *                    * deqps0*exp(-(CMt*(deqps0 + peeqOld))/y)*(dw 
     *                    + dwz0)** (1.d0/2)*(f0 - f00)*(f00 + exp(-(CMt
     *                    * (deqps0+ peeqOld))/y)*(f0 - f00) - 1))
     *                    / (b*y*(f00 + exp(-(CMt*(deqps0 + peeqOld))/y)
     *                    * (f0 -f00))**(2))

                      J22 = (3*CMt*db*deqps0*(1-exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0-f00)-f00)**(2.d0/3))/(b
     *                    * (Ck00+exp(-CMt*Cb*(deqps0+peeqOld))*(Ck0
     *                    - Ck00))*((dw+dwz0)*(f00 +exp(-(CMt*(deqps0 
     *                    + peeqOld))/y)*(f0-f00))-(dc + dcz0)*(f00 
     *                    + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 - f00)
     *                    - 1))**(1.d0/2)) - (CMt*Ckw*deqps0)
     *                    / ((CMt*deqps0)/(dt*yc))**(1.d0/Cnw)
     *                    - (3**(1.d0/2)*CMt*db*deqps0*(f00+exp(-(CMt
     *                    * (deqps0 + peeqOld))/y)*(f0 - f00) - 1))/(2*b
     *                    * (dw + dwz0)**(1.d0/2)*(f00+exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0-f00)))-1

                      J23 = (3*CMt*db*deqps0*(1-exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0-f00)-f00)**(5.d0/3))/(b
     *                    * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0
     *                    - f00))*(Ck00 + exp(-CMt*Cb*(deqps0+peeqOld))
     *                    * (Ck0-Ck00))*((dw+dwz0)*(f00+exp(-(CMt
     *                    * (deqps0+peeqOld))/y)*(f0 - f00))-(dc+dcz0)
     *                    * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0
     *                    - f00)-1))**(1.d0/2))

                      J31 = ((3**(1.d0/2)*CMt*da*(dw+dwz0)**(1.d0/2))
     *                    / (3*b)-(CMt*Ckc*(dc+dcz0))/((CMt*deqps0)/(dt
     *                    * yc))**(1.d0/Cnc)-(6*CMt*db*((dw + dwz0)*(f00
     *                    + exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00))-(dc
     *                    + dcz0)*(f00 + exp(-(CMt*(deqps0+peeqOld))/y)
     *                    * (f0-f00)-1))**(1.d0/2))/(b*(Ck00 + exp(-CMt
     *                    * Cb*(deqps0+peeqOld))*(Ck0-Ck00))*(1-exp(
     *                    - (CMt*(deqps0+peeqOld))/y)*(f0-f00)-f00)
     *                    ** (1.d0/3))- (3*CMt*db*deqps0*((CMt*exp(-(CMt
     *                    * (deqps0 + peeqOld))/y)*(dc+dcz0)*(f0 - f00))
     *                    / y - (CMt*exp(-(CMt*(deqps0+peeqOld))/y)*(dw 
     *                    + dwz0)*(f0 - f00))/y))/(b*(Ck00
     *                    + exp(-CMt*Cb*(deqps0 + peeqOld))*(Ck0 -Ck00))
     *                    * ((dw+dwz0)*(f00+exp(-(CMt*(deqps0+peeqOld))
     *                    / y)*(f0-f00))-(dc+dcz0)*(f00+exp(-(CMt
     *                    * (deqps0+peeqOld))/y)*(f0 - f00) - 1))
     *                    ** (1.d0/2)*(1 -exp(-(CMt*(deqps0 + peeqOld))
     *                    / y)*(f0 - f00) - f00)**(1.d0/3)) + (CMt**(2)
     *                    * Ckc*deqps0*(dc + dcz0))/(Cnc*dt*yc*((CMt
     *                    * deqps0)/(dt*yc))**(1.d0/Cnc+1))+(2*CMt**(2)
     *                    * db*deqps0*exp(-(CMt*(deqps0+peeqOld))
     *                    / y)*(f0 - f00)*((dw + dwz0)*(f00 + exp(-(CMt
     *                    * (deqps0 + peeqOld))/y)*(f0-f00))-(dc+dcz0)
     *                    * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00)
     *                    - 1))**(1.d0/2))/(b*y*(Ck00+exp(-CMt*Cb
     *                    * (deqps0+peeqOld))*(Ck0 - Ck00))*(1 - exp(
     *                    - (CMt*(deqps0+peeqOld))/y)*(f0-f00)-f00)
     *                    ** (4.d0/3))-(6*CMt**(2)*Cb*db*deqps0*exp(-CMt
     *                    * Cb*(deqps0 + peeqOld))*(Ck0 - Ck00)*((dw 
     *                    + dwz0)*(f00 + exp(-(CMt*(deqps0
     *                    + peeqOld))/y)*(f0 - f00)) - (dc + dcz0)*(f00
     *                    + exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00)-1))
     *                    ** (1.d0/2))/(b*(Ck00 + exp(-CMt*Cb*(deqps0
     *                    + peeqOld))*(Ck0 - Ck00))**(2)*(1 - exp(-(CMt
     *                    * (deqps0+peeqOld))/y)*(f0-f00)-f00)
     *                    ** (1.d0/3)))

                      J32 = (3**(1.d0/2)*CMt*da*deqps0)/(6*b*(dw + dwz0)
     *                    ** (1.d0/2))-(3*CMt*db*deqps0*(f00+exp(-(CMt
     *                    * (deqps0 + peeqOld))/y)*(f0-f00)))/(b*(Ck00
     *                    + exp(-CMt*Cb*(deqps0+peeqOld))*(Ck0-Ck00))
     *                    * ((dw+dwz0)*(f00+exp(-(CMt*(deqps0+peeqOld))
     *                    / y)*(f0-f00))-(dc+dcz0)*(f00+exp(-(CMt
     *                    * (deqps0+peeqOld))/y)*(f0-f00)-1))**(1.d0/2)
     *                    * (1-exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00)
     *                    - f00)**(1.d0/3))

                      J33 = -(CMt*Ckc*deqps0)/((CMt*deqps0)/(dt*yc))
     *                    ** (1.d0/Cnc)-(3*CMt*db*deqps0*(1-exp(-(CMt
     *                    * (deqps0+peeqOld))/y)*(f0 - f00)-f00)
     *                    ** (2.d0/3))/(b*(Ck00+exp(-CMt*Cb*(deqps0
     *                    + peeqOld))*(Ck0 - Ck00))*((dw+dwz0)*(f00 
     *                    + exp(-(CMt*(deqps0 +peeqOld))/y)*(f0-f00))
     *                    - (dc+dcz0)*(f00+exp(-(CMt*(deqps0+peeqOld))
     *                    / y)*(f0 - f00) - 1))**(1.d0/2)) - 1 

c                     Inverse of the Jacobian matrix
                      detJ = J11*J22*J33+J12*J23*J31+J13*J21*J32
     *                     - J31*J22*J13-J21*J12*J33-J11*J23*J32

                      idetJ  = 1.d0/detJ
                      invJ11 = idetJ*(J22*J33-J23*J32)
                      invJ12 = idetJ*(-1*(J12*J33-J32*J13))
                      invJ13 = idetJ*(J12*J23-J13*J22)
                      invJ21 = idetJ*(-1*(J21*J33-J31*J23))
                      invJ22 = idetJ*(J11*J33-J13*J31)
                      invJ23 = idetJ*(-1*(J11*J23-J21*J13))
                      invJ31 = idetJ*((J21*J32-J31*J22))
                      invJ32 = idetJ*(-1*(J11*J32-J31*J12))
                      invJ33 = idetJ*(J11*J22-J21*J12)
c                     InvJ
                  invJ = reshape([invJ11,invJ12,invJ13,invJ21,invJ22,
     *                invJ23,invJ31,invJ32,invJ33],shape(invJ))
                      invJ = transpose(invJ)
c                     Calculation of the target value matrix:mbz
                      mbz=czz-matmul(invJ,ff)
c                     Value of variables
                      deqps = mbz(1,1)
                      dwz   = mbz(2,1)
                      dcz   = mbz(3,1)
c                     Tolerance of variables
                      deqpswc=abs((deqps-deqps0)/deqps)
                      dwzwc=abs((dwz-dwz0)/dwz)
                      dczwc=abs((dcz-dcz0)/dcz)
c                     Maximum tolerance of variables:wcmax3
                      wcmax=reshape([deqpswc,dwzwc,dczwc],shape(wcmax))
                      wcmax2=transpose(wcmax)
                      wcmax3=maxval(wcmax2) 
                      if (wcmax3.gt.wc) then
                          if (deqps.lt.x0) then
                              deqps=x0
                          endif
                          deqps0=deqps
                          dwz0=dwz
                          dcz0=dcz
                      else
                          exit
                      endif
                  enddo
                  if (deqps.lt.x0) then
                      factor=1
                      deqps=0
                      peeq_rate=0


                      stateOld(k,1)=peeqOld
                      stateOld(k,2)=peeq_rate
                      stateOld(k,3)=sig_y
                      stateOld(k,4)=f
                      stateOld(k,5)=Pt
                      stateOld(k,6)=d
                      stateOld(k,7)=dc
                      stateOld(k,8)=dw
                      stateOld(k,9)=L
                      stateOld(k,10)=deqps


                      stateNew(k,1)=stateOld(k,1)
                      stateNew(k,2)=stateOld(k,2)
                      stateNew(k,3)=stateOld(k,3)
                      stateNew(k,4)=stateOld(k,4)
                      stateNew(k,5)=stateOld(k,5)
                      stateNew(k,6)=stateOld(k,6)
                      stateNew(k,7)=stateOld(k,7)
                      stateNew(k,8)=stateOld(k,8)
                      stateNew(k,9)=stateOld(k,9)
                      stateNew(k,10)=stateOld(k,10)
                      go to 2
                  endif
c                 Calculate dislocation wall volume fraction:f
                  f=f00+(f0-f00)*exp((-CMt*(peeqOld+deqps))/y)
c                 Calculate varible:CK
                  CK=Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld+deqps))
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
     *            ((CMt*deqps/(dt*yc))**(1.d0/Cm))
c                 Updating Equivalent Plastic Strain and Equivalent Plastic Strain Rate
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
                  stateOld(k,9)=L
                  stateOld(k,10)=deqps
c                 Update variable data (StateNew)
                  stateNew(k,1)=stateOld(k,1)
                  stateNew(k,2)=stateOld(k,2)
                  stateNew(k,3)=stateOld(k,3)
                  stateNew(k,4)=stateOld(k,4)
                  stateNew(k,5)=stateOld(k,5)
                  stateNew(k,6)=stateOld(k,6)
                  stateNew(k,7)=stateOld(k,7)
                  stateNew(k,8)=stateOld(k,8)
                  stateNew(k,9)=stateOld(k,9)
                  stateNew(k,10)=stateOld(k,10)
c             If no plastic flow occurs 
              else
                  factor=1
                  deqps=0
                  peeq_rate=0

                  stateOld(k,1)=peeqOld
                  stateOld(k,2)=peeq_rate
                  stateOld(k,3)=sig_y
                  stateOld(k,4)=f
                  stateOld(k,5)=Pt
                  stateOld(k,6)=d
                  stateOld(k,7)=dc
                  stateOld(k,8)=dw
                  stateOld(k,9)=L
                  stateOld(k,10)=deqps


                  stateNew(k,1)=stateOld(k,1)
                  stateNew(k,2)=stateOld(k,2)
                  stateNew(k,3)=stateOld(k,3)
                  stateNew(k,4)=stateOld(k,4)
                  stateNew(k,5)=stateOld(k,5)
                  stateNew(k,6)=stateOld(k,6)
                  stateNew(k,7)=stateOld(k,7)
                  stateNew(k,8)=stateOld(k,8)
                  stateNew(k,9)=stateOld(k,9)
                  stateNew(k,10)=stateOld(k,10)
              end if
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
              plasticWorkInc=sig_y*deqps
          enerInelasNew(k) =enerInelasOld(k)+stressWorkInc/density(k)
          enddo
      endif 
      return
      end subroutine




