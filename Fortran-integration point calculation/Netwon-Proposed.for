   
      program Netwon Proposed
      Implicit none
c     Date type
      real(8):: y,Cb,b,f,Cnc,Cnw,yc,db,CMt,dwz1,dwz2,dwz3,dwz,dcz1,
     *dcz2,dcz3,dcz,CK,Pt,d,deqps0,deqpsn,fdeqps0,fdeqps,ddeqps0,deqps,
     *e,xnu,S1,da,Ckc,Ckw,Cm,Ca,G,f0,f00,Ck0,Ck00,dc0,dw0,ds,Pt0,twomu,
     *alamda,vk,thremu,s11,s12,s13,s21,s22,s23,s31,s32,s33,sm,vmises,
     *peeqOld,peeq_rate,sig_y,dc,dw,sigdif,wc,x0,y0,dy0,deqps1,deqps2,
     *dwz0,dcz0,f1,f2,f3,J11,J12,J13,J21,J22,J23,J31,J32,J33,invJ11,
     *invJ12,invJ13,invJ21,invJ22,invJ23,invJ31,invJ32,invJ33,detJ,
     *idetJ,wcmax3,dt,a,deqpswc,dwzwc,dczwc
      integer::L,max_iter
      real(8):: J(3,3),invJ(3,3),cz(1,3),czz(3,1),hsjz(1,3),ff(3,1),
     *mbz(3,1),wcmax(1,3),wcmax2(3,1)
c     Dislocation density constitutive variable    
      da=1.54d-1
      db=7.8d-2
      Ckc=1.86d1
      Ckw=3.28d1
      Cnc=8.98d1
      Cnw=9.03d1
      S1=7.92d2
      CMt=3.06d0
      G=8.2d4
      b=2.48d-7
      ca=2.5d-1
      f00=6d-2
      f0=2.5d-1
      yc=1d7
      Cm=6.08d1
      thremu= 242307.690085744      
      y=3.2d0
      Ck0=100
      Ck00=1
      Cb=2.6d-1
c     Initial Data for Integration Points
      vmises  = 1008.86134912587     
      peeqOld = 1.803302802727558E-005 
         dw   = 553267328.000000     
         dc   = 142940032.000000     
         dt   = 4.2936843E-10
c     Solving a system of nonlinear equations
      deqps0=1.d-5
      dwz0=1d5  
      dcz0=1d5  
      x0=1.d-50
c     Tolerance
      wc = 1.d-3
c     Maximum number of iterations
      max_iter = 1000  

      do L = 1, max_iter
c         Initial value matrix
          cz=reshape([deqps0,dwz0,dcz0],shape(cz))
          czz=transpose(cz)
c         Nonlinear equations:f1
          f1 = S1+CMt*G*b*ca*((f00+(f0-f00)*exp((-CMt*(peeqOld+deqps0))
     *       / y))*sqrt(dw+dwz0)+(1-(f00+(f0-f00)*exp((-CMt*(peeqOld
     *       + deqps0))/y)))*sqrt(dc+dcz0))*((CMt*deqps0/(dt*yc))
     *       **(1.d0/Cm))+thremu*deqps0-vmises
   
c        Nonlinear equations:f2         
         f2 = (CMt*deqps0*6*db*((1-(f00+(f0-f00)
     *      * exp((-CMt*(peeqOld+deqps0))/y)))**(2.d0/3)))
     *      / (((Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld+deqps0)))
     *      / sqrt((f00+(f0-f00)*exp((-CMt*(peeqOld+deqps0))/y))
     *      * (dw+dwz0)+(1-(f00+(f0-f00)*exp((-CMt*(peeqOld+deqps0))
     *      / y)))*(dc+dcz0)))*b*(f00+(f0-f00)*exp((-CMt*(peeqOld
     *      + deqps0))/y)))+(sqrt(3.d0)*sqrt(dw+dwz0)*db*(1-(f00
     *      + (f0-f00)*exp((-CMt*(peeqOld+deqps0))/y)))
     *      * CMt*deqps0)/(b*(f00+(f0-f00)*exp((-CMt*(peeqOld+deqps0))
     *      / y)))-CMt*deqps0*(dw+dwz0)*Ckw*((CMt*deqps0/(dt*yc))
     *      ** (1.d0/(-Cnw)))-dwz0
c          Nonlinear equations:f3
           f3 = (da*CMt*deqps0*sqrt(dw+dwz0))/(b*sqrt(3.d0))
     *        - (db*CMt*deqps0*6)/(b*(1-(f00+(f0-f00)
     *        * exp((-CMt*(peeqOld+deqps0))/y)))**(1.d0/3)
     *        * ((Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld
     *        + deqps0)))/sqrt((f00+(f0-f00)*exp((-CMt
     *        * (peeqOld+deqps0))/y))*(dw+dwz0)+(1-(f00
     *        + (f0-f00)*exp((-CMt*(peeqOld+deqps0))/y)))*(dc
     *        + dcz0))))-(dcz0+dc)*CMt*deqps0*Ckc*((CMt
     *        * deqps0/(dt*yc))**(1.d0/(-Cnc)))-dcz0
      
c         Nonlinear equations:ff
          hsjz=reshape([f1,f2,f3],shape(hsjz))
          ff=transpose(hsjz)

c         Jacobian matrix       
c         J=[diff(fcdeqps,deqps) diff(fcdeqps,dwz) diff(fcdeqps,dcz) 
c            diff(fcdwz,  deqps) diff(fcdwz,  dwz) diff(fcdwz,  dcz)
c            diff(fcdcz,  deqps) diff(fcdcz,dwz)   diff(fcdcz,  dcz)      ]
          J11 = thremu+CMt*Ca*G*b*((CMt*exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(dc+dcz0)**(1.d0/2)*(f0-f00))/y
     *        - (CMt*exp(-(CMt*(deqps0+peeqOld))/y)
     *        * (dw+dwz0)**(1.d0/2)*(f0-f00))/y)*((CMt
     *        * deqps0)/(dt*yc))**(1.d0/Cm)-(CMt**(2)*Ca*G*b
     *        * ((dc+dcz0)**(1.d0/2)*(f00+exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(f0-f00)-1)-(dw+dwz0)**(1.d0/2)
     *        * (f00+exp(-(CMt*(deqps0+peeqOld))/y)
     *        * (f0-f00)))*((CMt*deqps0)/(dt*yc))
     *        ** ((1.d0/Cm)-1))/(Cm*dt*yc)

          J12 = (CMt*Ca*G*b*(f00+exp(-(CMt*(deqps0+peeqOld))
     *        / y)*(f0-f00))*((CMt*deqps0)/(dt*yc))
     *        ** (1.d0/Cm))/(2*(dw+dwz0)**(1.d0/2))

          J13 = -(CMt*Ca*G*b*(f00+exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(f0-f00)-1)*((CMt*deqps0)
     *        / (dt*yc))**(1.d0/Cm))/(2*sqrt(dc+dcz0))
      
          J21 = (6*CMt*db*((dw+dwz0)*(f00+exp(-(CMt*(deqps0+peeqOld))/y)
     *        * (f0 - f00)) - (dc + dcz0)*(f00 + exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(f0 - f00)-1))**(1.d0/2)*(1-exp(-(CMt
     *        * (deqps0 +peeqOld))/y)*(f0 - f00)-f00)**(2.d0/3))/(b
     *        * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0 - f00))*(Ck00 
     *        + exp(-CMt*Cb*(deqps0 + peeqOld))*(Ck0-Ck00)))
     *        - (3**(1.d0/2)*CMt*db*(dw + dwz0)**(1.d0/2)*(f00 
     *        + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 - f00) - 1))
     *        / (b*(f00+exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 - f00)))
     *        - (CMt*Ckw*(dw + dwz0))/((CMt*deqps0)/(dt*yc))**(1.d0/Cnw)
     *        + (CMt**(2)*Ckw*deqps0*(dw+dwz0))/(Cnw*dt*yc*((CMt*deqps0)
     *        / (dt*yc))**(1.d0/Cnw + 1))+(3*CMt*db*deqps0*((CMt*exp(
     *        - (CMt*(deqps0+peeqOld))/y)*(dc+dcz0)
     *        * (f0 -f00))/y - (CMt*exp(-(CMt*(deqps0 + peeqOld))/y)*(dw
     *        + dwz0)*(f0 - f00))/y)*(1 - exp(-(CMt*(deqps0 + peeqOld))
     *        / y)*(f0-f00)-f00)**(2.d0/3))/(b*(f00+exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(f0 -f00))*(Ck00 + exp(-CMt*Cb*(deqps0 
     *        + peeqOld))*(Ck0 - Ck00))*((dw + dwz0)*(f00 + exp(-(CMt
     *        * (deqps0 + peeqOld))/y)*(f0 - f00))-(dc + dcz0)*(f00 
     *        + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0-f00)-1))**(1.d0/2))
     *        + (3**(1.d0/2)*CMt**(2)*db*deqps0*exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(dw + dwz0)**(1.d0/2)*(f0-f00))
     *        / (b*y*(f00 + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0-f00)))
     *        + (4*CMt**(2)*db*deqps0*exp(-(CMt*(deqps0+peeqOld))/y)*(f0
     *        - f00)*((dw + dwz0)*(f00 + exp(-(CMt*(deqps0+peeqOld))/y)
     *        * (f0 - f00))-(dc + dcz0)*(f00 + exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(f0 - f00) - 1))**(1.d0/2))/(b*y*(f00 
     *        + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 - f00))*(Ck00 
     *        + exp(-CMt*Cb*(deqps0 + peeqOld))*(Ck0 - Ck00))
     *        * (1 - exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 - f00) -f00)
     *        ** (1.d0/3))+(6*CMt**(2)*db*deqps0*exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(f0-f00)*((dw + dwz0)*(f00 + exp(-(CMt
     *        * (deqps0 + peeqOld))/y)*(f0 - f00)) 
     *        - (dc + dcz0)*(f00 + exp(-(CMt*(deqps0 + peeqOld))/y)
     *        * (f0 - f00) -1))**(1.d0/2)*(1-exp(-(CMt*(deqps0+peeqOld))
     *        / y)*(f0 - f00) - f00)**(2.d0/3))/(b*y*(f00  
     *        + exp(-(CMt*(deqps0+peeqOld))/y)*(f0 - f00))**(2)*(Ck00 
     *        + exp(-CMt*Cb*(deqps0+peeqOld))*(Ck0 - Ck00)))
     *        + (6*CMt**(2)*Cb*db*deqps0*exp(-CMt*Cb*(deqps0+peeqOld))
     *        * (Ck0 - Ck00)*((dw + dwz0)*(f00 + exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(f0 - f00)) - (dc + dcz0)*(f00 
     *        + exp(-(CMt*(deqps0 + peeqOld))/y)
     *        * (f0 - f00) -1))**(1.d0/2)*(1-exp(-(CMt*(deqps0+peeqOld))
     *        / y)*(f0- f00) -f00)**(2.d0/3))/(b*(f00+exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(f0 - f00))*(Ck00 + exp(-CMt*Cb*(deqps0 
     *        + peeqOld))*(Ck0-Ck00))**(2)) - (3**(1.d0/2)*CMt**(2)*db
     *        * deqps0*exp(-(CMt*(deqps0 + peeqOld))/y)*(dw + dwz0)
     *        ** (1.d0/2)*(f0 - f00)*(f00 + exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(f0 - f00) - 1))
     *        / (b*y*(f00 + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 -f00))
     *        **(2))

          J22 = (3*CMt*db*deqps0*(1-exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(f0-f00)-f00)**(2.d0/3))/(b
     *        * (Ck00+exp(-CMt*Cb*(deqps0+peeqOld))*(Ck0
     *        - Ck00))*((dw+dwz0)*(f00 +exp(-(CMt*(deqps0 
     *        + peeqOld))/y)*(f0-f00))-(dc + dcz0)*(f00 
     *        + exp(-(CMt*(deqps0 + peeqOld))/y)*(f0 - f00)
     *        - 1))**(1.d0/2)) - (CMt*Ckw*deqps0)
     *        / ((CMt*deqps0)/(dt*yc))**(1.d0/Cnw)
     *        - (3**(1.d0/2)*CMt*db*deqps0*(f00+exp(-(CMt
     *        * (deqps0 + peeqOld))/y)*(f0 - f00) - 1))/(2*b
     *        * (dw + dwz0)**(1.d0/2)*(f00+exp(-(CMt*(deqps0
     *        + peeqOld))/y)*(f0-f00)))-1

              J23 = (3*CMt*db*deqps0*(1-exp(-(CMt*(deqps0
     *            + peeqOld))/y)*(f0-f00)-f00)**(5.d0/3))/(b
     *            * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0
     *            - f00))*(Ck00 + exp(-CMt*Cb*(deqps0+peeqOld))
     *            * (Ck0-Ck00))*((dw+dwz0)*(f00+exp(-(CMt
     *            * (deqps0+peeqOld))/y)*(f0 - f00))-(dc+dcz0)
     *            * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0
     *            - f00)-1))**(1.d0/2))

              J31 = ((3**(1.d0/2)*CMt*da*(dw+dwz0)**(1.d0/2))
     *            / (3*b)-(CMt*Ckc*(dc+dcz0))/((CMt*deqps0)/(dt
     *            * yc))**(1.d0/Cnc)-(6*CMt*db*((dw + dwz0)*(f00
     *            + exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00))-(dc
     *            + dcz0)*(f00 + exp(-(CMt*(deqps0+peeqOld))/y)
     *            * (f0-f00)-1))**(1.d0/2))/(b*(Ck00 + exp(-CMt
     *            * Cb*(deqps0+peeqOld))*(Ck0-Ck00))*(1-exp(
     *            - (CMt*(deqps0+peeqOld))/y)*(f0-f00)-f00)
     *            ** (1.d0/3))- (3*CMt*db*deqps0*((CMt*exp(-(CMt
     *            * (deqps0 + peeqOld))/y)*(dc+dcz0)*(f0 - f00))
     *            / y - (CMt*exp(-(CMt*(deqps0+peeqOld))/y)*(dw 
     *            + dwz0)*(f0 - f00))/y))/(b*(Ck00
     *            + exp(-CMt*Cb*(deqps0 + peeqOld))*(Ck0 -Ck00))
     *            * ((dw+dwz0)*(f00+exp(-(CMt*(deqps0+peeqOld))
     *            / y)*(f0-f00))-(dc+dcz0)*(f00+exp(-(CMt
     *            * (deqps0+peeqOld))/y)*(f0 - f00) - 1))
     *            ** (1.d0/2)*(1 -exp(-(CMt*(deqps0 + peeqOld))
     *            / y)*(f0 - f00) - f00)**(1.d0/3)) + (CMt**(2)
     *            * Ckc*deqps0*(dc + dcz0))/(Cnc*dt*yc*((CMt
     *            * deqps0)/(dt*yc))**(1.d0/Cnc+1))+(2*CMt**(2)
     *            * db*deqps0*exp(-(CMt*(deqps0+peeqOld))
     *            / y)*(f0 - f00)*((dw + dwz0)*(f00 + exp(-(CMt
     *            * (deqps0 + peeqOld))/y)*(f0-f00))-(dc+dcz0)
     *            * (f00+exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00)
     *            - 1))**(1.d0/2))/(b*y*(Ck00+exp(-CMt*Cb
     *            * (deqps0+peeqOld))*(Ck0 - Ck00))*(1 - exp(
     *            - (CMt*(deqps0+peeqOld))/y)*(f0-f00)-f00)
     *            ** (4.d0/3))-(6*CMt**(2)*Cb*db*deqps0*exp(-CMt
     *            * Cb*(deqps0 + peeqOld))*(Ck0 - Ck00)*((dw 
     *            + dwz0)*(f00 + exp(-(CMt*(deqps0
     *            + peeqOld))/y)*(f0 - f00)) - (dc + dcz0)*(f00
     *            + exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00)-1))
     *            ** (1.d0/2))/(b*(Ck00 + exp(-CMt*Cb*(deqps0
     *            + peeqOld))*(Ck0 - Ck00))**(2)*(1 - exp(-(CMt
     *            * (deqps0+peeqOld))/y)*(f0-f00)-f00)
     *            ** (1.d0/3)))

              J32 = (3**(1.d0/2)*CMt*da*deqps0)/(6*b*(dw + dwz0)
     *            ** (1.d0/2))-(3*CMt*db*deqps0*(f00+exp(-(CMt
     *            * (deqps0 + peeqOld))/y)*(f0-f00)))/(b*(Ck00
     *            + exp(-CMt*Cb*(deqps0+peeqOld))*(Ck0-Ck00))
     *            * ((dw+dwz0)*(f00+exp(-(CMt*(deqps0+peeqOld))
     *            / y)*(f0-f00))-(dc+dcz0)*(f00+exp(-(CMt
     *            * (deqps0+peeqOld))/y)*(f0-f00)-1))**(1.d0/2)
     *            * (1-exp(-(CMt*(deqps0+peeqOld))/y)*(f0-f00)
     *            - f00)**(1.d0/3))

              J33 = -(CMt*Ckc*deqps0)/((CMt*deqps0)/(dt*yc))
     *            ** (1.d0/Cnc)-(3*CMt*db*deqps0*(1-exp(-(CMt
     *            * (deqps0+peeqOld))/y)*(f0 - f00)-f00)
     *            ** (2.d0/3))/(b*(Ck00+exp(-CMt*Cb*(deqps0
     *            + peeqOld))*(Ck0 - Ck00))*((dw+dwz0)*(f00 
     *            + exp(-(CMt*(deqps0 +peeqOld))/y)*(f0-f00))
     *            - (dc+dcz0)*(f00+exp(-(CMt*(deqps0+peeqOld))
     *            / y)*(f0 - f00) - 1))**(1.d0/2)) - 1  
c            Inverse of the Jacobian matrix
             detJ = J11*J22*J33+J12*J23*J31+J13*J21*J32
     *            - J31*J22*J13-J21*J12*J33-J11*J23*J32
             
             invJ11=(1.d0/detJ)*(J22*J33-J23*J32)
             invJ12=(1.d0/detJ)*(-1*(J12*J33-J32*J13))
             invJ13=(1.d0/detJ)*(J12*J23-J13*J22)
             invJ21=(1.d0/detJ)*(-1*(J21*J33-J31*J23))
             invJ22=(1.d0/detJ)*(J11*J33-J13*J31)
             invJ23=(1.d0/detJ)*(-1*(J11*J23-J21*J13))
             invJ31=(1.d0/detJ)*((J21*J32-J31*J22))
             invJ32=(1.d0/detJ)*(-1*(J11*J32-J31*J12))
             invJ33=(1.d0/detJ)*(J11*J22-J21*J12)
             invJ=reshape([invJ11,invJ12,invJ13,invJ21,invJ22,invJ23,
     *       invJ31,invJ32,invJ33],shape(invJ))

             invJ=transpose(invJ)
c            Calculation of the target value matrix:mbz
             mbz=czz-matmul(invJ,ff)
c            Value of variables
             deqps=mbz(1,1)
             dwz=mbz(2,1)
             dcz=mbz(3,1)
c            Tolerance of variables
             deqpswc=abs((deqps-deqps0)/deqps)
             dwzwc=abs((dwz-dwz0)/dwz)
             dczwc=abs((dcz-dcz0)/dcz)
c            Maximum tolerance of variables:wcmax3
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
             end if
         enddo 
c        Result
         write(*,*)'Result:'
         write(*,*)'Number of iterations=',L
         write(*,*)'Plastic strain increment=',deqps
         write(*,*)'Cell wall dislocation density increment=',dwz
         write(*,*)'Intracellular dislocation density increment=',dcz
      end
           

   
   

     
               
      