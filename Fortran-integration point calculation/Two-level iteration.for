
      program two level iteration
      Implicit None
      real(8)::y,Cb,b,f,Cnc,Cnw,yc,db,CMt,dwz1,dwz2,dwz3,dwz,dcz1,
     *dcz2,dcz3,dcz,CK,Pt,d,deqps0,deqpsn,fdeqps0,fdeqps,ddeqps0,deqps,
     *e,xnu,S1,da,Ckc,Ckw,Cm,Ca,G,f0,f00,Ck0,Ck00,dc0,dw0,ds,Pt0,twomu,
     *alamda,vk,thremu,s11,s12,s13,s21,s22,s23,s31,s32,s33,sm,vmises,
     *peeqOld,peeq_rate,sig_y,dc,dw,sigdif,wc,x0,y0,dy0,deqps1,deqps2,
     *sigama,dt
      integer::max_iter,j,n,i
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
         f    = 0.249996721744537
         d    = 6.381899118423462E-003
         dw   = 553267328.000000     
         dc   = 142940032.000000     
         dt   = 4.2936843E-10 
c     Calculate plastic strain increment-1
      deqps0=1.d-6
      x0=1.d-50
c     Tolerance
      wc = 1.d-3
c     The discrepancy in plastic strain increments computed twice.
      sigama=1.d-3
c     Maximum number of iterations
      max_iter = 1000
c     Original function      
      y0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *   * ((CMt*x0/(dt*yc))**(1.d0/Cm))+thremu*x0-vmises
c     Derivative of the original function 
      dy0 = CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *    * ((CMt*x0/(dt*yc))**((1.d0/Cm)-1))*(CMt/(yc*dt))+thremu
      do j = 1, max_iter
          if (deqps0.Gt.x0) then
              fdeqps0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *                  * ((CMt*deqps0/(dt*yc))**(1.d0/Cm))+thremu
     *                  * deqps0-vmises

              ddeqps0 = CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *                * ((CMt*deqps0/(dt*yc))**((1.d0/Cm)-1))*(CMt
     *                / (yc*dt))+thremu
          else
              fdeqps0 = y0+dy0*(deqps0-x0)
              ddeqps0 = dy0
          endif
          deqps=deqps0-fdeqps0/ddeqps0
          if (abs((deqps-deqps0)/deqps).gt.wc) then
              deqps0=deqps
          else           
              exit
          endif
      enddo  
c     The plastic strain increment from the first computation:deqps1
      deqps1=deqps
c     
      n=0
1     continue
c     Calculate dislocation wall volume fraction:f
      f=f00+(f0-f00)*exp((-CMt*(peeqOld+deqps1))/y)
c     Calculate varible:CK
      CK=Ck00+(Ck0-Ck00)*exp(-Cb*CMt*(peeqOld+deqps1))
c     Calculate the increment of dislocation density in the interior of the cell
      dcz1=(da*CMt*deqps1*sqrt(dw))/(b*sqrt(3.d0))
      dcz2=(db*CMt*deqps1*6.d0)/(b*(1-f)**(1.d0/3.d0)*d)
      dcz3=dc*CMt*deqps1*Ckc
     *    *((CMt*deqps1/(dt*yc))**(1.d0/(-Cnc)))
      dcz=dcz1-dcz2-dcz3
c     Calculate the increment of dislocation density in the cell walls
      dwz1=(CMt*deqps1*6.d0*db*((1.d0-f)**(2.d0/3.d0)))
     *    /(d*b*f)
      dwz2=(sqrt(3.d0)*sqrt(dw)*db*(1.d0-f)*CMt*deqps1)
     *    /(b*f)
      dwz3=CMt*deqps1*dw*Ckw
     *    *((CMt*deqps1/(dt*yc))**(1.d0/(-Cnw)))
      dwz=dwz1+dwz2-dwz3
c     Update the dislocation density in the interior of the cell
      dc=dc+dcz
c     Update the dislocation density in the cell walls
      dw=dw+dwz
c     Update total dislocation density
      Pt=f*(dw)+(1-f)*(dc)
c     Update dislocation cell size
      d=CK/sqrt(Pt)
c     Calculate plastic strain increment-2
      deqps0=1.d-6
      x0=1.d-50
      y0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *   * ((CMt*x0/(dt*yc))**(1.d0/Cm))+thremu*x0-vmises
     
      dy0 = CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *    * ((CMt*x0/(dt*yc))**((1.d0/Cm)-1))*(CMt/(yc*dt))+thremu

      do i = 1, max_iter
          if (deqps0.GT.x0) then
              fdeqps0 = S1+CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))
     *                * ((CMt*deqps0/(dt*yc))**(1.d0/Cm))+thremu*deqps0
     *                - vmises

              ddeqps0 = CMt*G*b*Ca*(f*sqrt(dw)+(1-f)*sqrt(dc))*(1.d0/Cm)
     *                * ((CMt*deqps0/(dt*yc))**((1.d0/Cm)-1))*(CMt/(yc
     *                * dt))+thremu
          else
              fdeqps0=y0+dy0*(deqps0-x0)
              ddeqps0=dy0
          endif
          deqpsn=deqps0-fdeqps0/ddeqps0
          if (abs((deqpsn-deqps0)/deqpsn).gt.wc) then
              deqps0=deqpsn
          else
              exit
          end if
      enddo
c     The plastic strain increment from the second computation.
      deqps2=deqpsn
c     Calculate the error in plastic strain increment between the two computations.
      if (abs(deqps1-deqps2)/deqps1.GT.sigama) then
          n=n+1
          deqps1=deqps2      
          go to 1
      endif
      write(*,*)'Result'
      write(*,*)'Plastic strain increment=',deqps2
      write(*,*)'Cell wall dislocation density increment=',dwz
      write(*,*)'Intracellular dislocation density increment=',dcz
      write(*,*)'The number of iterations.=',n
      end




