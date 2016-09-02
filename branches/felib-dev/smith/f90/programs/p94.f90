program p94     
!------------------------------------------------------------------------------
!      program 9.4 plane strain consolidation of a Biot elastic
!      solid using 8-node solid quadrilateral elements
!      coupled to 4-node fluid elements  -  pcg version
!------------------------------------------------------------------------------
 use new_library       ; use geometry_lib   ;   implicit none
 integer::nels,nxe,nye,neq,nn,nr,nip,nodof=3,nod=8,nodf=4,nst=3,ndim=2,      &
          ndof,i,k,l,iel,ns,nstep,ntot,nodofs=2,cjiters,cjits,inc
 real::permx,permy,e,v,det,dtim,theta,x1,x2,time,up,alpha,beta,big,cjtol
 logical :: cj_converged   ;     character(len=15) :: element='quadrilateral'
!---------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: dee(:,:),points(:,:),coord(:,:),derivf(:,:),pmul(:),&
                         jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),   &
                         derf(:,:),funf(:), coordf(:,:), bee(:,:), km(:,:),  &
                         eld(:), sigma(:), kp(:,:), ke(:,:), g_coord(:,:),   &
                         kd(:,:),fun(:), c(:,:), width(:), depth(:),loads(:),&
                         vol(:), storke(:,:,:), ans(:) ,volf(:,:) ,          &
                         p(:),x(:),xnew(:),u(:),diag_precon(:),d(:),         &
                         utemp(:), storkd(:,:,:)  
 integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:) , g_g(:,:)            
!-------------------------input and initialisation-----------------------------
  open (10,file='p94.dat',status=    'old',action='read')
  open (11,file='p94.res',status='replace',action='write')           
  read (10,*) nels,nxe,nye,nn,nip,                                           &
              permx, permy, e,v, dtim, nstep, theta , cjits , cjtol
  ndof=nod*2; ntot=ndof+nodf                                                   
  allocate (dee(nst,nst),points(nip,ndim),coord(nod,ndim),derivf(ndim,nodf), &
            jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),     &
            derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),bee(nst,ndof),      &
            km(ndof,ndof),eld(ndof),sigma(nst),kp(nodf,nodf),g_g(ntot,nels), &
            ke(ntot,ntot),kd(ntot,ntot),fun(nod),c(ndof,nodf),width(nxe+1),  &
            depth(nye+1),vol(ndof),nf(nodof,nn), g(ntot), volf(ndof,nodf),   &
            g_coord(ndim,nn),g_num(nod,nels),num(nod),weights(nip),          &
            storke(ntot,ntot,nels),storkd(ntot,ntot,nels),                   &
            pmul(ntot),utemp(ntot))                                           
            kay=0.0; kay(1,1)=permx; kay(2,2)=permy
            read (10,*)width , depth                                           
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)
  call deemat (dee,e,v); call sample(element,points,weights)
!----------------- loop the elements to  set up global arrays------------------
 elements_1: do iel = 1 , nels
             call geometry_8qxv(iel,nxe,width,depth,coord,num)
             inc = 0
             do i=1,8; do k=1,2;inc=inc+1;g(inc)=nf(k,num(i)); end do; end do
             do i=1,7,2;inc=inc+1;g(inc)=nf(3,num(i)); end do
             g_num(:,iel)=num; g_coord(:,num)=transpose(coord); g_g(:,iel) = g
 end do elements_1                                                             
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                              "Element ",k,"        ",g_num(:,k); end do       
  write(11,'(a,i5,a)') "There are  ",neq, "  equations to be solved"     
  allocate(loads(0:neq),ans(0:neq),p(0:neq),x(0:neq),xnew(0:neq),u(0:neq),    &
           diag_precon(0:neq),d(0:neq))
    loads = .0 ; p = .0; xnew = .0; diag_precon = .0
!------------ element stiffness integration , storage and preconditioner ----- 
      elements_2:  do iel = 1 , nels 
               num = g_num(:,iel); coord=transpose(g_coord(:,num)) 
               g = g_g( :, iel ) ; coordf = coord(1 : 7 : 2, : )
               km = .0; c = .0; kp = .0
           gauss_points_1: do i = 1 , nip
              call shape_der(der,points,i);  jac = matmul(der,coord) 
              det = determinant(jac);call invert(jac); deriv = matmul(jac,der)
              call beemat(bee,deriv); vol(:)=bee(1,:)+bee(2,:)   
              km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
!--------------------------now the fluid contribution--------------------------
               call shape_fun(funf,points,i)
               call shape_der(derf,points,i)  ; derivf=matmul(jac,derf)
         kp=kp+matmul(matmul(transpose(derivf),kay),derivf)*det*weights(i)*dtim
               do l=1,nodf; volf(:,l)=vol(:)*funf(l); end do
               c= c+volf*det*weights(i)               
           end do gauss_points_1
          call fmkdke(km,kp,c,ke,kd,theta) 
          storke(: , : , iel) = ke   ;    storkd ( : , : , iel ) = kd
          do k=1,ndof;diag_precon(g(k))=diag_precon(g(k))+theta*km(k,k);end do
          do k=1 , nodf
           diag_precon(g(ndof+k))=diag_precon(g(ndof+k))-theta*theta*kp(k,k)
          end do
      end do elements_2    
          diag_precon(1:neq) = 1./ diag_precon(1:neq) ; diag_precon(0) = .0    
! ------------------------ enter the time-stepping loop------------------------
      time = .0     
  time_steps:  do ns = 1 , nstep
      ans = .0    ;   time=time+dtim;write(11,'(a,e12.4)')"The time is",time  
      elements_3 : do iel = 1 , nels
                     g = g_g( : , iel );  kd = storkd ( : , : , iel )
                     pmul = loads ( g ) ; utemp = matmul ( kd , pmul )
!dir$ ivdep
           do i = 1 , ntot
            ans( g ( i )) = ans( g ( i )) + utemp ( i )
           end do
      end do elements_3          ;  ans(0) = .0
!--------------------------   ramp loading  -----------------------------------
        x1=(.1*ns+.1*(theta-1.))/6.; x2=x1*4.
        if(ns>10) then
           ans(1)=ans(1)-1./6.; ans(3)=ans(3)-2./3.
           ans(4)=ans(4)-1./6.
        else if(ns<10) then
           ans(1)=ans(1)-x1;ans(3)=ans(3)-x2; ans(4)=ans(4)-x1
        end if             ;  
      d = diag_precon*ans    ;   p = d ; x = .0  ! depends on starting x = .0
!--------------------   solve the simultaneous equations by pcg ---------------
       cjiters = 0
      conjugate_gradients:  do
       cjiters = cjiters + 1 ; u = .0
      elements_4 : do iel = 1 , nels
                      g = g_g( : , iel ); ke = storke( : , : ,iel)
                      pmul = p(g); utemp=matmul(ke,pmul)
!dir$ ivdep
           do i = 1 , ntot
              u(g(i)) = u(g(i)) +  utemp(i)
           end do
      end do elements_4 
!----------------------------pcg process --------------------------------------
    up =dot_product(ans,d); alpha=up/dot_product(p,u)
    xnew = x + p* alpha; ans = ans - u*alpha; d = diag_precon*ans
    beta = dot_product(ans,d)/up; p = d + p * beta
    big = .0; cj_converged = .true.
    do i = 1,neq; if(abs(xnew(i))>big)big=abs(xnew(i)); end do
    do i = 1,neq; if(abs(xnew(i)-x(i))/big>cjtol)cj_converged=.false.;end do
    x = xnew
    if(cj_converged.or.cjiters==cjits) exit
      end do conjugate_gradients
!----------- end of pcg process----------------------------------------------
    write(11,'(a,i5,a)')                                                       &
             "Conjugate gradients took ",cjiters, "  iterations to converge"
    ans = xnew           ; ans(0) = .0        ; loads = ans

        write(11,'(a)') " The nodal displacements and porepressures are    :"
        do k=1,23,22; write(11,'(i5,a,3e12.4)')k,"    ",ans(nf(:,k)) ; end do   
!-------------------recover stresses at centroidal gauss-point------------
      elements_5 :  do iel = 1 , nels
               num = g_num(:,iel); coord=transpose(g_coord(:,num))
               g = g_g( : , iel ); eld = ans( g ( 1 : ndof ) )
         !    print*,"The Gauss Point effective stresses for element",iel,"are"
            gauss_pts_2: do i = 1,nip
              call shape_der (der,points,i);  jac= matmul(der,coord)
              call invert ( jac );    deriv= matmul(jac,der)
              bee= 0.;call beemat(bee,deriv);sigma= matmul(dee,matmul(bee,eld))
         !     print*,"Point    ",i       ;!  print*,sigma
           end do gauss_pts_2 
      end do elements_5
   end do time_steps
 end program p94
