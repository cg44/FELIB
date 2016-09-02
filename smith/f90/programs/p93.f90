program p93     
!------------------------------------------------------------------------------
!      program 9.4 plane strain consolidation of a Biot elasto-plastic
!      solid using 8-node solid quadrilateral elements
!      coupled to 4-node fluid elements : incremental version
!      Mohr-Coulomb failure criterion : viscoplastic strain method
!------------------------------------------------------------------------------
 use new_library       ; use geometry_lib    ;    implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=3,nod=8,nodf=4,nst=4,      &
          ndim=2,ndof, i,k,l,iel,ns,nstep,ntot, nodofs=2 ,iters,limit,inc
 real   ::permx,permy,e,v,det,dtim,theta,phi,snph,coh,psi,cons,p0,tol,      &
          dt,f,dsbar,dq1,dq2,dq3,lode_theta,sigm,pi,dpore,time
 logical::converged    ;  character(len=15) :: element = 'quadrilateral'
!--------------- dynamic arrays----------------------------
 real    ,allocatable :: dee(:,:), points(:,:), coord(:,:), derivf(:,:),    &
                         jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),  &
                         derf(:,:),funf(:), coordf(:,:), bee(:,:), km(:,:), &
                         eld(:), sigma(:), kp(:,:), ke(:,:), g_coord(:,:),  &
                         fun(:),c(:,:), width(:), depth(:), bk(:),stress(:),&
                         vol(:), loads(:), ans(:) ,volf(:,:),tensor(:,:,:), &
                         store_kp(:,:,:),phi0(:),phi1(:),bdylds(:),eps(:),  &
                         evpt(:,:,:),bload(:),eload(:),erate(:),            &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),  &
                         tempdis(:),newdis(:),displ(:)   
 integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:) , g_g(:,:)            
!-------------------------input and initialisation-----------------------------
  open (10,file='p93.dat',status=    'old',action='read')
  open (11,file='p93.res',status='replace',action='write')                    
  read (10,*) nels,nxe,nye,nn,nip,                                         &
              permx, permy, phi, coh, psi, e, v, dtim, nstep, theta ,      &
              cons , p0 , tol , limit
  ndof=nod*2; ntot=ndof+nodf                                             
  allocate (dee(nst,nst),points(nip,ndim),coord(nod,ndim),derivf(ndim,nodf), &
            jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),     &
            derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),bee(nst,ndof),      &
            km(ndof,ndof),eld(ndof),sigma(nst),kp(nodf,nodf),g_g(ntot,nels), &
            ke(ntot,ntot),fun(nod),c(ndof,nodf),width(nxe+1),phi0(nodf),     &
            depth(nye+1),vol(ndof),nf(nodof,nn), g(ntot), volf(ndof,nodf),   &
            g_coord(ndim,nn),g_num(nod,nels),num(nod),weights(nip),phi1(nodf),&
            store_kp(nodf,nodf,nels),tensor(nst+1,nip,nels),eps(nst),evp(nst),&
            evpt(nst,nip,nels),bload(ndof),eload(ndof),erate(nst),devp(nst),  &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst),stress(nst))    
        kay=0.0; kay(1,1)=permx; kay(2,2)=permy
        read (10,*)width , depth  ; pi = acos( -1.); snph=sin(phi*pi/180.)
        dt = 4.*(1.+v)*(1.-2.*v)/(e*(1.-2.*v*snph*snph))
  write(11,'(a,e12.4)') "The viscoplastic timestep is ", dt                
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)
  call deemat (dee,e,v); call sample(element,points,weights)                   
!--------- loop the elements to find nband and set up global arrays------------
  nband = 0                                                                     
 elements_1: do iel = 1 , nels
             call geometry_8qxv(iel,nxe,width,depth,coord,num)
             inc = 0
             do i=1,8; do k=1,2; inc=inc+1;g(inc)=nf(k,num(i));end do;end do
             do i=1,7,2 ; inc=inc+1 ; g(inc) = nf(3,num(i)); end do
             g_num(:,iel)=num; g_coord(:,num)=transpose(coord); g_g(:,iel) = g
             if(nband<bandwidth(g))nband=bandwidth(g)
 end do elements_1                                                            
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                 &
                          "Element ",k,"        ",g_num(:,k); end do        
  write(11,'(2(a,i5))')                                                       &
          "There are ",neq ,"  equations and the half-bandwidth is   ",nband
  allocate(bk(neq*(nband+1)),loads(0:neq),ans(0:neq),bdylds(0:neq),           &
            displ(0:neq),newdis(0:neq),tempdis(0:neq))
             bk = .0 ; loads = .0 ; displ = .0 ;  tensor = .0
!-------------- element stiffness integration and assembly--------------------- 
      elements_2:  do iel = 1 , nels 
               num = g_num(: , iel ); coord=transpose(g_coord(:,num)) 
               g = g_g( : , iel )   ; coordf = coord(1 : 7 : 2, : )
               km = .0; c = .0; kp = .0
           gauss_points_1: do i = 1 , nip
              call shape_der(der,points,i);  jac = matmul(der,coord) 
              det = determinant(jac)  ; call invert(jac)
              deriv = matmul(jac,der) ; tensor(1:2,i,iel) = cons
              tensor(4,i,iel) = cons; tensor(5,i,iel) = p0
              call beemat(bee,deriv); vol(:)=bee(1,:)+bee(2,:)                 
              km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
!----------------------now the fluid contribution------------------------------
               call shape_fun(funf,points,i)
               call shape_der(derf,points,i)  ; derivf=matmul(jac,derf)
         kp=kp+matmul(matmul(transpose(derivf),kay),derivf)*det*weights(i)*dtim
               do l=1,nodf; volf(:,l)=vol(:)*funf(l); end do
               c= c+volf*det*weights(i)               
           end do gauss_points_1
         store_kp( : , : , iel) = kp
         call formke(km,kp,c,ke,theta);call formkv(bk,ke,g,neq)
      end do elements_2
!------------------------reduce left hand side--------------------------------
    call banred(bk,neq)  ; bdylds = .0 ; evpt = .0      ; tempdis = .0        
! --------- enter the time-stepping (load increment) loop---------------------
  time = .0
  time_steps:  do ns = 1 , nstep
       time = time + dtim ; write(11,'(a,e12.4)') "The time is    ", time
                   ans = .0 ; bdylds = .0 ; evpt = .0 ;newdis = .0
      elements_3 : do iel = 1 , nels
                     g = g_g(: , iel ) ; kp = store_kp( : , : , iel)
                     phi0 = loads ( g (ndof + 1  : ))  ! gather
                     phi1 = matmul(kp,phi0)
                     ans(g(ndof+1:))=ans(g(ndof+1:))+ phi1;ans(0)=.0;! scatter
      end do elements_3
!---------------------------- constant loading --------------------------------
           ans(1)=ans(1)-1./24.; ans(3)=ans(3)-1./6. 
           ans(5)=ans(5)-1./12. ; ans(7)=ans(7)-1./6. ; ans(9) = ans(9)-1./24. 
!-------------------------   iteration loop ----------------------------------
        iters = 0
   iterations: do
    iters=iters+1;  loads = ans + bdylds ;   call bacsub(bk,loads)
!------------------------   check convergence  --------------------------------
      newdis = loads ; newdis(nf(3,:))=.0
      call checon(newdis,tempdis,tol,converged)
      if(iters==1)converged=.false. ;  if(converged.or.iters==limit)bdylds=.0
!--------------------- go round the Gauss Points-------------------------------
      elements_4: do iel = 1 , nels
       num = g_num( : , iel ) ; coord = transpose(g_coord( : , num ))
       g = g_g( : , iel )  ;  eld = loads ( g( 1 : ndof) )    ;  bload = .0     
       gauss_points_2 : do i = 1 , nip
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ;   call invert(jac)
          deriv = matmul(jac,der) ; call beemat (bee,deriv);eps=matmul(bee,eld)
          eps = eps -evpt( : , i , iel)    ;        sigma=matmul(dee,eps)
          stress = sigma + tensor( 1:4 , i , iel )
          call invar(stress,sigm,dsbar,lode_theta)                             
!------------------  check whether yield is violated --------------------------
call mocouf (phi, coh , sigm, dsbar , lode_theta , f )
         if(converged.or.iters==limit) then
         devp=stress 
           else
           if(f>=.0) then
           call mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)
           call formm(stress,m1,m2,m3)
           flow=f*(m1*dq1+m2*dq2+m3*dq3)     ;   erate=matmul(flow,stress)
           evp=erate*dt; evpt(:,i,iel)=evpt(:,i,iel)+evp; devp=matmul(dee,evp) 
         end if; end if
      if(f>=.0) then
        eload=matmul(transpose(bee),devp) ; bload=bload+eload*det*weights(i)
      end if
     if(converged.or.iters==limit) then
! ------------------- update stresses and porepressures -----------------------
       tensor (1:4 , i , iel ) = stress   ;   dpore = .0
       call shape_fun( funf , points , i )
       do k = 1 , nodf; dpore = dpore + funf(k)*loads(g(k+ndof)) ; end do
       tensor( 5 , i , iel ) = tensor ( 5 , i , iel ) + dpore
     end if
    end do gauss_points_2
!--------------------  compute the total bodyloads vector ---------------------
    bdylds(g( 1 : ndof)) = bdylds( g( 1 : ndof) ) + bload  ; bdylds(0) = .0
  end do elements_4             
  if(converged.or.iters==limit)exit
 end do iterations
 displ = displ + loads
 write(11,'(a,i5,a)') "It took",iters,"  iterations to converge"
 write(11,'(a,e12.4)') "The displacements are : ", displ(1)
 if(iters==limit)stop
end do time_steps
end program p93
