program p91     
!------------------------------------------------------------------------------
!      program 9.1 plane strain consolidation of a Biot elastic
!      solid using 8-node solid quadrilateral elements
!      coupled to 4-node fluid elements
!------------------------------------------------------------------------------
 use new_library       ; use geometry_lib   ;   implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=3,nod=8,nodf=4,nst=3,     &
          ndim=2,ndof, i,k,l,iel,ns,nstep,ntot, nodofs=2 ,inc
 real::permx,permy,e,v,det,dtim,theta,x1,x2,time
 character (len=15) :: element = 'quadrilateral'
!---------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: dee(:,:), points(:,:), coord(:,:), derivf(:,:),    &
                         jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),  &
                         derf(:,:),funf(:), coordf(:,:), bee(:,:), km(:,:), &
                         eld(:), sigma(:), kp(:,:), ke(:,:), g_coord(:,:),  &
                         kd(:,:),fun(:), c(:,:), width(:), depth(:), bk(:), &
                         vol(:), pb(:,:), loads(:), ans(:) ,volf(:,:)      
 integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:) , g_g(:,:)             
!--------------------------input and initialisation----------------------------
  open (10,file='p91.dat',status=    'old',action='read')
  open (11,file='p91.res',status='replace',action='write')                   
  read (10,*) nels,nxe,nye,nn,nip,                                          &
              permx, permy, e,v, dtim, nstep, theta 
  ndof=nod*2; ntot=ndof+nodf                                                   
  allocate (dee(nst,nst),points(nip,ndim),coord(nod,ndim),derivf(ndim,nodf), &
            jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),     &
            derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),bee(nst,ndof),      &
            km(ndof,ndof),eld(ndof),sigma(nst),kp(nodf,nodf),g_g(ntot,nels), &
            ke(ntot,ntot),kd(ntot,ntot),fun(nod),c(ndof,nodf),width(nxe+1),  &
            depth(nye+1),vol(ndof),nf(nodof,nn), g(ntot), volf(ndof,nodf),   &
            g_coord(ndim,nn),g_num(nod,nels),num(nod),weights(nip))          
            kay=0.0; kay(1,1)=permx; kay(2,2)=permy
            read (10,*)width , depth                                         
  nf=1; read(10,*) nr; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)
  call deemat (dee,e,v); call sample(element,points,weights)                  
!--------- loop the elements to find nband and set up global arrays------------
  nband = 0                                                                     
 elements_1: do iel = 1 , nels
             call geometry_8qxv(iel,nxe,width,depth,coord,num)
             inc=0
             do i=1,8; do k=1,2; inc=inc+1;g(inc)=nf(k,num(i));end do;end do
             do i=1,7,2; inc=inc+1;g(inc)=nf(3,num(i)); end do
             g_num(:,iel)=num; g_coord(:,num)=transpose(coord); g_g(:,iel)= g
             if(nband<bandwidth(g))nband=bandwidth(g)
 end do elements_1    
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do 
  write(11,'(2(a,i5))')                                                       &
          "There are ",neq,"  equations and the half-bandwidth is   ",nband
  allocate(bk(neq*(nband+1)),pb(neq,2*(nband+1)-1),loads(0:neq),ans(0:neq))
    pb = .0 ; bk = .0 ; loads = .0
!--------------- element stiffness integration and assembly---------------------

      elements_2:  do iel = 1 , nels 
               num = g_num(: , iel ); coord=transpose(g_coord(:,num)) 
               g = g_g( : , iel )   ; coordf = coord(1 : 7 : 2, : )
               km = .0; c = .0; kp = .0
           gauss_points_1: do i = 1 , nip
              call shape_der(der,points,i);  jac = matmul(der,coord) 
              det = determinant(jac ); call invert(jac);deriv = matmul(jac,der)
              call beemat(bee,deriv); vol(:)=bee(1,:)+bee(2,:)                 
              km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
!--------------------------now the fluid contribution--------------------------
               call shape_fun(funf,points,i)
               call shape_der(derf,points,i)  ; derivf=matmul(jac,derf)
         kp=kp+matmul(matmul(transpose(derivf),kay),derivf)*det*weights(i)*dtim
               do l=1,nodf; volf(:,l)=vol(:)*funf(l); end do
               c= c+volf*det*weights(i)               
           end do gauss_points_1
 call fmkdke(km,kp,c,ke,kd,theta);call formkv(bk,ke,g,neq);call formtb(pb,kd,g)
      end do elements_2
!------------------------factorise left hand side------------------------------
    call banred(bk,neq)                       
! ------------------------ enter the time-stepping loop------------------------
      time = .0
  time_steps:  do ns = 1 , nstep
      time=time+dtim;write(11,'(a,e12.4)')"The time is",time
     call bantmul(pb,loads,ans)
!   ramp loading
        x1=(.1*ns+.1*(theta-1.))/6.; x2=x1*4.
        if(ns>10) then
           ans(1)=ans(1)-1./6.; ans(3)=ans(3)-2./3.
           ans(4)=ans(4)-1./6.
        else if(ns<10) then
           ans(1)=ans(1)-x1;ans(3)=ans(3)-x2; ans(4)=ans(4)-x1
        end if
        call bacsub(bk,ans) ; loads=ans
        write(11,'(a)') " The nodal displacements and porepressures are    :"
        do k=1,23,22; write(11,'(i5,a,3e12.4)')k,"    ",ans(nf(:,k)) ; end do  
!-------------------recover stresses at  Gauss-points--------------------------
      elements_3 :  do iel = 1 , nels
               num = g_num(:,iel); coord=transpose(g_coord(:,num))
               g = g_g( : , iel );  eld = ans( g ( 1 : ndof ) )
         !    print*,"The Gauss Point effective stresses for element",iel,"are"
            gauss_points_2: do i = 1,nip
              call shape_der (der,points,i);  jac= matmul(der,coord)
              call invert ( jac );    deriv= matmul(jac,der)
              bee= 0.;call beemat(bee,deriv);sigma= matmul(dee,matmul(bee,eld))
         !     print*,"Point    ",i       ;!  print*,sigma
           end do gauss_points_2 
      end do elements_3
   end do time_steps
 end program p91
