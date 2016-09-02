program p90     
!-----------------------------------------------------------------------------
!      program 9.0 steady state Navier-Stokes equation
!      using 8-node velocity quadrilateral elements
!      coupled to 4-node pressure quadrilateral elements ; u-p-v order
!-----------------------------------------------------------------------------
 use new_library     ; use geometry_lib     ;   implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=3,nod=8,nodf=4,ndim=2,       &
          i,k,iel,ntot,limit ,fixed_nodes ,iters , inc
 real::visc, rho, det ,ubar, vbar , tol  ; logical :: converged
 character (len=15) :: element = 'quadrilateral'
!----------------------------- dynamic arrays----------------------------------
 real    ,allocatable :: points(:,:), coord(:,:),derivf(:,:),fun(:),work(:,:),&
                         jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:)  ,  &
                         derf(:,:),funf(:), coordf(:,:),ke(:,:), g_coord(:,:),&
                         width(:), depth(:),c11(:,:),c21(:,:),c12(:,:),val(:),&
                         c23(:,:),c32(:,:), pb(:,:), loads(:), oldlds(:) ,    &
                         funny(:,:),row1(:,:),row2(:,:),uvel(:),vvel(:)  ,    &
                         funnyf(:,:),rowf(:,:)
 integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:) , g_g(:,:) ,no(:),    &
                         sense(:), node(:)                                     
!----------------------------input and initialisation--------------------------
  open (10,file='p90.dat',status=    'old',action='read')
  open (11,file='p90.res',status='replace',action='write')                  
  read (10,*) nels,nxe,nye,nn,nip,visc,rho,tol,limit 
        ntot=nod+nodf+nod    
  allocate (points(nip,ndim),coord(nod,ndim),derivf(ndim,nodf), &
            jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),     &
            derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),funny(nod,1),       &
            g_g(ntot,nels),c11(nod,nod),c12(nod,nodf),c21(nodf,nod),g(ntot), &
            ke(ntot,ntot),fun(nod),width(nxe+1),depth(nye+1),nf(nodof,nn),   &
            g_coord(ndim,nn),g_num(nod,nels),num(nod),weights(nip),          &
            c32(nod,nodf),c23(nodf,nod),uvel(nod),vvel(nod),                 &
            row1(1,nod),row2(1,nod),funnyf(nodf,1),rowf(1,nodf))              
      read(10,*) width , depth
      uvel =.0; vvel =.0 ; kay=0.0; kay(1,1)=visc/rho; kay(2,2)=visc/rho      
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)  ;  call sample(element,points,weights)
!------- loop the elements to find nband and set up global arrays------------
  nband = 0                                                                     
 elements_1: do iel = 1 , nels
             call geometry_8qxv(iel,nxe,width,depth,coord,num)
             inc=0
             do i=1,8;inc=inc+1;g(inc)=nf(1,num(i));end do
             do i=1,7,2;inc=inc+1;g(inc)=nf(2,num(i));end do
             do i=1,8;inc=inc+1;g(inc)=nf(3,num(i));end do
             g_num(:,iel )=num; g_coord(:,num)=transpose(coord); g_g(:,iel)=g
             if(nband<bandwidth(g))nband=bandwidth(g)
 end do elements_1    
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do
  write(11,'(2(a,i5))')                                                       &
          "There are ",neq,"  equations and the half-bandwidth is   ",nband
  allocate(pb(neq,2*(nband+1)-1),loads(0:neq),oldlds(0:neq),work(nband+1,neq))
    loads = .0  ; oldlds =.0   ; iters = 0
    read(10,*) fixed_nodes  
       allocate(node(fixed_nodes),sense(fixed_nodes),val(fixed_nodes),        &
                no(fixed_nodes))
       read(10,*) (node(i),sense(i),val(i),i=1,fixed_nodes )
!-------------------iteration loop  -----------------------------------------
  iterations: do
               iters = iters + 1   ; converged = .false.
     pb = .0; work = .0; ke = .0
!------------ element stiffness integration and assembly---------------------

      elements_2:  do iel = 1 , nels 
               num = g_num(: , iel ); coord=transpose(g_coord(:,num)) 
               g = g_g( : , iel )   ; coordf = coord(1 : 7 : 2, : )
               uvel = (loads(g(1:nod))+oldlds(g(1:nod)))*.5
               do i = nod + nodf + 1 , ntot
                  vvel(i-nod-nodf) = (loads(g(i))+oldlds(g(i)))*.5
               end do
               c11 = .0; c12 = .0; c21 = .0; c23 = .0; c32 = .0 
           gauss_points_1: do i = 1 , nip
!--------------------- velocity contribution ----------------------------------
              call shape_fun(fun,points,i) ;funny(:,1) = fun
              ubar = dot_product(fun,uvel);vbar = dot_product(fun,vvel)
              if(iters==1) then; ubar = 1.; vbar = 0.; end if
              call shape_der(der,points,i);  jac = matmul(der,coord) 
              det = determinant(jac )     ; call invert(jac)
              deriv = matmul(jac,der);row1(1,:)=deriv(1,:);row2(1,:)=deriv(2,:) 
              c11 = c11 + matmul(matmul(transpose(deriv),kay),deriv) &
                     *det* weights(i) + &
                          matmul(funny,row1)*det*weights(i)*ubar + &
                          matmul(funny,row2)*det*weights(i)*vbar
!----------------------now the pressure contribution--------------------------
               call shape_fun(funf,points,i); funnyf(:,1)=funf
               call shape_der(derf,points,i)  ;jac=matmul(derf,coordf) 
               det=determinant(jac)      ;     call invert(jac)
               derivf=matmul(jac,derf)
               rowf(1,:) = derivf(1,:)
               c12 = c12 + matmul(funny,rowf)*det*weights(i)/rho
               rowf(1,:) = derivf(2,:)
               c32 = c32 + matmul(funny,rowf)*det*weights(i)/rho
               c21 = c21 + matmul(funnyf,row1)*det*weights(i)
               c23 = c23 + matmul(funnyf,row2)*det*weights(i)               
           end do gauss_points_1
         call formupv(ke,c11,c12,c21,c23,c32) ; call formtb(pb,ke,g)
      end do elements_2
!----------- prescribed values of velocity and pressure ----------------------
       loads = .0      
       do i=1, fixed_nodes; no(i) = nf(sense(i),node(i))  ; end do
         pb( no ,nband+1)=pb( no ,nband+1) + 1.e20
         loads(no) = pb(no,nband+1) * val 
!------------------------ solve the simultaneous equations -------------------
    call gauss_band(pb,work); call solve_band(pb,work,loads); loads(0) = .0   
    call checon(loads,oldlds,tol,converged);if(converged.or.iters==limit) exit
  end do iterations 
    write(11,'(a,i5,a)')"The solution took",iters,"  iterations to converge"
        write(11,'(a)') " The nodal velocities and porepressures are    :"
        write(11,'(a)')"   Node   u - velocity   pressure    v - velocity"
            do k=1,nn; write(11,'(i5,a,3e12.4)')k,"    ",loads(nf(:,k));end do  
 end program p90 
