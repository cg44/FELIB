 program p114       
!-----------------------------------------------------------------------------
!      program 11.4 forced vibration of a rectangular elastic
!      solid in plane strain using uniform 8-node quadrilateral elements
!      numbered in the y direction - lumped or consistent mass
!      implicit integration by theta method : pcg version
!------------------------------------------------------------------------------
 use new_library;     use  geometry_lib   ;      implicit none
 integer::nels,nye,neq,nn,nr,nip,nodof=2,nod=8,nst=3,ndof,                    &
          i,k,iel,ndim=2,nstep,npri,iters,limit
 real::aa,bb,e,v,det,rho,alpha,beta,omega,theta,period,pi,dtim,area,          &
       c1,c2,c3,c4,time,tol,big,up; character(len=15)::element='quadrilateral'
 logical :: consistent = .false.  , converged
!----------------------------- dynamic arrays---------------------------------
 real    ,allocatable :: loads(:),points(:,:),dee(:,:),coord(:,:),temp(:,:),  &
                         fun(:),jac(:,:), der(:,:),deriv(:,:), weights(:),    &
                         bee(:,:),km(:,:),g_coord(:,:),x1(:),d1x1(:),d2x1(:), &
                         emm(:,:),ecm(:,:),x0(:),d1x0(:),d2x0(:),             &
                         store_km(:,:,:),store_mm(:,:,:),u(:),p(:),d(:),      &
                         x(:),xnew(:),pmul(:),utemp(:),diag_precon(:)
 integer, allocatable :: nf(:,:), g(:) , num(:)  , g_num(:,:) , g_g(:,:)        
!------------------------input and initialisation------------------------------
  open (10,file='p114.dat',status=    'old',action='read')
  open (11,file='p114.res',status='replace',action='write')  
  read (10,*) nels,nye,nn,nip,aa,bb,rho,e,v,                                  &
              alpha,beta,nstep,npri,theta,omega,tol,limit 
  ndof=nod*nodof         
  allocate ( nf(nodof,nn), points(nip,ndim),g(ndof), g_coord(ndim,nn),        &
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),weights(nip),         &
            der(ndim,nod), deriv(ndim,nod), bee(nst,ndof), km(ndof,ndof),     &
            num(nod),g_num(nod,nels),g_g(ndof,nels),emm(ndof,ndof),           &
            ecm(ndof,ndof),fun(nod),store_km(ndof,ndof,nels),utemp(ndof),     &
            pmul(ndof),store_mm(ndof,ndof,nels),temp(ndof,ndof))  
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf (nf);neq=maxval(nf)                                              
  pi=acos(-1.)   ; period = 2.*pi/omega ; dtim =period/20.                     
  c1=(1.-theta)*dtim; c2=beta-c1; c3=alpha+1./(theta*dtim); c4=beta+theta*dtim
  call deemat (dee,e,v); call sample(element,points,weights)
!-------------- loop the elements to find neq and store globals ---------------
  elements_1: do iel = 1 , nels
              call geometry_8qy(iel,nye,aa,bb,coord,num)
              call num_to_g ( num , nf , g ) ; g_num(:,iel)=num
              g_coord(:,num)=transpose(coord);g_g(:,iel) = g
  end do elements_1                                                            
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do       
    write(11,'(a,i5,a)') "There are ",neq,"  equations to be solved"
  allocate(x0(0:neq),d1x0(0:neq),x1(0:neq),diag_precon(0:neq),u(0:neq),       &
           d2x0(0:neq),loads(0:neq),d1x1(0:neq),d2x1(0:neq),                  &
           d(0:neq),p(0:neq),x(0:neq),xnew(0:neq))  
    xnew=.0; p=.0; diag_precon=.0  ;    store_km = .0; store_mm = .0      
!------ element stiffness and mass integration ,storage and preconditioner ---- 
 elements_2: do iel = 1 , nels
             num = g_num( : , iel ); coord = transpose(g_coord(: , num )) 
             g = g_g( : , iel )    ; km=0.0   ; area = .0 ; emm = .0 
          gauss_points_1: do i = 1 , nip     
               call shape_der (der,points,i) ; jac = matmul(der,coord) 
               det = determinant(jac)  ; call invert(jac)
               deriv = matmul(jac,der) ; call beemat (bee,deriv) 
             km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               area = area + det*weights(i); call shape_fun(fun,points,i)
               if(consistent) then
                call ecmat(ecm,fun,ndof,nodof); ecm=ecm*det*weights(i)*rho
                emm = emm + ecm    
              end if
          end do gauss_points_1                  
       if(.not.consistent) then
         do i=1,ndof; emm(i,i)=area*rho*.2 ; end do
         do i=1,13,4 ; emm(i,i)=emm(3,3)*.25; end do
         do i=2,14,4 ; emm(i,i)=emm(3,3)*.25 ; end do
       end if
   store_km (: , : , iel ) = km ; store_mm( : , : , iel ) = emm
  do k=1,ndof;diag_precon(g(k))=diag_precon(g(k))+emm(k,k)*c3+km(k,k)*c4;end do
 end do elements_2 
 diag_precon(1:neq) = 1. / diag_precon(1:neq); diag_precon(0) = .0 
!-----------------------initial conditions -----------------------------------
            x0 = .0; d1x0 = .0; d2x0 = .0                                     
!----------------------- time stepping loop ----------------------------------
   time = .0
   write(11,'(a)') "    Time t  cos(omega*t) Displacement Iterations"
  timesteps: do i = 1 , nstep
              time = time + dtim   ; loads = .0                  
    u = .0
    elements_3 : do iel = 1 , nels    ! gather for rhs multiply
                    g = g_g( : , iel ) ; km = store_km(:, :, iel)
                    emm = store_mm(:,:,iel) ; pmul = x0(g); temp=km*c2+emm*c3 
                    utemp = matmul(temp , pmul)
                    u ( g ) = u ( g ) + utemp   ! scatter
                    pmul = d1x0(g) ! velocity bit
                    temp=emm/theta  ; utemp=matmul(temp,pmul);u(g)=u(g)+utemp
    end do elements_3                                                         
    loads(neq)=theta*dtim*cos(omega*time)+c1*cos(omega*(time-dtim))
    u(0) = .0; loads = u +   loads
!------------------ solve simultaneous equations by pcg ----------------------
       d = diag_precon*loads; p = d; x = .0
       iters = 0
     iterations  :      do 
             iters = iters + 1     ;    u = 0.
       elements_4 : do iel = 1, nels
                      g = g_g( : , iel ) ; km = store_km(:,:,iel)
                      emm=store_mm(:,:,iel); temp=emm*c3+km*c4 ;  pmul = p(g)
                      utemp = matmul(temp,pmul); u(g) = u(g)+ utemp 
       end do elements_4    ; u(0) = .0
!---------------------------pcg equation solution------------------------------
           up=dot_product(loads,d); alpha= up/ dot_product(p,u)
           xnew = x + p* alpha ; loads=loads - u*alpha;  d = diag_precon*loads
           beta=dot_product(loads,d)/up; p=d+p*beta
           big = .0; converged = .true.
           u = xnew ; where(u<.0) u = -u; big = maxval(u)
           u = (xnew - x)/big ;where(u<.0) u=-u ; big = maxval(u)
           if(big>tol) converged=.false.   ;     x=xnew    
           if(converged .or. iters==limit) exit
     end do iterations
               x1=xnew 
               d1x1=(x1-x0)/(theta*dtim)-d1x0*(1.-theta)/theta
               d2x1=(d1x1-d1x0)/(theta*dtim)-d2x0*(1.-theta)/theta
              if(i/npri*npri==i) then
                write(11,'(3e12.4,i10)')time,cos(omega*time),x1(neq),iters
              end if 
                x0 = x1; d1x0 = d1x1; d2x0 = d2x1                              
  end do timesteps
end program p114
