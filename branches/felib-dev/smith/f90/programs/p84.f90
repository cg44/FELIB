    program p84      
!------------------------------------------------------------------------------
!      program 8.4 conduction equation on rectangular 
!      area using 4-node quadrilateral elements : pcg version
!      implicit integration in time using 'theta' method
!------------------------------------------------------------------------------
 use new_library    ;  use geometry_lib  ;      implicit none
 integer::nels,nxe,neq,nn,nr,nip,nodof=1,nod=4,ndof,ndim=2,                  &
          i,j,k,iel,nstep,npri,nres,iters,limit
 real::aa,bb,permx,permy,det,theta,dtim,val0,time,tol,alpha,beta,up,big
 logical::converged      ; character(len=15) :: element='quadrilateral'
!------------------------- dynamic arrays--------------------------------------
 real ,allocatable :: loads(:),u(:),p(:),points(:,:),kay(:,:),coord(:,:),    &
                      fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),d(:),   &
                      kp(:,:), pm(:,:), funny(:,:),g_coord(:,:),             &
                      storka(:,:,:),storkb(:,:,:),x(:),xnew(:),pmul(:),      &
                      utemp(:),diag_precon(:)
 integer, allocatable :: nf(:,:), g(:) , num(:) , g_num(:,:) ,g_g(:,:)         
!----------------------input and initialisation--------------------------------
  open (10,file='p84.dat',status=    'old',action='read')
  open (11,file='p84.res',status='replace',action='write')                
  read (10,*) nels,nxe,nn,nip,aa,bb,permx,permy ,                            &
              dtim,nstep,theta,npri,nres,tol,limit 
  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),weights(nip),kay(ndim,ndim),     &
            coord(nod,ndim), fun(nod), jac(ndim,ndim),g_coord(ndim,nn),      &
            der(ndim,nod), deriv(ndim,nod), pm(ndof,ndof),g_num(nod,nels),   &
            kp(ndof,ndof), g(ndof),funny(1,nod),num(nod),g_g(ndof,nels),     &
            storka(ndof,ndof,nels),storkb(ndof,ndof,nels),utemp(ndof),       &
            pmul(ndof))
  kay=0.0 ; kay(1,1)=permx; kay(2,2)=permy                            
  call sample (element,points,weights)
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)                                             
!-------------loop the elements to  set up global arrays --------------------- 
    elements_1: do iel = 1 , nels
               call geometry_4qx(iel,nxe,aa,bb,coord,num)
               g_num(:,iel) = num; g_coord(: , num ) = transpose(coord)  
               call num_to_g (num,nf,g);    g_g( : , iel ) = g  
    end do elements_1     
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)')  "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do  
    allocate(loads(0:neq),diag_precon(0:neq),u(0:neq),d(0:neq),p(0:neq),      &
             x(0:neq),xnew(0:neq)) ; storka = .0; storkb = .0   
      write(11,'(a,i5,a)') "There are ",neq,"  equations to be solved"
      p = .0; diag_precon = .0; xnew = .0
!----------- element integration ,storage and build preconditioner ------------
    elements_2: do iel = 1 , nels
             num = g_num(:,iel) ; coord = transpose( g_coord( : , num )) 
             g = g_g( : , iel )     ;     kp=0.0 ; pm=0.0                
       gauss_pts:  do i =1 , nip
               call shape_der (der,points,i) ; call shape_fun(fun,points,i)
               funny(1,:)=fun(:) ; jac = matmul(der,coord)
               det=determinant(jac); call invert(jac); deriv = matmul(jac,der) 
               kp = kp + matmul(matmul(transpose(deriv),kay),deriv) &
                    *det*weights(i)
               pm  =  pm + matmul( transpose(funny),funny)*det*weights(i) 
       end do gauss_pts
       storka(:,:,iel)=pm+kp*theta*dtim;storkb(:,:,iel)=pm-kp*(1.-theta)*dtim
     do k=1,ndof; diag_precon(g(k))=diag_precon(g(k))+ storka(k,k,iel); end do
    end do elements_2
    diag_precon(1:neq) = 1./diag_precon(1:neq) ; diag_precon(0) = .0
!---------------------------initial conditions --------------------------------
      read(10,*) val0; loads=val0   ; loads(0) = .0
!----------------------time stepping recursion --------------------------------
   write(11,'(a,i5,a)') "  Time   Pressure at node",nres," Iterations"                   
  timesteps: do j=1,nstep
               time=j*dtim                                                
    u = .0
    elements_3 : do iel = 1 , nels    ! gather for rhs multiply
                    g = g_g( : , iel ) ; kp = storkb(:, :, iel)
                    pmul = loads(g) ; utemp = matmul(kp , pmul)
                    u ( g ) = u ( g ) + utemp   ! scatter
    end do elements_3
    u(0) = .0 ; loads = u
!------------------- solve simultaneous equations by pcg ----------------------
       d = diag_precon*loads; p = d; x = .0
       iters = 0
     iterations  :      do 
             iters = iters + 1     ;    u = 0.
       elements_4 : do iel = 1, nels
                      g = g_g( : , iel ) ; kp = storka(:,:,iel);  pmul = p(g)
                      utemp = matmul(kp,pmul); u(g) = u(g)+ utemp 
       end do elements_4    ; u(0) = .0
!--------------------------pcg equation solution-------------------------------
           up=dot_product(loads,d); alpha= up/ dot_product(p,u)
           xnew = x + p* alpha ; loads=loads - u*alpha;  d = diag_precon*loads
           beta=dot_product(loads,d)/up; p=d+p*beta
           big = .0; converged = .true.
           u = xnew ; where(u<.0) u = -u; big = maxval(u)
           u = (xnew - x)/big ;where(u<.0) u=-u ; big = maxval(u)
           if(big>tol) converged=.false.   ;     x=xnew    
           if(converged .or. iters==limit) exit
     end do iterations
               loads=xnew
               if(j/npri*npri==j) then 
                write(11,'(2e12.4,7x,i5)')time,loads(nf(:,nres)),iters
               end if
  end do timesteps
end program p84
