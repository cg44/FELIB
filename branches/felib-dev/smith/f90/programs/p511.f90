program p511         
!------------------------------------------------------------------------------
!      program 5.11 three dimensional analysis of an elastic
!      solid using 20-node brick elements
!      preconditioned conjugate gradient solver  ;  only integrate one element
!      diagonal preconditioner diag_precon   ;  vectorised version
!------------------------------------------------------------------------------
 use new_library  ;  use  geometry_lib    ;     implicit none
 integer::nxe,nze,neq,nn,nr,nip,nodof=3,nod=20,nst=6,ndof,loaded_nodes,     &
          i,k,ndim=3,iters,limit,iel,nels   
 real::aa,bb,cc,e,v,det,tol,up,alpha,beta,big
 logical :: converged  ; character(len=15) :: element = 'hexahedron'      
!--------------------------- dynamic arrays------------------------------------
 real    ,allocatable :: points(:,:),dee(:,:),coord(:,:), weights(:),        &
                         g_coord(:,:), jac(:,:), der(:,:), deriv(:,:),       &
                         bee(:,:), km(:,:),eld(:),eps(:),sigma(:),           &
                         diag_precon(:),p(:),r(:),x(:),xnew(:),              &
                         u(:),g_pmul(:,:),g_utemp(:,:),d(:)
 integer, allocatable :: nf(:,:), g(:), num(:), g_num(:,:) ,g_g( : , :)      
!--------------------------input and initialisation----------------------------
  open (10,file='p511.dat',status=    'old',action='read')
  open (11,file='p511.res',status='replace',action='write')                   
  read (10,*) nels,nxe,nze,nn,nip,aa,bb,cc,e,v,   tol,limit ;  ndof=nod*nodof   
  allocate ( nf(nodof,nn), points(nip,ndim),dee(nst,nst),coord(nod,ndim),    &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),                    &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),       &  
            g(ndof),g_pmul(ndof,nels),g_utemp(ndof,nels), g_coord(ndim,nn),  &
            g_num(nod,nels),weights(nip),num(nod),g_g(ndof,nels))      
   nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
         call formnf(nf);neq=maxval(nf)   
  allocate(p(0:neq),r(0:neq),x(0:neq),xnew(0:neq),u(0:neq),&
           diag_precon(0:neq),d(0:neq))   
       r=0.; p=0.; x=0.; xnew=0.  ; diag_precon=0.         
   call deemat(dee,e,v);   call sample(element,points,weights) 
! ----------------- single element stiffness integration ----------------------
             iel=1
             call geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
             km=0.0                                
      gauss_pts_1:  do i=1,nip
                call shape_der (der,points,i) ; jac = matmul(der,coord)
                det = determinant(jac)        ; call invert(jac)
                deriv = matmul(jac,der) ;call beemat (bee,deriv)  
                km=km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
      end do gauss_pts_1
! -------------- store global arrays and build the preconditioner -------------
     elements_1: do iel = 1,nels
                  call geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
                  g_num(:, iel) = num; g_coord(: ,num) = transpose(coord) 
                  call num_to_g(num,nf,g);  g_g( : , iel) = g
               do k=1,ndof;diag_precon(g(k))=diag_precon(g(k))+km(k,k);end do 
     end do elements_1   
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,3e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,27i3)')                                 &
                              "Element ",k,"    ",g_num(:,k); end do  
    write(11,'(a,i5)') "The number of equations is  ",neq                  
!--------------------invert the preconditioner and get starting r--------------
         read(10,*) loaded_nodes,(k,r(nf(:,k)),i=1,loaded_nodes)
         write(11,'(a,e12.4)') "The total load is", sum(r)  
        diag_precon(1:neq)=1./ diag_precon(1:neq)  ; diag_precon(0) = .0
                 d=diag_precon*r  ; p = d   
!----------------------preconditioned c. g. iterations-------------------------
       iters = 0
     iterations  :      do 
             iters = iters + 1     ;    u = 0.      
       elements_2 : do iel = 1 , nels                ! gather
                      g_pmul(: , iel) = p( g_g( : , iel))  
       end do elements_2
!------------------------- global matrix multiply ----------------------------
          g_utemp = matmul( km , g_pmul )
!dir$ ivdep
       elements_2a : do iel = 1 , nels                ! scatter
                      u(g_g(:,iel))=u(g_g(:,iel))+g_utemp(:,iel)
       end do elements_2a     ! let this vectorise by compiler directive
!----------------------------pcg equation solution-----------------------------
    up=dot_product(r,d); alpha= up/ dot_product(p,u)
    xnew = x + p* alpha ; r=r - u*alpha;  d = diag_precon*r
    beta=dot_product(r,d)/up; p=d+p*beta  ; p(0) = .0 
    big=0.; converged = .true.
    do i=1,neq; if(abs(xnew(i))>big) big= abs(xnew(i)) ; end do
    do i=1,neq; if(abs(xnew(i)-x(i))/big>tol)converged=.false.; end do; x=xnew
    if(converged .or. iters==limit) exit
                  end do iterations
       write(11,'(a,i5)') "The number of iterations to convergence was  ",iters 
       write(11,'(a)')    "The nodal displacements are   :"
   do k=1,22; write(11,'(i5,a,3e12.4)') k,"    ",xnew(nf(:,k)); end do
!-------------------recover stresses at centroidal gauss-point-----------------
  nip=1; deallocate(points,weights); allocate(points(nip,ndim),weights(nip))
  elements_3:do iel = 1, nels
                 num = g_num(: ,iel)  ; coord =transpose( g_coord(:,num))
                 g = g_g( : , iel) ; eld=xnew(g)
                 write(11,'(a,i5,a)')                                          &
                       "The Gauss point stresses for element",iel,"  are :"    
     gauss_pts_2: do i= 1 , nip 
       call shape_der(der,points,i); jac= matmul(der,coord)
       call invert (jac);   deriv= matmul(jac,der)
       bee= 0.;call beemat(bee,deriv); sigma = matmul (dee,matmul(bee,eld))
       write(11,'(a,i5)') "Point  ",i   ; write(11,'(6e12.4)') sigma
     end do gauss_pts_2 
  end do elements_3
 end program p511
