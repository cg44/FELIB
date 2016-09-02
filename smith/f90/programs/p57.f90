program p57
!----------------------------------------------------------------------------
! program 5.7 three-dimensional elastic analysis using 14-node brick elements  
!----------------------------------------------------------------------------
   use new_library  ;        use   geometry_lib   ;  implicit  none
   integer ::nels,neq,nn,nr,nip,nodof=3,nod=14,nst=6,ndof,fixed_nodes,        &
             iel,i,k,ii,jj,kk,ll,ndim=3 
   real    ::e,v,det  ; character(len=15) :: element = 'hexahedron'      
!----------------------  dynamic  arrays  -------------------------------------
   real    , allocatable :: dee(:,:),points(:,:),weights(:),                  &
                            coord(:,:),jac(:,:),der(:,:),deriv(:,:),          &
                            bee(:,:),km(:,:),eld(:),eps(:),sigma(:),          &
                            kv(:),loads(:),g_coord(:,:), value(:)
   integer  , allocatable ::g(:),nf(:,:),kdiag(:),num(:),g_num(:,:),g_g(:,:), &
                            no(:),sense(:),node(:)
! ---------------------   input and initialisation    -------------------------
   open(10,file='p57.dat',status=    'old',action='read')
   open(11,file='p57.res',status='replace',action='write') 
   read(10,*) nels,nn,nip,e,v         ;          ndof = nod * nodof
   allocate( nf(nodof,nn),dee(nst,nst),coord(nod,ndim),num(nod),              &
             jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g(ndof),            &
             bee(nst,ndof),km(ndof,ndof),eld(ndof),sigma(nst),eps(nst),       &
             g_g(ndof,nels),g_coord(ndim,nn),g_num(nod,nels))
     nf = 1; read(10,*) nr; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
     call formnf(nf); neq = maxval(nf)     ; call deemat(dee,e,v)
    allocate(loads(0:neq), kdiag(neq)) ; loads = .0 ; kdiag = 0
    read(10,*) g_num ; read(10,*) g_coord(:,1:16)
    do i=17,nn
    read(10,*)ii,jj,kk,ll
          g_coord(:,i)=.25*(g_coord(:,ii)+g_coord(:,jj)+g_coord(:,kk)+        &
          g_coord(:,ll))
    end do   
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,3e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,14i4)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do  
! --------  loop the elements to set up global g and find kdiag ---------------
    elements_1 :  do iel = 1 , nels
                    num = g_num(:,iel)  ; call num_to_g (num , nf , g )
                    call fkdiag(kdiag,g);   g_g( : , iel ) = g
    end do elements_1
    kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do
    write(11,'(2(a,i5))')                                                      &
         "There are",neq,"  equations and the skyline storage is :",kdiag(neq)
    allocate( kv(kdiag(neq))) ;  kv = .0
! -----------------  element stiffness integration and assembly ---------------
  allocate(weights(nip),points(nip,ndim));call sample(element, points, weights) 
 elements_2: do iel = 1 , nels       
             num = g_num(:,iel) ; g = g_g(:,iel)
             coord = transpose(g_coord(:,num)) ;        km=0.0      
    gauss_pts_1:  do i =1 , nip
               call shape_der(der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac);   call invert (jac)
               deriv = matmul(jac,der); call beemat (bee,deriv) 
               km=km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
    end do gauss_pts_1 
   call fsparv (kv,km,g,kdiag)
 end do elements_2    
 read(10,*) fixed_nodes
 if(fixed_nodes/=0) then
   allocate(no(fixed_nodes),node(fixed_nodes),                               &
            sense(fixed_nodes),value(fixed_nodes))
   read(10,*)(node(i),sense(i),value(i),i=1,fixed_nodes)
   do i=1,fixed_nodes; no(i) = nf(sense(i),node(i)); end do
   kv(kdiag(no)) = kv(kdiag(no)) + 1.e20; loads(no)=kv(kdiag(no)) * value
 end if  
!------------------------------equation solution-------------------------------
    call sparin(kv,kdiag) ;call spabac(kv,loads,kdiag)
    write(11,'(a)') "The nodal displacements are"
    do k=1,nn; write(11,'(i5,a,3e12.4)') k,"    ",loads(nf(:,k)); end do
!-------------------recover stresses at element centroids---------------------- 
 nip = 1; deallocate(points,weights); allocate(points(nip,ndim),weights(nip))
 elements_3 : do iel = 1 , nels
                  num = g_num(:,iel); coord = transpose(g_coord( : , num ))
                  g=g_g(:,iel)       ;        eld = loads( g )
                  write(11,'(a,i5,a)')                                        &
                           "The centroid stresses for element",iel,"  are"
    gauss_pts_2: do i = 1 , nip     
                    call shape_der(der,points,i);jac = matmul(der,coord)
                    call invert(jac);   deriv= matmul(jac,der)
                    call beemat(bee,deriv);sigma = matmul(dee,matmul(bee,eld))
                    write(11,'(a,i5)') "Point",i   ; write(11,'(6e12.4)') sigma
   end do gauss_pts_2 
 end do elements_3
end program p57 
