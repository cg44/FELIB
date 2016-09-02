    program p58     
!------------------------------------------------------------------------------
!      program 5.8 three dimensional analysis of an elastic
!      solid using uniform 20-node hexahedral brick elements
!------------------------------------------------------------------------------
 use new_library      ;  use geometry_lib  ;     implicit none
 integer::nels,nxe,nze,neq,nn,nr,nip,nodof=3,nod=20,nst=6,ndof,loaded_nodes, &
          i,k,iel,ndim=3
 real::aa,bb,cc,e,v,det    ;  character(len=15) :: element = 'hexahedron'      
!------------------------- dynamic arrays--------------------------------------
 real    ,allocatable :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:),     &
                         jac(:,:),weights(:), der(:,:), deriv(:,:),bee(:,:), &
                         km(:,:),eld(:),sigma(:),g_coord(:,:)
 integer, allocatable :: nf(:,:), g(:), kdiag(:) ,num(:) ,g_num(:,:),g_g(:,:)  
!-----------------------input and initialisation-------------------------------
  open (10,file='p58.dat',status=    'old',action='read')
  open (11,file='p58.res',status='replace',action='write')                    
  read (10,*) nels,nxe,nze,nn,nip,aa,bb,cc,e,v        ;  ndof=nod*nodof  
  allocate ( nf(nodof,nn), points(nip,ndim),dee(nst,nst),coord(nod,ndim),    &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g(ndof),            &
            bee(nst,ndof), km(ndof,ndof),eld(ndof),sigma(nst),g_g(ndof,nels),&
            g_coord(ndim,nn),g_num(nod,nels),weights(nip),num(nod))
  nf=1; read(10,*) nr ;if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
        call formnf(nf); neq=maxval(nf)
  allocate  ( loads(0:neq),kdiag(neq)  ) 
  call deemat (dee,e,v); call sample(element,points,weights) 
  kdiag=0      
!--------- loop the elements to set up global geometry and kdiag  -------------
  elements_1  :  do iel = 1 , nels
     call geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
     call num_to_g(num,nf,g);  call fkdiag(kdiag,g);g_num(:,iel)=num
     g_coord(:,num)=transpose(coord); g_g(:,iel)=g
  end do elements_1
  kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do
    allocate(kv(kdiag(neq)))  ; kv=0.0
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,3e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,20i3)')                                 &
                              "Element ",k,"    ",g_num(:,k); end do  
    write(11,'(2(a,i5))')                                                      &
          "There are",neq,"  equations and the skyline storage is",kdiag(neq)
    loads=0.0 ; read (10,*) loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
    write(11,'(a,e12.4)') "The total load is ",sum(loads)            
!--------------- element stiffness integration and assembly--------------------
 elements_2: do iel = 1 , nels
             num = g_num(:,iel) ; g = g_g(:,iel)
             coord = transpose(g_coord(:,num)) ;        km=0.0      
    gauss_pts_1:  do i =1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac) ;   call invert (jac)
               deriv = matmul(jac,der); call beemat (bee,deriv) 
               km=km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
    end do gauss_pts_1 
    call fsparv (kv,km,g,kdiag)
 end do elements_2                                             
!---------------------------equation solution----------------------------------
    call sparin(kv,kdiag) ;call spabac(kv,loads,kdiag)
    write(11,'(a)') "The nodal displacements are"
    do k=1,nn; write(11,'(i5,a,3e12.4)') k,"    ",loads(nf(:,k)); end do
!-------------------recover stresses at element centroids---------------------- 
  nip = 1; deallocate(points,weights); allocate(points(nip,ndim),weights(nip))  
 elements_3 : do iel = 1 , nels
                  num = g_num(:,iel); coord = transpose(g_coord( : , num ))
                  g=g_g(:,iel)       ;        eld = loads( g )
                  write(11,'(a,i5,a)')                                         &
                           "The centroid stresses for element",iel,"  are"
    gauss_pts_2: do i = 1 , nip     
                    call shape_der(der,points,i);jac = matmul(der,coord)
                    call invert(jac);   deriv= matmul(jac,der)
                    call beemat(bee,deriv);sigma = matmul(dee,matmul(bee,eld))
                    write(11,'(a,i5)') "Point",i   ; write(11,'(6e12.4)') sigma
   end do gauss_pts_2 
 end do elements_3
end program p58
