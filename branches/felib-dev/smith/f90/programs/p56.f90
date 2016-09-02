    program p56     
!------------------------------------------------------------------------
!      program 5.6 three dimensional analysis of an elastic
!      solid using 4-node tetrahedral elements
!------------------------------------------------------------------------
 use new_library      ;  use geometry_lib  ;     implicit none
 integer::nels,neq,nn,nr,nip,nodof=3,nod=4,nst=6,ndof,loaded_nodes,          &
          i,k,iel,ndim=3
 real:: e,v,det    ; character(len=15) :: element = 'tetrahedron'               
!--------------------------- dynamic arrays------------------------------------
 real    ,allocatable :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:),     &
                         jac(:,:),weights(:), der(:,:), deriv(:,:),bee(:,:), &
                         km(:,:),eld(:),sigma(:),g_coord(:,:)
 integer, allocatable :: nf(:,:), g(:), kdiag(:) ,num(:) ,g_num(:,:),g_g(:,:)  
!--------------------------input and initialisation----------------------------
  open (10,file='p56.dat',status=    'old',action='read')
  open (11,file='p56.res',status='replace',action='write')                    
  read (10,*) nels,nn,nip,e,v        ;     ndof=nod*nodof  
  allocate ( nf(nodof,nn), points(nip,ndim),dee(nst,nst),coord(nod,ndim),    &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g(ndof),            &
            bee(nst,ndof), km(ndof,ndof),eld(ndof),sigma(nst),g_g(ndof,nels),&
            g_coord(ndim,nn),g_num(nod,nels),weights(nip),num(nod))
  read (10, *) g_coord ; read (10, *)  g_num
  nf=1; read(10,*) nr ;if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
        call formnf(nf); neq=maxval(nf) ; allocate ( loads(0:neq),kdiag(neq) ) 
  call deemat (dee,e,v); call sample(element,points,weights) 
  kdiag=0      
!  ------------- loop the elements to set up  g_g  and  kdiag  ----------------
  elements_1  :  do iel = 1 , nels  
     num = g_num(:,iel); call num_to_g(num,nf,g)      
     g_g(:,iel) = g ;  call fkdiag(kdiag,g)
  end do elements_1
  kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do
    allocate(kv(kdiag(neq)))  ; kv=0.0
   write(11,'(a)') "Global Coordinates"
   do k=1,nn;write(11,'(a,i5,a,3e12.4)')"Node",k,"      ",g_coord(:,k);end do
   write(11,'(a)') "Global Node Numbers"
   do k=1,nels; write(11,'(a,i5,a,4i5)')                                      &
                         "Element",k,"        ",g_num(:,k); end do
   write(11,'(2(a,i5))')                                                      &
        "There are ",neq,"  equations and the skyline storage is",kdiag(neq)
    loads=0.0 ; read (10,*) loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
   write(11,'(a,e12.4)') " The total load is ", sum(loads)            
!--------------- element stiffness integration and assembly------------------- 
 elements_2: do iel = 1 , nels
             num = g_num(:,iel) ; g = g_g(:,iel)
             coord = transpose(g_coord(:,num)) ;        km=0.0      
    gauss_pts_1:  do i =1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac);   call invert (jac)
               deriv = matmul(jac,der); call beemat (bee,deriv) 
               km=km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
    end do gauss_pts_1 
    call fsparv (kv,km,g,kdiag)
 end do elements_2                                             
!---------------------------equation solution----------------------------------
    call sparin(kv,kdiag) ;call spabac(kv,loads,kdiag)
    write(11,'(a)') "The nodal displacements are   :"
    write(11,'(a)') "  Node            Displacement"
    do k=1,nn; write(11,'(i5,a,3e12.4)')k,"      ",loads(nf(:,k)); end do
!----------------------recover stresses at element Gauss-points---------------- 
 elements_3 : do iel = 1 , nels
                  num = g_num(:,iel); coord = transpose(g_coord( : , num ))
                  g=g_g(:,iel)       ;        eld = loads( g )
                  write(11,'(a,i5,a)')                                        &
                       "The Gauss point stresses for element",iel," are"
    gauss_pts_2: do i = 1,nip     
                    call shape_der(der,points,i);jac = matmul(der,coord)
                    call invert(jac);   deriv= matmul(jac,der)
                    call beemat(bee,deriv);sigma = matmul(dee,matmul(bee,eld))
                    write(11,'(a,i5)') "Point  ",i ;write(11,'(6e12.4)')sigma
   end do gauss_pts_2 
 end do elements_3
end program p56
