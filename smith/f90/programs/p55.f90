 program p55       
!------------------------------------------------------------------------
!      program 5.5 non-axisymmetric strain of an axisymmetric elastic solid
!      using uniform 8-node quadrilateral elements numbered in the x direction
!------------------------------------------------------------------------
 use new_library   ;   use geometry_lib      ;     implicit none
 integer::nels,nre,nde,neq,nn,nr,nip,nodof=3,nod=8,nst=6,ndof,loaded_nodes, &
          i,k,iel,ndim=2,lth,iflag
 real:: e,v,det,aa,bb,chi,pi,ca,sa,radius 
 character (len=15) :: element = 'quadrilateral'  
!----------------------------- dynamic arrays---------------------------------
 real  , allocatable :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:),     &
                        fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),     &
                        bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:),      &
                        value(:)
integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),kdiag(:),&
                        node(:),no(:),sense(:) 
!------------------------input and initialisation------------------------------
 open (10,file='p55.dat',status=    'old',action='read')
 open (11,file='p55.res',status='replace',action='write')                     
  read(10,*) nels,nre,nde,nn,nip,aa,bb,e,v
  read(10,*) lth,iflag,chi                ;    ndof=nod*nodof 
  allocate (nf(nodof,nn),points(nip,ndim),g(ndof),g_coord(ndim,nn),         &
            dee(nst,nst),coord(nod,ndim),fun(nod),jac(ndim,ndim),           &
            weights(nip),der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),       &
            km(ndof,ndof),eld(ndof),sigma(nst),num(nod),g_num(nod,nels),    &
            g_g(ndof,nels))
 nf=1; read(10,*) nr ;if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
       call formnf(nf); neq=maxval(nf) ;  allocate (kdiag(neq))
       call deemat(dee,e,v)    ; call sample( element, points, weights)
  pi=acos(-1.); chi=chi*pi/180.; ca=cos(chi); sa=sin(chi)
!-------- loop the elements to set up global geometry and kdiag -----------
  kdiag=0
  elements_1   : do iel =1,nels
                   call geometry_8qx(iel,nre,aa,bb,coord,num)
                   call num_to_g(num,nf,g);   g_num(:,iel)=num  
                   g_coord(:,num)=transpose(coord)
                   g_g( : , iel ) = g       ;    call fkdiag(kdiag,g) 
  end do elements_1     
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                              "Element ",k,"        ",g_num(:,k); end do  
        kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do
     write(11,'(2(a,i5))')                                                     &
          "There are",neq,"  equations and the skyline storage is",kdiag(neq)
        allocate(kv(kdiag(neq)),loads(0:neq)); kv=0.0
!--------------- element stiffness integration and assembly--------------------
 elements_2: do iel = 1 , nels
               num= g_num(:, iel);  g = g_g( : , iel )
               coord = transpose(g_coord(:,num)) ;  km = .0
          integrating_pts_1: do i = 1 , nip
               call shape_fun(fun,points,i); call shape_der(der,points,i)  
               jac = matmul(der,coord); det= determinant(jac) 
               call invert(jac); deriv= matmul(jac,der)  
               call bmat_nonaxi(bee,radius,coord,deriv,fun,iflag,lth)
               det=det*radius
               km= km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
          end do integrating_pts_1
          call fsparv (kv,km,g,kdiag)
 end do elements_2                                                            
loads=0.0;read(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
!------------------------equation solution--------------------------------
    call sparin(kv,kdiag) ;call spabac(kv,loads,kdiag)
    write(11,'(a)') "The Nodal Displacements Are :"
    write(11,'(a)') "Node         Displacement"
    do k=1,nn; write(11,'(i5,a,3e12.4)') k,"   ",loads(nf(:,k)); end do        
!-------------------recover stresses at  element centroids ----------------
 i = 1; points = .0
 elements_3:do iel = 1 , nels
             num = g_num(:, iel);  coord = transpose(g_coord(:,num)) 
             g = g_g(: ,iel )  ;    eld=loads(g)
    write(11,'(a,i5,a)') "The centre point stresses for element",iel,"  are :" 
        call shape_fun(fun,points,i); call shape_der(der,points,i)  
        jac = matmul(der,coord); call invert(jac); deriv= matmul(jac,der)  
        call bmat_nonaxi(bee,radius,coord,deriv,fun,iflag,lth)
        bee(1:4,:)=bee(1:4,:)*ca; bee(5:6,:)=bee(5:6,:)*sa
        sigma = matmul(dee,matmul(bee,eld)) 
       write(11,'(a,i5)') "Point",i    ;   write(11,'(6e12.4)') sigma
 end do elements_3
end program p55 
