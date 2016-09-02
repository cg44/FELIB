 program p52       
!------------------------------------------------------------------------
!      program 5.2 plane strain of an elastic solid using uniform 4-node
!      quadrilateral elements numbered in the y direction.
!      Analytical forms of km and bee matrices
!------------------------------------------------------------------------
 use new_library  ;  use  geometry_lib ;  use vlib  ;  implicit none
 integer::nels,nye,neq,nn,nr,nip,nodof=2,nod=4,nst=3,ndof,loaded_nodes,    &
          i,k,iel,ndim=2
 real:: e,v,det,aa,bb      ;   character(len=15) :: element='quadrilateral'   
!----------------------------- dynamic arrays---------------------------------
 real  , allocatable  :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:),   &
                         jac(:,:),der(:,:),deriv(:,:),weights(:),          &
                         bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:)
 integer, allocatable :: nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),kdiag(:)
!---------------------------input and initialisation---------------------------
   open (10,file='p52.dat',status=    'old',action='read')
   open (11,file='p52.res',status='replace',action='write')                   
  read (10,*) nels,nye,nn,nip,aa,bb,e,v 
  ndof=nod*nodof 
  allocate (nf(nodof,nn),points(nip,ndim),g(ndof),g_coord(ndim,nn),        & 
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),g_g(ndof,nels),    &
            weights(nip),der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),      &
            km(ndof,ndof),eld(ndof),sigma(nst),num(nod),g_num(nod,nels))
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)    ;  allocate (kdiag(neq))
  call deemat(dee,e,v)   ;  call sample(element,points,weights)               
!-------- loop the elements to set up global geometry and kdiag -------------
  kdiag=0
  elements_1  : do iel =1,nels
                  call geometry_4qy(iel,nye,aa,bb,coord,num)
                  call num_to_g( num , nf, g);  g_num(:,iel)=num
                  g_coord(:,num)=transpose(coord)
                  g_g(:,iel)=g   ;    call fkdiag(kdiag,g) 
  end do elements_1               
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                  &
                              "Element ",k,"        ",g_num(:,k); end do  
  kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do
  write(11,'(2(a,i5))')                                                        &
           "There are",neq,"  equations and the skyline storage is ",kdiag(neq)
  allocate(kv(kdiag(neq)),loads(0:neq)); kv=0.0   
!---------------- element stiffness integration and assembly-------------------
 elements_2: do iel = 1 , nels
               num= g_num(: ,iel); coord =transpose(g_coord(:,num));g=g_g(:,iel)
               call analy4(km,coord,e,v) ;  call fsparv (kv,km,g,kdiag)
 end do elements_2
 loads=0.0 ; read (10,*) loaded_nodes,(k,loads(nf(:,k)), i =1,loaded_nodes)
!---------------------------equation solution--------------------------------
    call sparin(kv,kdiag) ;call spabac(kv,loads,kdiag)
    write(11,'(a)') "The nodal displacements Are :"
    write(11,'(a)') "Node         Displacement"
    do k=1,nn; write(11,'(i5,a,2e12.4)') k,"   ",loads(nf(:,k)); end do
!-------------------recover stresses at element Gauss-points-----------------
 elements_3:do iel = 1 , nels
               num = g_num(:,iel);  coord =transpose( g_coord(: ,num)) 
               g = g_g(:,iel)     ;    eld=loads(g)
              write(11,'(a,i5,a)')                                             &
                       "The Gauss Point stresses for element",iel,"   are :" 
            integrating_pts_2: do i = 1 , nip
                 call bee4(coord,points,i,det,bee)
                 sigma = matmul (dee,matmul(bee,eld)) 
                 write(11,'(a,i5)') "Point",i    ;   write(11,'(3e12.4)') sigma
            end do integrating_pts_2 
 end do elements_3
end  program p52
