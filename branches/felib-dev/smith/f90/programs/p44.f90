program p44
!------------------------------------------------------------------------------
! program 4.4 equilibrium of rigid-jointed frames using beam/rod elements
!             in 1-, 2- or 3-dimensions
!------------------------------------------------------------------------------
use new_library  ;  use  geometry_lib   ;  use  vlib  ;   implicit none
integer::nels,neq,nn,nband,nr,nod=2,nodof,ndof,iel,i,k,ndim,loaded_nodes,  &
         fixed_nodes,nprops,np_types   
!----------------------------dynamic arrays------------------------------------
real,allocatable::km(:,:),eld(:),kv(:),loads(:),coord(:,:),                &
                  action(:),g_coord(:,:),value(:),prop(:,:),gamma(:) 
integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),no(:),g_g(:,:),        &
                     node(:),sense(:),etype(:)  
!-------------------------input and initialisation-----------------------------
open (10 , file = 'p44.dat' , status = 'old' ,    action ='read')
open (11 , file = 'p44.res' , status = 'replace', action='write')             
read(10,*)nels,nn,ndim,nprops,np_types
select case(ndim) 
  case(1); nodof=2; case(2); nodof=3; case(3); nodof=6 
  case default; write(11,'(a)')"Wrong number of dimensions input"
end select
ndof=nod*nodof
allocate(nf(nodof,nn),km(ndof,ndof),coord(nod,ndim),g_coord(ndim,nn),      &
         eld(ndof),action(ndof),g_num(nod,nels),num(nod),g(ndof),          &
         gamma(nels),g_g(ndof,nels),prop(nprops,np_types),etype(nels)) 
read(10,*)prop; etype=1; if(np_types>1)read(10,*)etype
if(ndim == 3)read(10,*)gamma
read(10,*)g_coord; read(10,*)g_num
read(10,*)nr 
nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf) 
!---------------loop the elements to find global array sizes------------------ 
nband=0
elements_1: do iel=1,nels    
              num=g_num(:,iel); call num_to_g ( num ,nf , g  )  
              g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g)
end do elements_1
allocate(kv(neq*(nband+1)),loads(0:neq)); kv=0.0  
write(11,'(a)')"Global coordinates"
do k=1,nn; write(11,'(a,i5,a,3e12.4)')                                     &
    "Node    ",k,"    ",g_coord(:,k); end do
write(11,'(a)')"Global node numbers"
do k=1,nels; write(11,'(a,i5,a,27i3)')                                     &
    "Element ",k,"      ",g_num(:,k); end do
write(11,'(2(a,i5),/)')                                                    &
    "There are ",neq,"  equations and the half-bandwidth is ",nband    
!---------------------global stiffness matrix assembly-------------------------
elements_2: do iel=1,nels
              num=g_num(:,iel); coord=transpose(g_coord(:,num))
              call rigid_jointed(km,prop,gamma,etype,iel,coord)
              g=g_g(:,iel)   ;     call formkv(kv,km,g,neq)
end do elements_2        
!--------------------read loads and/or displacements---------------------------
loads=0.0; read(10,*)loaded_nodes
if(loaded_nodes/=0)read(10,*)(k,loads(nf(:,k)),i =1,loaded_nodes)
read(10,*)fixed_nodes
if(fixed_nodes/=0)then
  allocate(node(fixed_nodes),no(fixed_nodes),                              &
           sense(fixed_nodes),value(fixed_nodes))
  read(10,*)(node(i),sense(i),value(i),i=1,fixed_nodes)
  do i=1,fixed_nodes; no(i)=nf(sense(i),node(i)); end do
  kv(no)=kv(no)+1.e20; loads(no)=kv(no)*value 
end if     
!-----------------------------equation solution -------------------------------
call banred(kv,neq); call bacsub(kv,loads)
write(11,'(a)')"The nodal displacements are:"
do k=1,nn; write(11,'(i5,6e12.4)')k,loads(nf(:,k)); end do   
!-----------------retrieve element end actions---------------------------------
write(11,'(a)')"The element 'actions' are:"
elements_3: do iel=1,nels
              num=g_num(:,iel); coord=transpose(g_coord(:,num)) 
              g=g_g(:,iel); eld=loads(g)
              call rigid_jointed(km,prop,gamma,etype,iel,coord) 
              action=matmul(km,eld)   
              if(ndim<3)then
                write(11,'(i5,6e12.4)')iel,action
              else
                write(11,'(i5,6e12.4)')iel,   action(1: 6:1)
                write(11,'(a,6e12.4)')"     ",action(7:12:1)      
              end if  
end do elements_3    
end program p44


