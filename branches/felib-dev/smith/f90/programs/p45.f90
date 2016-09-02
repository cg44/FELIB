program p45  
!------------------------------------------------------------------------------
! program 4.5 elasto-plastic analysis of rigid-jointed frames 
! using beam and beam/rod elements in 1-, 2- or 3-dimensions
!------------------------------------------------------------------------------
use new_library  ;  use  geometry_lib  ;   use  vlib    ;   implicit none  
integer::nels,neq,nn,nband,nr,nod=2,nodof,ndof,iel,i,k,ndim,loaded_nodes,  &
         incs,limit,iters,iy,nprops,np_types
real::tol,total_load      ;logical::converged  
!----------------------------dynamic arrays------------------------------------
real,allocatable::km(:,:),eld(:),kv(:),loads(:),coord(:,:),action(:),      &
                  g_coord(:,:),gamma(:),prop(:,:),bdylds(:),eldtot(:),     &
                  holdr(:,:),oldsps(:),react(:),val(:,:),dload(:)
integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),no(:),etype(:) 
!-----------------------input and initialisation-------------------------------
open (10 , file = 'p45.dat' , status = 'old' ,    action ='read')
open (11 , file = 'p45.res' , status = 'replace', action='write')              
read(10,*)nels,nn,ndim,nprops,np_types,limit,tol
select case(ndim) 
  case(1); nodof=2; case(2); nodof=3; case(3); nodof=6 
  case default; write(11,'(a)')"Wrong number of dimensions input"
end select
ndof=nod*nodof
allocate(nf(nodof,nn),km(ndof,ndof),coord(nod,ndim),g_coord(ndim,nn),      &
         eld(ndof),action(ndof),g_num(nod,nels),num(nod),g(ndof),          &
         gamma(nels),g_g(ndof,nels),holdr(ndof,nels),react(ndof),          &
         prop(nprops,np_types),etype(nels)) 
read(10,*)prop; etype=1; if(np_types>1) read(10,*)etype
if(ndim==3)read(10,*)gamma
read(10,*)g_coord; read(10,*)g_num
read(10,*)nr
nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf)  
!----------------loop the elements to find global array sizes------------------
nband=0
elements_1: do iel=1,nels    
              num=g_num(:,iel) ; call num_to_g (num , nf , g )
              g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g)
end do elements_1
allocate(kv(neq*(nband+1)),loads(0:neq),eldtot(0:neq),bdylds(0:neq),       &
           oldsps(0:neq)); kv=0.0; holdr=0.0
write(11,'(a)')"Global coordinates"
do k=1,nn; write(11,'(a,i5,a,3e12.4)')                                      & 
    "Node    ",k,"    ",g_coord(:,k); end do
write(11,'(a)')"Global node numbers"
do k=1,nels; write(11,'(a,i5,a,27i3)')                                      &
    "Element ",k,"      ",g_num(:,k); end do
write(11,'(2(a,i5),/)')                                                     &
    "There are ",neq,"  equations and the half-bandwidth is ",nband    
!--------------------global stiffness matrix assembly--------------------------
elements_2: do iel=1,nels
              num=g_num(:,iel); coord=transpose(g_coord(:,num))
              call rigid_jointed(km,prop,gamma,etype,iel,coord); g=g_g(:,iel)
              call formkv(kv,km,g,neq)
end do elements_2   
read(10,*)loaded_nodes; allocate(no(loaded_nodes),val(loaded_nodes,nodof))
read(10,*)(no(i),val(i,:),i=1,loaded_nodes)
read(10,*)incs; allocate(dload(incs))  ;   read(10,*)dload    
!-------------------------equation factorisation-------------------------------
call banred(kv, neq)
!--------------------------load increment loop---------------------------------
total_load=0.0 
load_increment: do iy=1,incs
   total_load=total_load+dload(iy)
   oldsps=0.0; iters=0
   iterations: do 
     iters=iters+1; loads=0.0
     do i=1,loaded_nodes
        loads(nf(:,no(i)))=dload(iy)*val(i,:)
     end do
     loads=loads+bdylds; bdylds=0.0    
!------------------------forward/back-substitution-----------------------------
     call bacsub(kv, loads)  
!-----------------------------check convergence--------------------------------
     call checon(loads,oldsps,tol,converged) 
!---------------------inspect moments in all elements--------------------------
     elements_3: do iel=1,nels
        num=g_num(:,iel); coord=transpose(g_coord(:,num)) 
        g=g_g(:,iel); eld=loads(g)
        call rigid_jointed(km,prop,gamma,etype,iel,coord) 
        action=matmul(km,eld); react=0.0   
!-------------if plastic moments exceeded generate correction vector-----------
        if(limit/=1)then
           call hinge(coord,holdr,action,react,prop,iel,etype,gamma)
           bdylds(g)=bdylds(g)-react; bdylds(0)=0.0
        end if   
!----------------at convergence update element reactions-----------------------
      if(iters==limit.or.converged)holdr(:,iel)=holdr(:,iel)+react(:)+action(:)
     end do elements_3                 
     if(iters==limit .or. converged)exit iterations
   end do iterations
   eldtot=loads+eldtot
   write(11,'(a,e12.4)')"Load factor   ",total_load
   write(11,'(a)')"The nodal displacements are:"
   do i=1,loaded_nodes
     write(11,'(i5,6e12.4)')no(i),eldtot(nf(:,no(i)))
   end do
   write(11,'(a,i5,a,/)')"Converged in  ",iters,"  iterations"
   if(iters==limit .and. limit/=1)exit load_increment
end do load_increment    
end program p45


