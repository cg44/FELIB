program p43
!------------------------------------------------------------------------------
! program 4.3 beam on an elastic foundation
!             numerically integrated beam and foundation stiffness
!------------------------------------------------------------------------------
use new_library  ;  use  geometry_lib   ;   implicit  none
integer::nels,neq,nn,nband,nr,nod=2,nodof=2,ndof=4,iel,i,k,l,ndim=1,       &
         loaded_nodes,fixed_nodes,nip,np_types
real::fs,fs0,fs1,x,samp_pt ; character(len=15) :: element = 'line'    
!------------------------------dynamic arrays----------------------------------
real,allocatable::km(:,:),mm(:,:),eld(:),kv(:),loads(:),coord(:,:),        &
                  action(:),g_coord(:,:),value(:),ftf(:,:),dtd(:,:),       &
                  der2(:),fun(:),mom(:),store_km(:,:,:),points(:,:),       &
                  weights(:),ei(:),ell(:)
integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),no(:),g_g(:,:),        &
                       node(:),sense(:),etype(:)   
!------------------------input and initialisation------------------------------
open (10 , file = 'p43.dat' , status = 'old' ,    action ='read')
open (11 , file = 'p43.res' , status = 'replace', action='write')    
read(10,*)nels,nip,np_types; nn=nels+1
allocate(nf(nodof,nn),km(ndof,ndof),coord(nod,ndim),g_coord(ndim,nn),      &
         eld(ndof),action(ndof),g_num(nod,nels),num(nod),g(ndof),          &
         g_g(ndof,nels),mm(ndof,ndof),ftf(ndof,ndof),ei(np_types),         &
         ell(nels),dtd(ndof,ndof),store_km(ndof,ndof,nels),der2(ndof),     &
         fun(ndof),mom(nn),points(nip,ndim),weights(nip),etype(nels)) 
read(10,*)fs0,fs1,ei; etype=1; if(np_types>1)read(10,*)etype
read(10,*)ell,nr
nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf)  
!--------------loop the elements to find global array sizes--------------------
nband=0
elements_1: do iel=1,nels    
              call geometry_2l(iel,ell(iel),coord,num);call num_to_g(num,nf,g) 
              g_num(:,iel)=num; g_coord(:,num)=transpose(coord) 
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
!--------numerical integration of beam and foundation stiffness----------------
!----------------global stiffness matrix assembly------------------------------
call sample(element,points,weights)
x=0.0
elements_2: do iel=1,nels
    km=0.0; mm=0.0; g=g_g(:,iel)
    integrating_pts: do i=1,nip
        samp_pt=x+ell(iel)*0.5*(points(i,1)+1.0)
        fs=samp_pt/(ell(iel)*nels)*(fs1-fs0)+fs0
        call fmbeam(der2,fun,points,i,ell(iel))
        do k=1,ndof; do l=1,ndof
          ftf(k,l)=fun(k)*fun(l)*weights(i)*0.5*ell(iel)*fs
          dtd(k,l)=der2(k)*der2(l)*weights(i)*8.0*ei(etype(iel))/(ell(iel)**3)
        end do; end do
        mm=mm+ftf; km=km+dtd
    end do integrating_pts
    km=km+mm; store_km(:,:,iel)=km(:,:); x=x+ell(iel)
    call formkv(kv,km,g,neq)
end do elements_2      
!-----------------------------read loads---------------------------------------
loads=0.0; read(10,*)loaded_nodes
if(loaded_nodes/=0)read(10,*)(k,loads(nf(:,k)),i=1,loaded_nodes)
read (10,*)fixed_nodes
if(fixed_nodes /=0)then
  allocate( node(fixed_nodes),no(fixed_nodes),&
  sense(fixed_nodes),value(fixed_nodes))
  read(10,*) (node(i),sense(i),value(i),i=1,fixed_nodes)
  do i=1,fixed_nodes; no(i)=nf(sense(i),node(i)); end do
  kv(no)=kv(no)+1.e20; loads(no)=kv(no)*value 
end if    
!-----------------------------equation solution -------------------------------
call banred(kv,neq); call bacsub(kv,loads)    
!-----------------------retrieve element end actions---------------------------
elements_3: do iel=1,nels
              km(:,:)=store_km(:,:,iel) 
              g=g_g(:,iel); eld=loads( g ); action=matmul(km,eld)   
              mom(iel)=action(2)
              if(iel==nels)mom(iel+1)=-action(4)
end do elements_3                 
write(11,'(a)')"   Node   Displacement Moment "
nodes: do i=1,nn
         write(11,'(i5,a,2e12.4)')i,"   ",loads(2*i-1),mom(i)
end do nodes    
end program p43


