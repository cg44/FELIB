 program p46
!------------------------------------------------------------------------------
! program 4.6 stability analysis of rigid-jointed frames using 
!             2-d beam/rod elements
!------------------------------------------------------------------------------
use new_library  ;  use  geometry_lib  ;  use vlib   ;implicit none
integer::nels,neq,nn,nband,nr,nod=2,nodof=3,ndof=6,iel,i,k,ndim=2,         &
         loaded_nodes,iy,iters,limit,ksc,incs,nprops,np_types
real::total_load,tol,det  ;  logical::converged     
!-------------------------dynamic arrays---------------------------------------
real,allocatable::km(:,:),kp(:,:),eld(:),kv(:),loads(:),coord(:,:),        &
                    action(:),g_coord(:,:),val(:,:),prop(:,:),             &
                    axif(:),axip(:),eldtot(:),dload(:),gamma(:),           &
                    oldsps(:),kcop(:),disps(:),local(:)
integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),no(:),g_g(:,:),etype(:) 
!-----------------------input and initialisation-------------------------------
open (10 , file = 'p46.dat' , status = 'old' ,    action ='read')
open (11 , file = 'p46.res' , status = 'replace', action='write')             
read(10,*)nels,nn,nprops,np_types,limit,tol
allocate(nf(nodof,nn),km(ndof,ndof),coord(nod,ndim),g_coord(ndim,nn),      &
          eld(ndof),action(ndof),g_num(nod,nels),num(nod),g(ndof),         &
          g_g(ndof,nels),prop(nprops,np_types),                            &
          axif(nels),axip(nels),local(ndof),kp(ndof,ndof),etype(nels)) 
read(10,*)prop; etype=1; if(np_types>1)read(10,*)etype
read(10,*)g_coord; read(10,*)g_num
read(10,*)nr
nf = 1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf) 
!---------------loop the elements to find global array sizes-------------------
nband=0
elements_1: do iel=1,nels    
              num=g_num(:,iel)  ;  call num_to_g (num , nf , g)
              g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g)
end do elements_1
allocate(kv(neq*(nband+1)),kcop(neq*(nband+1)),loads(0:neq),disps(0:neq),  &
            eldtot(0:neq),oldsps(0:neq))    
write(11,'(a)')"Global coordinates"
do k=1,nn; write(11,'(a,i5,a,3e12.4)')                                      &
    "Node    ",k,"    ",g_coord(:,k); end do
write(11,'(a)')"Global node numbers"
do k=1,nels; write(11,'(a,i5,a,27i3)')                                      &
    "Element ",k,"      ",g_num(:,k); end do
write(11,'(2(a,i5),/)')                                                     &
    "There are ",neq,"  equations and the half-bandwidth is ",nband   
axif=0.0; axip=0.0; eldtot=0.0
read(10,*)loaded_nodes; allocate(no(loaded_nodes),val(loaded_nodes,nodof))
read(10,*)(no(i),val(i,:),i=1,loaded_nodes)
read(10,*)incs; allocate(dload(incs)) ;    read(10,*)dload    
!----------------------load increment loop-------------------------------------
total_load=0.0 
load_increments: do iy=1,incs
   total_load=total_load+dload(iy)
   loads=0.0
   do i = 1,loaded_nodes; loads(nf(:,no(i)))=dload(iy)*val(i,:); end do
   oldsps=0.0; iters=0
   iterations: do 
      iters=iters+1; kv = 0.0    
!-------------------global stiffness matrix assembly---------------------------
      elements_2 : do iel = 1, nels
         num=g_num(:,iel)  ; coord=transpose(g_coord(:,num))
         call rigid_jointed(km,prop,gamma,etype,iel,coord) 
         call beam_kp(kp,coord,axif(iel))
         km=km+kp; g=g_g(:,iel)
         call formkv(kv,km,g,neq)
      end do elements_2
      kcop=kv; call kvdet(kcop,neq,nband,det,ksc); disps=loads   
!-------------------------equation solution -----------------------------------
      call banred(kv, neq); call bacsub(kv, disps)   
!---------------------------check convergence----------------------------------
      call checon(disps,oldsps,tol,converged)
      elements_3: do iel=1,nels
         num=g_num(:,iel); coord=transpose(g_coord(:,num)) 
         g=g_g(:,iel); eld=disps(g)
         call rigid_jointed(km,prop,gamma,etype,iel,coord) 
         call beam_kp(kp,coord,axif(iel))
         km=km+kp   ;            action=matmul(km,eld)   
         call glob_to_loc(local,action,gamma(iel),coord)
         axif(iel)=axip(iel)+local(4)
      end do elements_3                 
      if(iters==limit .or. converged)exit iterations
   end do iterations   
!-------------at convergence update displacements and axial forces-------------
   axip=axif; eldtot=eldtot+disps
   write(11,'(a,e12.4)')"Load factor   ",total_load
   write(11,'(a)')"The nodal displacements are:"
   do i=1,loaded_nodes
     write(11,'(i5,6e12.4)')no(i),eldtot(nf(:,no(i)))
   end do
   write(11,'(a,e12.4)')"The determinant is  ",det
   write(11,'(a,i5,a,/)')"Converged in  ",iters,"  iterations"
   if(iters == limit)exit 
end do load_increments     
end program p46


