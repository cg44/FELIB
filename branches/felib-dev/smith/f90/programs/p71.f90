program p71
!------------------------------------------------------------------------
!      program 7.1 solution of Laplace's equation
!      for plane free-surface flow using 4-node quadrilaterals
!------------------------------------------------------------------------
 use new_library   ;   use geometry_lib ; use vlib   ;  implicit none 
 integer::nels,nxe,neq,nn,nr,nodof=1,nod=4,ndof=4,np_types,i,k,iel,ndim=2,&
          fixed_up,fixed_down,fixed_seep,iters,limit,nband
 real::tol,upstream,downstream
 logical::converged
!--------------------------- dynamic arrays------------------------------------
 real,allocatable::kv(:),kvh(:),loads(:),coord(:,:),kp(:,:),g_coord(:,:), &
                   kay(:,:),disps(:),oldpot(:),width(:),surf(:),prop(:,:),&
                   angs(:)
 integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),node_up(:), &
                      no_up(:),etype(:),node_down(:),no_down(:),          &
                      node_seep(:),no_seep(:)
!-----------------------input and initialisation------------------------------
 open (10,file='p71.dat',status='old'    , action ='read' )
 open (11,file='p71.res',status='replace', action = 'write')
 read(10,*) nels,nxe,tol,limit,np_types; nn=(nxe+1)*(nels/nxe+1)
 allocate(nf(nodof,nn),g(ndof),g_coord(ndim,nn),coord(nod,ndim),          &
          width(nxe+1),surf(nxe+1),angs(nxe+1),kp(ndof,ndof),num(nod),    &
          g_num(nod,nels),prop(ndim,np_types),g_g(ndof,nels),             &
          kay(ndim,ndim),etype(nels))
 read(10,*)prop; etype=1; if(np_types>1)read(10,*)etype
 read(10,*)width; read(10,*)angs; read(10,*)surf
 read(10,*)nr  
 nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf)
 read(10,*)upstream,fixed_up ; allocate (node_up(fixed_up),no_up(fixed_up))
 read(10,*) node_up          ; read (10,*)downstream,fixed_down 
 allocate (node_down(fixed_down),no_down(fixed_down))
 read(10,*)node_down      ;  fixed_seep=nels/nxe-fixed_down
 allocate (node_seep(fixed_seep),no_seep(fixed_seep))
 do i=1,fixed_seep; node_seep(i)=i*(nxe+1)+1; end do      
!--------- loop the elements to find nband and set up global arrays ------
 nband=0
 elements_1: do iel=1,nels    
      call geometry_freesurf(iel,nxe,fixed_seep,fixed_down,downstream,&
                             width,angs,surf,coord,num)  
      call num_to_g ( num , nf , g )
      g_coord(:,num)=transpose(coord); g_num(:,iel)=num
      if(nband<bandwidth(g))nband=bandwidth(g)
 end do elements_1    
 write(11,'(a)')"Initial global coordinates"
 do k=1,nn 
   write(11,'(a,i5,a,3e12.4)')"Node    ",k,"    ",g_coord(:,k); end do  
 write(11,'(2(a,i5),/)')                                                   &
     "There are ",neq,"  equations and the half bandwidth is ",nband  
 allocate(kv(neq*(nband+1)),kvh(neq*(nband+1)),loads(0:neq),disps(0:neq), &
          oldpot(0:neq)); oldpot=0.0
!--------------- element conductivity integration and assembly-----------------
 iters=0
 iterations: do
    iters=iters+1; kv=0.0
    elements_2: do iel=1,nels
         kay=0.0; do i=1,ndim; kay(i,i)=prop(i,etype(iel)); end do
         call geometry_freesurf(iel,nxe,fixed_seep,fixed_down,downstream, &
                                width,angs,surf,coord,num)  
         call num_to_g ( num , nf , g )
         g_coord(:,num)=transpose(coord)
         call seep4(coord,kay,kp)  ;   call formkv(kv,kp,g,neq)
    end do elements_2
    kvh=kv
!----------   specify fixed potentials and factorise equations     
    loads=0.0 
    do i=1,fixed_up; no_up(i)=nf(1,node_up(i)); end do
    kv(no_up)=kv(no_up)+1.e20; loads(no_up)=kv(no_up)*upstream 
    do i=1,fixed_down; no_down(i)=nf(1,node_down(i)); end do
    kv(no_down)=kv(no_down)+1.e20; loads(no_down)=kv(no_down)*downstream 
    do i=1,fixed_seep 
      no_seep(i)=nf(1,node_seep(i))
      kv(no_seep(i))=kv(no_seep(i))+1.e20  
      loads(no_seep(i))=kv(no_seep(i))*                                   &
      (downstream+(surf(1)-downstream)*(fixed_seep+1-i)/(fixed_seep+1))
    end do     
!------------------------------equation solution-------------------------------
    call banred(kv,neq);call bacsub(kv,loads)
    surf(1:nxe)=loads(1:nxe)    
!-------------------------------check convergence-----------------------------
    call checon(loads,oldpot,tol,converged)
    if(converged.or.iters==limit)exit
 end do iterations
! -------------------- write out the results --------------------------------
 write(11,'(a)')"Final global coordinates for plotting"
 do k=1,nn
    write(11,'(a,i5,a,3e12.4)')"Node    ",k,"    ",g_coord(:,k); end do
 write(11,'(a)')"Global node numbers"
 do k=1,nels 
   write(11,'(a,i5,a,27i3)')"Element ",k,"      ",g_num(:,k); end do
 call linmul(kvh,loads,disps)
 write(11,'(a)')"The nodal values are:"
 write(11,'(a)')"          Potentials  Flow rate"
 do k=1,nn
   write(11,'(i5,a,2e12.4)')k,"   ",loads(nf(1,k)),disps(nf(1,k)); end do
 write(11,'(a)')"  Inflow      Outflow"
 write(11,'(2e12.4)')sum(disps,mask=disps<0.),sum(disps,mask=disps>0.)
 write(11,'(a,i5)')"Number of iterations =",iters   
end program p71


