 program p73
!------------------------------------------------------------------------
!      program 7.3 general program for two- or three-dimensional
!                  analysis of Laplace's equation
!------------------------------------------------------------------------
 use new_library   ;    use  geometry_lib    ;  implicit none
 integer::nels,neq,nband,nn,nr,nip,nodof,nod,ndof,i,k,iel,ndim,           &
          loaded_nodes,fixed_nodes,np_types         
 real::det    ;     character (len=15) :: element    
!----------------------------- dynamic arrays----------------------------------
 real,allocatable::kv(:),kvh(:),loads(:),disps(:),points(:,:),            &
                   coord(:,:),jac(:,:),der(:,:),deriv(:,:),weights(:),    &
                   prop(:,:),kp(:,:),g_coord(:,:),value(:),kay(:,:)  
 integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),g_g( :, :),no(:),    &
                      node(:),etype(:)          
!-----------------------input and initialisation------------------------------
 open (10, file = 'p73.dat' , status = 'old' , action = 'read')
 open (11, file = 'p73.res' , status='replace',action = 'write')                                                
 read(10,*)element,nels,nn,nip,nodof,nod,ndim,np_types; ndof=nod*nodof
 allocate(nf(nodof,nn),points(nip,ndim),g_coord(ndim,nn),coord(nod,ndim), &
             etype(nels),jac(ndim,ndim),weights(nip),num(nod),            &
             g_num(nod,nels),der(ndim,nod),deriv(ndim,nod),kp(ndof,ndof), &
             g(ndof),g_g(ndof,nels),kay(ndim,ndim),prop(ndim,np_types))
 read(10,*)prop
 etype=1; if(np_types>1)read(10,*)etype
 read(10,*)g_coord; read(10,*)g_num     
 nf=1; read(10,*)nr; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
 call formnf(nf); neq=maxval(nf); call sample(element,points,weights)   
!------------- loop the elements to find nband and store steering vectors ---
 nband=0
 elements_1: do iel =1,nels
         num=g_num(:,iel) ; call num_to_g(num,nf,g); g_g(:,iel)=g   
         if(nband<bandwidth(g))nband=bandwidth(g) 
 end do elements_1
 write(11,'(a)')"Global coordinates"
 do k=1,nn 
   write(11,'(a,i5,a,3e12.4)')"Node    ",k,"    ",g_coord(:,k); end do
 write(11,'(a)')"Global node numbers"
 do k=1,nels 
   write(11,'(a,i5,a,20i4)')"Element ",k,"      ",g_num(:,k); end do
 write(11,'(2(a,i5),/)')                                                   &
   "There are ",neq,"  equations and the half bandwidth is ",nband
 allocate( kv(neq*(nband+1)),kvh(neq*(nband+1)),loads(0:neq),disps(0:neq)) 
 kv=0.0; loads =0.0 
!------------- element stiffness integration and assembly-------------------
 elements_2: do iel=1,nels
            kay=0.0; do i=1,ndim; kay(i,i)=prop(i,etype(iel)); end do 
            num=g_num(:,iel); coord=transpose(g_coord(:,num))
            g=g_g(:,iel); kp=0.0   
            integrating_pts_1: do i=1,nip
                  call shape_der(der,points,i); jac=matmul(der,coord)
                  det=determinant(jac); call invert(jac)
                  deriv=matmul(jac,der)
                  kp= kp+matmul(matmul(transpose(deriv),kay),deriv)       &
                      *det*weights(i)
            end do integrating_pts_1    
            call formkv(kv,kp,g,neq)  
 end do elements_2
 kvh=kv     
 read(10,*)loaded_nodes
 if(loaded_nodes/=0)read(10,*)(k,loads(nf(:,k)),i=1,loaded_nodes)   
 read(10,*)fixed_nodes 
 if(fixed_nodes/=0)then
        allocate(node(fixed_nodes),no(fixed_nodes),value(fixed_nodes))
        read(10,*)(node(i),value(i),i=1,fixed_nodes)
        do i=1,fixed_nodes; no(i)=nf(1,node(i)); end do
        kv(no)=kv(no)+1.e20; loads(no)=kv(no)*value
 end if
!------------------------equation solution------------------------------------   
 call banred(kv,neq); call bacsub(kv,loads)

!------------------------retrieve flow rates-----------------------------------
 call linmul(kvh,loads,disps)
 write(11,'(a)')"The nodal values are:"
 write(11,'(a)')"          Potentials  Flow rates"
 do k=1,nn
    write(11,'(i5,a,2f12.2)')k,"   ",loads(nf(1,k)),disps(nf(1,k)); end do
 write(11,'(a)')"      Inflow      Outflow"
 write(11,'(2f12.2)')sum(disps,mask=disps>0.),sum(disps,mask=disps<0.)
end program p73


