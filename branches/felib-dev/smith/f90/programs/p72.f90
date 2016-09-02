 program p72       
!------------------------------------------------------------------------
!      program 7.2 solution of Laplace's equation
!      over an axisymmetric region using 4-node quadrilaterals
!------------------------------------------------------------------------
 use new_library   ;    use   geometry_lib ;   implicit none
 integer::nels,nxe,neq,nn,nr,nip,nodof=1,nod=4,ndof=4,lOaded_nodes,i,k,   &
               iel,ndim=2,fixed_nodes
 real::det,aa,bb,radius  ;  character (len=15) :: element='quadrilateral'   
!--------------------------- dynamic arrays---------------------------------
 real,allocatable::kv(:),kvh(:),loads(:),points(:,:),coord(:,:),jac(:,:), &
                   der(:,:),deriv(:,:),weights(:),kp(:,:),g_coord(:,:),   &
                   value(:),kay(:,:),disps(:),perms(:),fun(:)
 integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),kdiag(:),   &
                      node(:),no(:)
!------------------------input and initialisation--------------------------
 open(10,file='p72.dat',status='old',    action='read')
 open(11,file='p72.res',status='replace',action='write')
 read (10,*)nels,nxe,nip,aa,bb; nn=(nxe+1)*(nels/nxe+1)
 allocate(nf(nodof,nn),points(nip,ndim),g(ndof),g_coord(ndim,nn),         &
          coord(nod,ndim),jac(ndim,ndim),fun(ndof),weights(nip),          &
          der(ndim,nod),deriv(ndim,nod),kp(ndof,ndof),num(nod),           &
          g_num(nod,nels),g_g(ndof,nels),kay(ndim,ndim),perms(ndim))
 read(10,*)perms
 kay=0.0; do i=1,ndim; kay(i,i)=perms(i); end do
 read(10,*)nr  
 nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf)
 allocate (kdiag(neq))
!------- loop the elements to set up global geometry and kdiag ----------------
 kdiag=0
 elements_1: do iel=1,nels
               call geometry_4qx(iel,nxe,aa,bb,coord,num)
               g_num(:,iel)=num; g_coord(:,num)=transpose(coord)
               call num_to_g(num,nf,g);   g_g(:,iel)=g  ; call fkdiag(kdiag,g) 
 end do elements_1
 write(11,'(a)')"Global coordinates"
 do k=1,nn 
   write(11,'(a,i5,a,3e12.4)')"Node    ",k,"    ",g_coord(:,k); end do
 write(11,'(a)')"Global node numbers"
 do k=1,nels 
   write(11,'(a,i5,a,27i3)')"Element ",k,"      ",g_num(:,k); end do    
 kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do
 write(11,'(2(a,i5),/)')                                                   &
   "There are ",neq,"  equations and the skyline storage is ",kdiag(neq)
 allocate(kv(kdiag(neq)),kvh(kdiag(neq)),loads(0:neq),disps(0:neq)) 
 kv=0.0; loads=0.0   ; call sample(element,points,weights)
!--------------- element conductivity integration and assembly----------------
 elements_2: do iel=1,nels
        kp=0.0
        num=g_num(:,iel); coord=transpose(g_coord(:,num)); g=g_g(:,iel)
        integrating_pts_1: do i=1,nip
           call shape_der(der,points,i); jac=matmul(der,coord) 
           det=determinant(jac); call invert(jac)
           deriv=matmul(jac,der); call shape_fun(fun,points,i)
           radius=sum(fun(:)*coord(:,1))
           kp=kp+matmul(matmul(transpose(deriv),kay),deriv) &
              *radius*det*weights(i)
        end do integrating_pts_1
        call fsparv(kv,kp,g,kdiag)
 end do elements_2
 kvh=kv    
 read (10,*)loaded_nodes
 if(loaded_nodes/=0)read(10,*)(k,loads(nf(1,k)),i=1,loaded_nodes)
 read (10,*) fixed_nodes
 if(fixed_nodes/=0)then
     allocate( node(fixed_nodes),no(fixed_nodes),value(fixed_nodes))
     read(10,*) (node(i),value(i),i=1,fixed_nodes)
     do i=1,fixed_nodes; no(i)=nf(1,node(i)); end do
     kv(kdiag(no))=kv(kdiag(no))+1.e20; loads(no)=kv(kdiag(no))*value 
 end if     
!------------------------equation solution-------------------------------------
 call sparin(kv,kdiag); call spabac(kv,loads,kdiag)      
!------------------------ retrieve flow rate ---------------------------------
 call linmul_sky(kvh,loads,disps,kdiag)
 write(11,'(a)')"The nodal values are:"
 write(11,'(a)')"             Potentials   Flow rate"
 do k=1,nn
   write(11,'(i5,a,2f12.2)')k,"   ",loads(nf(1,k)),disps(nf(1,k)); end do
 write(11,'(a)')"      Inflow     Outflow" 
 write(11,'(2f12.2)')sum(disps,mask=disps<.0),sum(disps,mask=disps>.0) 
end program p72                         


