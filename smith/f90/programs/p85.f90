    program p85      
!------------------------------------------------------------------------------
!      program 8.5 conduction equation on rectangular area using 4-node
!      quadrilateral elements and a simple explicit algorithm
!------------------------------------------------------------------------------
 use new_library     ;  use  geometry_lib      ;     implicit none
 integer::nels,nxe,neq,nn,nr,nip,nodof=1,nod=4,ndof,ndim=2,                  &
          i,j,k,iel,nstep,npri,nres
 real::aa,bb,permx,permy,det,dtim,val0,time
 character (len=15) :: element = 'quadrilateral'
!------------------------- dynamic arrays--------------------------------------
 real ,allocatable :: loads(:),points(:,:),kay(:,:),coord(:,:),mass(:),       &
                      jac(:,:),der(:,:),deriv(:,:),weights(:),kp(:,:),        &
                      pm(:,:), funny(:,:),g_coord(:,:),globma(:),fun(:),      &
                      store_pm(:,:,:), newlo(:)
 integer, allocatable :: nf(:,:), g(:) , num(:) , g_num(:,:) ,g_g(:,:)         
!-----------------------input and initialisation-------------------------------
  open (10,file='p85.dat',status=    'old',action='read')
  open (11,file='p85.res',status='replace',action='write')                   
  read (10,*) nels,nxe,nn,nip,aa,bb,permx,permy ,                             &
              dtim,nstep,npri,nres 
  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),weights(nip),kay(ndim,ndim),      &
            coord(nod,ndim), fun(nod), jac(ndim,ndim),g_coord(ndim,nn),       &
            der(ndim,nod), deriv(ndim,nod), pm(ndof,ndof),g_num(nod,nels),    &
            kp(ndof,ndof), g(ndof),funny(1,nod),num(nod), g_g(ndof,nels),     &
            globma(0:nn),store_pm(ndof,ndof,nels),mass(ndof))
  kay=0.0 ; kay(1,1)=permx; kay(2,2)=permy     
  globma = .0   ;  do i = 1,nn; nf(nodof,i)=i; end do
  call sample(element,points,weights)                                    
!------------ loop the elements for integration and to store globals  --------  
 elements_1: do iel = 1 , nels
             call geometry_4qx(iel,nxe,aa,bb,coord,num)
             g_num(:,iel) = num;  g_coord(:,num) =transpose( coord )
             call num_to_g(num,nf,g);     kp=0.0 ; pm=0.0                     
       gauss_pts:  do i =1 , nip
               call shape_der (der,points,i) ; call shape_fun(fun,points,i)
               funny(1,:)=fun(:) ; jac = matmul(der,coord)
               det=determinant(jac); call invert(jac); deriv = matmul(jac,der) 
               kp=kp+matmul(matmul(transpose(deriv),kay),deriv)*det*weights(i)
               pm  =  pm + matmul( transpose(funny),funny)*det*weights(i) 
       end do gauss_pts
      do i=1,ndof; mass(i) = sum(pm(i,:)); end do
      pm = .0 ; do i = 1 , ndof; pm(i,i) = mass(i); end do
      store_pm(:,:,iel) = pm - kp*dtim  ; globma(g) = globma(g) + mass
 end do elements_1   
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do      
!--------------take account of initial and boundary conditions----------------
   nf=1; read(10,*) nr ;if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
   call formnf(nf);neq=maxval(nf)
   write(11,'(a,i5)') "The number of equations is : ",   neq
   allocate(loads(0:neq), newlo(0:neq))
   j = 0         ;     globma(0) = .0      
   do i=1,nn; if(nf(1,i)/=0)then;j=j+1;globma(j)=1./globma(i);end if ; end do
   read(10,*) val0; loads=val0   ; loads(0) = .0                  
!---------------  go round the elements  for revised g -----------------------
    elements_2 : do iel = 1 , nels   ! g is different 
                    call geometry_4qx(iel,nxe,aa,bb,coord,num)
                    call num_to_g(num,nf,g);   g_g ( : , iel) = g
    end do elements_2                                                        
!-------------------time stepping recursion-----------------------------------
   write(11,'(a,i5)') "    Time     Pressure at node",nres 
 timesteps: do j=1,nstep
               time=j*dtim        ; newlo = .0      
!---------------  go round the elements  --------------------------------------
    elements_3 : do iel = 1 , nels   ! g is new one 
                    g = g_g ( : , iel) ;   pm = store_pm(: , : , iel) 
                    loads(0) = .0 ; mass = loads(g) ; fun = matmul(pm , mass)
                    newlo ( g ) = newlo ( g ) + fun 
    end do elements_3 
    newlo(0) = .0;     loads = newlo * globma    
  if(j/npri*npri==j)write(11,'(2e12.4)')time,loads(nf(:,nres))
 end do timesteps
end program p85
