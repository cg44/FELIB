    program p86      
!------------------------------------------------------------------------------
!      program 8.6 conduction equation on rectangular area using 4-node
!      quadrilateral elements and an ebe product algorithm
!------------------------------------------------------------------------------
 use new_library     ;  use  geometry_lib      ;     implicit none
 integer::nels,nxe,neq,nn,nr,nip,nodof=1,nod=4,ndof,ndim=2,                  &
          i,j,k,iel,nstep,npri,nres
 real::aa,bb,permx,permy,det,dtim,theta,val0,time
 character (len=15) :: element = 'quadrilateral' 
!---------------------------- dynamic arrays-----------------------------------
 real ,allocatable :: loads(:),points(:,:),kay(:,:),coord(:,:),mass(:),fun(:),&
                      jac(:,:),der(:,:),deriv(:,:),weights(:),kp(:,:),        &
                      pm(:,:), funny(:,:),g_coord(:,:),globma(:),store_kp(:,:,:)
 integer, allocatable :: nf(:,:), g(:) , num(:) , g_num(:,:) ,g_g(:,:)         
!-------------------------input and initialisation-----------------------------
  open (10,file='p86.dat',status=    'old',action='read')
  open (11,file='p86.res',status='replace',action='write')                
  read (10,*) nels,nxe,nn,nip,aa,bb,permx,permy ,                             &
              dtim,nstep,theta,npri,nres 
  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),weights(nip),kay(ndim,ndim),      &
            coord(nod,ndim), fun(nod), jac(ndim,ndim),g_coord(ndim,nn),       &
            der(ndim,nod), deriv(ndim,nod), pm(ndof,ndof),g_num(nod,nels),    &
            kp(ndof,ndof), g(ndof),funny(1,nod),num(nod),g_g(ndof,nels),      &
            globma(0:nn),store_kp(ndof,ndof,nels),mass(ndof))
  kay=0.0 ; kay(1,1)=permx; kay(2,2)=permy     
  globma = .0   ;  do i = 1,nn; nf(nodof,i)=i; end do
  call sample(element,points,weights)                                        
!----------- loop the elements for integration and to store globals  ----------
 elements_1: do iel = 1 , nels
             call geometry_4qx(iel,nxe,aa,bb,coord,num)
             g_num(:,iel) = num;  g_coord(:,num) =transpose( coord)
             call num_to_g(num,nf,g) ;g_g( : , iel ) = g ;  kp=0.0 ; pm=0.0    
       gauss_pts:  do i =1 , nip
               call shape_der (der,points,i) ; call shape_fun(fun,points,i)
               funny(1,:)=fun(:) ; jac = matmul(der,coord)
               det=determinant(jac); call invert(jac); deriv = matmul(jac,der) 
               kp=kp+matmul(matmul(transpose(deriv),kay),deriv)*det*weights(i)
               pm  =  pm + matmul( transpose(funny),funny)*det*weights(i) 
       end do gauss_pts
      store_kp(:,:,iel) = kp 
      do i=1,ndof ; globma(g(i))=globma(g(i))+sum(pm(i,:)); end do
      globma ( 0 ) = .0
 end do elements_1   
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do   
!-------------- recover  element A and B matrices  ---------------------------
 elements_2 : do iel = 1 , nels        ;           g = g_g( : , iel )
               kp = - store_kp(:,:,iel) * (1. - theta) * dtim * .5
               pm =  store_kp(:,:,iel) * theta * dtim * .5
               do i = 1,ndof
                pm (i,i) = pm (i,i) + globma ( g(i))
                kp (i,i) = kp (i,i) + globma ( g(i))
               end do     
               call invert ( pm )  ; pm = matmul( pm , kp)
               store_kp ( : , : , iel) =  pm
 end do elements_2
!--------------take account of initial and boundary conditions----------------
   nf=1; read(10,*) nr ;if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
   call formnf(nf);neq=maxval(nf)
   write(11,'(a,i5)') "The number of equations is : ", neq   
   allocate(loads(0:neq))
   read(10,*) val0; loads=val0   ; loads(0) = .0
!-------------------time stepping recursion-----------------------------------
   write(11,'(a,i5)') "    Time     Pressure at node"   ,nres 
 timesteps: do j=1,nstep
               time=j*dtim 
!---------------  first pass 1 to nels ----------------------------------------
    elements_3 : do iel = 1 , nels   ! g is different 
                    call geometry_4qx(iel,nxe,aa,bb,coord,num)
                    call num_to_g(num,nf,g) ; g_g( : , iel ) = g
                    pm = store_kp(: , : , iel) ; mass = loads(g)
                    fun = matmul(pm , mass); loads(g) = fun; loads(0) = .0
    end do elements_3     
!---------------  second pass nels to 1 ----------------------------------------
    elements_4 : do iel = nels , 1 , -1 ;    g = g_g( : , iel )
                    pm = store_kp(: , : , iel) ; mass = loads(g)
                    fun = matmul(pm , mass); loads(g) = fun; loads(0) = .0
    end do elements_4
  if(j/npri*npri==j)write(11,'(2e12.4)')time,loads(nf(:,nres))
 end do timesteps
end program p86   
