    program p80      
!------------------------------------------------------------------------------
!      program 8.0 conduction equation on rectangular 
!      area using 4-node quadrilateral elements
!      implicit integration in time using 'theta' method
!------------------------------------------------------------------------------
 use new_library    ;  use geometry_lib  ;      implicit none
 integer::nels,nxe,neq,nband,nn,nr,nip,nodof=1,nod=4,ndof,ndim=2,            &
          i,j,k,iel,nstep,npri,nres
 real::aa,bb,permx,permy,det,theta,dtim,val0,time
 character(len=15) :: element = 'quadrilateral'
!------------------------------ dynamic arrays--------------------------------
 real ,allocatable :: bp(:),bk(:),loads(:),points(:,:),kay(:,:),coord(:,:),  &
                      fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),        &
                      kp(:,:), pm(:,:), newlo(:) ,funny(:,:),g_coord(:,:)
 integer, allocatable :: nf(:,:), g(:) , num(:) , g_num(:,:) ,g_g(:,:)         
!------------------------input and initialisation-----------------------------
  open (10,file='p80.dat',status=    'old',action='read')
  open (11,file='p80.res',status='replace',action='write')                    
  read (10,*) nels,nxe,nn,nip,aa,bb,permx,permy ,                            &
              dtim,nstep,theta,npri,nres 
  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),weights(nip),kay(ndim,ndim),     &
            coord(nod,ndim), fun(nod), jac(ndim,ndim),g_coord(ndim,nn),      &
            der(ndim,nod), deriv(ndim,nod), pm(ndof,ndof),g_num(nod,nels),   &
            kp(ndof,ndof), g(ndof),funny(1,nod),num(nod),g_g(ndof,nels))
  kay=0.0 ; kay(1,1)=permx; kay(2,2)=permy                            
  call sample (element,points,weights)
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf)  ;    neq=maxval(nf)                                        
!-------------loop the elements to find nband and set up global arrays--------- 
   nband = 0
   elements_1: do iel = 1 , nels
               call geometry_4qx(iel,nxe,aa,bb,coord,num)
               g_num(:,iel) = num; g_coord(: , num ) = transpose(coord)  
               call num_to_g (num,nf,g);   g_g( : , iel ) = g
               if(nband<bandwidth(g)) nband = bandwidth(g)
   end do elements_1     
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                  &
                              "Element ",k,"        ",g_num(:,k); end do  
    allocate(bp(neq*(nband+1)),bk(neq*(nband+1)),loads(0:neq),newlo(0:neq))
      bp = 0.;  bk = 0.
      write(11,'(2(a,i5))')                                                    &
              "There are ",neq,"  equations and the half-bandwidth is ",nband
!-------------------- element integration and assembly------------------------- 
 elements_2: do iel = 1 , nels
             num = g_num(:,iel) ; coord = transpose( g_coord( : , num )) 
             g = g_g( : , iel )     ;     kp=0.0 ; pm=0.0                
       gauss_pts:  do i =1 , nip
               call shape_der (der,points,i) ; call shape_fun(fun,points,i)
               funny(1,:)=fun(:) ; jac = matmul(der,coord)
               det=determinant(jac); call invert(jac);deriv = matmul(jac,der)  
               kp = kp + matmul(matmul(transpose(deriv),kay),deriv) &
                    *det*weights(i)*theta*dtim
               pm  =  pm + matmul( transpose(funny),funny)*det*weights(i) 
       end do gauss_pts                                                      
   call formkv (bk,kp,g,neq) ; call formkv(bp,pm,g,neq)
 end do elements_2
!------------------------factorise left hand side-----------------------------
           bp=bp+bk; bk=bp-bk/theta ; call banred(bp,neq)                    
    read(10,*) val0; loads=val0   
!-------------------time stepping recursion-----------------------------------
  write(11,'(a,i5)') "    Time     Pressure at node" , nres    
 timesteps: do j=1,nstep
               time=j*dtim ; call linmul(bk,loads,newlo)
               call bacsub(bp,newlo) ;  loads=newlo
               if(j/npri*npri==j)write(11,'(2e12.4)')time,loads(nf(:,nres))
            end do timesteps
end program p80                          
