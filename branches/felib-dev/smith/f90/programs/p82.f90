    program p82      
!------------------------------------------------------------------------------
!      program 8.2 diffusion - convection equation on rectangular 
!      area using 4-node quadrilateral elements
!      self-adjoint transformation
!      implicit integration in time using 'theta' method
!------------------------------------------------------------------------------
 use new_library   ;  use geometry_lib    ;    implicit none
 integer::nels,nxe,neq,nband,nn,nr,nip,nodof=1,nod=4,ndof,ndim=2,            &
          i,j,k,iel,nstep,npri
 real::aa,bb,permx,permy,det,theta,dtim,ux,uy,time,f1,f2
 character (len=15) :: element = 'quadrilateral' 
!---------------------------- dynamic arrays---------------------------------- 
real ,allocatable ::kb(:,:),pb(:,:),loads(:),points(:,:),kay(:,:),coord(:,:),& 
                    fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),          &
                    kp(:,:), pm(:,:), ans(:) ,funny(:,:),g_coord(:,:)
 integer, allocatable :: nf(:,:), g(:) , num(:) , g_num(:,:) ,g_g(:,:)         
!-------------------------input and initialisation-----------------------------
  open (10,file='p82.dat',status=    'old',action='read')
  open (11,file='p82.res',status='replace',action='write')                    
  read (10,*) nels,nxe,nn,nip,aa,bb,permx,permy,ux,uy,                       &
              dtim,nstep,theta,npri 
  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),weights(nip),kay(ndim,ndim),     &
            coord(nod,ndim), fun(nod), jac(ndim,ndim),g_coord(ndim,nn),      &
            der(ndim,nod), deriv(ndim,nod), pm(ndof,ndof),g_num(nod,nels),   &
            kp(ndof,ndof), g(ndof),funny(1,nod),num(nod),g_g(ndof,nels))
  kay=0.0 ; kay(1,1)=permx; kay(2,2)=permy                            
  call sample(element,points,weights)
  nf=1;  read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)                                           
!----------loop the elements to find nband and set up global arrays -----------
  nband = 0
   elements_1: do iel = 1 , nels
               call geometry_4qx(iel,nxe,aa,bb,coord,num)
               g_num(:,iel) = num; g_coord(:,num) = transpose(coord)  
               call num_to_g(num,nf,g);  g_g( : , iel ) = g
               if(nband<bandwidth(g)) nband = bandwidth(g)
   end do elements_1  
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                  &
                          "Element ",k,"        ",g_num(:,k); end do   
    allocate(kb(neq,nband+1),pb(neq,nband+1),loads(0:neq),ans(0:neq))
      kb = 0.;  pb = 0. ; loads = .0
      write(11,'(2(a,i5))')                                                    &
              "There are ",neq," equations and the half-bandwidth is ",nband
!------- element integration and assembly------------------------------------- 
 elements_2: do iel = 1 , nels
             num = g_num(:,iel) ; coord =  transpose(g_coord(: , num )) 
             g = g_g( : , iel )    ;    kp=0.0 ; pm=0.0                 
       gauss_pts:  do i =1 , nip
               call shape_der (der,points,i) ; call shape_fun(fun,points,i)
               funny(1,:)=fun(:) ; jac = matmul(der,coord)
               det=determinant(jac); call invert(jac); deriv = matmul(jac,der) 
               kp = kp + matmul(matmul(transpose(deriv),kay),deriv) &
                    *det*weights(i)
               pm  =  pm + matmul( transpose(funny),funny)*det*weights(i) 
       end do gauss_pts
               kp = kp + pm*(ux*ux/permx+uy*uy/permy)*.25
               pm = pm/(theta*dtim)
!------------------- derivative boundary conditions ---------------------------
       if(iel==1) then
          kp(2,2)=kp(2,2)+uy*aa/6.; kp(2,3)=kp(2,3)+uy*aa/12.
          kp(3,2)=kp(3,2)+uy*aa/12.; kp(3,3)=kp(3,3)+uy*aa/6.
       else if(iel==nels) then
         kp(1,1)=kp(1,1)+uy*aa/6.; kp(1,4)=kp(1,4)+uy*aa/12.
         kp(4,1)=kp(4,1)+uy*aa/12.; kp(4,4)=kp(4,4)+uy*aa/6.
       end if       
    call formkb (kb,kp,g) ; call formkb(pb,pm,g)
 end do elements_2
!------------------------factorise left hand side----------------------------
           f1=uy*aa/(2.*theta); f2 = f1
           pb = pb + kb; kb = pb - kb/theta ; call cholin(pb) 
!-------------------time stepping recursion-----------------------------------
   write(11,'(a,i5)') "    Time     Concentration at node "  , nn 
 timesteps: do j=1,nstep
               time=j*dtim ; call banmul(kb,loads,ans)
               ans(neq)=ans(neq)+f1; ans(neq-1) = ans(neq-1)+f2
               call chobac(pb,ans) ;  loads=ans        
               if(j/npri*npri==j)write(11,'(2e12.4)')  &
                         time,loads(nf(:,nn))*exp(ux/2./permx)*exp(uy/2./permy)
            end do timesteps
end program p82
