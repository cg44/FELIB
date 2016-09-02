    program p83      
!------------------------------------------------------------------------------
!      program 8.3 diffusion - convection equation on rectangular 
!      area using 4-node quadrilateral elements
!      untransformed solution by Galerkin's method
!      implicit integration in time using 'theta' method
!------------------------------------------------------------------------------
 use new_library        ;  use  geometry_lib        ;    implicit none
 integer::nels,nxe,neq,nband,nn,nr,nip,nodof=1,nod=4,ndof,ndim=2,           &
          i,j,k,l,iel,nstep,npri,nfix
 real::aa,bb,permx,permy,det,theta,dtim,ux,uy,time,part1,part2
 character (len=15) :: element = 'quadrilateral'
!-------------------------- dynamic arrays-------------------------------------
real ,allocatable ::kb(:,:),pb(:,:),loads(:),points(:,:),kay(:,:),coord(:,:),& 
                    fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),          &
                    kp(:,:), pm(:,:), ans(:) ,funny(:,:),g_coord(:,:),       &
                    storpb(:),work(:,:),copy(:,:) , dtkd(:,:)
integer, allocatable :: nf(:,:), g(:) , num(:) , g_num(:,:) ,g_g(:,:) , no(:) 
!-------------------------input and initialisation-----------------------------
  open (10,file='p83.dat',status=    'old',action='read')
  open (11,file='p83.res',status='replace',action='write')                    
  read (10,*) nels,nxe,nn,nip,aa,bb,permx,permy,ux,uy,                       &
              dtim,nstep,theta,npri,nfix 
  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),weights(nip),kay(ndim,ndim),    &
             coord(nod,ndim), fun(nod),jac(ndim,ndim),g_coord(ndim,nn)     ,&  
             der(ndim,nod), deriv(ndim,nod), pm(ndof,ndof),g_num(nod,nels), &
             kp(ndof,ndof), g(ndof),funny(1,nod),num(nod),g_g(ndof,nels),   &
             no(nfix),storpb(nfix),dtkd(ndof,ndof))
  kay=0.0 ; kay(1,1)=permx; kay(2,2)=permy                             
  call sample(element,points,weights)
  nf=1;  read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf);neq=maxval(nf)                                           
!-------------loop the elements to find nband and set up global arrays--------- 
   nband = 0
   elements_1: do iel = 1 , nels
               call geometry_4qx(iel,nxe,aa,bb,coord,num)
               g_num( : , iel ) = num; g_coord(:,num) = transpose(coord)  
               call num_to_g(num,nf,g)  ;   g_g( : , iel ) = g
               if(nband<bandwidth(g)) nband = bandwidth(g)
   end do elements_1    
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                 &
                          "Element ",k,"        ",g_num(:,k); end do 
    allocate(kb(neq,2*nband+1),pb(neq,2*nband+1),loads(0:neq),ans(0:neq),&
             work(nband+1,neq),copy(nband+1,neq))
      kb = 0.;  pb = 0. ; work = .0 ; loads = .0
      write(11,'(2(a,i5))')                                                    &
              "There are ",neq,"  equations and the half-bandwidth is",nband
!------------------ element integration and assembly------------------------   
 elements_2: do iel = 1 , nels
             num = g_num(: , iel ) ; coord =  transpose(g_coord( : , num )) 
             g = g_g( : , iel )    ;    kp=0.0 ; pm=0.0                        
       gauss_pts:  do i =1 , nip
               call shape_der (der,points,i) ; call shape_fun(fun,points,i)
               funny(1,:)=fun(:) ; jac = matmul(der,coord)
               det=determinant(jac); call invert(jac); deriv = matmul(jac,der)
               do k=1,nod;do l=1,nod
                 part1=permx*deriv(1,k)*deriv(1,l)+permy*deriv(2,k)*deriv(2,l)
                 part2=ux*fun(k)*deriv(1,l)+uy*fun(k)*deriv(2,l)
                 dtkd(k,l)=(part1-part2)*det*weights(i)
               end do; end do
               kp = kp + dtkd   
               pm  =  pm + matmul( transpose(funny),funny)*det*weights(i) 
       end do gauss_pts
               pm = pm/(theta*dtim)
    call formtb (kb,kp,g) ; call formtb(pb,pm,g)
 end do elements_2
!------------------------specify fixed nodal values --------------------------
           pb = pb + kb; kb = pb - kb / theta   ; read(10,*) no  
           pb(no,nband+1) = pb(no,nband+1) + 1.e20 ; storpb = pb(no,nband+1)
!------------------------factorise left hand side-----------------------------
           call gauss_band(pb,work) 
!-------------------time stepping recursion-----------------------------------
   write(11,'(a)') "    Time     Concentration"
 timesteps: do j=1,nstep
               time=j*dtim ; copy = work ; call bantmul(kb,loads,ans); ans(0)=.0
               if(time<=.2) then ;ans(no)=storpb; else; ans(no) = .0; end if
               call solve_band(pb,copy,ans) ;  ans(0)=.0; loads=ans 
               if(j/npri*npri==j)write(11,'(2e12.4)') time,loads(nf(:,3))
            end do timesteps
end program p83

