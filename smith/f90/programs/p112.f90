 program p112       
!------------------------------------------------------------------------
!      program 11.2 forced vibration of a rectangular elastic
!      solid in plane strain using uniform 8-node quadrilateral elements
!      numbered in the y direction - lumped or consistent mass
!      implicit integration by Wilson method
!------------------------------------------------------------------------
 use new_library;     use  geometry_lib   ;      implicit none
 integer::nels,nye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=3,ndof,              &
          i,k,iel,ndim=2,nstep,npri
 real   ::aa,bb,e,v,det,rho,alpha,beta,omega,theta,period,pi,dtim,area,       &
          c1,c2,c3,c4,c5,c6,c7,c8,time
 logical:: consistent = .false.; character(len=15) :: element='quadrilateral' 
!----------------------------- dynamic arrays---------------------------------
 real    ,allocatable :: kv(:),mv(:),loads(:),points(:,:),dee(:,:),coord(:,:),&
                         fun(:),jac(:,:), der(:,:),deriv(:,:), weights(:),    &
                         bee(:,:),km(:,:),g_coord(:,:),                       &
                         emm(:,:),ecm(:,:),f1(:),x0(:),d1x0(:),d2x0(:),       &
                         x1(:),d1x1(:),d2x1(:)
 integer, allocatable :: nf(:,:), g(:) , num(:)  , g_num(:,:) , g_g(:,:)       
!-------------------------input and initialisation----------------------------
  open (10,file='p112.dat',status=    'old',action='read')
  open (11,file='p112.res',status='replace',action='write')               
  read (10,*) nels,nye,nn,nip,aa,bb,rho,e,v,                                  &
              alpha,beta,nstep,npri,theta,omega 
  ndof=nod*nodof         
  allocate ( nf(nodof,nn), points(nip,ndim),g(ndof), g_coord(ndim,nn),        &
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),weights(nip),         &
            der(ndim,nod), deriv(ndim,nod), bee(nst,ndof), km(ndof,ndof),     &
            num(nod),g_num(nod,nels),g_g(ndof,nels),                          &
            emm(ndof,ndof),ecm(ndof,ndof),fun(nod))  
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf (nf);neq=maxval(nf)
  nband = 0      ; pi=acos(-1.)   ; period = 2.*pi/omega ; dtim =period/20.
  call deemat (dee,e,v); call sample(element,points,weights)
!----------------loop the elements to find bandwidth and neq-------------------
  elements_1: do iel = 1 , nels
              call geometry_8qy(iel,nye,aa,bb,coord,num)
              call num_to_g ( num , nf , g )
              g_num(:,iel)=num; g_coord(:,num)=transpose(coord);g_g(:,iel) = g
              if(nband<bandwidth(g))nband=bandwidth(g)
             end do elements_1  
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do       
    write(11,'(2(a,i5))')                                                     &
            "There are ",neq,"  equations and the half-bandwidth is", nband
 allocate(kv(neq*(nband+1)),mv(neq*(nband+1)),x0(0:neq),d1x0(0:neq),x1(0:neq),&
        d2x0(0:neq),loads(0:neq),d1x1(0:neq),d2x1(0:neq),f1(neq*(nband+1)))  
   kv=.0; mv=.0; x0=.0; d1x0=.0; d2x0=.0       
!------------ element stiffness and mass integration and assembly-------------- 
 elements_2: do iel = 1 , nels
             num = g_num( : , iel ); coord = transpose(g_coord(: , num )) 
             g = g_g( : , iel )    ; km=0.0   ; area = .0 ; emm = .0 
          gauss_points_1: do i = 1 , nip     
               call shape_der (der,points,i) ; jac = matmul(der,coord) 
               det = determinant(jac)  ; call invert(jac)
               deriv = matmul(jac,der) ; call beemat (bee,deriv) 
             km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               area = area + det*weights(i); call shape_fun(fun,points,i)
               if(consistent) then
                call ecmat(ecm,fun,ndof,nodof); ecm=ecm*det*weights(i)*rho
                emm = emm + ecm    
              end if
          end do gauss_points_1                  
       if(.not.consistent) then
         do i=1,ndof; emm(i,i)=area*rho*.2 ; end do
         do i=1,13,4 ; emm(i,i)=emm(3,3)*.25; end do
         do i=2,14,4 ; emm(i,i)=emm(3,3)*.25 ; end do
       end if
   call formkv (kv,km,g,neq); call formkv(mv,emm,g,neq)
 end do elements_2 
!-----------------------------factorisation------------------------------------
   c1=6./(theta*dtim)**2; c2=6./(theta*dtim); c3=dtim**2/6.; c4 = 2.
   c5=3.*alpha/(theta*dtim); c6=3.*beta/(theta*dtim); c7=.5*alpha*theta*dtim
   c8=.5*beta*theta*dtim;   f1 = (c1+c5)*mv +(1.+c6)*kv  ; call banred(f1,neq)
!------------------------- time stepping loop ---------------------------------
   time = .0
   write(11,'(a)') "  Time t    cos(omega*t) Displacement"
  timesteps: do i = 1 , nstep
              time = time + dtim  ; loads = .0
              x1 = (c1+c5)*x0 + (c2+2.*alpha)*d1x0 + (2. + c7)*d2x0 
              loads(neq)=theta*cos(omega*time)+(1.-theta)*cos(omega*(time-dtim))
              call linmul(mv,x1,d1x1)  ; d1x1 = loads + d1x1 
              loads = c6*x0 + 2.*beta*d1x0 + c8 * d2x0 
              call linmul(kv,loads,x1); x1 = x1 + d1x1
              call bacsub(f1,x1) 
              d2x1=(x1-x0)*c1-d1x0*c2-d2x0*c4 ;   d2x1=d2x0+(d2x1-d2x0)/theta
              d1x1=d1x0+.5*dtim*(d2x1+d2x0); x1=x0+dtim*d1x0+2.*c3*d2x0+d2x1*c3
              if(i/npri*npri==i)write(11,'(3e12.4)')time,cos(omega*time),x1(neq)
              x0 = x1; d1x0 = d1x1; d2x0 = d2x1
  end do timesteps
end program p112                                                              
