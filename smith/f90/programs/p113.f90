 program p113       
!------------------------------------------------------------------------
!      program 11.3 forced vibration of a rectangular elastic
!      solid in plane strain using uniform 4-node quadrilateral elements
!      numbered in the y direction - lumped or consistent mass
!      mixed explicit/implicit integration 
!------------------------------------------------------------------------
 use new_library;     use  geometry_lib   ;      implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=2,nod=4,nst=3,ndof,          &
          i,k,iel,ndim=2,nstep
 real   ::aa,bb,e,v,det,rho,beta,gamma,dtim,area,c1,c2,time
 character (len=15) :: element = 'quadrilateral'
!----------------------------- dynamic arrays----------------------------------
 real    ,allocatable :: kv(:),mv(:),loads(:),points(:,:),dee(:,:),coord(:,:),&
                         fun(:),jac(:,:), der(:,:),deriv(:,:), weights(:),    &
                         bee(:,:),km(:,:),g_coord(:,:),                       &
                         emm(:,:),ecm(:,:),x0(:),d1x0(:),d2x0(:),             &
                         x1(:),d1x1(:),d2x1(:)
 integer, allocatable :: nf(:,:), g(:),num(:),g_num(:,:),g_g(:,:),kdiag(:)
 logical, allocatable :: type(:)     
!-----------------------input and initialisation-----------------------------
  open (10,file='p113.dat',status=    'old',action='read')
  open (11,file='p113.res',status='replace',action='write')                  
  read (10,*) nxe,nye,nn,nip,aa,bb,rho,e,v,                                  &
              gamma,beta,nstep,dtim 
  nels = nxe*nye ; ndof=nod*nodof ;c1=1./dtim/dtim/beta;c2=gamma/dtim/beta    
  allocate ( nf(nodof,nn), points(nip,ndim),g(ndof), g_coord(ndim,nn),        &
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),weights(nip),         &
            der(ndim,nod), deriv(ndim,nod), bee(nst,ndof), km(ndof,ndof),     &
            num(nod),g_num(nod,nels),g_g(ndof,nels),type(nels),               &
            emm(ndof,ndof),ecm(ndof,ndof),fun(nod))  
  read(10,*) type
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
  call formnf (nf);neq=maxval(nf) ; allocate(kdiag(neq))
  call deemat (dee,e,v);call sample(element,points,weights); kdiag= 0;nband = 0
!------------loop the elements to set globals and find nband ------------------
  elements_1: do iel = 1 , nels
              call geometry_4qy(iel,nye,aa,bb,coord,num)
              call num_to_g (num , nf , g )
              g_num(:,iel)=num; g_coord(:,num)=transpose(coord);g_g(:,iel) = g
              if(nband<bandwidth(g)) nband = bandwidth(g)                      
              if (.not.type(iel)) call fkdiag(kdiag,g) 
  end do elements_1                                                           
      where(kdiag==0)kdiag=1   ; kdiag(1) = 1
      do i=2,neq; kdiag(i)=kdiag(i) + kdiag(i-1); end do                       
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do       
    write(11,'(3(a,i5))')                                                     &
           "There are ",neq," equations,nband=",nband," Sky store=",kdiag(neq)
  allocate(kv(kdiag(neq)),mv(neq*(nband+1)),x0(0:neq),d1x0(0:neq),x1(0:neq),  &
        d2x0(0:neq),loads(0:neq),d1x1(0:neq),d2x1(0:neq))  
        kv= .0 ; mv=.0         
!-------------- element stiffness and mass integration and assembly------------ 
 elements_2: do iel = 1 , nels
             num = g_num( : , iel ); coord = transpose(g_coord(: , num )) 
             g = g_g( : , iel )    ; km=0.0   ; area = .0 ; emm = .0 
             if(type(iel)) then ; do i=1,ndof; emm(i,i) = 1.; end do; end if
          gauss_points_1: do i = 1 , nip     
               call shape_der (der,points,i) ; jac = matmul(der,coord) 
               det = determinant(jac)  ; call invert(jac)
               deriv = matmul(jac,der) ; call beemat (bee,deriv) 
             km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               area = area + det*weights(i); call shape_fun(fun,points,i)
               if(.not. type(iel)) then
                call ecmat(ecm,fun,ndof,nodof); ecm=ecm*det*weights(i)*rho*c1
                emm = emm + ecm    
              end if
          end do gauss_points_1  
         area = area/nod * rho                 
       if( type(iel) ) then
         do i=1,ndof; emm(i,i)=emm(i,i)*area*c1 ; end do
         call fsparv(kv,emm,g,kdiag)   ; km = - km
         call formkv(mv,km,g,neq); call formkv(mv,emm,g,neq)
       else
         call fsparv(kv,km,g,kdiag) ; call fsparv(kv,emm,g,kdiag)
         call formkv(mv,emm,g,neq)
       end if             
 end do elements_2             
!-----------------   initial conditions and factorisation of l.h.s. ----------
  x0 = .0 ;  d1x0 = 1. ;  d2x0 = .0;   call sparin(kv,kdiag)
!----------------------- time stepping loop ----------------------------------
   time = .0
   write(11,'(a)') "   Time  displacement velocity acceleration"
  timesteps: do i = 1 , nstep
              time = time + dtim  ; loads = .0
              d1x1=x0+d1x0*dtim+d2x0*.5*dtim*dtim*(1.-2.*beta)
              call linmul(mv,d1x1,x1); x1=loads+x1 ; call spabac(kv,x1,kdiag)
              d2x1=(x1-d1x1)/dtim/dtim/beta
              d1x1=d1x0+d2x0*dtim*(1.-gamma)+d2x1*dtim*gamma
              write(11,'(f8.5,3e12.4)')time,x1(neq),d1x1(neq),d2x1(neq)
              x0 = x1; d1x0 = d1x1; d2x0 = d2x1
  end do timesteps
end program p113
