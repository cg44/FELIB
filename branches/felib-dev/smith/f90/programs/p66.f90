    program p66      
!-----------------------------------------------------------------------
!      program 6.6 axisymmetric 'undrained' strain of an elastic-plastic
!      (Mohr-Coulomb) solid
!      using 8-node quadrilateral elements; viscoplastic strain method
!------------------------------------------------------------------------
 use new_library    ;   use  geometry_lib  ;      implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,          &
          i,j,k,iel,iters,limit,incs,iy,ndim=2,loaded_nodes
 logical::converged      ; character (len=15) :: element='quadrilateral'
 real::e,v,det,phi,c,psi,dt,f,dsbar,dq1,dq2,dq3,lode_theta,                   &
           sigm,pi,snph,bulk,cons,presc,ptot,radius,tol                        
!----------------------------- dynamic arrays----------------------------------
 real    ,allocatable :: kv(:),loads(:),points(:,:),bdylds(:),totd(:),        & 
                         evpt(:,:,:),oldis(:),width(:),depth(:),stress(:),    &
                         dee(:,:),coord(:,:),jac(:,:),weights(:),storkv(:),   &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),    &
                         tensor(:,:,:),etensor(:,:,:),pore(:,:),fun(:)
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)   
!--------------------------input and initialisation----------------------------
  open (10,file='p66.dat',status=    'old',action='read')
  open (11,file='p66.res',status='replace',action='write')                    
  read (10,*) phi,c,psi,e,v,bulk,cons,     nels,nxe,nye,nn,nip
  ndof=nod*nodof 
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            width(nxe+1),depth(nye+1),num(nod),evpt(nst,nip,nels),            &
            coord(nod,ndim),g_g(ndof,nels),tensor(nst,nip,nels),fun(nod),     &
            etensor(nst,nip,nels),dee(nst,nst),pore(nip,nels),stress(nst),    & 
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst))              
  nf=1; read (10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf);read(10,*) width , depth
!------------ loop the elements to find nband and set up global arrays --------
     nband = 0
       elements_1:   do iel = 1 , nels
                       call geometry_8qyv(iel,nye,width,depth,coord,num)
                       call num_to_g(num,nf,g) ;    g_num(:,iel)=num
                       g_coord(: , num )=transpose(coord); g_g( : , iel ) = g
                       if (nband<bandwidth(g)) nband = bandwidth(g)
       end do elements_1 
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do   
   write(11,'(a,i5,a,i5)')                                                     &
           "The system has ",neq,"  equations and the half-bandwidth is",nband
allocate(kv(neq*(nband+1)),loads(0:neq),bdylds(0:neq),oldis(0:neq),totd(0:neq)) 
           kv=0.0; oldis=0.0; totd=0.0 ; tensor = 0.0; etensor = 0.0  
  call deemat(dee,e,v); call sample(element,points,weights)
!------------------     fluid bulk modulus is "bulk" --------------------------
  do i=1,nst; do j=1,nst;if(i/=3.and.j/=3)dee(i,j)=dee(i,j)+bulk; end do; end do
  pi = acos( -1. ); snph = sin(phi*pi/180.)
  dt = 4.*(1.+ v)*(1.-2.*v)/(e*(1.-2.*v+snph*snph))
!---------- element stiffness integration and assembly & initial conditions----
 elements_2: do iel = 1 , nels 
                num = g_num(: ,iel ) ; coord = transpose (g_coord(: ,num )) 
                g = g_g( : ,iel )    ;      km=0.0     
             gauss_pts_1:  do i =1 , nip    ; call shape_fun(fun,points,i)
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv=matmul(jac,der);call bmataxi(bee,radius,coord,deriv,fun) 
             km=km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)*radius
               tensor(1:2,i,iel)=cons; tensor(4,i,iel)=cons 
             end do gauss_pts_1  
   call formkv (kv,km,g,neq)
 end do elements_2    
!--------------- prescribe displacements and factorise l.h.s. ---------------
      read(10,*) loaded_nodes ; allocate(no(loaded_nodes),storkv(loaded_nodes))
          read(10,*)no , presc  , incs , tol , limit
          do i=1,loaded_nodes    
             kv(nf(2,no(i)))=kv(nf(2,no(i)))+1.e20 ; storkv(i)=kv(nf(2,no(i)))
          end do               ;      call banred(kv,neq)
!-------------------displacement increment loop--------------------------------
    call deemat(dee,e,v)
   load_increments: do iy=1,incs
    ptot = presc * iy
    write(11,'(/,a,i5)') 'Load increment',iy ;  iters=0;  bdylds=.0;  evpt=.0
!---------------------------   iteration loop  --------------------------------
   iterations: do
    iters=iters+1;  loads = .0
     do i=1,loaded_nodes;loads(nf(2,no(i)))=storkv(i)*presc; end do 
     loads = loads + bdylds  ;  call bacsub(kv,loads)
!------------------------   check convergence ---------------------------------
      call checon(loads,oldis,tol,converged)
      if(iters==1)converged=.false. 
!----------------------- go round the Gauss Points ----------------------------
      elements_3: do iel = 1 , nels
       bload=.0
       num = g_num( : , iel ) ; coord = transpose( g_coord( : , num ))
       g = g_g( : , iel )     ; eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          call shape_fun(fun,points,i)
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ;  call invert(jac)
          deriv = matmul(jac,der) ;  call bmataxi (bee,radius,coord,deriv,fun)
          eps=matmul(bee,eld); det = det * radius;  eps=eps-evpt(:,i,iel)
          sigma=matmul(dee,eps)  ;    stress=sigma+tensor(: , i , iel)
          call invar(stress,sigm,dsbar,lode_theta) 
!--------------------  check whether yield is violated ------------------------
          call mocouf (phi, c , sigm, dsbar , lode_theta , f )
          if(f>=.0) then
           call mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)
           call formm(stress,m1,m2,m3)
           flow=f*(m1*dq1+m2*dq2+m3*dq3)     ;   erate=matmul(flow,stress)
           evp=erate*dt; evpt(:,i,iel)=evpt(:,i,iel)+evp;devp=matmul(dee,evp) 
           eload=matmul(devp,bee)    ; bload=bload+eload*det*weights(i)
          end if
        if(converged.or.iters==limit)  then
!----------------   update stresses and calculate porepressures  --------------
          tensor(:,i,iel)=stress
          etensor(:,i,iel)=etensor(:,i,iel)+eps+evpt(:,i,iel)
          pore(i,iel)=(etensor(1,i,iel)+etensor(2,i,iel)+etensor(4,i,iel))*bulk
        end if
    end do gauss_points_2
!         compute the total bodyloads vector
    bdylds( g ) = bdylds( g ) + bload      ; bdylds(0) = .0
  end do elements_3             
  if(converged.or.iters==limit)exit
 end do iterations
  totd = totd + loads
 write(11,'(a,e12.4)') "     Displacement" , ptot
 write(11,'(a,3e12.4)')                                                        &
         "     Effective stresses",tensor(1,1,1),tensor(2,1,1),tensor(4,1,1)
 write(11,'(a,2e12.4)') "     Deviator stress and porepressure",dsbar,pore(1,1) 
 write(11,'(a,i5,a)') "It took",iters,"  iterations to converge"
 if(iters==limit)stop
end do load_increments
end program p66
