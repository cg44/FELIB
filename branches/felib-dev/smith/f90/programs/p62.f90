    program p62      
!-----------------------------------------------------------------------
!      program 6.2 plane strain of an elastic-plastic(Mohr-Coulomb) solid
!      using 8-node quadrilateral elements; viscoplastic strain method
!------------------------------------------------------------------------
 use new_library   ;   use geometry_lib      ; implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,          &
          i,k,iel,iters,limit,incs,iy,ndim=2
 logical::converged       ;    character(len=15) :: element='quadrilateral'
 real::e,v,det,phi,c,psi,gama,dt,f,dsbar,dq1,dq2,dq3,lode_theta,              &
           sigm,pi,tnph,phif,snph,cf,tol                                       
!---------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: kb(:,:),loads(:),points(:,:),bdylds(:),bot(:),fos(:),&
                         evpt(:,:,:),oldis(:),top(:),depth(:),gravlo(:),      &
                         dee(:,:),coord(:,:),fun(:),jac(:,:),weights(:),      &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:)  
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)    
!-------------------------input and initialisation-----------------------------
  open (10,file='p62.dat',status=    'old',action='read')
  open (11,file='p62.res',status='replace',action='write')                    
  read (10,*) phi,c,psi,gama,e,v,    nels,nxe,nye,nn,nip,tol,limit
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            top(nxe+1),depth(nye+1),num(nod),dee(nst,nst),evpt(nst,nip,nels), &
            bot(nxe+1),coord(nod,ndim),fun(nod),g_g(ndof,nels),               &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst))                  
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)    ;  read(10,*) top , bot , depth           
!---------- loop the elements to find nband and set up global arrays ----------
     nband = 0
       elements_1:   do iel = 1 , nels
                       call slope_geometry(iel,nye,top,bot,depth,coord,num)
                       call num_to_g(num,nf,g);      g_num(:,iel)=num
                       g_coord(:,num)=transpose(coord);   g_g(:,iel) = g
                       if (nband<bandwidth(g)) nband = bandwidth(g)
      end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do    
   write(11,'(a,i5,a,i5)')                                                     &
           "The system has",neq,"  equations and the half-bandwidth is ",nband
allocate(kb(neq,nband+1),loads(0:neq),bdylds(0:neq),oldis(0:neq),gravlo(0:neq)) 
           kb=0.0; oldis=0.0; gravlo=0.0 
  call deemat(dee,e,v);      call sample(element,points,weights)
  pi = acos( -1. ); tnph = tan(phi*pi/180.)
!----------------- element stiffness integration and assembly------------------
 elements_2: do iel = 1 , nels
                num = g_num(:,iel) ; coord = transpose(g_coord( : , num )) 
                g = g_g( : , iel )  ;      km=0.0   ;  eld = .0
             gauss_pts_1:  do i =1 , nip    ; call shape_fun(fun,points,i)
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               do k=2,ndof,2;eld(k)=eld(k)+fun(k/2)*det*weights(i);end do
            end do gauss_pts_1   
   call formkb (kb,km,g)
   gravlo ( g ) = gravlo ( g ) - eld * gama ; gravlo(0) = .0
 end do elements_2                                                            
!------------------------ factorise left hand side----------------------------- 
          call cholin(kb)                                                  
!------------------------trial factor of safety loop---------------------------
    read(10,*) incs;   allocate ( fos (incs ))  ;    read(10,*) fos
    load_increments: do iy=1,incs
    phif = atan(tnph/fos(iy))*180./pi; snph = sin(phif*pi/180.)
    dt = 4.*(1.+v)*(1.-2.*v)/(e*(1.-2.*v+snph**2)) ; cf = c/fos(iy)
    write(11,'(a,i5)') "Load increment",iy ;  iters=0;  bdylds=.0;  evpt=.0
!----------------------------   iteration loop  -------------------------------
   iterations: do    
    iters=iters+1;  loads = gravlo + bdylds   ;  call chobac(kb,loads)
!-------------------------   check convergence --------------------------------
      call checon(loads,oldis,tol,converged)
      if(iters==1)converged=.false. ;  if(converged.or.iters==limit)bdylds=.0
!----------------------- go round the Gauss Points  ---------------------------
      elements_3: do iel = 1 , nels
       bload=.0
       num = g_num( : , iel ) ; coord =transpose( g_coord( : ,num ))
       g = g_g( : , iel )  ;       eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac); call invert(jac) ; deriv = matmul(jac,der) 
          call beemat (bee,deriv);   eps=matmul(bee,eld)
          eps = eps -evpt( : , i , iel)    ;        sigma=matmul(dee,eps)
          call invar(sigma,sigm,dsbar,lode_theta)                             
!------------------  check whether yield is violated  -------------------------
         call mocouf (phif, cf , sigm, dsbar , lode_theta , f )
         if(converged.or.iters==limit) then
         devp=sigma 
           else
           if(f>=.0) then
           call mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)
           call formm(sigma,m1,m2,m3)   ;    flow=f*(m1*dq1+m2*dq2+m3*dq3) 
           erate=matmul(flow,sigma)     ;    evp=erate*dt
           evpt(:,i,iel)=evpt(:,i,iel)+evp;  devp=matmul(dee,evp) 
         end if; end if
      if(f>=.0) then
        eload=matmul(devp,bee) ; bload=bload+eload*det*weights(i)
      end if
    end do gauss_points_2
!-------------------  compute the total bodyloads vector ----------------------
    bdylds( g ) = bdylds( g ) + bload      ; bdylds(0) = .0
  end do elements_3             
  if(converged.or.iters==limit)exit
 end do iterations
 write(11,'(a)') "    fos    max displacement"
 write(11,'(2e12.4)')fos(iy),maxval(abs(loads)) 
 write(11,'(a,i5,a)') "It took",iters,"   iterations to converge"
 if(iters==limit)stop
end do load_increments
end program p62
