    program p68      
!-----------------------------------------------------------------------
!      program 6.8 general strain of an elastic-plastic(Mohr-Coulomb)
!      solid --- viscoplastic strain method
!------------------------------------------------------------------------
 use new_library  ;  use  geometry_lib       ; implicit none
 integer::nels,neq,nn,nr,nip,nodof,nod,nst,ndof,nprops,np_types,              &
          i,k,iel,iters,limit,incs,iy,ndim,loaded_nodes,fixed_nodes
 logical::converged   ;  character (len=15) :: element
 real::e,v,phi,c,psi,                                                         &
       det,dt,ddt,f,dsbar,dq1,dq2,dq3,lode_theta,sigm,pi,snph,ptot,tol         
!---------------------------- dynamic arrays-----------------------------------
 real   ,allocatable ::  kv(:),loads(:),points(:,:),bdylds(:),oldis(:),       &
                         dee(:,:),coord(:,:),fun(:),jac(:,:),weights(:),      &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),    &
                         val(:),storkv(:),tensor(:,:,:),stress(:),totd(:),    &
                         evpt(:,:,:),value(:),load_store(:),prop(:,:) 
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:), &
                         kdiag(:), sense(:), node(:) , etype(:)        
!---------------------------input and initialisation---------------------------
  open (10,file='p68.dat',status=    'old',action='read')
  open (11,file='p68.res',status='replace',action='write')               
  read (10,*) element,nels,nn,nip,nodof,nod,nst,ndim,incs,tol,limit
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            num(nod),dee(nst,nst),evpt(nst,nip,nels),tensor(nst,nip,nels),    &
            coord(nod,ndim),g_g(ndof,nels),stress(nst),etype(nels),           &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst))
   read(10,*) nprops , np_types 
    allocate(prop(nprops,np_types)) ; read(10,*) prop 
    etype = 1 ; if(np_types>1) read(10,*) etype
!---------------------- read geometry and connectivity ------------------------
  read(10,*) g_coord; read(10,*) g_num   
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)    ;  allocate(kdiag(neq))     ; kdiag = 0
!--------- loop the elements to set up global arrays and kdiag ----------------
      elements_1:   do iel = 1 , nels
                       num = g_num(:,iel); coord = transpose(g_coord(:,num))
                       call num_to_g( num , nf , g )
                       g_g( : , iel ) = g ;   call fkdiag(kdiag,g)
      end do elements_1  
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,3e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,27i3)')                                 &
                             "Element ",k,"   ",g_num(:,k); end do
      kdiag(1)=1; do i=2,neq; kdiag(i)=kdiag(i)+kdiag(i-1); end do             
 write(11,'(a,i5,a,i5,a)')                                                     &
         "The skyline storage is", kdiag(neq),"and there are",neq," equations"
 allocate(kv(kdiag(neq)),loads(0:neq),bdylds(0:neq),oldis(0:neq),totd(0:neq), &
          load_store(0:neq)) 
           kv=0.0; oldis=0.0; totd=0.0 
  call sample(element,points,weights) ;  pi = acos( -1. )
  dt = 100.
  do iel = 1 , nels
   e = prop(1,etype(iel)); v = prop(2,etype(iel))
   snph = sin(prop(3,etype(iel))*pi/180.)
   ddt=4.*(1.+v)*(1.-2.*v)/(e*(1.-2.*v+snph*snph))
   if(ddt<dt)dt=ddt
  end do
 write(11,'(a,e12.4)') "The critical timestep is ",dt
!---------- element stiffness integration and assembly & set stresses-------- 
 elements_2: do iel = 1 , nels
                num = g_num(: , iel) ; coord = transpose(g_coord(: , num )) 
                g = g_g( : , iel )    ;  km=0.0   ; tensor = .0     
             gauss_pts_1:  do i =1 , nip    
               tensor(1:3 , i , iel) = prop( 6 , etype(iel) )
               e=prop(1,etype(iel)); v=prop(2,etype(iel)); call deemat(dee,e,v)
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;   call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
             end do gauss_pts_1   
             call fsparv (kv,km,g,kdiag)
 end do elements_2                                                             
!------------ read prescribed loads/displacements and factorise l.h.s. ------- 
    read(10,*) loaded_nodes
    if(loaded_nodes/=0) then
     read(10,*)(k,loads(nf(:,k)),i=1,loaded_nodes); load_store = loads
    end if                                                                    
    read(10,*) fixed_nodes     
    if(fixed_nodes /=0) then
     allocate(node(fixed_nodes),sense(fixed_nodes),value(fixed_nodes),        &
              no(fixed_nodes),storkv(fixed_nodes))
     read(10,*) (node(i), sense(i), value(i),i=1,fixed_nodes)
     do i=1,fixed_nodes; no(i)=nf(sense(i),node(i)); end do 
     kv(kdiag(no)) = kv(kdiag(no)) + 1.e20  ; storkv = kv(kdiag(no)) 
    end if 
         call sparin (kv,kdiag)                                                 
!-------------------displacement increment loop-------------------------------
    load_increments: do iy=1,incs
    ptot = value(1) * iy
    write(11,'(a,i5)')"Load increment",iy ;  iters=0;  bdylds=.0;  evpt=.0
!--------------------------   iteration loop  --------------------------------
   iterations: do
    iters=iters+1; loads =.0 
      if(loaded_nodes/=0) loads = load_store 
      if(fixed_nodes/=0) loads(no) = storkv * value 
      loads=loads+bdylds; call spabac(kv,loads,kdiag)
!--------------------------   check convergence -------------------------------
      call checon(loads,oldis,tol,converged)
      if(iters==1)converged=.false. ;  if(converged.or.iters==limit)bdylds=.0
!------------------------ go round the Gauss Points ---------------------------
      elements_3: do iel = 1 , nels
       bload=.0
       num = g_num( : , iel ) ; coord = transpose(g_coord( : , num ))
       g = g_g( : , iel )     ; eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ; call invert(jac)
          deriv = matmul(jac,der) ; call beemat (bee,deriv);eps=matmul(bee,eld)
          eps=eps-evpt(:,i,iel);  sigma=matmul(dee,eps) 
          stress = sigma + tensor(:,i,iel)
          call invar(stress,sigm,dsbar,lode_theta)                            
!-------------------  check whether yield is violated -------------------------
         phi = prop(3,etype(iel));c=prop(4,etype(iel)); psi=prop(5,etype(iel))
         call mocouf (phi, c , sigm, dsbar , lode_theta , f )
         if(converged.or.iters==limit) then
         devp=stress 
           else
           if(f>=.0) then
           call mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)
           call formm(stress,m1,m2,m3)
           flow=f*(m1*dq1+m2*dq2+m3*dq3)     ;   erate=matmul(flow,stress)
           evp=erate*dt; evpt(:,i,iel)=evpt(:,i,iel)+evp; devp=matmul(dee,evp) 
         end if; end if
      if(f>=.0) then
        eload=matmul(devp,bee) ; bload=bload+eload*det*weights(i)
      end if
    if(converged.or.iters==limit) then
!----------------------- update the Gauss point stresses ----------------------
                 tensor ( : , i , iel) = stress
    end if
    end do gauss_points_2
!----------------------- compute the total bodyloads vector -------------------
    bdylds( g ) = bdylds( g ) + bload      ; bdylds(0) = .0
  end do elements_3             
  if(converged.or.iters==limit)exit
 end do iterations
 totd = totd + loads
 write(11,'(a,e12.4)') "The displacement is ",totd(1) 
 write(11,'(a)')"   The  stresses  are  "
 write(11,'(6e12.4)') tensor(: , 1 , 1)
 write(11,'(a,i5,a)') "It took",iters,"  iterations to converge"
 if(iters==limit)stop
end do load_increments
end program p68
