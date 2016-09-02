    program p61      
!---------------------------------------------------------------------------
!      program 6.1 plane strain of an elastic-plastic(Von Mises) solid
!      using 8-node quadrilateral elements; viscoplastic strain method 
!      mesh - free method using preconditioned conjugate gradients 
!---------------------------------------------------------------------------
 use new_library   ;    use  geometry_lib  ;      ; implicit none
 integer::nels,nxe,nye,neq,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,loaded_nodes,   &
          i,k,iel,plasiters,plasits,cjiters,cjits,incs,iy,ndim=2,cjtot
 logical:: plastic_converged, cj_converged
 real::e,v,det,cu,dt,ptot,f,dsbar,dq1,dq2,dq3,lode_theta,sigm,                &
     up,alpha,beta,big,plastol,cjtol;character(len=15)::element='quadrilateral' 
!---------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: loads(:),points(:,:),totd(:),bdylds(:),pmul(:),      &
                         evpt(:,:,:),oldis(:),width(:),depth(:),qinc(:),      &
                         tensor(:,:,:),val(:,:),stress(:),storkm(:,:,:),      &
                         dee(:,:),coord(:,:),jac(:,:),weights(:),             &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),    &
                         p(:),x(:),xnew(:),u(:),diag_precon(:),d(:),utemp(:)  
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)    
!-----------------------input and initialisation-------------------------------
  open (10,file='p61.dat',status=    'old',action='read')
  open (11,file='p61.res',status='replace',action='write')                
  read (10,*) cu,e,v,  nels,nxe,nye,nn,nip,    plasits,cjits,plastol,cjtol
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            width(nxe+1),depth(nye+1),num(nod),dee(nst,nst),coord(nod,ndim),  &
            evpt(nst,nip,nels), tensor(nst,nip,nels),pmul(ndof),utemp(ndof),  &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst),stress(nst),    &
            g_g(ndof,nels),storkm(ndof,ndof,nels))                              
      nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
      call formnf(nf); neq=maxval(nf)   ;      read(10,*) width, depth        
! ---------------loop the elements to set up global arrays --------------------
        elements_1:   do iel = 1 , nels
                        call geometry_8qyv(iel,nye,width,depth,coord,num)
                        call num_to_g( num , nf , g );  g_num(:,iel)=num
                        g_coord(:,num)=transpose(coord);   g_g(:,iel) = g      
        end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do  
    write(11,'(a,i5)') "The number of equations to be solved is",neq          
 allocate(loads(0:neq),bdylds(0:neq),oldis(0:neq),totd(0:neq),p(0:neq),       &
          x(0:neq),xnew(0:neq),u(0:neq),diag_precon(0:neq),d(0:neq)) 
   oldis=0.0; totd=0.0 ; tensor = 0.0
   p = .0;  xnew = .0; diag_precon = .0
  call deemat(dee,e,v); call sample(element,points,weights)
  dt=4.*(1.+v)/(3.*e)
!---------- element stiffness integration,storage and preconditioner -------- 
 elements_2: do iel = 1 , nels
                num = g_num(:,iel) ; coord = transpose(g_coord(:  , num ))
                g = g_g(:,iel)     ;          km=0.0
             gauss_pts_1:  do i =1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
            end do gauss_pts_1   
       storkm(:,:,iel) = km
       do k=1,ndof; diag_precon(g(k))=diag_precon(g(k))+km(k,k); end do
 end do elements_2  
       diag_precon(1:neq)=1./diag_precon(1:neq) ; diag_precon(0)=.0           
!-----------------------read load weightings ----------------------------------
    read(10,*)loaded_nodes  ; allocate(no(loaded_nodes),val(loaded_nodes,ndim))
    read(10,*)(no(i),val(i,:),i=1,loaded_nodes)              
!---------------------------load increment loop--------------------------------
    read(10,*)incs  ; allocate(qinc(incs));  read(10,*)qinc
    ptot=.0
   load_increments: do iy=1,incs
    write (11,'(a,i5)') "Load increment" ,iy 
    ptot=ptot+qinc(iy) ; plasiters=0;bdylds=.0;evpt=.0  ;  cjtot = 0
!--------------------------   plastic iteration loop   ------------------------
   plastic_iterations: do
    plasiters=plasiters+1;  loads=.0
      do i=1,loaded_nodes ; loads(nf(:,no(i)))=val(i,:)*qinc(iy) ;  end do
      loads=loads+bdylds      ; d=diag_precon*loads ; p = d   ; x = .0
!-------------------   solve the simultaneous equations by pcg ----------------
       cjiters = 0
      conjugate_gradients:  do
       cjiters = cjiters + 1 ; u = .0
      elements_3 : do iel = 1 , nels
                      g = g_g( : , iel ); km = storkm( : , : ,iel)
                      pmul = p(g); utemp=matmul(km,pmul)
!dir$ ivdep
           do i = 1 , ndof
              u(g(i)) = u(g(i)) +  utemp(i)
           end do
      end do elements_3 
!-------------------------------pcg process -----------------------------------
    up =dot_product(loads,d); alpha=up/dot_product(p,u)
    xnew = x + p* alpha; loads = loads - u*alpha; d = diag_precon*loads
    beta = dot_product(loads,d)/up; p = d + p * beta
    big = .0; cj_converged = .true.
    do i = 1,neq; if(abs(xnew(i))>big)big=abs(xnew(i)); end do
    do i = 1,neq; if(abs(xnew(i)-x(i))/big>cjtol)cj_converged=.false.;end do
    x = xnew
    if(cj_converged.or.cjiters==cjits) exit
      end do conjugate_gradients
      cjtot = cjtot + cjiters  
!---------------------------- end of pcg process ------------------------------
    loads = xnew           ; loads(0) = .0
!-----------------------   check plastic convergence  -------------------------
      call checon(loads,oldis,plastol,plastic_converged)
      if(plasiters==1)plastic_converged=.false. 
      if(plastic_converged.or.plasiters==plasits)bdylds=.0
!------------------------ go round the Gauss Points ---------------------------
      elements_4: do iel = 1 , nels
       bload=.0
       num = g_num( : ,iel) ; coord =transpose( g_coord( : , num ))
       g = g_g( : , iel )    ; eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac);   call invert(jac) ; deriv = matmul(jac,der)
          call beemat (bee,deriv);eps=matmul(bee,eld)
          eps = eps - evpt(: ,i ,iel) ;sigma=matmul(dee,eps)
          stress = sigma+tensor(: , i, iel)
          call invar(stress,sigm,dsbar,lode_theta)                            
!----------------------  check whether yield is violated  ---------------------
         f=dsbar-sqrt(3.)*cu
         if(plastic_converged.or.plasiters==plasits) then
         devp=stress 
           else
           if(f>=.0) then
           dq1=.0; dq2=1.5/dsbar; dq3=.0     ;   call formm(stress,m1,m2,m3)
           flow=f*(m1*dq1+m2*dq2+m3*dq3)     ;   erate=matmul(flow,stress)
           evp=erate*dt; evpt(:,i,iel)=evpt(:,i,iel)+evp; devp=matmul(dee,evp) 
         end if; end if
      if(f>=.0) then
        eload=matmul(transpose(bee),devp) ; bload=bload+eload*det*weights(i)
      end if
      if(plastic_converged.or.plasiters==plasits)then                          
!---------------------- update the Gauss Point stresses  ----------------------
          tensor(: , i , iel) = stress          
      end if
    end do gauss_points_2
!-------------compute the total bodyloads vector ; dependency if vectorised----
!dir$ ivdep
    do i=1,ndof
         bdylds( g(i) ) = bdylds( g(i) ) + bload(i)
    end do  ;    bdylds(0) = .0
  end do elements_4             
  if(plastic_converged.or.plasiters==plasits)exit
 end do plastic_iterations
 totd=totd+loads
 write(11,'(a,e12.4)')"The total load is  ",ptot 
 write(11,'(a,10e12.4)')"Displacements are",(totd(nf(2,no(i))),i=1,loaded_nodes)
 write(11,'(a,i12)')"The total number of cj iterations was    ",cjtot 
 write(11,'(a,i12)')"The number of plastic iterations was     ",plasiters 
 write(11,'(a,f11.2)')"cj iterations per plastic iteration were   ", &
  & real(cjtot)/real(plasiters)
 if(plasiters==plasits)stop
end do load_increments
end program p61    

