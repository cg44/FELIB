    program p65                   
!-----------------------------------------------------------------------
!      program 6.5 plane strain of an elastic-plastic(Von Mises)
!      solid using 8-node quadrilateral elements; initial stress method
!      consistent return algorithm for problem of p60 ; pcg version
!------------------------------------------------------------------------
 use new_library  ; use geometry_lib    ;   implicit none
 integer::nels,nxe,nye,neq,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,loaded_nodes,  &
          i,k,iel,plasiters,plasits,cjiters,cjits,incs,iy,ndim=2,cjtot
 logical::plastic_converged,cj_converged
 character(len=15)::element='quadrilateral'
 real::e,v,det,cu,ptot,fnew,ff,fstiff,dlam,dslam,dsbar,lode_theta,sigm,top,  &
       bot,tload,tloads,residual,up,alpha,beta,big,plastol,cjtol,fftol,ltol
!------------------------------- dynamic arrays--------------------------------
 real    ,allocatable :: loads(:),points(:,:),totd(:),bdylds(:),              &
                         width(:),depth(:),tensor(:,:,:),pmul(:),utemp(:),    &
                         dee(:,:),coord(:,:),fun(:),jac(:,:),weights(:),      &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),elso(:),g_coord(:,:),     &
                         oldis(:),val(:,:),stress(:),qinc(:),ddylds(:),       &
                         dl(:,:),dload(:),vmfl(:),caflow(:),dsigma(:),        &
                         ress(:),rmat(:,:),acat(:,:),acatc(:,:),qmat(:,:),    &
                         qinva(:),daatd(:,:),temp(:,:),vmflq(:),storkm(:,:,:),&
                         vmfla(:),qinvr(:),vmtemp(:,:),                       &
                         p(:),x(:),xnew(:),u(:),diag_precon(:),d(:)
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)    
!------------------------input and initialisation------------------------------
  open (10,file='p65.dat',status=    'old',action='read')
  open (11,file='p65.res',status='replace',action='write')                
  read (10,*) cu,e,v, nels,nxe,nye,nn,nip,                                   &
              plasits,cjits,plastol,cjtol,fftol,ltol
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            width(nxe+1),depth(nye+1),num(nod),dee(nst,nst),pmul(ndof),       &
            tensor(nst,nip,nels),g_g(ndof,nels),storkm(ndof,ndof,nels),       &
            coord(nod,ndim),stress(nst),temp(nst,nst),dl(nip,nels),           &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),elso(nst),g(ndof),vmfl(nst),qinvr(nst),   &
            dload(ndof),caflow(nst),dsigma(nst),ress(nst),rmat(nst,nst),      &
            acat(nst,nst),acatc(nst,nst),qmat(nst,nst),qinva(nst),            &
            daatd(nst,nst),vmflq(nst),vmfla(nst),vmtemp(1,nst),utemp(ndof))    
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)   ;  read(10,*) width, depth
  temp = .0; temp(1,1)=1.;temp(2,2)=1.;temp(3,3)=3.;temp(4,4)=1.
  temp(1,2)=-.5;temp(2,1)=-.5;temp(1,4)=-.5
  temp(2,4)=-.5;temp(4,1)=-.5;temp(4,2)=-.5
!---------------- loop the elements to set up global arrays -------------------
       elements_1:   do iel = 1 , nels
                        call geometry_8qyv(iel,nye,width,depth,coord,num)
                        call num_to_g(num,nf,g);          g_num(:,iel)=num
                        g_coord(:,num)=transpose(coord);  g_g( : , iel ) = g
       end do elements_1  
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do        
  allocate(loads(0:neq),bdylds(0:neq),oldis(0:neq),totd(0:neq),ddylds(0:neq), &
           p(0:neq),x(0:neq),xnew(0:neq),u(0:neq),diag_precon(0:neq),d(0:neq)) 
  oldis=0.0; totd=0.0 ; p = .0; xnew = .0; diag_precon = .0
  call deemat(dee,e,v); call sample(element,points,weights); tensor= .0; dl=.0 
!---------- starting element stiffness integration,storage,preconditioner------
 elements_2: do iel = 1 , nels
                num = g_num(: , iel ) ; coord = transpose(g_coord(: , num )) 
                g = g_g( : , iel )    ;     km=0.0
             gauss_pts_1:  do i =1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
            end do gauss_pts_1   
   storkm(:,:,iel) = km
   do k=1,ndof; diag_precon(g(k))=diag_precon(g(k)) + km(k,k); end do
 end do elements_2
   diag_precon(1:neq) = 1./diag_precon(1:neq) ; diag_precon(0) = .0
!---------------read load weightings ----------------------------------------- 
  read(10,*) loaded_nodes ; allocate (no(loaded_nodes),val(loaded_nodes,ndim))
  read(10,*)(no(i),val(i,:),i=1,loaded_nodes)   
!----------------- load increment loop----------------------------------------
      read(10,*) incs ; allocate(qinc(incs));  read(10,*)qinc ;   ptot = .0
     load_increments: do iy=1,incs
         write(11,'(/,a,i5)') ' Load increment  ',iy
         ptot=ptot + qinc(iy) ; plasiters = 0;bdylds=.0 ;loads=.0; cjtot=0
      do i=1,loaded_nodes ; loads(nf(:,no(i)))=val(i,:)*qinc(iy) ;  end do
!----------------------   plastic iteration loop  -----------------------------
   plastic_iterations: do
    plasiters=plasiters+1;  if(plasiters/=1)loads=.0   ;loads = loads + bdylds
    if(abs(sum(loads))<1.e-5) then; plasiters=plasiters-1;exit ; end if 
       write(11,'(a,i5)') "Plastic iteration number",plasiters  
        bdylds = .0    ; ddylds = .0 ; d=diag_precon*loads ; p = d   ; x = .0
!---------------   solve the simultaneous equations by pcg  -------------------
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
!----------------------------pcg process --------------------------------------
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
!------------------------- end of pcg process ---------------------------------
    loads = xnew           ; loads(0) = .0    ; diag_precon = .0
!----------------------- go round the elements---------------------------------
      elements_4: do iel = 1 , nels
       bload=.0   ;   dload = .0
       num = g_num( : , iel ) ; coord = transpose(g_coord( : , num ))
       g = g_g( : , iel )     ;    km = .0; eld = loads(g)   
!------------------------ go round the Gauss points --------------------------- 
       gauss_points_2 : do i = 1 , nip
          elso = .0 ; call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ;   call invert(jac)
          deriv = matmul(jac,der) ;   call beemat (bee,deriv)
          eps = matmul(bee,eld)   ;   call deemat(dee,e,v)
          stress = tensor(: , i , iel)
          call invar(stress,sigm,dsbar,lode_theta) ; ff = dsbar-sqrt(3.)*cu
           if(ff>fftol) then
             dlam = dl(i,iel) ; call vmflow(stress,dsbar,vmfl)
             call fmrmat(vmfl,dsbar,dlam,dee,temp,rmat)
             caflow = matmul(rmat,vmfl); bot=dot_product(vmfl,caflow)
             call formaa(vmfl,rmat,daatd); dee = rmat - daatd/bot
           end if
          sigma = matmul(dee,eps)
          stress = sigma + tensor( : , i , iel)
          call invar(stress,sigm,dsbar,lode_theta)                            
!---------------------  check whether yield is violated  ----------------------
          fnew = dsbar - sqrt(3.)*cu ; fstiff = fnew
          if (fnew>=.0) then
              call deemat(dee,e,v) ; call vmflow(stress,dsbar,vmfl)
              caflow = matmul(dee,vmfl); bot=dot_product(vmfl,caflow)
              dlam = fnew/bot; elso = caflow*dlam 
              stress = tensor( : , i , iel) + sigma - elso
              call invar(stress,sigm,dsbar,lode_theta);fnew=dsbar-sqrt(3.)*cu
             iterate_on_fnew : do
               call vmflow(stress,dsbar,vmfl); caflow = matmul(dee,vmfl)*dlam
               ress = stress - (tensor(: , i , iel) +sigma - caflow)
               call fmacat(vmfl,temp,acat); acat = acat / dsbar
               acatc = matmul(dee,acat); qmat = acatc*dlam
               do k=1,4; qmat(k,k)=qmat(k,k)+1.; end do; call invert(qmat)
               vmtemp(1,:)=vmfl; vmtemp = matmul(vmtemp,qmat);vmflq=vmtemp(1,:)
               top = dot_product(vmflq,ress)
               vmtemp = matmul(vmtemp,dee);vmfla=vmtemp(1,:) 
               bot = dot_product(vmfla,vmfl) ; dslam = (fnew - top)/bot
               qinvr = matmul(qmat,ress); qinva=matmul(matmul(qmat,dee),vmfl)
               dsigma = -qinvr - qinva*dslam; stress = stress + dsigma
               call invar(stress,sigm,dsbar,lode_theta)
               fnew = dsbar - sqrt(3.)*cu;  dlam = dlam + dslam
               if (fnew<plastol) exit
             end do iterate_on_fnew
             dl(i,iel) = dlam
             elso = tensor( : , i , iel) + sigma - stress
             eload=matmul(elso,bee);bload=bload+eload*det*weights(i)
             call vmflow(stress,dsbar,vmfl)
             call fmrmat(vmfl,dsbar,dlam,dee,temp,rmat)
             caflow=matmul(rmat,vmfl);bot = dot_product(vmfl,caflow)
             call formaa(vmfl,rmat,daatd)
             dee = rmat - daatd/bot
          end if
          if(fstiff<.0) call deemat(dee,e,v)
              km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
!---------------------  update the Gauss Point stresses -----------------------
          tensor( : , i , iel) = tensor( : , i , iel) + sigma - elso
          stress = tensor ( : , i , iel)
          eload=matmul(stress,bee); dload=dload+eload*det*weights(i)
    end do gauss_points_2
!         compute the total bodyloads vector
    bdylds( g ) = bdylds( g ) + bload  ; bdylds(0) = .0   
    ddylds( g ) = ddylds( g ) + dload  ; ddylds(0) = .0
    storkm(:,:,iel) = km
   do k =1,ndof; diag_precon(g(k))=diag_precon(g(k))+km(k,k); end do 
  end do elements_4 
    diag_precon(1:neq)=1./diag_precon(1:neq) ; diag_precon(0) = .0
          tload = sum(ddylds)  ; tloads = sum(bdylds)
   if(plasiters==1)plastic_converged=.false.
   if(plasiters/=1.and.tloads<ltol)plastic_converged=.true.
   residual = (2.*ptot+tload-tloads)/(2.*ptot)
   write(11,'(a,10e12.4)')"tloads,tload,residual are",tloads,tload,residual
   totd = totd + loads             
  if(plastic_converged.or.plasiters==plasits)exit
 end do plastic_iterations
 totd=totd+loads  
 write(11,'(a,e12.4)')"The total load is   ",ptot
 write(11,'(a,10e12.4)')"Displacements are",totd(nf(2,no))
 write(11,'(a,i12)')"The number of cj iterations was    ",cjtot
 write(11,'(a,i12)')"The number of plastic iterations was  ",plasiters
 write(11,'(a,f11.2)')"cj iterations per plastic iteration were   ",          &
       &  real(cjtot)/real(plasiters)
 if(plasiters==plasits)stop
end do load_increments
end program p65
