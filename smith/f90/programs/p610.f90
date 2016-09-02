   program p610      
!-----------------------------------------------------------------------
!      program 6.10 plane strain of an elastic-plastic(Mohr-Coulomb) solid
!      using 8-node quadrilateral elements; viscoplastic strain method
!      construction of an excavation in stages introducing air elements
!------------------------------------------------------------------------
 use new_library      ; use geometry_lib   ;     implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,          &
          i,k,iel,iters,limit,incs,iy,ndim=2,layers,ii,iq,nofix,noexe
 logical::converged          ; character (len=15) :: element='quadrilateral'
 real::es,ea,det,phi,psi,c,gama,                                              &
       dt,f,dsbar,dq1,dq2,dq3,lode_theta,sigm,pi,snph,e,v,epk0,tol             
!------------------------------ dynamic arrays---------------------------------
 real    ,allocatable :: kb(:,:),loads(:),points(:,:),bdylds(:),disps(:,:),   &
                         evpt(:,:,:),oldis(:),width(:),depth(:),dloads(:),    &
                         dee(:,:),coord(:,:),fun(:),jac(:,:),weights(:),      &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),    &
                         totd(:),tensor(:,:,:),gc(:),s(:),lod(:,:)   
 integer, allocatable :: nf(:,:) , g(:),  num(:), g_num(:,:) ,                &
                         prop(:) , lnf(:,:) , fixnod(:) , exele(:)             
!------------------------input and initialisation------------------------------
  open (10,file='p610.dat',status=    'old',action='read')
  open (11,file='p610.res',status='replace',action='write')                  
  read (10,*) nxe,nye,nn,nip,incs,limit,tol,layers,epk0,                      &
              ea,es,v,c,phi,psi,gama 
  ndof=nod*nodof   ; pi = acos( -1. )    ; snph = sin(phi*pi/180.)  
 ! calculate the total number of elements
   nels = nxe*nye;   write(11,'(a,i5)')"The total number of elements is ",nels
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            depth(nye+1),num(nod),dee(nst,nst),evpt(nst,nip,nels),            &
            width(nxe+1),coord(nod,ndim),fun(nod),prop(nels),                 &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst),s(nst),         &
            disps(nodof,nn),gc(ndim),tensor(nst,nip,nels),lnf(nodof,nn))
!------ nf is an index array of 1s and 0s : lnf is the local nf  --------------
            nf=1; read(10,*) nr; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
            lnf = nf;   call formnf(lnf); neq=maxval(lnf)
            write(11,'(a,i5)')"The total possible number of equations is:",neq
  allocate(totd(0:neq))         ;  read(10,*)  width , depth    
!--------------------- set the element type -----------------------------------
   prop = 1  
!-------- set up the global node numbers and global nodal coordinates ---------
   call fmglob(nxe,nye,g_num)
   call fmcoco(g_num,g_coord,width,depth,nxe,nye)
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do       
    disps = .0;     call sample(element,points,weights)
!--------------- loop the elements to set starting stresses -------------------
        elements_1:   do iel = 1 , nels     
                       num = g_num(:,iel); coord=transpose(g_coord(:,num))
                gauss_points_1: do i=1,nip
                       call shape_fun(fun,points,i); gc = matmul(fun,coord)  
                       tensor(2,i,iel)=gc(2)*gama
                       tensor(1,i,iel)=epk0*tensor(2,i,iel)
                       tensor(4,i,iel)=tensor(1,i,iel) ;  tensor(3,i,iel)=.0
               end do gauss_points_1        
        end do elements_1        
! -----------------------  excavate another layer -----------------------------
  layer_number : do ii = 1 , layers ;write(11,'(a,i5)') "Layer no",ii
!--------- recalculate the number of freedoms neq and half-bandwidth nband ----
          read(10,*)nofix ;   allocate(fixnod(nofix))   
             read(10,*)fixnod ; nf(:,fixnod) = 0 ; lnf = nf 
             call formnf(lnf) ;   neq = maxval(lnf)
             nband = 0
             elements_1a : do iel = 1 , nels
                          num = g_num(:,iel);call num_to_g(num,lnf,g)
                          if(nband<bandwidth(g))nband = bandwidth(g)
             end do elements_1a
      write(11,'(/,3(a,i5))')                                                  &
              "There are ",neq, " freedoms and nband is",nband," in step",ii
allocate(kb(neq,nband+1),loads(0:neq),bdylds(0:neq),oldis(0:neq),dloads(0:neq)) 
           kb=0.0;   loads = .0
!-------------------  specify the elements to be removed ----------------------
    read(10,*)noexe  ; allocate(exele(noexe),lod(ndof,noexe))  
    read(10,*)exele; prop(exele)=0                                             
!----------------------  calculate excavation load ----------------------------
    s = .0;  lod = .0
    elements_2 : do iel = 1 , noexe
       iq = exele(iel) ; bload = .0; eld = .0 ;  num = g_num(:,iq)
       call num_to_g(num,lnf,g);coord = transpose(g_coord(:,num))
       gauss_points_2 : do i = 1 , nip 
          call shape_fun(fun,points,i) 
          call shape_der (der,points,i);  jac = matmul(der,coord) 
          det = determinant(jac)  ;   call invert(jac)
          deriv = matmul(jac,der) ;   call beemat (bee,deriv)           
          s = tensor(:,i,iq) ; eload = matmul(s,bee)
          bload = bload + eload * det * weights(i)
          do k=2,ndof,2;eld(k)=eld(k)+fun(k/2)*det*weights(i)*gama;end do
       end do gauss_points_2
          lod(:,iel) = eld + bload; loads(g)=loads(g)+lod(:,iel);loads(0)=.0
    end do elements_2
!------------------ element stiffness integration and assembly-----------------
 elements_3: do iel = 1 , nels
              if(prop(iel)==0) e = ea
              if(prop(iel)==1) e = es
              km = .0; eld = .0; call deemat(dee,e,v)
              num=g_num(:,iel);call num_to_g(num,lnf,g)
              coord = transpose(g_coord(: , num))
              gauss_points_3: do i = 1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;   call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
              end do gauss_points_3   
   call formkb (kb,km,g)
 end do elements_3 
!--------------- factorise l.h.s.  ------------------------------------------- 
          call cholin(kb)                                                  
!------------------ factor excavation load by incs----------------------------
     write(11,'(a,i5,a,e12.4)')                                                &
             "The total gravity load in lift", ii, " is" , sum(loads)
       loads = loads / incs
!------------ apply excavation loads incrementally ----------------------------
    load_increments: do iy=1,incs
              iters=0; oldis =.0; bdylds=.0; evpt=.0
!--------------------------   iteration loop ----------------------------------
   iterations: do
    iters=iters+1;  dloads = .0; dloads= loads+bdylds ;  call chobac(kb,dloads)
!--------------------------   check convergence -------------------------------
      call checon(dloads,oldis,tol,converged)
      if(iters==1)converged=.false. 
      if(converged.or.iters==limit) then
        bdylds=.0  
          do iq=1,nn;do i=1,nodof     
           if(lnf(i,iq)/=0)disps(i,iq) = disps(i,iq) + dloads(lnf(i,iq))
          end do;end do
      end if
!------------------------ go round the Gauss Points ---------------------------
 elements_4: do iel = 1 , nels
              if(prop(iel)==0) then;e = ea;dt=1.e10; end if
              if(prop(iel)==1) then
                e = es;dt=(4.*(1.+v)*(1.-2.*v))/(e*(1.-2.*v+snph*snph))
              end if
              bload = .0; call deemat(dee,e,v)
              num=g_num(:,iel);call num_to_g(num,lnf,g)
              coord = transpose(g_coord(: , num)) ;   eld = dloads(g)
              gauss_points_4: do i = 1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               eps = matmul(bee,eld); eps=eps-evpt(:,i,iel)
               s = tensor(:,i,iel) + matmul(dee,eps)
!---------------------- air element stresses are zero -------------------------
               if(prop(iel)==0) s = .0
               call invar(s,sigm,dsbar,lode_theta)
!--------------------  check whether yield is violated ------------------------
         call mocouf (phi, c , sigm, dsbar , lode_theta , f )
         if(converged.or.iters==limit) then
         devp=s 
           else
           if(f>=.0) then
           call mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3);call formm(s,m1,m2,m3)
           flow=f*(m1*dq1+m2*dq2+m3*dq3)     ;   erate=matmul(flow,s)
           evp=erate*dt; evpt(:,i,iel)=evpt(:,i,iel)+evp; devp=matmul(dee,evp) 
         end if; end if
      if(f>=.0) then
        eload=matmul(devp,bee)   ; bload=bload+eload*det*weights(i)
      end if
!--------------- if appropriate update the Gauss point stresses ---------------
      if(converged.or.iters==limit) tensor(:,i,iel) = s
    end do gauss_points_4
!-------------------   compute the total bodyloads vector ---------------------
    bdylds( g ) = bdylds( g ) + bload      ; bdylds(0) = .0
  end do elements_4             
  if(converged.or.iters==limit)exit
 end do iterations
 write(11,'(2(a,i5))') "Lift number",ii," gravity load increment",iy
 write(11,'(a,i5,a)') "It took ",iters, "  iterations to converge"
 if(iy==incs.or.iters==limit)then
     write(11,'(a)') "The displacements are :"
     write(11,'(5e12.4)')disps(1,61),disps(2,61),disps(1,37),disps(2,37) 
 end if
 if(iters==limit)stop
end do load_increments
deallocate(kb,loads,bdylds,oldis,dloads,fixnod,exele,lod)
end do layer_number
end program p610
