    program p63      
!-----------------------------------------------------------------------
!      program 6.3 plane strain of an elastic-plastic(Mohr-Coulomb) solid
!      using 8-node quadrilateral elements; initial stress method
!------------------------------------------------------------------------
 use new_library ;   use geometry_lib  ; implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,          &
          loaded_nodes, i,k,iel,iters,limit,incs,iy,ndim=2
 logical::converged       ;  character (len=15) :: element='quadrilateral'
 real::e,v,det,c,phi,psi,epk0,ptot,fac,f,fnew,dsbar,lode_theta,sigm,presc,    &
       pav,gama,tol                                                            
!----------------------------- dynamic arrays----------------------------------
 real    ,allocatable :: kb(:,:),loads(:),points(:,:),totd(:),bdylds(:),      &
                         width(:),depth(:),tensor(:,:,:),                     &
                         dee(:,:),coord(:,:),fun(:),jac(:,:),weights(:),      &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),elso(:),g_coord(:,:),     &
                         oldis(:),storkb(:),stress(:),pl(:,:),gc(:)
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)   
!------------------------input and initialisation------------------------------
  open (10,file='p63.dat',status=    'old',action='read')
  open (11,file='p63.res',status='replace',action='write')                    
  read (10,*) phi,c,psi,gama,epk0,e,v,    nels,nxe,nye,nn,nip
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            width(nxe+1),depth(nye+1),num(nod),dee(nst,nst),fun(nod),gc(ndim),&
            tensor(nst,nip,nels),g_g(ndof,nels),coord(nod,ndim),stress(nst),  &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),pl(nst,nst),elso(nst),g(ndof))              
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)     ;  read(10,*) width, depth               
!----------- loop the elements to find nband and set up global arrays ---------
     nband = 0
       elements_1:   do iel = 1 , nels
                        call geometry_8qyv(iel,nye,width,depth,coord,num)
                        call num_to_g(num,nf,g);          g_num(:,iel)=num
                        g_coord(:,num)=transpose(coord);  g_g(:,iel) = g
                        if (nband<bandwidth(g)) nband = bandwidth(g)
      end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do  
   write(11,'(a,i5,a,i5)')                                                     &
           "The system has ",neq,"  equations and the half-bandwidth is ",nband
  allocate(kb(neq,nband+1),loads(0:neq),bdylds(0:neq),oldis(0:neq),totd(0:neq)) 
  kb=0.0; oldis=0.0; totd=0.0 ;bdylds=.0;  tensor = 0.0
  call deemat(dee,e,v);    call sample(element,points,weights)                
!---------------- element stiffness integration and assembly-------------------
 elements_2: do iel = 1 , nels
                num = g_num(:, iel) ; coord =transpose( g_coord(: ,num))
                g = g_g( : , iel )  ;        km=0.0
             gauss_pts_1:  do i =1 , nip
               call shape_fun(fun,points,i);  gc = matmul(fun,coord)
!--------------------starting stress state  -----------------------------------
          tensor(2,i,iel)=gc(2)*gama;tensor(1,i,iel)=gc(2)*gama*epk0
          tensor(4,i,iel)=tensor(1,i,iel) 
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
            end do gauss_pts_1   
   call formkb (kb,km,g)
 end do elements_2
!---------------read prescribed displacements and factorise l.h.s.-------------
    read(10,*) loaded_nodes ; allocate (no(loaded_nodes),storkb(loaded_nodes)) 
    read(10,*)(no(i),i=1,loaded_nodes),presc,incs,tol,limit             
    do i=1,loaded_nodes
            kb(nf(1,no(i)),nband+1)=kb(nf(1,no(i)),nband+1)+1.e20 
            storkb(i) = kb(nf(1,no(i)),nband+1)
    end do;         call cholin(kb)                    
!-------------------displacement increment loop--------------------------------
      write(11,'(a)') "Displacement    Force         Iterations"
     displacement_increments: do iy=1,incs
            ptot=presc*iy; iters=0      
!-------------------------   iteration loop  ----------------------------------
   iterations: do
    iters=iters+1;  loads=.0   ;loads = loads + bdylds    
      do i=1,loaded_nodes ; loads(nf(1,no(i)))=storkb(i)*presc ;  end do
        call chobac(kb,loads)  ; bdylds = .0
!-------------------------   check convergence  -------------------------------
      call checon(loads,oldis,tol,converged)
      if(iters==1)converged=.false. 
!------------------------- go round the Gauss Points --------------------------
      elements_3: do iel = 1 , nels
       bload=.0
       num = g_num( : ,iel ) ; coord = transpose(g_coord( : , num ))
       g = g_g( : , iel )    ;    eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          elso = .0 ; call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ;   call invert(jac)
          deriv = matmul(jac,der) ; call beemat (bee,deriv)
          eps = matmul(bee,eld);    sigma = matmul(dee,eps)
          stress = sigma + tensor(:,i,iel)
          call invar(stress,sigm,dsbar,lode_theta)                            
!-------------------  check whether yield is violated -------------------------
         call mocouf(phi,c,sigm,dsbar,lode_theta,fnew)
         if (fnew>.0) then
           stress = tensor(:,i,iel) ;  call invar(stress,sigm,dsbar,lode_theta)
           call mocouf(phi,c,sigm,dsbar,lode_theta,f) 
           fac = fnew / (fnew - f);    stress = (1.-fac)*sigma+tensor(:,i,iel)
           call mocopl(phi,psi,e,v,stress,pl); pl = fac * pl
           elso =  matmul(pl,eps) ;  eload = matmul(elso,bee)
           bload = bload + eload * det * weights(i)
         end if
       if(converged.or.iters==limit)then
!-------------------    update the Gauss Point stresses  ----------------------
          tensor(:,i,iel) = tensor(:,i,iel) + sigma - elso
       end if
    end do gauss_points_2
!---------------------- compute the total bodyloads vector --------------------
    bdylds( g ) = bdylds( g ) + bload  ; bdylds(0) = .0    
  end do elements_3             
  if(converged.or.iters==limit)exit
 end do iterations
 totd=totd+loads  
 pav = .5*((depth(1)-depth(2))*(tensor(1,1,1)+tensor(1,3,1)) &
          +(depth(2)-depth(3))*(tensor(1,1,2)+tensor(1,3,2)) &
          +(depth(3)-depth(4))*(tensor(1,1,3)+tensor(1,3,3)) &
          +(depth(4)-depth(5))*(tensor(1,1,4)+tensor(1,3,4)))    
 write(11,'(2e12.4,i12)') ptot,pav,iters 
 if(iters==limit)stop
end do displacement_increments
end program p63 
