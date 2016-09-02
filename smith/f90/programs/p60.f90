    program p60      
!-----------------------------------------------------------------------
!      program 6.0 plane strain of an elastic-plastic(Von Mises) solid
!      using 8-node quadrilateral elements; viscoplastic strain method
!------------------------------------------------------------------------
 use new_library      ; use geometry_lib  ;     implicit none
 integer::nels,nxe,nye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,         &
          loaded_nodes,i,k,iel,iters,limit,incs,iy,ndim=2
 character (len = 15) :: element = 'quadrilateral';   logical::converged
 real::e,v,det,cu,dt,ptot,f,dsbar,dq1,dq2,dq3,lode_theta,sigm,tol
!--------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: kb(:,:),loads(:),points(:,:),totd(:),bdylds(:),      &
                         evpt(:,:,:),oldis(:),width(:),depth(:),              &
                         tensor(:,:,:),val(:,:),stress(:),qinc(:),            &
                         dee(:,:),coord(:,:),jac(:,:),weights(:),             &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:)
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)    
!----------------input and initialisation----------------------
  open (10,file='p60.dat',status=    'old',action='read')
  open (11,file='p60.res',status='replace',action='write')                     
  read (10,*) cu,e,v,nels,nxe,nye,nn,nip,tol,limit
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            width(nxe+1),depth(nye+1),num(nod),dee(nst,nst),g_g(ndof,nels),   &
            evpt(nst,nip,nels), tensor(nst,nip,nels),coord(nod,ndim),         &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst),stress(nst))     
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
   call formnf(nf); neq=maxval(nf)     ;    read(10,*) width, depth           
!---------- loop the elements to find nband and set up global arrays ----------
     nband = 0
       elements_1:   do iel = 1 , nels
                        call geometry_8qyv(iel,nye,width,depth,coord,num)
                        call num_to_g(num,nf,g); g_num(:,iel)=num
                        g_coord(:,num)=transpose(coord);  g_g(:,iel) = g
                        if (nband<bandwidth(g)) nband = bandwidth(g)
      end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k=1,nels; write(11,'(a,i5,a,8i5)')                                      &
                         "Element ",k,"        ",g_num(:,k); end do  
    write(11,'(a,i5,a,i5)')                                                    &
            "There are ",neq ,"   equations and the half-bandwidth is ",nband
  allocate(kb(neq,nband+1),loads(0:neq),bdylds(0:neq),oldis(0:neq),totd(0:neq)) 
  kb=0.0; oldis=0.0; totd=0.0 ; tensor = 0.0
  call deemat(dee,e,v); call sample(element,points,weights)
  dt=4.*(1.+v)/(3.*e)
!------------------ element stiffness integration and assembly-----------------
 elements_2: do iel = 1 , nels
                num = g_num(:,iel) ; coord = transpose(g_coord(: ,num)) 
                g = g_g(: ,iel ) ;          km=0.0
             gauss_pts_1:  do i =1 , nip
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
            end do gauss_pts_1   
   call formkb (kb,km,g)
 end do elements_2                                                            
!------------------read load weightings and factorise equations---------------- 
    read(10,*)loaded_nodes; allocate(no(loaded_nodes),val(loaded_nodes,ndim))
    read(10,*)(no(i),val(i,:),i=1,loaded_nodes)               
    call cholin(kb)                                                  
!-------------------load increment loop---------------------------------------
    read(10,*)incs;  allocate(qinc(incs));  read(10,*) qinc
    ptot=.0
   load_increments: do iy=1,incs
   write(11,'(a,i5)') 'increment',iy 
    ptot=ptot+qinc(iy) ; iters=0;  bdylds=.0;  evpt=.0
!--------------------   iteration loop  ---------------------------------------
   iterations: do
    iters=iters+1;  loads=.0
      do i=1,loaded_nodes ; loads(nf(:,no(i)))=val(i,:)*qinc(iy) ;  end do
      loads=loads+bdylds     ;    call chobac(kb,loads)
!----------------------   check convergence  ----------------------------------
      call checon(loads,oldis,tol,converged)
      if(iters==1)converged=.false. ;  if(converged.or.iters==limit)bdylds=.0
!------------------- go round the Gauss Points --------------------------------
      elements_3: do iel = 1 , nels
       bload=.0
       num = g_num( : ,iel ) ; coord = transpose (g_coord( : , num)) 
       g = g_g( : , iel ) ;        eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ;   call invert(jac)
          deriv = matmul(jac,der) ; call beemat (bee,deriv);eps=matmul(bee,eld)
          eps = eps - evpt( : ,i , iel) ;   sigma=matmul(dee,eps)
          stress = sigma + tensor( : ,i , iel) 
          call invar(stress,sigm,dsbar,lode_theta)                             
!-----------------  check whether yield is violated  --------------------------
         f=dsbar-sqrt(3.)*cu
         if(converged.or.iters==limit) then
         devp=stress 
           else
           if(f>=.0) then
           dq1=.0; dq2=1.5/dsbar; dq3=.0     ;   call formm(stress,m1,m2,m3)
           flow=f*(m1*dq1+m2*dq2+m3*dq3)     ;   erate=matmul(flow,stress)
           evp = erate*dt; evpt(:,i,iel)=evpt(:,i,iel) + evp
           devp=matmul(dee,evp) 
         end if; end if
      if(f>=.0) then
        eload=matmul(transpose(bee),devp) ; bload=bload+eload*det*weights(i)
      end if
      if(converged.or.iters==limit)then                                       
!-------------------------- update the Gauss Point stresses -------------------
          tensor( : , i , iel) = stress
      end if
    end do gauss_points_2
!         compute the total bodyloads vector
    bdylds( g ) = bdylds( g ) + bload      ; bdylds(0) = .0
  end do elements_3             
  if(converged.or.iters==limit)exit
 end do iterations
 totd=totd+loads
 write(11,'(a,e12.4)') "The total load is ", ptot 
 write(11,'(a,10e12.4)')"Displacements are",                                  &
                        (totd(nf(2,no(i))),i=1,loaded_nodes)
 write(11,'(a,i5,a)')"It took ",iters ,"   iterations to converge"
 if(iters==limit)stop
end do load_increments
end program p60 
