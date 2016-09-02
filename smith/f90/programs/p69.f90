    program p69      
!-----------------------------------------------------------------------
!      program 6.9 plane strain of an elastic-plastic(Mohr-Coulomb) solid
!      using 8-node quadrilateral elements; viscoplastic strain method
!      construction of an embankment in layers on a foundation
!------------------------------------------------------------------------
 use new_library      ; use geometry_lib   ; implicit none
 integer::nels,lnxe,lnye,neq,nband,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,       &
          i,k,iel,iters,limit,incs,iy,ndim=2,fnxe,fnye,lifts,oldele,newele,  &
          oldnn,lnn,ii,itype
 logical::converged   ; character (len=15) :: element = 'quadrilateral'
 real::ef,es,vf,vs,det,phif,phis,cf,cs,psif,psis,gamaf,gamas,dt,f,dsbar,     &
       dq1,dq2,dq3,lode_theta,sigm,pi,snph,e,v,c,phi,psi,gama,epk0,tol         
!---------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: kb(:,:),loads(:),points(:,:),bdylds(:),              &
                         evpt(:,:,:),oldis(:),width(:),depth(:),gravlo(:),    &
                         dee(:,:),coord(:,:),fun(:),jac(:,:),weights(:),      &
                         der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),  &
                         sigma(:),bload(:),eload(:),erate(:),g_coord(:,:),    &
                         evp(:),devp(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),    &
                         totd(:),fwidth(:),fdepth(:),tensor(:,:,:),gc(:),s(:)
 integer, allocatable :: nf(:,:) , g(:),  num(:), g_num(:,:) ,g_g(:,:),       &
                         prop(:) , lnf(:,:)                                    
!------------------------input and initialisation------------------------------
  open (10,file='p69.dat',status=    'old',action='read')
  open (11,file='p69.res',status='replace',action='write')               
  read (10,*) fnxe,fnye,nn,nip,incs,limit,tol,                                &
              lifts,lnxe,lnye,itype,epk0,                                     &
              ef,vf,cf,phif,psif,gamaf,                                       &
              es,vs,cs,phis,psis,gamas 
  ndof=nod*nodof   ; pi = acos( -1. )      
!------------------ calculate the total number of elements --------------------
   k=0;do i=1,lnye-1;k=i+k;end do ; nels=fnxe*fnye+(lnxe*lnye-k)
         write(11,'(a,i5)') "The total number of elements is ",nels
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            depth(lnye+1),num(nod),dee(nst,nst),evpt(nst,nip,nels),           &
            width(lnxe+1),coord(nod,ndim),fun(nod),prop(nels),g_g(ndof,nels), &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),km(ndof,ndof),eld(ndof),eps(nst),sigma(nst),        &
            bload(ndof),eload(ndof),erate(nst),evp(nst),devp(nst),g(ndof),    &
            m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst),s(nst),         &
            fwidth(fnxe+1),fdepth(fnye+1),gc(ndim),tensor(nst,nip,nels))       
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)
         write(11,'(a,i5)') "The final number of equations is:",neq
  allocate(totd(0:neq))
  read(10,*) fwidth , fdepth , width , depth
!------------------------- set the element type  ------------------------------
   prop(1:fnxe*fnye)=1 ; prop(fnxe*fnye+1:nels)=2 
!----------- set up the global node numbers and element nodal coordinates -----
   call fmglem(fnxe,fnye,lnxe,1,g_num,lifts)
   call fmcoem(g_num,g_coord,fwidth,fdepth,width,depth,      &
               lnxe,lifts,fnxe,fnye,itype)   
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                  &
                             "Element ",k,"        ",g_num(:,k); end do   
    tensor = .0; totd = .0; call sample(element,points,weights)
!------------- loop the elements to find the global g -------------------------
        elements_1:   do iel = 1 , nels    ; num = g_num(:,iel)
                       call num_to_g (num,nf,g) ;   g_g(:,iel) = g
        end do elements_1
! -------------------  construct another lift ---------------------------------
  lift_number : do ii = 1 , lifts
! -------------      calculate how many elements there are --------------------
    if (ii<=lifts) then
        if(ii==1) then
           newele=fnxe*fnye; oldele = newele
        else
           newele = lnxe - (ii -2); oldele = oldele + newele
        end if
!--------- go round the elements and get nband from the g vectors -------------
       nband = 0
       elements_2 : do iel = 1 , oldele
                       g=g_g( : , iel )
                       if(nband<bandwidth(g)) nband = bandwidth( g )
       end do elements_2
! --------------      calculate how many nodes there are   --------------------
        if(ii==1) then
           lnn=(fnxe*2+1)*(fnye+1)+(fnxe+1)*fnye ; oldnn = lnn
        end if
        if(ii>1) then
           lnn=oldnn+(lnxe-(ii-2))*2+1+(lnxe-(ii-2)+1) ; oldnn = lnn
        end if
!------------------- now get the new node freedom array -----------------------
       allocate(lnf(nodof,lnn))     ; lnf = nf(:,1:lnn)
!-----------------  recalculate the number of freedoms neq --------------------
          neq = maxval(lnf)
      write(11,'(/,3(a,i5))')                                                  &
              "There are",neq," freedoms and",lnn," nodes in lift",ii
      write(11,'(a,i5,a,i5,a)')                                                &
              "There are ",oldele," elements and",newele," were added"
    end if
allocate(kb(neq,nband+1),loads(0:neq),bdylds(0:neq),oldis(0:neq),gravlo(0:neq)) 
           kb=0.0;  gravlo=0.0  ; loads = .0
!----------------- element stiffness integration and assembly----------------- 
 elements_3: do iel = 1 , oldele
              if(prop(iel)==1)then
                 gama = gamaf; e = ef ; v = vf
              else
                 gama = gamas; e = es ; v = vs
              end if
              if(iel<=(oldele-newele)) gama = .0
                num = g_num(: , iel) ; coord = transpose(g_coord(:,num )) 
                g = g_g(:,iel);  km=0.0   ;  call deemat(dee,e,v); eld = .0
             gauss_pts_1:  do i =1 , nip    
               call shape_fun(fun,points,i)  ;  gc = matmul ( fun , coord )
!--------------------- initial stress in foundation ---------------------------
              if(ii==1) then
               tensor(2,i,iel)=-1.*(fdepth(fnye+1)-gc(2))*gama
               tensor(1,i,iel)=epk0*tensor(2,i,iel)
               tensor(4,i,iel)=tensor(1,i,iel);tensor(3,i,iel)=.0
              end if 
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;   call beemat (bee,deriv)           
               km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
               do k=2,ndof,2;eld(k)=eld(k)+fun(k/2)*det*weights(i);end do
            end do gauss_pts_1   
   call formkb (kb,km,g)
   if(ii<=lifts) gravlo ( g ) = gravlo ( g ) - eld * gama ; gravlo(0) = .0
 end do elements_3 
!------------------------- factorise equations--------------------------------- 
          call cholin(kb)                                                  
!------------------ factor gravlo by incs-------------------------------------
     write(11,'(a,i5,a,e12.4)')                                                &
             "The total gravity load in lift", ii, " is" , sum(gravlo)
       gravlo = gravlo / incs
!-------------------- apply gravity loads incrementally -----------------------
    load_increments: do iy=1,incs
  write(11,'(a,i5)')                                                           &
          "Increment",iy ; iters=0; oldis =.0; bdylds=.0; evpt(:,:,1:oldele)=.0
!--------------------------   iteration loop  ---------------------------------
   iterations: do
    iters=iters+1;  loads = .0; loads= gravlo+bdylds ;  call chobac(kb,loads)
!-------------------------   check convergence --------------------------------
      call checon(loads,oldis,tol,converged)
      if(iters==1)converged=.false. ;  if(converged.or.iters==limit)bdylds=.0
!----------------------- go round the Gauss Points ----------------------------
      elements_4: do iel = 1 , oldele
                   if(prop(iel)==1)then
                      phi = phif; c = cf;  e = ef ; v = vf; psi = psif
                   else
                      phi = phis; c = cs;  e = es;  v = vs; psi = psis
                   end if
             snph=sin(phi*pi/180.);dt=4.*(1.+v)*(1.-2.*v)/(e*(1.-2.*v+snph**2))
             call deemat(dee,e,v);    bload=.0
             num = g_num( : , iel ) ; coord = transpose(g_coord( : , num ))
             g = g_g( : , iel )     ; eld = loads ( g )         
       gauss_points_2 : do i = 1 , nip
          call shape_der ( der,points,i); jac=matmul(der,coord)
          det = determinant(jac)  ;   call invert(jac)
          deriv = matmul(jac,der) ; call beemat (bee,deriv);eps=matmul(bee,eld)
          eps = eps -evpt( : , i , iel)    ;        sigma=matmul(dee,eps)
          if(ii==1)then;s=tensor(:,i,iel);else;s=tensor(:,i,iel)+sigma;end if
          call invar(s,sigm,dsbar,lode_theta)                                 
!------------------  check whether yield is violated --------------------------
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
        eload=matmul(devp,bee) ; bload=bload+eload*det*weights(i)
      end if
!---------------- if appropriate update the Gauss point stresses --------------
      if(converged.or.iters==limit) then
        if(ii/=1) tensor(:,i,iel) = s
      end if
    end do gauss_points_2
!---------------- compute the total bodyloads vector --------------------------
    bdylds( g ) = bdylds( g ) + bload      ; bdylds(0) = .0
  end do elements_4             
  if(converged.or.iters==limit)exit
 end do iterations
 if(ii/=1) totd(:neq) = totd(:neq) + loads(:neq)
 write(11,'(2(a,i5))') "Lift number",ii," gravity load increment",iy
 write(11,'(a,i5,a)') "It took ",iters, " iterations to converge"
 if(iy==incs.or.iters==limit)write(11,'(a,e12.4)')                             &
                                     "Max displacement is",maxval(abs(loads)) 
 if(iters==limit)stop
end do load_increments
deallocate(lnf,kb,loads,bdylds,oldis,gravlo)
end do lift_number
end program p69
