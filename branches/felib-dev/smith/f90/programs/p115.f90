    program p115      
!-----------------------------------------------------------------------
!      program 11.5 forced vibration of an elastic-plastic(Von Mises) solid
!      using 8-node quadrilateral elements; viscoplastic strain method
!      rectangular mesh : lumped mass , explicit integration
!------------------------------------------------------------------------
 use new_library      ; use geometry_lib     ;      implicit none
 integer::nels,nxe,neq,nn,nr,nip,nodof=2,nod=8,nst=4,ndof,loaded_nodes,      &
          i,k,iel,ndim=2,jj,nstep ,npri
 real   ::aa,bb,rho,dtim,time,e,v,det,sbary,pload,sigm,f,fnew,fac,           &
          area,sbar,dsbar,lode_theta
 character (len = 15) :: element = 'quadrilateral' 
!---------------------------- dynamic arrays-----------------------------------
 real    ,allocatable :: points(:,:),bdylds(:),x1(:),d1x1(:),stress(:),      &
                         pl(:,:),emm(:),d2x1(:),tensor(:,:,:),etensor(:,:,:),&
                         val(:,:),mm(:),dee(:,:),coord(:,:),jac(:,:),        &
                         weights(:), der(:,:),deriv(:,:),bee(:,:),eld(:),    &
                         eps(:),sigma(:),bload(:),eload(:),g_coord(:,:)
 integer, allocatable :: nf(:,:) , g(:), no(:) ,num(:), g_num(:,:) ,g_g(:,:)   
!-----------------------input and initialisation-------------------------------
  open (10,file='p115.dat',status=    'old',action='read')
  open (11,file='p115.res',status='replace',action='write')                  
  read (10,*) aa,bb,sbary,e,v,rho,pload,   &
              nels,nxe,nn,nip,loaded_nodes,dtim,nstep,npri
  ndof=nod*nodof   
  allocate (nf(nodof,nn), points(nip,ndim),weights(nip),g_coord(ndim,nn),     &
            num(nod),dee(nst,nst),tensor(nst,nip,nels),no(loaded_nodes),      &
            coord(nod,ndim), pl(nst,nst), etensor(nst,nip,nels),              &
            jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),g_num(nod,nels),     &
            bee(nst,ndof),eld(ndof),eps(nst),sigma(nst),emm(ndof),            &
            bload(ndof),eload(ndof),g(ndof), stress(nst),                     &
            val(loaded_nodes,ndim),g_g(ndof,nels))                     
  nf=1; read(10,*) nr ; if(nr>0) read(10,*)(k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)                                            
    read(10,*)(no(i),val(i,:),i=1,loaded_nodes)                                
  ! loop the elements to set up global arrays
       elements_1:   do iel = 1 , nels
                        call geometry_8qx(iel,nxe,aa,bb,coord,num)
                        call num_to_g (num , nf , g); g_num(:,iel)=num
                        g_coord(:,num)=transpose(coord); g_g( : , iel ) = g
       end do elements_1                                                      
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,8i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do      
    write(11,'(a,i5,a)') "There are ",neq," equations to be solved"
  allocate(bdylds(0:neq),x1(0:neq),d1x1(0:neq),d2x1(0:neq),mm(0:neq)) 
  tensor = .0;  etensor = .0
  x1=0.0; d1x1=0.0; d2x1=0.0; mm=0.0  ;  call sample(element,points,weights)
!--------------------- explicit integration loop -----------------------------  
 write(11,'(a)') "   Time     Displacement  Velocity   Acceleration"   
 time = .0
 write(11,'(4e12.4)')time,x1(neq),d1x1(neq),d2x1(neq)
 time_steps : do jj = 1 , nstep
!-------------------------  apply the load ------------------------------------
 time = time + dtim
 x1 = x1 +(d1x1+d2x1*dtim*.5)*dtim  ;  bdylds = .0 
!-------------------- element stress-strain relationship ----------------------
 elements_2: do iel = 1 , nels           
                num = g_num(:,iel) ; coord = transpose(g_coord(: , num )) 
                g = g_g( : , iel ) ; area = 0.0 ; bload = .0 ; eld = x1 ( g )
             gauss_pts_1:  do i =1 , nip  ;  dee = .0 ; call deemat(dee,e,v)
               call shape_der (der,points,i);  jac = matmul(der,coord) 
               det = determinant(jac)  ;   call invert(jac)
               deriv = matmul(jac,der) ;  call beemat (bee,deriv)           
               area = area + det * weights(i)*rho; eps = matmul ( bee , eld )   
               eps = eps - etensor( : , i , iel )
               sigma= matmul ( dee , eps ); stress = sigma+tensor (: , i, iel)
               call invar(stress,sigm,dsbar,lode_theta); fnew = dsbar - sbary
!---------------------- check whether yield is violated -----------------------
           if(fnew>=.0) then   
               stress= tensor(:,i,iel); call invar(stress,sigm,sbar,lode_theta)
               f = sbar - sbary; fac = fnew/(fnew - f) 
               stress = tensor ( : , i , iel )+(1.-fac) * sigma
               call vmpl(e,v,stress,pl); dee = dee - fac * pl
           end if
               sigma = matmul(dee ,eps) ;sigma = sigma + tensor( : , i , iel )
               eload=matmul(sigma,bee); bload= bload+ eload * det * weights(i)
!-----------------------update Gauss point stresses and strains ---------------
               tensor( : , i , iel) = sigma 
               etensor( : , i , iel ) = etensor( : , i , iel ) + eps
            end do gauss_pts_1   
            bdylds ( g ) = bdylds ( g ) - bload  ; bdylds(0) = .0
            if( jj == 1) then
              emm = .2 * area; emm(1:13:4)=.05*area; emm(2:14:4)=.05*area
              mm ( g ) = mm ( g ) + emm
            end if   
 end do elements_2
      do i=1,loaded_nodes 
         bdylds(nf(:,no(i)))=bdylds(nf(:,no(i)))+val(i,:)*pload 
      end do
      bdylds(1:neq) = bdylds(1:neq) / mm(1:neq)
      d1x1=d1x1+(d2x1+bdylds)*.5*dtim  ;  d2x1 = bdylds
      if(jj==jj/npri*npri)write(11,'(4e12.4)')time,x1(neq),d1x1(neq),d2x1(neq)
end do time_steps
end program p115
