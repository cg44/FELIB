program p104    
!------------------------------------------------------------------------------
!      program 10.4 eigenvalues and eigenvectors of a
!      rectangular elastic solid in plane strain using
!      uniform 4-node quadrilateral elements : kind = 1  precision 
!      for lumped mass this is done element by element
!------------------------------------------------------------------------------
 use libks  ;  use new_library ; use geometry_lib ;  implicit none
 integer::nels,nye,neq,nn,nr,nip,nodof=2,nod=4,nst=3,ndof,                  &
          i,k,iel,ndim=2,nmodes,jflag,iflag=-1,itape=1,lp=6 ,               &
          lalfa=500,leig=20, lx=80, lz=500  ,iters  ,neig = 0
 real::aa,bb,rho,e,v,det  , el,er,  acc = 1.e-6 
 character (len=15) :: element = 'quadrilateral'
!--------------------------- dynamic arrays------------------------------------
 real,allocatable :: points(:,:),dee(:,:),coord(:,:),vdiag(:),    &
                         fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),  & 
                         bee(:,:),km(:,:),emm(:,:),ecm(:,:),g_coord(:,:), &
                         ua(:),va(:),eig(:),x(:),del(:), udiag(:),diag(:),&
                         alfa(:),beta(:),w1(:),y(:,:),z(:,:),pmul(:),utemp(:) 
 integer, allocatable :: nf(:,:), g(:)  , num(:)  , g_num(:,:) , g_g (:,:),&
                         nu(:),jeig(:,:)                            
!----------------------input and initialisation--------------------------------
  open (10,file='p104.dat',status=    'old',action='read')
  open (11,file='p104.res',status='replace',action='write')
  open ( 1,file='p104.tem',form='unformatted')                             
  read (10,*) nels,nye,nn,nip,aa,bb,rho,e,v,nmodes,el,er  
  ndof=nod*nodof   
  allocate ( nf(nodof,nn), points(nip,ndim),dee(nst,nst), g_coord(ndim,nn),  &
            coord(nod,ndim),fun(nod),jac(ndim,ndim), weights(nip),           &
            g_num(nod,nels),der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),     &
            num(nod),km(ndof,ndof),g(ndof),g_g(ndof,nels),emm(ndof,ndof),    &
            ecm(ndof,ndof),eig(leig),x(lx),del(lx),nu(lx),jeig(2,leig),      &
            alfa(lalfa),beta(lalfa),z(lz,leig),pmul(ndof),utemp(ndof))        
  nf=1; read(10,*) nr ; if(nr>0) read(10,*) (k,nf(:,k),i=1,nr)
        call formnf(nf); neq=maxval(nf)                                        
!------------------ loop the elements to set up global arrays  ----------------
 elements_1   : do iel =1,nels
                 call geometry_4qy(iel,nye,aa,bb,coord,num)
                 call num_to_g ( num , nf , g )
                 g_num(:,iel)=num;g_coord(:,num)=transpose(coord);g_g(:,iel)=g 
                end do elements_1                                             
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do      
       write(11,'(a,i5,a)') "There are ",neq,"  equations to be solved"
   allocate  ( ua(0:neq),va(0:neq),vdiag(0:neq),                              &
               diag(0:neq),udiag(0:neq),w1(0:neq), y(0:neq,leig)) 
    ua = .0 ; va = .0  ; eig = .0
    jeig = 0;  x=.0; del=.0; nu=0; alfa=.0; beta=.0   
    diag = .0 ; udiag = .0 ; w1 = .0 ; y=.0; z=.0          
    call sample( element, points, weights); call deemat(dee,e,v)
!--------------- element stiffness integration and assembly--------------------
 elements_2: do iel=1,nels
                num = g_num(:,iel); coord =transpose( g_coord(:, num ))
                g = g_g( : , iel );     km=0.0    ; emm=0.0    
                integrating_pts_1:  do i=1,nip
                  call shape_fun(fun,points,i)
                  call shape_der(der,points,i); jac=matmul(der,coord)
                  det= determinant(jac) ; call invert(jac) 
                  deriv = matmul(jac,der);call beemat(bee,deriv)
                  km= km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
                  call ecmat(ecm,fun,ndof,nodof);emm=emm+ecm*det*weights(i)*rho
                end do integrating_pts_1
                do i=1,ndof; diag(g(i))=diag(g(i))+sum(emm(i,:));end do
  end do elements_2
!------------------------------find eigenvalues--------------------------------
  diag = 1. / sqrt(diag) ; diag(0) = .0 ! diag holds l**(-1/2)
    do iters = 1 , lalfa
       call lancz1(neq,el,er,acc,leig,lx,lalfa,lp,itape,iflag,ua,va,          &
                     eig,jeig,neig,x,del,nu,alfa,beta)
       if(iflag==0) exit
       if(iflag>1) then
          write(11,'(a,i5)')                                                  &
                   " Lancz1 is signalling failure, with iflag = ",   iflag
          stop
       end if           
!----- iflag = 1 therefore form u + a * v  ( done element by element )---------
      vdiag = va ;    vdiag = vdiag * diag ! vdiag is l**(-1/2).va
      udiag = .0  ; vdiag(0)=.0
      elements_3 : do iel = 1 , nels
                      g = g_g( : , iel )
                      pmul = vdiag (g); utemp = matmul(km,pmul)
                      udiag(g) = udiag(g) + utemp     ! udiag is A.l**(-1/2).va
      end do elements_3
      udiag = udiag *  diag  ; ua = ua + udiag
    end do
!-------------- iflag = 0 therefore write out the spectrum --------------------
      write(11,'(2(a,e12.4))') "The range is",el,"  to ",er 
      write(11,'(a,i5,a)') "It took ",iters,"  iterations"
      write(11,'(a)') "The eigenvalues are   :"
      write(11,'(6e12.4)') eig(1:neig)       
!  calculate the eigenvectors
   if(neig>10)neig = 10
   call lancz2(neq,lalfa,lp,itape,eig,jeig,neig,alfa,beta,lz,jflag,y,w1,z)     
!------------------if jflag is zero  calculate the eigenvectors ---------------
   if (jflag==0) then
     write(11,'(a)') "The eigenvectors are  :"  
     do i = 1 , nmodes
        udiag(:) = y(:,i)  ; udiag = udiag * diag 
        write(11,'("Eigenvector number  ",i4," is: ")')  i
        write(11,'(6e12.4)') udiag(1:)
     end do
   else
! lancz2 fails
     write(11,'(a,i5)')" Lancz2 is signalling failure with jflag = ",  jflag
   end if
end program p104
