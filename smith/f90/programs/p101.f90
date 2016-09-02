program p101    
!------------------------------------------------------------------------
!      program 10.1 eigenvalues and eigenvectors of a rectangular
!      elastic solid in plane strain using uniform 4-node
!      quadrilateral elements  :     lumped mass
!------------------------------------------------------------------------
 use new_library  ;   use geometry_lib  ;  use vlib   ;  implicit none
 integer::nels,nye,neq,nband,nn,nr,nodof=2,nod=4,ndof,                      &
          i,j,k,iel,ndim=2,ifail,icount,nmodes  
 real::aa,bb,e,v,rho,tol=1.e-30
!------------------------------ dynamic arrays---------------------------------
 real    ,allocatable :: ku(:,:),loads(:),coord(:,:),km(:,:),g_coord(:,:),  &
                         diag(:),udiag(:),emm(:,:),kv(:),kh(:),rrmass(:)
 integer, allocatable :: nf(:,:), g(:)  , num(:)  , g_num(:,:) , g_g (:,:)     
!------------------------input and initialisation------------------------------
  open (10,file='p101.dat',status=    'old',action='read')
  open ( 6,file='p101.res',status='replace',action='write')                    
  read (10,*) nels,nye,nn,aa,bb,rho,e,v,nmodes  
  ndof=nod*nodof
  allocate ( nf(nodof,nn), g_coord(ndim,nn),coord(nod,ndim),emm(ndof,ndof),  &
            g_num(nod,nels),num(nod),km(ndof,ndof),g(ndof),g_g(ndof,nels))     
  nf=1; read(10,*) nr ; if(nr>0) read(10,*) (k,nf(:,k),i=1,nr)
  call formnf(nf); neq=maxval(nf)                                             
!------- loop the elements to find nband and set up global arrays ------------
      nband=0
 elements_1   : do iel =1,nels
                 call geometry_4qy(iel,nye,aa,bb,coord,num)
                 call num_to_g ( num , nf , g )
                 g_num(:,iel)=num;g_coord(:,num)=transpose(coord);g_g(:,iel)=g
                 if(nband<bandwidth(g))nband=bandwidth(g) 
 end do elements_1
        write(6,'(a)') "Global coordinates"
        do k=1,nn;write(6,'(a,i5,a,2e12.4)')"Node",k,"    ",g_coord(:,k);end do
        write(6,'(a)') "Global node numbers"
        do k=1,nels;write(6,'(a,i5,a,4i5)')                                   &
                         "Element",k,"      ",g_num(:,k) ; end do
     write(6,'(2(a,i5))')                                                     &
             "There are ",neq,"  equations and the half-bandwidth is", nband
   allocate( ku(neq,nband+1),loads(0:neq),diag(0:neq),udiag(0:neq),           &
             kv(neq*(nband+1)),kh(neq*(nband+1)),rrmass(0:neq))
        emm = .0; diag = .0; ku = .0
!--------element mass matrix is lumped----------------------------------------
        do i=1,ndof; emm(i,i)=.25*aa*bb*rho; end do
!------- element stiffness and mass integration and assembly------------------
 elements_2: do iel=1,nels
                num = g_num(:,iel); coord =transpose(g_coord(:,num))
                g = g_g( : , iel );   call analy4(km,coord,e,v)
                call formku (ku,km,g)   ;  call formlump(diag,emm,g)
  end do elements_2
  write(6,'(a)') "The global mass diagonal is:"
  write(6,'(6e12.4)') diag(1:neq)
!----------------------reduce to standard eigenvalue problem-----------------  
  rrmass(1:neq) = 1./sqrt(diag(1:neq))
  do i=1,neq
     if(i<=neq-nband)then;k=nband+1;else;k=neq-i+1;end if
     do j=1,k; ku(i,j)=ku(i,j)*rrmass(i)*rrmass(i+j-1); end do
  end do
  icount=0
  do j=1,nband+1;do i=1,neq; icount=icount+1;kh(icount)=ku(i,j);end do; end do
!----------------------extract the eigenvalues--------------------------------
  call bandred(ku,diag,udiag,loads);ifail=1; call bisect(diag,udiag,tol,ifail)
  write(6,'(a)') "The eigenvalues are:"    ;  write(6,'(6e12.4)') diag(1:neq)
 do i = 1 , nmodes
    kv = kh; kv(:neq)=kv(:neq)-diag(i);kv(1)=kv(1)+1.e20
    udiag=0.0; udiag(1)=kv(1)
    call banred(kv,neq);call bacsub(kv,udiag);udiag=rrmass*udiag
    write(6,'("Eigenvector number ",i3," is: ")')i
    write(6,'(6e12.4)')udiag(1:)/maxval(abs(udiag))
 end do
end program p101                                                            
