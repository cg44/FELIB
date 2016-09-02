program p100
!------------------------------------------------------------------------------
! program 10.0 eigenvalues and eigenvectors of a string of beam elements
!------------------------------------------------------------------------------
use new_library      ;    use geometry_lib ;        implicit none 
integer::nels,neq,nn,nband,nr,nod=2,nodof=2,ndof=4,iel,i,j,k,ndim=1,      &
         np_types,ifail,icount,nmodes
real::tol=1.e-30, el_ei , el_ell
!--------------------------dynamic arrays--------------------------------------
real,allocatable::km(:,:),ku(:,:),loads(:),coord(:,:),g_coord(:,:),ei(:), &
                  rhoa(:),ell(:),diag(:),udiag(:),emm(:,:),kv(:),kh(:),   &
                  rrmass(:)
integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),etype(:)
!--------------------input and initialisation----------------------------------
open(10,file='p100.dat'); open(11,file='p100.res')
read(10,*)nels,np_types,nmodes; nn=nels+1
allocate(nf(nodof,nn),km(ndof,ndof),coord(nod,ndim),g_coord(ndim,nn),     &
         g_num(nod,nels),num(nod),g(ndof),emm(ndof,ndof),ei(np_types),    &
         rhoa(np_types),ell(nels),g_g(ndof,nels),etype(nels)) 
read(10,*)(ei(i),rhoa(i),i=1,np_types)
etype=1; if(np_types>1)read(10,*)etype
read(10,*)ell,nr
nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf)
!-------------loop the elements to find global array sizes--------------------
nband=0
elements_1: do iel=1,nels    
              el_ell = ell(iel) ; call geometry_2l(iel,el_ell,coord,num)
              call num_to_g ( num , nf , g )
              g_num(:,iel)=num; g_coord(:,num)=transpose(coord) 
              g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g)
end do elements_1
allocate(ku(neq,nband+1),kv(neq*(nband+1)),kh(neq*(nband+1)),              &
         loads(0:neq),diag(0:neq),udiag(0:neq),rrmass(0:neq))
write(11,'(a)')"Global coordinates"
do k=1,nn; write(11,'(a,i5,a,3e12.4)')                                     &
    "Node    ",k,"    ",g_coord(:,k); end do
write(11,'(a)')"Global node numbers"
do k=1,nels; write(11,'(a,i5,a,27i3)')                                     &
    "Element ",k,"      ",g_num(:,k); end do
write(11,'(2(a,i5),/)')                                                    &
    "There are ",neq,"  equations and the half-bandwidth is ",nband
!-------------global stiffness and (lumped) mass matrix assembly---------------
diag=0.0; ku=0.0
elements_2: do iel=1,nels
              emm=0.0;   el_ei = ei(etype(iel)); el_ell = ell(iel)
              emm(1,1)=0.5*rhoa(etype(iel))*el_ell; emm(3,3)=emm(1,1)
              emm(2,2)=emm(1,1)*el_ell**2/12.; emm(4,4)=emm(2,2)
              call beam_km(km,el_ei,el_ell); g=g_g(:,iel)
              call formku(ku,km,g); call formlump(diag,emm,g)
end do elements_2
write(11,*)"The global mass diagonal is:"; write(11,'(6e12.4)')diag(1:)
!----------------------reduce to standard eigenvalue problem-----------------
rrmass(1:)=1./sqrt(diag(1:))
do i=1,neq
  if(i<=neq-nband)then; k=nband+1; else; k=neq-i+1; end if
  do j=1,k; ku(i,j)=ku(i,j)*rrmass(i)*rrmass(i+j-1); end do
end do
icount=0 
do j=1,nband+1; do i=1,neq 
icount=icount+1; kh(icount)=ku(i,j)
end do; end do
!-------------------------extract the eigenvalues------------------------------
call bandred(ku,diag,udiag,loads);ifail=1; call bisect(diag,udiag,tol,ifail)
write(11,*)"The eigenvalues are:"; write(11,'(6e12.4)')diag(1:)
!----------------------extract the eigenvectors-----------------------------
do i=1,nmodes
  kv=kh; kv(:neq)=kv(:neq)-diag(i); kv(1)=kv(1)+1.e20
  udiag=0.0; udiag(1)=kv(1)
  call banred(kv,neq); call bacsub(kv,udiag); udiag=rrmass*udiag
  write(11,'("Eigenvector number",i3," is:")')i
  write(11,'(6e12.4)')udiag(1:)/maxval(abs(udiag))
end do       
end program p100                                                    
 
