program p59    
!------------------------------------------------------------------------------
!      program 5.9 general analysis of elastic solids
!------------------------------------------------------------------------------
 use new_library      ;     use geometry_lib   ;   implicit none
 integer::nels,neq,nband,nn,nr,nip,nodof,nod,nst,ndof, &
          i,k,iel,ndim,loaded_nodes,fixed_nodes,nprops,np_types         
 real:: det            ;  character (len=15) :: element
!-------------------------- dynamic arrays-------------------------------------
 real    ,allocatable :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:),   &
                         jac(:,:),der(:,:),deriv(:,:),weights(:),prop(:,:),&
                         bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:),    &
                         value(:)  
 integer, allocatable :: nf(:,:), g(:) ,num(:), g_num(:,:) , g_g( :, :),   &
                         no(:),sense(:),node(:) , etype(:)                     
!-----------------------input and initialisation-------------------------------
  open (10,file='p59.dat',status=    'old',action='read')
  open (11,file='p59.res',status='replace',action='write')                     
  read (10,*) element,nels,nn,nip,nodof,nod,nst,ndim     ;  ndof=nod*nodof
  allocate ( nf(nodof,nn), points(nip,ndim),dee(nst,nst), g_coord(ndim,nn),   &
            coord(nod,ndim),etype(nels),jac(ndim,ndim),weights(nip),num(nod), &
            g_num(nod,nels),der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),      &
            km(ndof,ndof),eld(ndof),sigma(nst),g(ndof),g_g(ndof,nels))
         read(10,*) nprops , np_types
         allocate(prop(nprops,np_types))   ;  read (10,*) prop
         etype=1 ; if(np_types>1) read(10,*) etype
         read(10,*) g_coord ; read(10,*) g_num     
        nf=1; read(10,*) nr ; if(nr>0) read(10,*) (k,nf(:,k),i=1,nr)
        call formnf(nf); neq=maxval(nf) ; call sample(element,points,weights)  
!------------- loop the elements to find nband and store steering vectors  ----
      nband=0
   elements_1   : do iel =1,nels
                   num=g_num(:,iel) ; call num_to_g(num,nf,g);
                   g_g( : ,iel ) = g; if(nband<bandwidth(g))nband=bandwidth(g) 
   end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,3e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,27i3)')                                 &
                              "Element ",k,"    ",g_num(:,k); end do  
        write(11,'(2(a,i5))')                                                  &
                 "There are",neq,"  equations and the half-bandwidth is",nband
   allocate( kv(neq*(nband+1)),loads(0:neq)); kv= .0  ; loads =.0 
!--------------- element stiffness integration and assembly--------------------
  elements_2: do iel=1,nels
                call deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)))
                num = g_num(:,iel); coord = transpose(g_coord(: , num))
                g=g_g(:,iel)    ;   km = .0  
                integrating_pts_1:  do i=1,nip
                  call shape_der(der,points,i); jac=matmul(der,coord)
                  det= determinant(jac) ; call invert(jac)
                  deriv = matmul(jac,der);call beemat(bee,deriv)
                  km= km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
                end do integrating_pts_1    
               call formkv (kv,km,g,neq) 
  end do elements_2   
 read(10,*) loaded_nodes
 if(loaded_nodes/=0)read (10,*)(k,loads(nf(:,k)),i=1,loaded_nodes) 
  read(10,*) fixed_nodes 
 if(fixed_nodes/=0)then
          allocate(node(fixed_nodes),sense(fixed_nodes),                      &
                   value(fixed_nodes),no(fixed_nodes))
          read(10,*)(node(i),sense(i),value(i),i=1,fixed_nodes)
          do i=1,fixed_nodes; no(i)=nf(sense(i),node(i)); end do
          kv(no)=kv(no) + 1.e20; loads(no) = kv(no) * value
 end if
!-----------------------------equation solution-------------------------------- 
    call banred(kv,neq) ;call bacsub(kv,loads)
    write(11,'(a)') "The nodal displacements are:"
    do k=1,nn; write(11,'(i5,a,3e12.4)') k,"   ",loads(nf(:,k)); end do
!---------------------recover stresses at element Gauss-points-----------------
 elements_3:do iel=1,nels
             num = g_num(:,iel) ; coord = transpose(g_coord(:,num))
             g = g_g( : , iel )  ;  eld = loads(g)
             call deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)))
             write(11,'(a,i5,a)')                                             &
                      "The Gauss point stresses for element",iel,"  are :" 
             integrating_pts_2: do i=1,nip
                 call shape_der(der,points,i); jac=matmul(der,coord) 
                 call invert(jac); deriv=matmul(jac,der)
                 call beemat(bee,deriv); sigma=matmul(dee,matmul(bee,eld))
                 write(11,'(a,i5)') "Point",i  ;  write(11,'(6e12.4)') sigma
             end do integrating_pts_2 
 end do elements_3
end program p59         
