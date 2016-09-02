program p54    
!------------------------------------------------------------------------
!      program 5.4 axisymmetric strain of a rectangular section elastic
!      solid using variable 4-node quadrilateral elements numbered in y 
!      and variable material properties
!------------------------------------------------------------------------
 use new_library   ;  use  geometry_lib   ;      implicit none
 integer::nels,nre,nde,neq,nband,nn,nr,nip,nodof=2,nod=4,nst=4,ndof,    &
          i,k,iel,ndim=2,loaded_nodes                  
 real:: e,v,det ,radius    ;  character(len=15) :: element = 'quadrilateral'    
!-------------------------- dynamic arrays--------------------------------
 real    ,allocatable :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:), &
                         fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:), &
                         bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:),  &
                         width(:),depth(:) ,prop(:,:)
 integer, allocatable :: nf(:,:), g(:)  , num(:)  , g_num(:,:) , g_g(:,:)      
!-------------------------input and initialisation-----------------------------
  open (10,file='p54.dat',status=    'old',action='read')
  open (11,file='p54.res',status='replace',action='write')                     
  read (10,*) nels,nre,nde,nn,nip        ;    ndof=nod*nodof     
  allocate ( nf(nodof,nn), points(nip,ndim),dee(nst,nst), g_coord(ndim,nn),  &
            coord(nod,ndim),fun(nod),jac(ndim,ndim), weights(nip),           &
            g_num(nod,nels),der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),     &
            num(nod),km(ndof,ndof), eld(ndof),  sigma(nst), g(ndof),         &
            width(nre+1),depth(nde+1),prop(2,nels), g_g(ndof,nels))
  read(10,*) width ; read(10,*) depth     
  read(10,*)(prop(k,:),k=1,2)   
  nf=1; read(10,*) nr ;if(nr>0) read(10,*) (k,nf(:,k),i=1,nr)
      call formnf(nf); neq=maxval(nf)  ; call sample(element,points,weights)   
!-------- loop the elements to find nband and set up global arrays ----------
      nband=0
     elements_1   : do iel =1,nels     
                 call geometry_4qyv(iel,nde,width,depth,coord,num)
                 call num_to_g ( num , nf , g );  g_num(:,iel)=num
                 g_coord(:,num)=transpose(coord);g_g(:,iel)=g
                 if(nband<bandwidth(g))nband=bandwidth(g) 
     end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,4i5)')                                  &
                              "Element ",k,"        ",g_num(:,k); end do  
        write(11,'(2(a,i5))')                                                  &
                 "There are",neq,"  equations and the half-bandwidth is",nband
       allocate( kv(neq*(nband+1)),loads(0:neq)); kv=0.0 
!--------------- element stiffness integration and assembly--------------------
  elements_2: do iel=1,nels
                num=g_num(:,iel)  ; coord =transpose( g_coord(:,num)) 
                g = g_g(: ,iel)   ;     km=0.0
                e = prop(1 , iel);  v = prop(2 , iel); call deemat(dee,e,v)
                integrating_pts_1:  do i=1,nip
                  call shape_fun(fun,points,i); call shape_der(der,points,i)
                  jac=matmul(der,coord) ; det= determinant(jac) 
                  call invert(jac);   deriv = matmul(jac,der)
                  call bmataxi(bee,radius,coord,deriv,fun);det =det*radius
                  km= km+matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
                end do integrating_pts_1                                    
               call formkv (kv,km,g,neq)
  end do elements_2                                                            
  loads=0.0 ; read (10,*) loaded_nodes,(k,loads(nf(:,k)), i =1,loaded_nodes)
!------------------------equation solution--------------------------------      
    call banred(kv,neq) ;call bacsub(kv,loads)
    write(11,'(a)') "The nodal displacements are:"
    do k=1,nn; write(11,'(i5,a,2e12.4)') k,"   ",loads(nf(:,k)); end do
!-------------------recover stresses at element centroids-----------------
  nip = 1; deallocate(points,weights); allocate(points(nip,ndim),weights(nip))
       call sample (element , points , weights )   
 elements_3:do iel=1,nels
             num = g_num(:,iel) ; coord = transpose(g_coord(:,num))  
             g = g_g( : , iel)  ;     eld=loads(g)
             write(11,'(a,i5,a)')                                             &
                      "The centroidal stresses for element",iel,"  are :" 
              e = prop(1 ,iel);  v = prop(2 , iel); call deemat(dee,e,v)
              integrating_pts_2: do i = 1 , nip
                 call shape_fun(fun,points,i); call shape_der(der,points,i)
                 jac=matmul(der,coord);call invert(jac); deriv=matmul(jac,der)
                 call bmataxi(bee,radius,coord,deriv,fun)
                 sigma=matmul(dee,matmul(bee,eld))
                 write(11,'(a,i5)') "Point",i   ; write(11,'(4e12.4)') sigma
             end do integrating_pts_2 
 end do elements_3
end program p54
