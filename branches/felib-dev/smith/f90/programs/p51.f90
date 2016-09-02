 program p51       
!------------------------------------------------------------------------
!      program 5.1 plane strain of a rectangular elastic
!      solid using variable-sized 15-node triangular 
!      elements numbered in the y direction
!------------------------------------------------------------------------
 use new_library  ;  use geometry_lib    ;    implicit none
 integer::nels,nce,nye,neq,nband,nn,nr,nip,nodof=2,nod=15,nst=3,ndof,       &
          loaded_nodes,i,k,iel,ndim=2 , nde
 real:: e,v,det    ;      character(len=15):: element = 'triangle'
!--------------------- dynamic arrays-------------------------------------
 real    ,allocatable :: kb(:,:),loads(:),points(:,:),dee(:,:),coord(:,:),   &
                         jac(:,:), der(:,:),deriv(:,:), weights(:),          &
                         bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:),      &
                         width(:), depth(:)
 integer, allocatable :: nf(:,:), g(:) , num(:)  , g_num(:,:)  , g_g(: , :)     
!----------------input and initialisation---------------------------------
  open (10,file='p51.dat',status=    'old',action='read')
  open (11,file='p51.res',status='replace',action='write')                    
  read (10,*) nels,nce,nye,nn,nip,e,v 
  ndof=nod*nodof     ;  nde = (nye + 2) / 2  
  allocate ( nf(nodof,nn), points(nip,ndim),g(ndof), g_coord(ndim,nn),        &
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),weights(nip),         &
            der(ndim,nod), deriv(ndim,nod), bee(nst,ndof), km(ndof,ndof),     &
            eld(ndof),sigma(nst),num(nod),g_num(nod,nels),width(nce+1),       &
            depth(nde),g_g(ndof , nels))
  read(10,*) width , depth    
  nf=1; read(10,*) nr ; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
        call formnf (nf);neq=maxval(nf)
  call deemat (dee,e,v); call sample(element,points,weights)
!----------------loop the elements to find bandwidth and neq-------------------
 nband = 0
  elements_1: do iel = 1 , nels
              call geometry_15tyv(iel,nye,width,depth,coord,num)
              call num_to_g(num,nf,g);  g_num(:,iel) = num
              g_coord(:,num)=transpose( coord );g_g(:,iel)=g 
              if(nband<bandwidth(g))nband=bandwidth(g)
  end do elements_1
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,15i3)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do      
    write(11,'(2(a,i5))')                                                      &
             "There are ",neq,"  equations and the half-bandwidth is", nband
             allocate(kb(neq,nband+1),loads(0:neq)); kb=.0        
!------------- element stiffness integration and assembly--------------------  
 elements_2: do iel = 1 , nels
             num = g_num(:,iel);    g = g_g( : , iel )
             coord = transpose(g_coord(: ,num))  ;        km=0.0 
          gauss_points_1: do i = 1 , nip
               call shape_der(der,points,i) ; jac = matmul(der,coord) 
               det = determinant(jac); call invert(jac)
               deriv = matmul(jac,der) ; call beemat (bee,deriv) 
             km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
          end do gauss_points_1                   
   call formkb (kb,km,g)
 end do elements_2                                                             
 loads=.0; read(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)    
!------------------------equation solution--------------------------------
    call cholin(kb) ;call chobac(kb,loads)
    write(11,'(a)') "The nodal displacements Are :"
    write(11,'(a)') "Node         Displacement"
    do k=1,nn; write(11,'(i5,a,2e12.4)') k,"   ",loads(nf(:,k)); end do
!-------------------recover stresses at centroidal point------------------
  nip = 1; deallocate(points,weights);allocate(points(nip,ndim),weights(nip))
    call sample( element ,points , weights )
        write(11,'(a)') "The central  point stresses are :"
 elements_3:do iel = 1 , nels
         write(11,'(a,i5)') "Element no.  ",iel 
      num = g_num(:, iel); coord = transpose(g_coord(:,num))
      g = g_g( : ,iel); eld=loads(g)
    gauss_points_2 : do i = 1 , nip
       call shape_der (der,points,i); jac= matmul(der,coord)
       call invert(jac) ;    deriv= matmul(jac,der)
       call beemat(bee,deriv) ; sigma = matmul (dee,matmul(bee,eld)) 
       write(11,'(a,i5)') "Point  ",i   ;  write(11,'(3e12.4)') sigma
    end do gauss_points_2 
 end do elements_3
end program p51
