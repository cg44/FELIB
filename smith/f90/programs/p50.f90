 program p50       
!------------------------------------------------------------------------
!      program 5.0 plane stress of an elastic
!      solid using uniform 3-node triangular 
!      elements numbered in the x direction
!------------------------------------------------------------------------
 use new_library   ;   use  geometry_lib  ;         implicit none
 integer::nels,nce,neq,nband,nn,nr,nip,nodof=2,nod=3,nst=3,ndof,            &
          loaded_nodes,i,k,iel,ndim=2
 real:: e,v,det,aa,bb    ; character(len=15):: element = 'triangle'            
!--------------------- dynamic arrays----------------------------
 real    ,allocatable :: kv(:),loads(:),points(:,:),dee(:,:),coord(:,:),     &
                         jac(:,:), der(:,:),deriv(:,:), weights(:),          &
                         bee(:,:),km(:,:),eld(:),sigma(:),g_coord(:,:)
integer, allocatable :: nf(:,:), g(:) , num(:)  , g_num(:,:) , g_g(:,:)
      
!----------------input and initialisation----------------------
  open (10,file='p50.dat',status='old',    action='read')
  open (11,file='p50.res',status='replace',action='write')

  read (10,*) nels,nce,nn,nip,aa,bb,e,v 
  ndof=nod*nodof 
  allocate ( nf(nodof,nn), points(nip,ndim),g(ndof), g_coord(ndim,nn),        &
            dee(nst,nst),coord(nod,ndim),jac(ndim,ndim),weights(nip),         &
            der(ndim,nod), deriv(ndim,nod), bee(nst,ndof), km(ndof,ndof),     &
            eld(ndof),sigma(nst),num(nod),g_num(nod,nels),g_g(ndof,nels))
  
    nf=1; read(10,*) nr;if(nr>0)read(10,*)(k,nf(:,k),i=1,nr)
    call formnf (nf);neq=maxval(nf)
    nband = 0  
! this is a plane stress analysis 
  dee=.0; dee(1,1)=e/(1.-v*v);dee(2,2)=dee(1,1);dee(3,3)=.5*e/(1.+v)
  dee(1,2)=v*dee(1,1);dee(2,1)=dee(1,2); call sample(element,points,weights)
!--------loop the elements to find bandwidth and neq-------------------
  elements_1: do iel = 1 , nels
              call geometry_3tx(iel,nce,aa,bb,coord,num);call num_to_g(num,nf,g)
              g_num(:,iel)=num;g_coord(:,num)=transpose(coord);g_g(:,iel) = g
              if(nband<bandwidth(g))nband=bandwidth(g)
  end do elements_1                                                           
    write(11,'(a)') "Global coordinates "
    do k=1,nn;write(11,'(a,i5,a,2e12.4)')"Node",k,"       ",g_coord(:,k);end do
    write(11,'(a)') "Global node numbers "
    do k = 1 , nels; write(11,'(a,i5,a,3i5)')                                 &
                              "Element ",k,"        ",g_num(:,k); end do  
    write(11,'(2(a,i5))')                                                     &
             "There are ",neq,"  equations and the half-bandwidth is ", nband
             allocate(kv(neq*(nband+1)),loads(0:neq)); kv=.0        
!------- element stiffness integration and assembly-------------------- 

 elements_2: do iel = 1 , nels
             num = g_num(:, iel);     g = g_g( : , iel )
             coord = transpose(g_coord(:, num)) ;       km=0.0  
          gauss_pts_1: do i = 1 , nip
               call shape_der(der,points,i) ; jac = matmul(der,coord) 
               det = determinant(jac); call invert(jac)
               deriv = matmul(jac,der) ; call beemat (bee,deriv) 
             km = km + matmul(matmul(transpose(bee),dee),bee) *det* weights(i)
          end do gauss_pts_1                                                  
   call formkv (kv,km,g,neq)
 end do elements_2   
                                                         
 loads=.0; read(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)  
!------------------------equation solution--------------------------------
    call banred(kv,neq) ;call bacsub(kv,loads)
    write(11,'(a)') "The nodal displacements Are :"
    write(11,'(a)')  "Node         Displacement"
    do k=1,nn; write(11,'(i5,a,2e12.4)')   k,"   ",loads(nf(:,k)); end do
!-------------------recover stresses at centroidal gauss-point------------
    nip = 1; deallocate(points,weights);allocate(points(nip,ndim),weights(nip))
    call sample ( element , points , weights)
        write(11,'(a)') "The central point stresses are :"
 elements_3:do iel = 1 , nels
         write(11,'(a,i5)') "Element No.  ",iel 
    num = g_num(: , iel);    coord =transpose( g_coord(: ,num) )
    g = g_g( : ,iel )    ;     eld=loads(g)
    gauss_pts_2: do i = 1 , nip
       call shape_der (der,points,i); jac= matmul(der,coord)
       call invert(jac) ;    deriv= matmul(jac,der)
       call beemat(bee,deriv);  sigma = matmul (dee,matmul(bee,eld)) 
       write(11,'(a,i5)') "Point  ",i   ;  write(11,'(3e12.4)') sigma
    end do gauss_pts_2 
 end do elements_3
end program p50        

