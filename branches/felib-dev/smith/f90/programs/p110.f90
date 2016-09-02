program p110
!----------------------------------------------------------
! program 11.0 forced vibration analysis of beams
!----------------------------------------------------------
 use new_library     ;  use  geometry_lib    ;  implicit none 
 integer::nels,neq,nn,nband,nr,nod=2,nodof=2,ndof=4,iel,i,j,k,ndim=1,      &
  np_types,nln,lnode,lsense,np,ntp,non
 real::beta,gamma,fm,fk,dt,f1,f2,el_ei,el_ell    
!--------------------------dynamic arrays-------------------------------------
 real,allocatable::km(:,:),mm(:,:),kb(:,:),mb(:,:),ek(:,:),em(:,:),cb(:,:),&
  kp(:,:),a(:),d(:),v(:),a1(:),b1(:),vc(:),kd(:),coord(:,:),g_coord(:,:),  &
  ei(:),rhoa(:),ell(:),rt(:),rl(:),al(:,:),dis(:,:),vel(:,:),acc(:,:)
 integer,allocatable::nf(:,:),g(:),num(:),g_num(:,:),g_g(:,:),etype(:),    &
  lp(:),lf(:)
!------------------------input and initialisation------------------------------
 open(10,file='p110.dat'); open(11,file='p110.res')
 read(10,*)nels,np_types; nn=nels+1
 allocate(nf(nodof,nn),km(ndof,ndof),mm(ndof,ndof),coord(nod,ndim),        &
  g_coord(ndim,nn),g_num(nod,nels),num(nod),g(ndof),&
  ei(np_types),rhoa(np_types),ell(nels),g_g(ndof,nels),etype(nels)) 
 read(10,*)(ei(i),rhoa(i),i=1,np_types) 
 etype=1; if(np_types>1)read(10,*)etype
 read(10,*)ell,nr
 nf=1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf)
!--------------loop the elements to find global array sizes--------------------
 nband=0
 elements_1: do iel=1,nels 
              el_ell = ell(iel)   
              call geometry_2l(iel,el_ell,coord,num); call num_to_g(num,nf,g) 
              g_num(:,iel)=num; g_coord(:,num)=transpose(coord) 
              g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g)
 end do elements_1
 allocate(kb(neq,nband+1),mb(neq,nband+1),ek(neq,nband+1),em(neq,nband+1), &
  cb(neq,nband+1),kp(neq,nband+1),a1(0:neq),b1(0:neq),                     &
  vc(0:neq),kd(0:neq),a(0:neq),d(0:neq),v(0:neq)) 
 kb=0.0; mb=0.0
!--------global stiffness and mass matrix assembly-------------------------
 elements_2: do iel=1, nels   
              el_ei = ei(etype(iel)) ; el_ell = ell(iel)
              call beam_km(km,el_ei,el_ell); g=g_g(:,iel)
              call beam_mm(mm,rhoa(etype(iel)),ell(iel))
              call formkb(kb,km,g); call formkb(mb,mm,g)
 end do elements_2
 ek=kb; em=mb
!--------read newmark time-stepping parameters-----------------------------
!--------rayleigh damping parameters and initial conditions----------------
 read(10,*)beta,gamma,dt,fm,fk; d=0.0; v=0.0
! read initial conditions if /= 0.0 d(1:),v(1:)
!--------read number of loaded freedoms
!--------for each of these read number of points in load/time history
!--------and coordinates of load/time history
 read(10,*)nln; allocate(lf(nln))    
 do k=1,nln
   read(10,*)lnode,lsense,np;  allocate(rt(np),rl(np))
   lf(k)=nf(lsense,lnode)
   do j=1,np; read(10,*)rt(j),rl(j); end do
   if(k==1)then
     ntp=nint((rt(np)-rt(1))/dt)+1; allocate(al(ntp,nln))
   end if
   call interp(k,dt,rt,rl,al,ntp)
 end do
 f1=beta*dt**2; f2=beta*dt
 cb=fm*mb+fk*kb; kp=mb/f1+gamma*cb/f2+kb
 call cholin(em); call cholin(kp); a=0.0
 a(lf(:))=al(1,:)
 call banmul(cb,v,vc); call banmul(kb,d,kd)
 a=a-vc-kd; call chobac(em,a)
 read(10,*)non
 allocate(lp(non),dis(ntp,non),vel(ntp,non),acc(ntp,non))
 do k=1,non; read(10,*)lnode,lsense; lp(k)=nf(lsense,lnode); end do
 dis(1,:)=d(lp); vel(1,:)=v(lp); acc(1,:)=a(lp) 
!---------Time stepping----------------------------------------------------
 do j=2,ntp
   a1=d/f1+v/f2+a*(0.5/beta-1.)
   b1=gamma*d/f2-v*(1.-gamma/beta)-dt*a*(1.-0.5*gamma/beta)
   call banmul(mb,a1,vc); call banmul(cb,b1,kd)
   d=vc+kd; d(lf(:))=d(lf(:))+al(j,: )
   call chobac(kp,d); v=gamma*d/f2-b1; a=d/f1-a1
   dis(j,:)=d(lp); vel(j,:)=v(lp); acc(j,:)=a(lp)
 end do
 write(11,'(//,a)')"Computed time histories  "
 do i=1,non
!   write(11,'(/,a,i4)')"      Freedom",lp(i)
!   write(11,'(a)')"      Time    Disp        Velo       Accel"
   do j=1,ntp
     write(11,'(f10.2,4e12.4)')(j-1)*dt,dis(j,i),vel(j,i),acc(j,i)
   end do
 end do   
end program p110
