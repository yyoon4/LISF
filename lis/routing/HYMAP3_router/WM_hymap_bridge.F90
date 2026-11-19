! Aggregate runtime-selected cells by member/phase/rank. No MPI or RNG calls.
! HYMAP3 SWOT surface-elevation diagnostics, currently for a single nest.
! MARK values: river/flood/runoff-reservoir/baseflow-reservoir storage [m3],
! sum absolute stored WSE [m], outlet river+flood storage [m3], wrapper outlet
! flow [m3/s, already includes flood flow], outlet count (global sum must be1).
! VECTOR values: sum, sum(abs), max(abs), predicted post-profile storage [m3],
! dry-cell count, outlet value. For increment phases, values4/5 instead count
! abs(increment)>=10m and >5m. Predicted storage is not an observed increment.
! CORE values: river/flood storage [m3], runoff, river/flood inflow, river/flood
! outflow [m3/s], dt [s], outlet combined outflow, outlet count, predicted flux
! clipping [m3, use PRE_FLUX], predicted negative-runoff limit [m3, use POST_FLUX].
! Retain row order for repeated core substeps sharing the same model tick.
! Pre/post storage differences are direct; predictions and routing residuals
! require separate attribution (including evaporation/drainage outside hooks).
module WM_hymap_bridge_util
 use iso_fortran_env,only:int64,real64
 use LIS_coreMod,only:LIS_rc,LIS_localpet,LIS_ews_halo_ind,LIS_nss_halo_ind,&
  LIS_routing_gdeltas,LIS_routing_goffsets,LIS_npes
 use HYMAP3_routingMod,only:HYMAP3_routing_struc
 use WM_mass_diag_mod
 implicit none
 integer,allocatable,save::owner(:)
 logical,save::owner_ready=.false.
contains
 function wm_clock()result(t)
  integer(int64)::t
  t=int(LIS_rc%yr,int64)*100000000_int64+int(LIS_rc%mo,int64)*1000000_int64+&
    int(LIS_rc%da,int64)*10000_int64+int(LIS_rc%hr,int64)*100_int64+LIS_rc%mn
 end function
 integer function wm_cell(n,i)result(c)
  integer,intent(in)::n,i
  integer::x,y,l,j,ix,gc
  if(.not.owner_ready)then
   if(LIS_rc%nnest/=1)error stop 'WM diagnostics currently require a single nest'
   allocate(owner(LIS_rc%gnc(n)*LIS_rc%gnr(n)));owner=-1
   ! Match HYMAP3_gather_tiles' last-rank assignment if halo copies exist.
   ! Build once from its actual global sequence, without a new collective.
   do l=1,LIS_npes
    do j=1,LIS_routing_gdeltas(n,l-1)
     ix=j+LIS_routing_goffsets(n,l-1)
     x=HYMAP3_routing_struc(n)%seqx_glb(ix)+LIS_ews_halo_ind(n,l)-1
     y=HYMAP3_routing_struc(n)%seqy_glb(ix)+LIS_nss_halo_ind(n,l)-1
     gc=(y-1)*LIS_rc%gnc(n)+x
     if(gc<1.or.gc>size(owner))error stop 'WM owner mapping error'
     owner(gc)=l-1
    enddo
   enddo
   owner_ready=.true.
  endif
  x=HYMAP3_routing_struc(n)%seqx(i)+LIS_ews_halo_ind(n,LIS_localpet+1)-1
  y=HYMAP3_routing_struc(n)%seqy(i)+LIS_nss_halo_ind(n,LIS_localpet+1)-1
  c=(y-1)*LIS_rc%gnc(n)+x
  if(c<1.or.c>size(owner))error stop 'WM local/global mapping error'
  if(owner(c)/=LIS_localpet)c=0
 end function
end module

subroutine wm_hymap_mark(phase,n,member)
 use WM_hymap_bridge_util
 implicit none
 character(*),intent(in)::phase
 integer,intent(in)::n,member
 real(real64)::v(8)
 integer::m,i,c,k,m1,m2
 if(.not.wm_active(wm_clock()))return
 if(trim(LIS_rc%routingmodel)/='HYMAP3 router')return
 m1=1;m2=LIS_rc%nensem(n)
 if(member>0)then
  m1=member;m2=member
 endif
 do m=m1,m2
  v=0.;k=0
  do i=1,HYMAP3_routing_struc(n)%nseqall
   c=wm_cell(n,i)
   if(c<=0)cycle
   if(.not.wm_contains(c))cycle
   k=k+1
   v(1)=v(1)+dble(HYMAP3_routing_struc(n)%rivsto(i,m))
   v(2)=v(2)+dble(HYMAP3_routing_struc(n)%fldsto(i,m))
   v(3)=v(3)+dble(HYMAP3_routing_struc(n)%rnfsto(i,m))
   v(4)=v(4)+dble(HYMAP3_routing_struc(n)%bsfsto(i,m))
   v(5)=v(5)+dble(HYMAP3_routing_struc(n)%sfcelv(i,m))
   if(c==wm_outlet_cell)then
    v(6)=dble(HYMAP3_routing_struc(n)%rivsto(i,m))+dble(HYMAP3_routing_struc(n)%fldsto(i,m))
    v(7)=dble(HYMAP3_routing_struc(n)%rivout(i,m)) ! wrapper already includes flood.
    v(8)=1.
   endif
  enddo
  if(k>0)call wm_emit(LIS_localpet,wm_clock(),phase,m,k,v)
 enddo
end subroutine

subroutine wm_hymap_vector(phase,n,nval,vec,is_absolute)
 use WM_hymap_bridge_util
 use HYMAP3_modelMod,only:HYMAP3_get_volume_profile
 implicit none
 character(*),intent(in)::phase
 integer,intent(in)::n,nval,is_absolute
 real,intent(in)::vec(nval)
 real(real64)::v(6),h,vol
 integer::i,c,m,t,k
 if(.not.wm_active(wm_clock()))return
 if(trim(LIS_rc%routingmodel)/='HYMAP3 router')return
 if(nval/=HYMAP3_routing_struc(n)%nseqall*LIS_rc%nensem(n))error stop 'WM vector mapping mismatch'
 do m=1,LIS_rc%nensem(n)
  v=0.;k=0
  do i=1,HYMAP3_routing_struc(n)%nseqall
   c=wm_cell(n,i)
   if(c<=0)cycle
   if(.not.wm_contains(c))cycle
   k=k+1;t=(i-1)*LIS_rc%nensem(n)+m
   v(1)=v(1)+dble(vec(t));v(2)=v(2)+abs(dble(vec(t)))
   v(3)=max(v(3),abs(dble(vec(t))))
   if(is_absolute==2)then
    if(abs(vec(t))>=10.)v(4)=v(4)+1
    if(abs(vec(t))>5.)v(5)=v(5)+1
   else
    if(is_absolute==1)then
     h=dble(vec(t))
    else
     ! Match setSWOT's default-REAL addition before the REAL*8 profile call.
     h=dble(vec(t)+HYMAP3_routing_struc(n)%rivelv(i))
    endif
    call HYMAP3_get_volume_profile(HYMAP3_routing_struc(n)%nz,&
     dble(HYMAP3_routing_struc(n)%elevtn(i)),dble(HYMAP3_routing_struc(n)%fldhgt(i,:)),&
     dble(HYMAP3_routing_struc(n)%fldstomax(i,:,m)),dble(HYMAP3_routing_struc(n)%grarea(i)),&
     dble(HYMAP3_routing_struc(n)%rivstomax(i,m)),dble(HYMAP3_routing_struc(n)%rivelv(i)),&
     dble(HYMAP3_routing_struc(n)%rivlen(i)),dble(HYMAP3_routing_struc(n)%rivwth(i,m)),h,vol)
    v(4)=v(4)+dble(real(vol))
    if(h<=dble(HYMAP3_routing_struc(n)%rivelv(i)))v(5)=v(5)+1
   endif
   if(c==wm_outlet_cell)v(6)=dble(vec(t))
  enddo
  if(k>0)call wm_emit(LIS_localpet,wm_clock(),phase,m,k,v)
 enddo
end subroutine

subroutine wm_hymap_core(phase,n,nseq,rs,fs,runoff,ri,fi,ro,fo,dt)
 use WM_hymap_bridge_util
 implicit none
 character(*),intent(in)::phase
 integer,intent(in)::n,nseq
 real,intent(in)::rs(nseq),fs(nseq),runoff(nseq),ri(nseq),fi(nseq),ro(nseq),fo(nseq),dt
 real(real64)::v(12)
 integer::i,c,k
 if(.not.wm_active(wm_clock()))return
 v=0.;k=0
 do i=1,nseq
  c=wm_cell(n,i)
  if(c<=0)cycle
  if(.not.wm_contains(c))cycle
  k=k+1
  v(1)=v(1)+dble(rs(i));v(2)=v(2)+dble(fs(i));v(3)=v(3)+dble(runoff(i))
  v(4)=v(4)+dble(ri(i));v(5)=v(5)+dble(fi(i))
  v(6)=v(6)+dble(ro(i));v(7)=v(7)+dble(fo(i));v(8)=dble(dt)
  v(11)=v(11)+dble(max(0.,-(rs(i)+(ri(i)-ro(i))*dt)))+&
   dble(max(0.,-(fs(i)+(fi(i)-fo(i))*dt)))
  if(runoff(i)<0.)v(12)=v(12)+dble(max(0.,-runoff(i)*dt-rs(i)-fs(i)))
  if(c==wm_outlet_cell)then
   v(9)=dble(ro(i))+dble(fo(i)) ! core components, before wrapper adds.
   v(10)=1.
  endif
 enddo
 if(k>0)call wm_emit(LIS_localpet,wm_clock(),phase,wm_member,k,v)
end subroutine

subroutine wm_hymap_banner(marker,n)
 use WM_hymap_bridge_util
 implicit none
 character(*),intent(in)::marker
 integer,intent(in)::n
 logical,save::printed=.false.
 if(.not.wm_active(wm_clock()))return
 if(LIS_localpet==0.and..not.printed)then
  write(*,*)trim(marker),' MPI=',LIS_npes,' members=',LIS_rc%nensem(n),' selected=',count(selected)
  printed=.true.
 endif
end subroutine
