c $Id$
      subroutine wfat_et_det(handle,n,detval)
c     DESCRIPTION:
c     This subroutine is derived from the et_det subroutine where the
c     argument 'op', which gives the option whether to transpose or
c     not, has been removed.
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
#include "stdio.fh"
      integer handle
      integer n,i,j,info
      integer l_ipvt,k_ipvt,l_a,k_a,l_work,k_work
      double precision detval,det(2),ten
c
      ten = 10.0d0
      if (.not.ma_push_get(mt_int,n,'ipvt',
     1  l_ipvt,k_ipvt)) call errquit
     2  ('et_det: failed to allocate ipvt',0,0)
      if (.not.ma_push_get(mt_dbl,n*n,'a',
     1  l_a,k_a)) call errquit
     2  ('et_det: failed to allocate a',0,0)
      if (.not.ma_push_get(mt_dbl,n,'work',
     1  l_work,k_work)) call errquit
     2  ('et_det: failed to allocate work',0,0)
c
      call ga_get(handle, 1, n, 1, n, dbl_mb(k_a), n)
      
      call dgefa(dbl_mb(k_a),n,n,int_mb(k_ipvt),info)
      call dgedi(dbl_mb(k_a),n,n,int_mb(k_ipvt),det,dbl_mb(k_work),10)
      detval =  det(1) * ten**det(2)
c
      if (.not.ma_pop_stack(l_work)) call errquit
     1  ('et_det: failed to deallocate work',0,0)
      if (.not.ma_pop_stack(l_a)) call errquit
     1  ('et_det: failed to deallocate a',0,0)
      if (.not.ma_pop_stack(l_ipvt)) call errquit
     1  ('et_det: failed to deallocate ipvt',0,0)
      call ga_sync()
      return
      end

