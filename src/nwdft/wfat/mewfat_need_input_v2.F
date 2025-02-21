cdebugc=======================================================
cdebugc=======================================================
cdebug      subroutine mewfat_need_input()
cdebug      implicit none
cdebug#include "errquit.fh"
cdebug
cdebug      call errquit('stopping in mewfat_need_input', 52, UERR)
cdebug      end
cdebugc=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_dummy()
      implicit none
#include "errquit.fh"

      call errquit('stopping in dummy', 52, UERR)
      end
c=======================================================      
