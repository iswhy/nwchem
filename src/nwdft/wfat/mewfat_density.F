c=================================================
c=================================================
      subroutine mewfat_calc_rdm0(nb, n0, ocmo0, rdm0)

      implicit none
      integer nb                                 ! [input], no. of basis functions.
      integer n0(2)                              ! [input], no. of electrons in the neutral.
      double precision ocmo0(nb, n0(1)+n0(2))    ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision rdm0(nb,nb,2)             ! [output], neutral RDM.

      integer is, j, i, io, mo
      
      
c     ==== Calculate RDM's ====
      do is = 1, 2
         do j = 1, nb
            do i = 1, j
c              ==== Neutral RDM ====
               rdm0(i,j,is) = 0.0d0
               do io = 1, n0(is)
                  mo = (is-1)*n0(1) + io
                  rdm0(i,j,is) = rdm0(i,j,is) + ocmo0(i,mo)*ocmo0(j,mo)
               enddo
               rdm0(j,i,is) = rdm0(i,j,is)
            enddo
         enddo
      enddo

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_rdm(nb, n0, n1, ocmo0, ocmo1, rdm0, rdm1)

      implicit none
      integer nb                                 ! [input], no. of basis functions.
      integer n0(2)                              ! [input], no. of electrons in the neutral.
      integer n1(2)                              ! [input], no. of electrons in the cation.
      double precision ocmo0(nb, n0(1)+n0(2))    ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))    ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision rdm0(nb,nb,2)             ! [output], neutral RDM.
      double precision rdm1(nb,nb,2)             ! [output], cation RDM.

      integer is, j, i, io, mo
      
      
c     ==== Calculate RDM's ====
      do is = 1, 2
         do j = 1, nb
            do i = 1, j
c              ==== Neutral RDM ====
               rdm0(i,j,is) = 0.0d0
               do io = 1, n0(is)
                  mo = (is-1)*n0(1) + io
                  rdm0(i,j,is) = rdm0(i,j,is) + ocmo0(i,mo)*ocmo0(j,mo)
               enddo
               rdm0(j,i,is) = rdm0(i,j,is)

c              ==== Cation RDM ====
               rdm1(i,j,is) = 0.0d0
               do io = 1, n1(is)
                  mo = (is-1)*n1(1) + io
                  rdm1(i,j,is) = rdm1(i,j,is) + ocmo1(i,mo)*ocmo1(j,mo)
               enddo
               rdm1(j,i,is) = rdm1(i,j,is)
            enddo
         enddo
      enddo

      end
c=================================================
      
