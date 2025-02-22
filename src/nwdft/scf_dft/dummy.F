C                                                                                                                                  
C     Read in header to get file name, check, then read in movecs.  Note                                                            
C     that the g_movecs_pad is ns_ao x ns_ao, with the last few columns                                                            
C     possibly 0 (if lindep), which is how the SCF code does it, which                                                              
C     is unlike my way, which has ns_mo x ns_mo.                                                                                    
C                                                                                                                                  
C                                                                                                                                                                            
        call ga_zero (g_movecs_pad)

        call rt_tddft_movecs_fname (params, rt_movecs_fname)

        if (.not. movecs_read_header (rt_movecs_fname, mo_title,
     $       mo_basis_name, mo_scftype, mo_nbf, mo_nsets,
     $       mo_ns_mo, nsets))
     $       call errquit (pname//"Failed to read movecs header", 0, 0)


C                                                                                                                                                                            
C     Check that movecs are legit.                                                                                                                                          
C                                                                                                                                                                            
        if (mo_scftype .ne. "dft")
     $       call errquit (pname//
     $       'Initial movecs should have scftype "dft"', 0, 0)

        if (mo_nbf .ne. params%ns_ao)
     $       call errquit (pname//
     $       'Initial movecs wrong size: mo_nbf /= ns_ao', mo_nbf, 0)

        do is = 1, nsets
           if (mo_ns_mo(is) .ne. params%ns_mo)
     $          call errquit (pname//
     $          'Initial movecs wrong size: mo_ns_mo /= ns_mo',
     $          mo_ns_mo(is), 0)
        enddo

        if (mo_nsets .ne. nsets)
     $       call errquit (pname//"Wrong number of initial movecs,",
     $       mo_nsets, 0)


        call rt_tddft_movecs_print_header (params, rt_movecs_fname,
     $       mo_title, mo_basis_name, mo_scftype, mo_nbf,
     $       mo_nsets, mo_ns_mo)



C                                                                                                                                                                            
C     Allocate buffers and movecs (ns_ao x ns_ao padded with zero *not* ns_mo x ns_mo).                                                                                      
C                                                                                                                                                                            
        if (.not.ma_push_get(mt_dbl, params%ns_ao, 'occ', locc, iocc))
     &       call errquit(pname//'cannot allocate occ',0, MA_ERR)

        if (.not.ma_push_get(mt_dbl, params%ns_ao, 'evals',
     $       levals, ievals))
     &       call errquit(pname//'cannot allocate evals',0, MA_ERR)


C                                                                                                                                                                            
C     Read in movecs (note ispin).                                                                                                                                          
C                                                                                                                                                                            
        if (.not. movecs_read (rt_movecs_fname, 1, dbl_mb(iocc), !note 1!!!
     $       dbl_mb(ievals), g_movecs_pad))
     $       call errquit (pname//"Failed to read movecs data", 0, 0)

        call rt_tddft_movecs_print_evals (params,
     $       dbl_mb(iocc), dbl_mb(ievals))
