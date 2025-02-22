
      subroutine wfat_load_intgl()

      implicit none


c     ==== Other quantities to be printed ====
      atom coords.
      basis
      convergence crit.
      functional
      mulliken analysis
      orbital property
      

c     ==== Inputs to be retained ====
      gridqual   ! optional
      paraqn     ! optional
      movec_id   ! optional
      lmax       ! optional


c     ==== Inputs to be deleted ====
      fieldmag     ! mandatory when "print iony tmd tmd_atom" 
      beta         ! mandatory when "print srfc iony"
      gamma        ! mandatory when "print srfc iony"
      beta_tmd     ! mandatory when "print tmd tmd_atom"
      gamma_tmd    ! mandatory when "print tmd tmd_atom"
      print        ! optional (default is "iony")
      nktmd        ! mandatory when "print tmd tmd_atom"


      end
