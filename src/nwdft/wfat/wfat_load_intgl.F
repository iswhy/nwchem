
cOLD      subroutine wfat_load_intgl()
cOLD
cOLD      implicit none
cOLD
cOLD
cOLDc     ==== Other quantities to be printed ====
cOLD      atom coords.
cOLD      basis
cOLD      convergence crit.
cOLD      functional
cOLD      mulliken analysis
cOLD      orbital property
cOLD      
cOLD
cOLDc     ==== Inputs to be retained ====
cOLD      gridqual   ! optional
cOLD      paraqn     ! optional
cOLD      movec_id   ! optional
cOLD      lmax       ! optional
cOLD
cOLD
cOLDc     ==== Inputs to be deleted ====
cOLD      fieldmag     ! mandatory when "print iony tmd tmd_atom" 
cOLD      beta         ! mandatory when "print srfc iony"
cOLD      gamma        ! mandatory when "print srfc iony"
cOLD      beta_tmd     ! mandatory when "print tmd tmd_atom"
cOLD      gamma_tmd    ! mandatory when "print tmd tmd_atom"
cOLD      print        ! optional (default is "iony")
cOLD      nktmd        ! mandatory when "print tmd tmd_atom"
cOLD
cOLD
cOLD      end
