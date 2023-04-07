subroutine potrefewald(nG,Gval,Gvec,alpha,potref) 

! This routine computes the Fourier components of the 3D potential by combining the isotropic single-site contributions

implicit none
integer:: iG
real*8:: Gx,Gy,Gz,alpha,fourierpot,G2
integer:: nG
real*8,dimension(nG):: Gval
real*8,dimension(3,nG):: Gvec
real*8:: potref

!f2py intent(in) nG,Gval,Gvec,alpha 
!f2py intent(out) potref
!f2py depend(nG) Gval,Gvec

potref = 0.d0
!$OMP PARALLEL DEFAULT(private) &
!$OMP SHARED(nG,Gval,Gvec,alpha,potref)
!$OMP DO SCHEDULE(dynamic)
do iG=2,nG
   Gx = Gvec(1,iG)
   Gy = Gvec(2,iG)
   Gz = Gvec(3,iG)
   G2 = Gval(iG)**2
   fourierpot = dexp(-G2/(4.d0*alpha))/G2
   potref = potref + 2.0 * fourierpot
end do
!$OMP END DO
!$OMP END PARALLEL

return
end

