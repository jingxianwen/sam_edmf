! ==================================================================
      SUBROUTINE get_abcd(iin,jin,betap,betam,s,sumMs,tkin,a,b,c,d, massflux, kinflx)

      use grid
      use vars, only : rho, rhow 
      use sgs, only : sgs_field_sumM
      implicit none


!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! we need a,b,c,d - which all are vectors of order n 
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector

!-------------------------------------------------------------------

        ! input
        REAL, DIMENSION(nzm),    INTENT(in)  :: s, tkin    ! S and eddy viscosity at cell center
        REAL, DIMENSION(nz    ), INTENT(in)  :: sumMs      ! sum(MiSi) on faces
        LOGICAL,                 INTENT(in)  :: massflux   ! If true, then sumM is used otherwise sumM=0;
                                                           ! note that sumMs needs to be zero then (and this is 
                                                           ! ensured in edmf.f90)
        REAL, INTENT(in)                     :: kinflx,betap,betam
        INTEGER, INTENT(in) :: iin,jin
        ! output
        REAL, DIMENSION(nzm),    INTENT(out) :: a,b,c,d
        

        ! local
        INTEGER :: k
        REAL, DIMENSION(nz) :: tkf, sumMloc                   ! tk on faces
 
        ! these rhow values are actually never used
        tkf(1)  = 0.
        tkf(nz) = 0.

        if (massflux) then
          sumMloc = sgs_field_sumM(iin,jin,:,1)
        else
          sumMloc = 0.0
        end if


        DO k=1,nzm
           ! compute tk on faces
           IF (k.lt.nzm) tkf(k+1)  = 0.5*(tkin(k)+tkin(k+1))
           
           IF (k.gt.1.and.k.lt.nzm) THEN
              a(k) = rhow(k)/rho(k)/adz(k)/dz * betap *           &
                    (tkf(k)/adzw(k)/dz)
              !a(k) = rhow(k)/rho(k)/adz(k)/dz *                   &
              !      (betap*tkf(k)/adzw(k)/dz - 0.5 * sumMloc(k))
              b(k) = -1./dtn + betap* rhow(k+1)/rho(k)/adz(k)/dz * & 
                    (-tkf(k+1)/adzw(k+1)/dz ) - & 
                    betap*rhow(k)/rho(k)/adz(k)/dz*               &
                    (tkf(k)/adzw(k)/dz + sumMloc(k))
              !b(k) = -1./dtn +  rhow(k+1)/rho(k)/adz(k)/dz * & 
              !      (-betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1) ) - & 
              !       rhow(k)/rho(k)/adz(k)/dz*               &
              !      (betap*tkf(k)/adzw(k)/dz + 0.5 * sumMloc(k))
              c(k) = rhow(k+1)/rho(k)/adz(k)/dz * betap *         &
                     (tkf(k+1)/adzw(k+1)/dz + sumMloc(k+1))
              !c(k) = rhow(k+1)/rho(k)/adz(k)/dz *              &
              !       (betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1))
              d(k) = -s(k)/dtn - 1./rho(k)/adz(k)/dz * (                       &
                    tkf(k+1)*rhow(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
                   -tkf(k)*rhow(k)/adzw(k)/dz*betam*(s(k)-s(k-1))         &
                   +(rhow(k)*(sumMs(k) - betam*s(k)*sumMloc(k)))&
                   -(rhow(k+1)*(sumMs(k+1) - betam*s(k+1)*sumMloc(k+1))) )
              !d(k) = -s(k)/dtn - 1./rho(k)/adz(k)/dz * (                       &
              !      tkf(k+1)*rhow(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
              !     -tkf(k)*rhow(k)/adzw(k)/dz*betam*(s(k)-s(k-1))         &
              !     +(rhow(k)* sumMs(k))&
              !     -(rhow(k+1)*sumMs(k+1)) )
           ELSEIF (k.eq.1) THEN
              a(k) = 0.0
              b(k) = -1./dtn + betap* rhow(k+1)/rho(k)/adz(k)/dz * &
                    (-tkf(k+1)/adzw(k+1)/dz ) 
              !b(k) = -1./dtn +  rhow(k+1)/rho(k)/adz(k)/dz * &
              !      (-betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1) ) 
              c(k) = rhow(k+1)/rho(k)/adz(k)/dz * betap *         &
                     (tkf(k+1)/adzw(k+1)/dz + sumMloc(k+1))
              !c(k) = rhow(k+1)/rho(k)/adz(k)/dz *                  &
              !       (betap*tkf(k+1)/adzw(k+1)/dz + 0.5 * sumMloc(k+1))
              d(k) = -s(k)/dtn - 1./rho(k)/adz(k)/dz * (                       &
                    tkf(k+1)*rhow(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
                   -rhow(k+1)*(sumMs(k+1) - betam*s(k+1)*sumMloc(k+1)))
              !d(k) = -s(k)/dtn - 1./rho(k)/adz(k)/dz * (                       &
              !      tkf(k+1)*rhow(k+1)/adzw(k+1)/dz*betam*(s(k+1)-s(k))   &
              !     -rhow(k+1)*sumMs(k+1) )
              d(k) = d(k) - kinflx/adz(k)/dz
           ELSE
              a(k) = rhow(k)/rho(k)/adz(k)/dz * betap *                  &
                    (tkf(k)/adzw(k)/dz)
              !a(k) = rhow(k)/rho(k)/adz(k)/dz *                           &
              !      (betap*tkf(k)/adzw(k)/dz - 0.5 * sumMloc(k))
              b(k) = -1./dtn - betap* rhow(k)/rho(k)/adz(k)/dz *          &
                    (tkf(k)/adzw(k)/dz + sumMloc(k) )
              !b(k) = -1./dtn - rhow(k)/rho(k)/adz(k)/dz *          &
              !      (betap * tkf(k)/adzw(k)/dz + 0.5 * sumMloc(k) )
              c(k) = 0.0
              d(k) = -s(k)/dtn + 1./rho(k)/adz(k)/dz * (              &
                    tkf(k)*rhow(k)/adzw(k)/dz*betam*(s(k)-s(k-1))&
                   -rhow(k)*(sumMs(k) - betam*s(k)*sumMloc(k)))
              !d(k) = -s(k)/dtn + 1./rho(k)/adz(k)/dz * (              &
              !      tkf(k)*rhow(k)/adzw(k)/dz*betam*(s(k)-s(k-1))&
              !     -rhow(k)*sumMs(k) )
           END IF
        ENDDO


      END SUBROUTINE get_abcd

