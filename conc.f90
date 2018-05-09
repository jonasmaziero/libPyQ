real(8) function concurrence_2qb(rho)  ! Returns the entanglement measure concurrence, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(4,4)  ! Density matrix we want to compute the concurrence
complex(8) :: R(4,4), rho_tilde(4,4), s2_kp_s2(4,4)  ! Auxiliary matrices
complex(8) :: egv(4) ! Eigenvalues of R = rho*rho^tilde
real(8) :: egv_max  ! The greater eigenvalue of R
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

!f2py intent(in) :: rho

call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, s2_kp_s2)
rho_tilde = matmul( matmul(s2_kp_s2,conjg(rho)) , s2_kp_s2 )
R = matmul(rho,rho_tilde) ;   call lapack_zgeev('N', 4, R, egv)

egv_max = max( real(egv(1)), real(egv(2)), real(egv(3)), real(egv(4)) )
 concurrence_2qb = max( 0.d0, (2.d0*sqrt(egv_max)-sqrt(real(egv(1)))-sqrt(real(egv(2)))-sqrt(real(egv(3)))-sqrt(real(egv(4)))))

end

subroutine kronecker_product_c(M1, nr1, nc1, M2, nr2, nc2, M1_kp_M2)  ! Returns the tensor product of two general complex matrices
implicit none
integer :: nr1, nc1, nr2, nc2  ! Number of rows and columns of the two matrices
complex(8) :: M1(1:nr1,1:nc1), M2(1:nr2,1:nc2)  ! Matrices to take the tensor product of
complex(8) :: M1_kp_M2(1:nr1*nr2,1:nc1*nc2)  ! Matrix containing the tensor product of M1 and M2
integer :: i, j  ! Auxiliary variables for counters

M1_kp_M2 = 0.d0 ;   forall ( i = 1:nr1 , j = 1:nc1 ) M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2

end

subroutine lapack_zgeev(JOBVR, N, A, Wc)  ! Calls LAPACK's eigensolver for GENERAL complex matrices
! ZGEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
!character(1) :: JOBVL  !  JOBVL is CHARACTER*1; = 'N': left eigenvectors of A are not computed; = 'V': left eigenvectors of are computed.
character(1) :: JOBVR  !  JOBVR is CHARACTER*1; = 'N': right eigenvectors of A are not computed;  = 'V': right eigenvectors of A are computed.
character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA = N  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(1:N,1:N)  ! A is COMPLEX*16 array, dimension (LDA,N); On entry, the N-by-N matrix A. On exit, A has been overwritten.
complex(8) :: Wc(1:N)  ! Wc is COMPLEX*16 array, dimension (N). Wc contains the computed eigenvalues.
!integer :: LDVL = N  !LDVL is INTEGER. The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
complex(8) :: VL(1:N,1:N)  ! VL is COMPLEX*16 array, dimension (LDVL,N); If JOBVL = 'V', the left eigenvectors u(j) are stored one
                        ! after another in the columns of VL, in the same order as their eigenvalues.
                        ! If JOBVL = 'N', VL is not referenced. u(j) = VL(:,j), the j-th column of VL.
!integer :: LDVR = N  ! LDVR is INTEGER. The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
complex(8) :: VR(1:N,1:N)  ! VR is COMPLEX*16 array, dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one
                                     ! after another in the columns of VR, in the same order their eigenvalues.
                                     ! If JOBVR = 'N', VR is not referenced. v(j) = VR(:,j), the j-th column of VR.

!integer :: LWORK = 2*N  ! LWORK is INTEGER; The dimension of the array WORK.  LWORK >= max(1,2*N). For good performance, LWORK must generally be larger.
                  ! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                  ! this value as the first entry of the WORK array, and no error related to LWORK is issued by XERBLA.
complex(8) :: WORK(1:2*N)  ! WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
real(8) :: RWORK(1:2*N)  ! RWORK is DOUBLE PRECISION array, dimension (2*N)
integer :: INFO   ! INFO is INTEGER
                  ! = 0:  successful exit
                  ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                  ! > 0:  if INFO = i, the QR algorithm failed to compute all the
                  ! eigenvalues, and no eigenvectors have been computed; elements and i+1:N of W contain eigenvalues which have converged.

 call zgeev ('N',   JOBVR, N, A, N,   Wc, VL, N,    VR, N,    WORK, 2*N,   RWORK, INFO)
!call zgeev (JOBVL, JOBVR, N, A, LDA, Wc, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)

end


subroutine pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)  ! Defines the three Pauli's matrices and the identity matrix
implicit none
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

sigma_0 = 0.d0 ;   sigma_0(1,1) = 1.d0 ;   sigma_0(2,2) = 1.d0
sigma_1 = 0.d0 ;   sigma_1(1,2) = 1.d0 ;   sigma_1(2,1) = 1.d0
sigma_2 = 0.d0 ;   sigma_2(1,2) = -(0.d0,1.d0) ;   sigma_2(2,1) = (0.d0,1.d0)
sigma_3 = 0.d0 ;   sigma_3(1,1) = 1.d0 ;   sigma_3(2,2) = -1.d0

end
