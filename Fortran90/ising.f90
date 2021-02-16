module ising_framework
    implicit none
contains
    subroutine ising(L, T, next_I, prev_I, MCS, m, chi, C)
    implicit none
    INTEGER (kind=4),intent(in) :: L, MCS, next_I(L), prev_I(L)
    real (kind=8), INTENT(IN) :: T
    real (kind=8), intent(out) :: m, chi, C
    integer (kind=4):: lattice(L,L), i,j,k, delta_E, no_of_probes
    real (kind=8) :: boltzmann_factor(5), E=0.0d0, m_k=0.0d0, magnetizations = 0.0d0, magnetizations2 = 0.0d0, energies=0.0d0,&
             energies2=0.0d0, r

    no_of_probes=(MCS-30000)/100
    lattice = 1
    
    do i=-8,8,4
        boltzmann_factor(i/4 + 3) = exp(-i / T)
    end do

    magnetizations = 0.0d0
    magnetizations2 = 0.0d0
    energies = 0.0d0
    energies2 = 0.0d0
    do k=1, MCS
        do j=1, L
            do i=1, L
                delta_E = 2 * lattice(i,j) * (lattice(next_I(i),j) + lattice(prev_I(i),j) + &
                                                lattice(i,next_I(j)) + lattice(i, prev_I(j)))
                if (delta_E <= 0) then
                    lattice(i,j) = -lattice(i,j)
                else
                    call random_number(r) 
                    if (r <= boltzmann_factor(delta_E/4 + 3)) then
                        lattice(i,j) = -lattice(i,j)
                    end if
                end if
            end do
        end do

        if ((30000 < k) .and. (mod(k,100)==0)) then
            m_k = 0.0d0
            E = 0.0d0
            do j=1,L
                do i=1,L
                    m_k = m_k + lattice(i,j)
                    E = E + lattice(i,j) * (lattice(next_I(i),j) + lattice(i,next_I(j)))
                end do
            end do
            m_k = abs(m_k / L**2)
            magnetizations = magnetizations + m_k
            magnetizations2 = magnetizations2 + m_k**2
            energies = energies + E
            energies2 = energies2 + E**2
        end if
    end do
    m = magnetizations / no_of_probes
    chi = L**2 / T * (magnetizations2 - magnetizations**2 / no_of_probes)/no_of_probes
    C = (energies2 - energies**2 / no_of_probes)/no_of_probes / L**2 / T**2    
    end subroutine ising
end module ising_framework

program ising2d
    use ising_framework
    implicit none
    INTEGER (kind=4), PARAMETER :: MCS=230000, L=20, no_of_probes=(MCS-30000)/100
    REAL(kind=8), PARAMETER :: T1=0.5d0, T2=3.5d0, dT=0.1d0
    INTEGER(kind=4) :: next_I(L), prev_I(L), i
    real(kind=8) :: T, m, chi, C

    
    do i=1,L
        next_I(i) = i+1
        prev_I(i) = i-1
    end do
    next_I(L) = 1
    prev_I(1) = L


    OPEN(1, file="m.dat")
    OPEN(2, file="chi.dat")
    OPEN(3, file="C.dat")

    !trzeba zmienić taktykę z while na for, żeby wrzucić do openmp
    T=T1
    do
        if (T>T2) exit
        T = T + dT
        call ising(L, T, next_I, prev_I, MCS, m, chi, C)
        WRITE(1,*) T, m
        WRITE(2,*) T, chi
        WRITE(3,*) T, C
    end do

end program ising2d