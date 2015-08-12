!**********************************************************************
!**********************************************************************
! Ce module encapsule toutes les fonctions de gestion de MPI liées au
! solveur d'équations de Poisson.
!**********************************************************************
!**********************************************************************

module gestion_mpi

    use mpi
    implicit none
    
    integer, parameter :: requested_mpi_version = 3, requested_mpi_subversion = 0
    integer, parameter :: requested_thread_support = MPI_THREAD_FUNNELED
    integer, parameter :: root_process_world_rank = 0
    
    integer :: return_code, world_rank, world_size, thread_support
    
contains
    
    ! Initialisation du module
    subroutine initialisation_gestion_mpi
        ! Initialization de MPI
        call MPI_INIT_THREAD(requested_thread_support, thread_support, return_code)
        
        ! Récupérations du rang du processus dans le communicateur WORLD, et de la taille de ce dernier
        call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
        
        ! Vérification des préconditions sur l'implémentation MPI
        call ensure_preconditions_OK
    end subroutine initialisation_gestion_mpi
    
    ! Finalisation du module
    subroutine finalisation_gestion_mpi
        ! Finalisation de MPI
        call MPI_FINALIZE(return_code)
    end subroutine finalisation_gestion_mpi
    
    ! Vérification des préconditions du programme (version MPI, nombre de processus...)
    subroutine ensure_preconditions_OK
        if (.NOT.mpi_implementation_ok()) then
            if (world_rank == root_process_world_rank) then
                print *, "ERREUR : Le test des préconditions a échoué, le programme ne peut continuer"
                call MPI_ABORT(MPI_COMM_WORLD, 42, return_code)
            else
                call MPI_BARRIER(MPI_COMM_WORLD, return_code) ! On ne peut pas utiliser STOP ici, car MPI l'interprète comme un crash
            endif
        end if
    end subroutine ensure_preconditions_OK
    
    ! Test des préconditions sur l'implémentation de MPI. Si échec, affichage d'un message d'erreur explicatif.
    function mpi_implementation_ok() result(check_success)
        logical :: check_success
        
        logical :: display_errors
        
        display_errors = (world_rank == root_process_world_rank)
        check_success = .FALSE.
        if (MPI_VERSION < requested_mpi_version) then
            if (display_errors) print *, "ERREUR : L'implémentation MPI suit une norme de version trop ancienne"
        elseif ((MPI_VERSION == requested_mpi_version).AND.(MPI_SUBVERSION < requested_mpi_subversion)) then
            if (display_errors) print *, "ERREUR : L'implémentation MPI suit une norme de sous-version trop ancienne"
        elseif (thread_support < requested_thread_support) then
            if (display_errors) print *, "ERREUR : L'implémentation MPI ne supporte pas assez bien les threads"
        else
            check_success = .TRUE.
        end if
    end function mpi_implementation_ok
    
end module gestion_mpi
