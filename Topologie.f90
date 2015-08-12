!**********************************************************************
!**********************************************************************
! Jouons un peu avec les topologies ! Cet outil très pratique est
! fourni par MPI pour gérer simplement et optimiser la performance
! d'un système de processus ayant des "voisins" avec lesquels ils
! communiquent naturellement.
!
! Ce programme suppose l'existence d'un nombre pair de processus.
!**********************************************************************
!**********************************************************************

program topologie

    use mpi
    implicit none
    
    integer, parameter :: requested_mpi_version = 3, requested_mpi_subversion = 0
    integer, parameter :: root_process_world_rank = 0
    
    integer :: return_code, world_rank, mpi_world_size
    integer :: cartesian_line_length
    integer :: cartesian_communicator, cartesian_line_communicator
    integer :: cartesian_line_rank
    
    integer, dimension(2) :: cartesian_dims, my_cartesian_coords
    logical, dimension(2) :: cartesian_periodicity
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations du rang du processus dans le communicateur WORLD, et de la taille de ce dernier
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_world_size, return_code)
    
    ! Vérification des préconditions du programme (version de MPI, nombre de processus...)
    call ensure_preconditions_OK
    
    ! Création d'une topologie cartésienne non périodique à deux lignes
    cartesian_line_length = mpi_world_size/2
    call MPI_CART_CREATE(MPI_COMM_WORLD, &
                         2, &
                         (/ 2, cartesian_line_length /), &
                         (/ .false., .false. /), &
                         .true., &
                         cartesian_communicator, &
                         return_code)
        
    ! Identification de la la position du processus dans la topologie (et autres infos inutiles)
    call MPI_CART_GET(cartesian_communicator, &
                      2, &
                      cartesian_dims, &
                      cartesian_periodicity, &
                      my_cartesian_coords, &
                      return_code)
    
    ! Affichage d'informations sur le communicateur cartésien créé
    if (world_rank == root_process_world_rank) then
        print *, "On a créé un communicateur cartésien de dimensions ", cartesian_dims
        print *, "...de périodicité ", cartesian_periodicity
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, return_code)
    print *, "...où les coordonnées du processus ", world_rank, " sont ", my_cartesian_coords
    
    ! Subdivision en deux communicateurs avec MPI_COMM_SPLIT
    call MPI_COMM_SPLIT(cartesian_communicator, &
                        my_cartesian_coords(1), &
                        my_cartesian_coords(2), &
                        cartesian_line_communicator, &
                        return_code)
                        
    ! Identification du processus dans les communicateurs créés
    call MPI_COMM_RANK(cartesian_line_communicator, cartesian_line_rank, return_code)
    
    ! Affichage d'informations sur le statut de chaque processus dans les communicateurs
    if (world_rank == root_process_world_rank) then
        print *, "On a ensuite divisé le communicateur en deux, qui suivent les lignes de la topologie"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, return_code)
    print *, "...où le processus ", world_rank, &
             " est de rang ", cartesian_line_rank
        
    ! Finalization de MPI
    call MPI_FINALIZE(return_code)

contains

    ! Test des préconditions sur la version de MPI. Si échec, affichage d'un message d'erreur explicatif.
    ! Cette fonction utilise world_rank
    function mpi_version_preconditions_ok() result(check_success)
        logical :: check_success
        
        logical :: display_errors
        
        display_errors = (world_rank == root_process_world_rank)
        check_success = .FALSE.
        if (MPI_VERSION < requested_mpi_version) then
            if (display_errors) print *, "ERREUR : L'implémentation MPI suit une norme de version trop ancienne"
        elseif ((MPI_VERSION == requested_mpi_version).AND.(MPI_SUBVERSION < requested_mpi_subversion)) then
            if (display_errors) print *, "ERREUR : L'implémentation MPI suit une norme de sous-version trop ancienne"
        else
            check_success = .TRUE.
        end if
    end function mpi_version_preconditions_ok
    
    ! Test des préconditions sur le nombre de processus. Si échec, affichage d'un message d'erreur explicatif.
    ! Cette fonction utilise world_rank et mpi_world_size
    function process_amount_preconditions_ok() result(check_success)
        logical :: check_success
        
        logical :: display_errors
        
        display_errors = (world_rank == root_process_world_rank)
        check_success = .FALSE.
        if (modulo(mpi_world_size, 2) /= 0) then
            if (display_errors) print *, "ERREUR : Le nombre de processus doit être pair pour ce programme"
        else
            check_success = .TRUE.
        end if
    end function process_amount_preconditions_ok

    ! Vérification des préconditions du programme (version MPI, nombre de processus...)
    ! Cette fonction utilise world_rank, et peut altérer return_code
    subroutine ensure_preconditions_OK
        if (.NOT.(mpi_version_preconditions_ok().AND.process_amount_preconditions_ok())) then
            if (world_rank == root_process_world_rank) then
                print *, "ERREUR : Le test des préconditions a échoué, le programme ne peut continuer"
                call MPI_ABORT(MPI_COMM_WORLD, 42, return_code)
            else
                call MPI_BARRIER(MPI_COMM_WORLD, return_code) ! On ne peut pas utiliser STOP ici, car MPI l'interprète comme un crash
            endif
        end if
    end subroutine ensure_preconditions_OK

end program ! Fin du programme principal "topologie"
