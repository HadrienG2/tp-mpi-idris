!**********************************************************************
!**********************************************************************
! Ce programme n'est qu'une annexe aux programmes voisins de lecture
! parallèle de fichiers. Il sert à créer un fichier comportant 484
! valeurs entières.
!
! On suppose une exécution par un seul processus.
!**********************************************************************
!**********************************************************************

program creer_fichier

    use mpi
    implicit none
    
    integer, parameter :: requested_mpi_version = 3, requested_mpi_subversion = 0
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: nb_valeurs = 484
    
    integer :: i
    
    integer, dimension(nb_valeurs), parameter :: valeurs_a_ecrire = (/ (i, i=1,nb_valeurs) /)
    
    character(len=*), parameter :: nom_fichier = "donnees.dat"
    
    integer :: return_code, world_rank, world_size
    integer :: file_mode, file_descriptor
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations du rang du processus dans le communicateur WORLD, et de la taille de ce dernier
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Vérification des préconditions du programme (version de MPI, nombre de processus...)
    call ensure_preconditions_OK
    
    ! Ouverture du fichier cible
    file_mode = MPI_MODE_WRONLY + &
                MPI_MODE_CREATE + &
                MPI_MODE_UNIQUE_OPEN + &
                MPI_MODE_SEQUENTIAL
    call MPI_FILE_OPEN(MPI_COMM_WORLD, &
                       nom_fichier, &
                       file_mode, &
                       MPI_INFO_NULL, &
                       file_descriptor, &
                       return_code)
                       
    ! Ecriture des valeurs souhaitées
    call MPI_FILE_WRITE_ORDERED(file_descriptor, &
                                valeurs_a_ecrire(1), &
                                nb_valeurs, &
                                MPI_INTEGER, &
                                MPI_STATUS_IGNORE, &
                                return_code)
    
    ! Fermeture du fichier
    call MPI_FILE_CLOSE(file_descriptor, &
                        return_code)
    
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
    ! Cette fonction utilise world_rank et world_size
    function process_amount_preconditions_ok() result(check_success)
        logical :: check_success
        
        logical :: display_errors
        
        display_errors = (world_rank == root_process_world_rank)
        check_success = .FALSE.
        if (world_size /= 1) then
            if (display_errors) print *, "ERREUR : On doit n'avoir qu'un seul processus pour ce programme"
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
    
end program ! Fin du programme principal "creer_fichier"
