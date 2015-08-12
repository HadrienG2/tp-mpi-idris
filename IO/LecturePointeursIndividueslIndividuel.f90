!**********************************************************************
!**********************************************************************
! Ce programme accède aux 484 valeurs de donnees.dat en parallèle,
! par une méthode utilisant des pointeurs individuels dans le fichier, en
! mode individuel, et en extrait 121 valeurs
!
! Puis il écrit les valeurs lues dans un fichier dont le nom reflète
! la méthode utilisée et le nom du processus
!
! On suppose une exécution par 4 processus.
!**********************************************************************
!**********************************************************************

program lecture_pointeurs_individuels_individuel

    use mpi
    implicit none
    
    integer, parameter :: requested_mpi_version = 3, requested_mpi_subversion = 0
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: nb_valeurs_par_processus = 121
    
    integer, dimension(nb_valeurs_par_processus) :: valeurs_lues
    
    character(len=*), parameter :: nom_fichier_source = "donnees.dat"
    character(len=14), parameter :: base_nom_fichier_cible = "fichier_ptinin"
    character(len=4), parameter :: ext_fichier_cible = ".dat"
    character(len=*), parameter :: format_nom_fichier_cible = "(A14,I0,A4)"
    
    integer :: return_code, world_rank, world_size
    integer :: file_mode, file_descriptor
    integer :: octets_par_entier
    
    integer(kind=MPI_OFFSET_KIND) :: offset_vue
    
    character(len=len(base_nom_fichier_cible)+len(ext_fichier_cible)+1) nom_fichier_cible
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations du rang du processus dans le communicateur WORLD, et de la taille de ce dernier
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Vérification des préconditions du programme (version de MPI, nombre de processus...)
    call ensure_preconditions_OK
    
    ! Ouverture du fichier source
    file_mode = MPI_MODE_RDONLY
    call MPI_FILE_OPEN(MPI_COMM_WORLD, &
                       nom_fichier_source, &
                       file_mode, &
                       MPI_INFO_NULL, &
                       file_descriptor, &
                       return_code)
                       
    ! Lecture via des pointeurs individuels en mode individuel (avec une vue appropriée)
    call MPI_TYPE_SIZE(MPI_INTEGER, octets_par_entier, return_code)
    offset_vue = octets_par_entier*nb_valeurs_par_processus*world_rank
    call MPI_FILE_SET_VIEW(file_descriptor, &
                           offset_vue, &
                           MPI_INTEGER, &
                           MPI_INTEGER, &
                           "native", &
                           MPI_INFO_NULL, &
                           return_code)
    call MPI_FILE_READ(file_descriptor, &
                       valeurs_lues(1), &
                       nb_valeurs_par_processus, &
                       MPI_INTEGER, &
                       MPI_STATUS_IGNORE, &
                       return_code)
    
    ! Fermeture du fichier source
    call MPI_FILE_CLOSE(file_descriptor, &
                        return_code)
                        
    ! Ouverture du fichier cible
    write (nom_fichier_cible, format_nom_fichier_cible) base_nom_fichier_cible, world_rank, ext_fichier_cible
    print *, "Pour le processus ", world_rank, ", le nom de sortie est ", nom_fichier_cible
    file_mode = MPI_MODE_WRONLY + &
                MPI_MODE_CREATE + &
                MPI_MODE_UNIQUE_OPEN + &
                MPI_MODE_SEQUENTIAL
    call MPI_FILE_OPEN(MPI_COMM_SELF, &
                       nom_fichier_cible, &
                       file_mode, &
                       MPI_INFO_NULL, &
                       file_descriptor, &
                       return_code)
    
    ! Ecriture des données lues précédemment
    call MPI_FILE_WRITE_ORDERED(file_descriptor, &
                                valeurs_lues(1), &
                                nb_valeurs_par_processus, &
                                MPI_INTEGER, &
                                MPI_STATUS_IGNORE, &
                                return_code)
    
    ! Fermeture du fichier cible
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
        if (world_size /= 4) then
            if (display_errors) print *, "ERREUR : On doit avoir exactement 4 processus pour ce programme"
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
    
end program ! Fin du programme principal
