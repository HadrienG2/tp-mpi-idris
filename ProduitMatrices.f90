!**********************************************************************
!**********************************************************************
! Comment répartir un produit matriciel sur plusieurs processus ? Une
! méthode est d'envoyer une ligne de la matrice de gauche et une
! colonne de la matrice de droite à chaque processus, puis de leur
! faire échanger la colonne qu'ils ont reçue selon leurs besoins.
!
! Ce programme suppose que les matrices sont carrées, que leur ordre
! est multiple du nombre de processus actifs, et bien sûr que l'ordre
! des deux matrices est égal (sinon le produit n'a pas de sens)
!
! Pour le compiler avec des matrics 1024x1024, gfortran a besoin de
! l'option -fmax-array-constructor=1048576
!
! Sur ma machine dual-core, en -O3, les temps de calcul sont...
!   - 1.5s avec un seul processus (exécution séquentielle)
!   - 0.9s avec deux processus (optimal pour cette machine, 1.7x)
!   - 1.4s avec 4 processus (trop de processus, 1.1x)
!
! On voit qu'avec MPI, il est très coûteux de se planter dans son
! nombre de processus, contrairement au cas des threads où ça reste
! peu pénalisant avec OpenMP.
!**********************************************************************
!**********************************************************************

program produit_matrices

    use mpi
    implicit none
    
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: matrix_order = 1024
    integer, parameter :: message_tag = 42
    
    integer :: return_code, world_rank, world_size
    integer :: i, results_offset, initial_results_offset
    integer :: next_process_in_world, previous_process_in_world
    integer :: matrix_slice_width
    integer :: matrix_slice_type
    
    real, dimension(matrix_order, matrix_order), parameter :: left_matrix = &
        reshape((/ (i, i=1,matrix_order**2) /), (/ matrix_order, matrix_order /))
    real, dimension(matrix_order, matrix_order), parameter :: right_matrix = &
        reshape((/ (i, i=matrix_order**2,1,-1) /), (/ matrix_order, matrix_order /))

    real, dimension(matrix_order, matrix_order) :: transposed_left_matrix, results_matrix
    
    real, dimension(:, :), allocatable :: transposed_left_matrix_slice, right_matrix_slice, results_matrix_slice
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations du rang du processus dans le communicateur WORLD, et de la taille de ce dernier
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Vérification des préconditions sur les matrices utilisées
    call enforce_matrix_preconditions
    
    ! Transposition de la matrice gauche pour simplifier les communications collectives
    if (world_rank == root_process_world_rank) then
        transposed_left_matrix = transpose(left_matrix)
    end if
    
    ! Allocation des tableaux et définition des types MPI pour
    !   - Les tranches de matrice envoyées/reçues par le processus 0 (et échangées entre processus)
    !   - Le tampon où chaque processus esclave stocke ses résultats
    matrix_slice_width = matrix_order/world_size
    allocate(transposed_left_matrix_slice(matrix_order, matrix_slice_width), &
             right_matrix_slice(matrix_order, matrix_slice_width), &
             results_matrix_slice(matrix_order, matrix_slice_width))
    call MPI_TYPE_CONTIGUOUS(matrix_order*matrix_slice_width, &
                             MPI_REAL, &
                             matrix_slice_type, &
                             return_code)
    call MPI_TYPE_COMMIT(matrix_slice_type, return_code)
                      
    
    ! Répartition des colonnes des matrices gauche-transp et droite sur les différents processus
    call MPI_SCATTER(transposed_left_matrix(1,1), &
                     1, &
                     matrix_slice_type, &
                     transposed_left_matrix_slice(1,1), &
                     1, &
                     matrix_slice_type, &
                     root_process_world_rank, &
                     MPI_COMM_WORLD, &
                     return_code)
    call MPI_SCATTER(right_matrix(1,1), &
                     1, &
                     matrix_slice_type, &
                     right_matrix_slice(1,1), &
                     1, &
                     matrix_slice_type, &
                     root_process_world_rank, &
                     MPI_COMM_WORLD, &
                     return_code)
                     
    ! Produit réparti de matrices, avec transmission du bloc de matrice gauche de proche en proche au fil du calcul
    next_process_in_world = modulo(world_rank+1, world_size)
    previous_process_in_world = modulo(world_rank-1, world_size)
    initial_results_offset = world_rank*matrix_slice_width
    results_offset = initial_results_offset
    do
        ! Calcul de la portion de la matrice résultat correspondant aux données actuellement disponibles
        call matmul_transpleft(transposed_left_matrix_slice, &
                               right_matrix_slice, &
                               results_matrix_slice(results_offset+1:results_offset+matrix_slice_width, :))
        
        ! Passage à la région de la tranche de matrice résultat suivante
        results_offset = modulo(results_offset+matrix_slice_width, matrix_order)
        
        ! Est-ce qu'on a fini de calculer la tranche de matrice résultat pour ce processus ?
        if(results_offset == initial_results_offset) exit
        
        ! Sinon, passage de la tranche de matrice gauche courante au prochain voisin qui en aura besoin
        ! (processus de rang précédent dans le communicateur WORLD), et remplacement par la prochaine tranche gauche du
        ! voisin qui la possède actuellement (processus de rang suivant dans le communicateur WORLD)
        call MPI_SENDRECV_REPLACE(transposed_left_matrix_slice(1,1), &
                                  1, &
                                  matrix_slice_type, &
                                  previous_process_in_world, &
                                  message_tag, &
                                  next_process_in_world, &
                                  message_tag, &
                                  MPI_COMM_WORLD, &
                                  MPI_STATUS_IGNORE, &
                                  return_code)
    end do
    
    ! Collecte des résultats du produit de matrices sur le processus 0
    call MPI_GATHER(results_matrix_slice(1,1), &
                    1, &
                    matrix_slice_type, &
                    results_matrix(1,1), &
                    1, &
                    matrix_slice_type, &
                    root_process_world_rank, &
                    MPI_COMM_WORLD, &
                    return_code)
    
    ! Annonce de la fin du calcul par le processus 0
    if (world_rank == root_process_world_rank) then
        print *, "Produit de matrice terminé !"
    end if
    
    ! Libération des tableaux dynamiques
    deallocate(transposed_left_matrix_slice, right_matrix_slice, results_matrix_slice)
    
    ! Finalization de MPI
    call MPI_FINALIZE(return_code)

contains

    ! Version de matmul pour le cas où on a la transposée de la matrice de gauche
    ! Cette fonction suppose que les bornes inférieures de toutes les matrices sont à 1
    pure subroutine matmul_transpleft(transposed_left_matrix, right_matrix, results_matrix)
        real, dimension(:,:), intent(in) :: transposed_left_matrix, right_matrix
        real, dimension(:,:), intent(out) :: results_matrix
        
        integer :: i, j
        
        forall (i = 1:size(results_matrix,1), j = 1:size(results_matrix,2))
            results_matrix(i,j) = dot_product(transposed_left_matrix(:,i), right_matrix(:,j))
        end forall
    end subroutine matmul_transpleft

    ! Affichage d'une matrice de réels sur la sortie standard, avec un saut de ligne final
    subroutine print_real_matrix(matrix)
        real, dimension(:,:), intent(in) :: matrix
        
        do i = lbound(matrix, 1), ubound(matrix, 1)
            print *, matrix(i, :)
        end do
        print *, ""
    end subroutine print_real_matrix    

    ! Test des préconditions sur les matrices. Si échec, affichage d'un message d'erreur explicatif.
    ! Cette fonction utilise world_rank et world_size, et peut altérer return_code
    function matrix_preconditions_ok() result(check_success)
        logical :: check_success
        
        logical :: display_errors
        
        ! Vérification du fait que les matrices sont carrées et d'ordre égal
        display_errors = (world_rank == root_process_world_rank)
        check_success = .FALSE.
        if ((size(left_matrix, 1) /= matrix_order).OR.(size(left_matrix, 2) /= matrix_order)) then
            if (display_errors) print *, "ERREUR : Les matrices d'entrées doivent être carrées, d'ordre matrix_order"
        elseif (all(shape(left_matrix) /= shape(right_matrix))) then
            if (display_errors) print *, "ERREUR : Les matrices d'entrées doivent être de même forme"
        elseif (modulo(matrix_order, world_size) /= 0) then
            if (display_errors) print *, "ERREUR : L'ordre des matrices doit être multiple du nombre de processus"
        else
            check_success = .TRUE.
        end if
    end function matrix_preconditions_ok

    ! Vérification des préconditions sur les matrices, en cas d'échec suicide orchestré par le processus 0.
    ! Cette fonction utilise world_rank et world_size, et peut altérer return_code
    subroutine enforce_matrix_preconditions    
        if (.NOT.matrix_preconditions_ok()) then
            if (world_rank == root_process_world_rank) then
                print *, "ERREUR : Le test des préconditions sur les matrices a échoué, le programme ne peut continuer"
                call MPI_ABORT(MPI_COMM_WORLD, 42, return_code)
            else
                call MPI_BARRIER(MPI_COMM_WORLD, return_code) ! On ne peut pas utiliser STOP ici, car MPI l'interprète comme un crash
            endif
        end if
    end subroutine enforce_matrix_preconditions

end program ! Fin du programme principal "produit_matrices"
