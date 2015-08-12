!**********************************************************************
!**********************************************************************
! Pour démarrer en communications collectives, on essaye de tirer des
! nombres à pile ou face sur tous les processus jusqu'à ce que l'une
! de ces deux conditions soient atteintes :
!       - Tous les programmes ont tiré la même face de la pièce
!       - Le nombre d'itérations maximal a été atteint
!**********************************************************************
!**********************************************************************

program pile_ou_face

    use mpi
    implicit none
    
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: max_attempts = 1000
    
    real :: random_value
    
    logical :: coin_is_head, coin_is_tail, all_coins_are_head, all_coins_are_tail, all_coins_are_equal
    
    integer :: return_code, world_rank, iteration
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations du rang du processus dans le communicateur WORLD
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    
    ! Initialisation du générateur de nombres aléatoires
    call init_random_number_generator(world_rank)
    
    ! Lancer des pièces jusqu'à accord entre les processus ou dépassement du nombre maximum d'essais
    do iteration = 1, max_attempts
        ! Lancer une pièce
        call random_number(random_value)
        coin_is_head = (random_value < 0.5)
        coin_is_tail = .not.coin_is_head
        
        ! Afficher le résultat du lancer
        if (coin_is_head) then
            print *, "Pile, processus ", world_rank
        else
            print *, "Face, processus ", world_rank
        end if
        
        ! Utiliser des réductions réparties pour comparer les lancers de pièces
        ! (deux réductions sont nécessaires car un booléen ne contient pas assez d'information
        ! pour stocker les deux cas dans lesquels toutes les pièces sont identiques)
        call MPI_ALLREDUCE(coin_is_head, &
                           all_coins_are_head, &
                           1, &
                           MPI_LOGICAL, &
                           MPI_LAND, &
                           MPI_COMM_WORLD, &
                           return_code)
        call MPI_ALLREDUCE(coin_is_tail, &
                           all_coins_are_tail, &
                           1, &
                           MPI_LOGICAL, &
                           MPI_LAND, &
                           MPI_COMM_WORLD, &
                           return_code)
        all_coins_are_equal = all_coins_are_head.or.all_coins_are_tail
        
        ! Signaler qu'une réduction a eu lieu
        if (world_rank == root_process_world_rank) then
            print *, "------------REDUCTION----------"
        end if
        
        ! Si les résultats sont identiques, terminer le programme
        if(all_coins_are_equal) exit
    end do
    
    ! Faire afficher par le processus maître le résultat
    if (world_rank == root_process_world_rank) then
        if (all_coins_are_equal) then
            print *, "Les processus sont tombés d'accord après ", iteration, " lancers"
        else
            print *, "Les processus ne sont pas tombés d'accords après ", max_attempts, " lancers"
        end if
    end if
    
    ! Finalization de MPI
    call MPI_FINALIZE(return_code)
    
contains

    ! Initialisation du générateur de nombres aléatoires, via la date/heure et le rang du processus
    subroutine init_random_number_generator(world_rank)
        integer, intent(in) :: world_rank
        
        integer :: min_seed_size, seed_size, seed_index, seed_material_index
        integer, dimension(:), allocatable :: seed
        integer, parameter :: date_and_time_size = 8
        integer, parameter :: seed_material_size = date_and_time_size+1
        integer, dimension(1:seed_material_size) :: seed_material
        
        ! Préparer le tableau de graine, pour qu'elle vérifie les contraintes de RANDOM_SEED et soit
        ! assez grande pour accueillir tout le matériau de départ
        call random_seed(size=min_seed_size)
        seed_size = max(min_seed_size, seed_material_size)
        allocate(seed(1:seed_size))
        
        ! Préparer le matériau utilisé pour la génération de la graine
        seed_material(1) = world_rank
        call date_and_time(values=seed_material(2:date_and_time_size+1))
        
        ! Generer la graine depuis le matériau, en le périodisant au besoin
        do seed_index = 1, seed_size
            seed_material_index = mod(seed_index, seed_material_size)+1
            seed(seed_index) = seed_material(seed_material_index)
        end do
        
        ! Utiliser cette graine pour initialiser le générateur de nombres aléatoires
        call random_seed(put=seed)
        
        ! Libérer le tableau utilisé pour stocker la graine
        deallocate(seed)
    end subroutine init_random_number_generator

end program ! Fin du programme principal "pile_ou_face"
