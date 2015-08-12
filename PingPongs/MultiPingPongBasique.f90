!**********************************************************************
!**********************************************************************
! Ping-pong d'entier à répétition entre les processus pairs et impairs.
! Le programme suppose qu'il existe un nombre pair de processus.
!
! Ceci est une implémentation "naïve" à base de SEND et RECV. Pas
! facile à maintenir vu la duplication, et pas forcément performante
! non plus (env. 45 secondes pour 20 processus à 1000 itérations)
!**********************************************************************
!**********************************************************************

program ping_pong

    use mpi
    implicit none
    
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: message_tag = 42
    integer, parameter :: sent_message = 42
    integer, parameter :: number_of_iterations = 1000
    
    integer :: return_code, world_rank, world_size
    integer :: buddy_world_rank, received_message, iteration
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations des propriétés du processus dans le communicateur WORLD
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Si le nombre de processus n'est pas pair, commettre un suicide de masse
    call ensure_even_amount_of_processes
    
    ! Les processus pair envoient un message à leurs voisins impairs, qui accusent réception et renvoient un message
    if (modulo(world_rank, 2) == 0) then
        buddy_world_rank = world_rank+1
        print *, "Sending ", number_of_iterations, " pings - ", world_rank, " -> ", buddy_world_rank
        do iteration = 1, number_of_iterations
            call MPI_SEND(sent_message, &
                          1, &
                          MPI_INTEGER, &
                          buddy_world_rank, &
                          message_tag, &
                          MPI_COMM_WORLD, &
                          return_code)
            call MPI_RECV(received_message, &
                          1, &
                          MPI_INTEGER, &
                          buddy_world_rank, &
                          message_tag, &
                          MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE, &
                          return_code)
        end do
    else
        buddy_world_rank = world_rank-1
        do iteration = 1, number_of_iterations
            call MPI_RECV(received_message, &
                          1, &
                          MPI_INTEGER, &
                          buddy_world_rank, &
                          message_tag, &
                          MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE, &
                          return_code)
            call MPI_SEND(sent_message, &
                          1, &
                          MPI_INTEGER, &
                          buddy_world_rank, &
                          message_tag, &
                          MPI_COMM_WORLD, &
                          return_code)
        end do
        print *, "Received ", number_of_iterations, " pongs - ", world_rank, " <- ", buddy_world_rank
    end if
    
    ! Finalization de MPI
    call MPI_FINALIZE(return_code)

contains

    ! Vérification du fait que le nombre de processus est pair, sinon suicide orchestré par le processus 0.
    ! Cette fonction utilise world_rank et world_size, et peut altérer return_code
    subroutine ensure_even_amount_of_processes    
        if (modulo(world_size, 2) /= 0) then
            if (world_rank == root_process_world_rank) then
                print *, "ERREUR : Il faut un nombre pair de processus pour cet algorithme"
                call MPI_ABORT(MPI_COMM_WORLD, 42, return_code)
            else
                call MPI_BARRIER(MPI_COMM_WORLD, return_code) ! On ne peut pas utiliser STOP ici, car MPI l'interprète comme un crash
            endif
        end if
    end subroutine ensure_even_amount_of_processes

end program ! Fin du programme principal "ping_pong"
