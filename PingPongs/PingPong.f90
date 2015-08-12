!**********************************************************************
!**********************************************************************
! Ping-pong entre les processus pairs et les processus impairs. Le
! programme suppose qu'il existe un nombre pair de processus.
!**********************************************************************
!**********************************************************************

program ping_pong

    use mpi
    implicit none
    
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: message_tag = 42
    
    integer :: return_code, world_rank, world_size
    integer :: message
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations des propriétés du processus dans le communicateur WORLD
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Si le nombre de processus n'est pas pair, commettre un suicide de masse
    call ensure_even_amount_of_processes
    
    ! Les processus pair envoient un message à leurs voisins impairs, qui accusent réception
    if (mod(world_rank, 2) == 0) then
        message = 42
        print *, "Ping ! ", world_rank, " -> ", world_rank+1
        call MPI_SEND(message, &
                      1, &
                      MPI_INTEGER, &
                      world_rank+1, &
                      message_tag, &
                      MPI_COMM_WORLD, &
                      return_code)
    else
        call MPI_RECV(message, &
                      1, &
                      MPI_INTEGER, &
                      world_rank-1, &
                      message_tag, &
                      MPI_COMM_WORLD, &
                      MPI_STATUS_IGNORE, &
                      return_code)
        print *, "Pong ! ", world_rank, " <- ", world_rank-1
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
