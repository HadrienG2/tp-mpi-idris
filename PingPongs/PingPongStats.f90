!**********************************************************************
!**********************************************************************
! Test d'échange de messages toujours plus gros entre processus pairs
! et impairs, avec mesure du temps d'envoi. Le programme termine à une
! taille de message maximale supportée par la machine.
!
! Le programme suppose, comme toute la famille PingPong, qu'il y a un
! nombre de processus pair. La sortie est plus facile pour 2 processus.
!
! On voit, encore une fois, que les échanges de petites tailles sont
! horriblement coûteux, en proportion de la vitesse à laquelle on peut
! faire les gros échanges.
!
! Comme précédemment, il y aurait probablement de la perf à gagner en
! utilisant un SENDRECV au lieu d'un couple (Send, Recv). Le code
! en sortirait aussi sans doute plus simple et plus clair.
!**********************************************************************
!**********************************************************************

program ping_pong

    use mpi
    implicit none
    
    integer, parameter :: root_process_world_rank = 0
    integer, parameter :: message_tag = 42
    integer, parameter :: exchange_repetitions = 40
    integer, parameter :: max_message_size = 100000000
    integer, parameter :: size_iterations_amount = nint(log10(real(max_message_size)))+1
    
    integer :: return_code, world_rank, world_size
    integer :: buddy_world_rank, received_message, message_size
    integer :: size_iteration, repetition_iteration
    
    double precision :: timer, average_time
    
    real, dimension(max_message_size) :: message
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations des propriétés du processus dans le communicateur WORLD
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Si le nombre de processus n'est pas pair, commettre un suicide de masse
    call ensure_even_amount_of_processes
    
    ! On travaille à taille croissante, multipliant la taille du message par 10 à chaque fois.
    ! Les processus pairs envoient un message à leurs voisins impairs, qui renvoient, et mesurent le temps d'échange.
    if (modulo(world_rank, 2) == 0) then
        buddy_world_rank = world_rank+1
        
        print *, "Crunching random numbers at process ", world_rank
        call random_number(message)
        
        message_size = 1
        do size_iteration = 1, size_iterations_amount
            print *, "Exchanging ", exchange_repetitions, " times ", message_size, " reals, ", &
                     world_rank, " -> ", buddy_world_rank
            timer = MPI_WTIME()
            do repetition_iteration = 1, exchange_repetitions
                call MPI_SEND(message, &
                              message_size, &
                              MPI_REAL, &
                              buddy_world_rank, &
                              message_tag, &
                              MPI_COMM_WORLD, &
                              return_code)
                call MPI_RECV(message, &
                              message_size, &
                              MPI_REAL, &
                              buddy_world_rank, &
                              message_tag, &
                              MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, &
                              return_code)
            end do
            timer = MPI_WTIME() - timer
            average_time = timer/exchange_repetitions
            print *, "Exchanged ", exchange_repetitions, " times ", message_size, " reals, average time ", average_time
            message_size = message_size * 10
        end do
    else
        buddy_world_rank = world_rank-1
        
        message_size = 1
        do size_iteration = 1, size_iterations_amount
            do repetition_iteration = 1, exchange_repetitions
                call MPI_RECV(message, &
                              message_size, &
                              MPI_REAL, &
                              buddy_world_rank, &
                              message_tag, &
                              MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, &
                              return_code)
                call MPI_SEND(message, &
                              message_size, &
                              MPI_REAL, &
                              buddy_world_rank, &
                              message_tag, &
                              MPI_COMM_WORLD, &
                              return_code)
            end do
            message_size = message_size * 10
        end do
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
