!**********************************************************************
!**********************************************************************
! Un petit hello world en MPI, parce qu'on en est tous passés par là ;)
!**********************************************************************
!**********************************************************************

program hello_mpi

    use mpi
    implicit none
    
    integer :: return_code, world_rank, world_size
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations des propriétés du processus dans le communicateur WORLD
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Affichage d'un message de statut
    print *, "Salutations ! Je suis le processus ",world_rank," sur un total de ",world_size
    
    ! Finalization de MPI
    call MPI_FINALIZE(return_code)
    
end  ! Fin du programme principal "hello_mpi"
