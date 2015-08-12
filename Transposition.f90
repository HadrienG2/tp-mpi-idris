!**********************************************************************
!**********************************************************************
! Une propriété amusante des communications MPI est qu'elles n'opèrent
! pas une vérification rigoureuse de la compatibilité entre les types
! de données envoyés et reçus. En effet, quand on cherche une perf
! maximale, l'information de type peut être trop coûteuse à échanger.
!
! Tant que les quantités de données échangées sont équivalentes, il
! est donc possible de convertir une donnée d'un type en une donnée
! d'un autre type au cours d'une communication. Cette possibilité est
! utilisée ici pour transposer une matrice.
!
! Ce programme est conçu pour être exécuté avec deux processus.
!**********************************************************************
!**********************************************************************

program transposition

    use mpi
    implicit none
    
    integer, parameter :: sender_world_rank = 0, receiver_world_rank = 1
    integer, parameter :: message_tag = 42
    integer, parameter :: dim_longue_matrice = 17, dim_courte_matrice = 4
    
    integer :: return_code, world_rank, world_size
    integer :: type_colonne_matrice_verticale, type_ligne_matrice_horizontale
    integer :: i
    
    integer, dimension(dim_longue_matrice, dim_courte_matrice), parameter :: matrice_verticale = &
        reshape((/ (i, i=1,dim_longue_matrice*dim_courte_matrice) /), (/ dim_longue_matrice, dim_courte_matrice /))

    integer, dimension(dim_courte_matrice, dim_longue_matrice) :: matrice_horizontale
    
    ! Initialization de MPI
    call MPI_INIT(return_code)
    
    ! Récupérations du rang du processus dans le communicateur WORLD, et de la taille de ce dernier
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, return_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, return_code)
    
    ! Vérification de la précondition sur le nombre de processus
    call ensure_there_are_two_processes
            
    ! Séparation entre activités d'envoi et réception
    if (world_rank == sender_world_rank) then
        ! Affichage de la matrice verticale avant envoi
        print *, "Envoi de la matrice verticale :"
        do i = 1, ubound(matrice_verticale, 1)
            print *, matrice_verticale(i, :)
        end do
        
        ! Définition du type matrice verticale
        call MPI_TYPE_CONTIGUOUS(dim_longue_matrice, &
                                 MPI_INTEGER, &
                                 type_colonne_matrice_verticale, &
                                 return_code)
        call MPI_TYPE_COMMIT(type_colonne_matrice_verticale, return_code)
        
        ! Envoi de la matrice verticale
        do i = 1, ubound(matrice_verticale, 2)
            call MPI_SEND(matrice_verticale(1, i), &
                          1, &
                          type_colonne_matrice_verticale, &
                          receiver_world_rank, &
                          message_tag, &
                          MPI_COMM_WORLD, &
                          return_code)
        end do
    else
        ! Définition du type matrice horizontale
        call MPI_TYPE_VECTOR(dim_longue_matrice, &
                             1, &
                             dim_courte_matrice, &
                             MPI_INTEGER, &
                             type_ligne_matrice_horizontale, &
                             return_code)
        call MPI_TYPE_COMMIT(type_ligne_matrice_horizontale, return_code)
        
        ! Réception de la matrice horizontale
        do i = 1, ubound(matrice_horizontale, 1)
            call MPI_RECV(matrice_horizontale(i, 1), &
                          1, &
                          type_ligne_matrice_horizontale, &
                          sender_world_rank, &
                          message_tag, &
                          MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE, &
                          return_code)
        end do
        
        ! Affichage de la matrice reçue
        print *, "Réception de la matrice horizontale :"
        do i = 1, ubound(matrice_horizontale, 1)
            print *, matrice_horizontale(i, :)
        end do
    end if
    
    ! Finalization de MPI
    call MPI_FINALIZE(return_code)

contains

    ! Vérification qu'on a deux processus, sinon suicide orchestré par le processus 0.
    ! Cette fonction utilise world_rank et world_size, et peut altérer return_code
    subroutine ensure_there_are_two_processes    
        if (world_size /= 2) then
            if (world_rank == 0) then
                print *, "ERREUR : Ce programme est conçu pour s'exécuter avec deux processus"
                call MPI_ABORT(MPI_COMM_WORLD, 42, return_code)
            else
                call MPI_BARRIER(MPI_COMM_WORLD, return_code) ! On ne peut pas utiliser STOP ici, car MPI l'interprète comme un crash
            endif
        end if
    end subroutine ensure_there_are_two_processes 

end program ! Fin du programme principal "transposition"
