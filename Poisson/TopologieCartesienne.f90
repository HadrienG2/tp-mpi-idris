!**********************************************************************
!**********************************************************************
! Ce module encapsule toutes les fonctions de gestion de MPI liées à
! l'utilisation de la topologie cartésienne pour la résolution
! d'équations de Poisson.
!**********************************************************************
!**********************************************************************

module topologie_cartesienne

    use equation_specifique
    use gestion_mpi
    use mpi
    implicit none
    
    ! Pour la dimension et la position des voisins, on utilise la convention MPI
    !      - La coordonnée 1 est x, la coordonnée 2 est y, etc
    !      - Les coordonnées x et y démarrent à 0
    !      - Les coordonnees x croissent vers la droite
    !      - Les coordonnees y croissent vers le haut
    
    integer, dimension(nb_dimensions) :: dimensions_topologie, position_topologie
    
    type rangs_voisinage
        integer :: inferieur, superieur
    end type
    type(rangs_voisinage), dimension(nb_dimensions) :: rang_voisins
    
    integer :: comm_topologie, rang_topologie
    
contains

    ! Initialisation complète du paquet de gestion de topologie cartésienne
    subroutine initialisation_topologie
        integer :: dimension
    
        ! Construire une topologie adaptée au problème
        call construire_topologie
        
        ! Se repérer dedans, ainsi que ses voisins
        call MPI_COMM_RANK(comm_topologie, rang_topologie, return_code)
        call MPI_CART_COORDS(comm_topologie, rang_topologie, nb_dimensions, position_topologie, return_code)
        do dimension = 1, nb_dimensions
            call MPI_CART_SHIFT(comm_topologie,                     &
                                dimension-1,                        &
                                1,                                  &
                                rang_voisins(dimension)%inferieur,  &
                                rang_voisins(dimension)%superieur,  &
                                return_code)
        end do
    end subroutine initialisation_topologie

    ! Construction d'une topologie cartésienne adaptée au problème
    subroutine construire_topologie   
        ! Laisser MPI calculer les dimensions de la topologie
        dimensions_topologie(:) = 0
        call MPI_DIMS_CREATE(world_size,            &
                             nb_dimensions,         &
                             dimensions_topologie,  &
                             return_code)
        
        ! Créer la topologie avec des propriétés données par le problème
        call MPI_CART_CREATE(MPI_COMM_WORLD,        &
                             nb_dimensions,         &
                             dimensions_topologie,  &
                             periodicite_u,         &
                             .true.,                &
                             comm_topologie,        &
                             return_code)
    end subroutine construire_topologie

end module topologie_cartesienne
