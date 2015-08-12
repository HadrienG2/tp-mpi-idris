!**********************************************************************
!**********************************************************************
! Comme dernier exercice de MPI, on va essayer de résoudre une équation
! de Poisson 2D en parallèle par décomposition de domaine.
!
! L'équation spécifique résolue ici est détaillée dans le module
! EquationSpecifique. La méthode utilisé est un solveur Jacobi sur une
! discrétisation de type différences finies.
!
! Et comme je suis d'humeur coquine, on va même essayer de le faire en
! programmation hybride OpenMP/MPI, ce qui permettra de tester la
! performance comparée des threads et des messages sur ma machine.
! Cette dernière est dual-core, donc l'optimum est 2 process/threads
!
! Performance en optimisation -O3, sur une grille 2000x2000, à 1e-3
! (avec sortie texte désactivée, bien entendu)
!   - 37.771s en séquentiel
!   - 21.326s avec 2 threads OpenMP (1.77x)
!   - 19.741s avec 2 process MPI (1.91x)
!   - 22.870s avec 2 threads et 2 process (1.65x)
!   - 50.424s avec 4 process (0.75x)
!   - 21.236s avec 4 threads (1.78x)
!
! Performance en optimisation -O3, sur une grille 4000x4000, à 5e-3
!   - 30.505s en séquentiel
!   - 15.921s avec 2 threads OpenMP (1.92x)
!   - 15.965s avec 2 process MPI (1.91x)
!   - 17.244s avec 2 threads et 2 process (1.76x)
!   - 24.858s avec 4 process (1.23x)
!   - 16.498s avec 4 threads (1.85x)
!
! On constate qu'OpenMP n'aime pas du tout les grilles trop petites (on
! sent l'overhead de lancement/synchro des threads !), et que MPI
! n'aime pas du tout qu'on lance trop de processus (problème d'affinité
! processeur ou communication excessive ?)
!**********************************************************************
!**********************************************************************

program poisson
    
    use equation_specifique
    use gestion_mpi
    use mpi
    use resolution_2D
    use topologie_cartesienne
    implicit none
    
    character(len=*), parameter :: nom_fichier = "resultat.dat"
    
    integer, dimension(nb_dimensions), parameter :: taille_active_tableau_u = (/ 8000, 8000 /)
    integer, dimension(nb_dimensions), parameter :: taille_tableau_u = taille_active_tableau_u + nb_valeurs_limite_u
    
    real(kind=u_kind), parameter :: erreur_globale_desiree = 1.e-1
    
    type(donnees_resolution_poisson) :: donnees_resolution
    
    real(kind=u_kind) :: erreur_locale, erreur_globale
    
    ! Initialization de la gestion MPI
    call initialisation_gestion_mpi
    
    ! Initialiser la topologie et se localiser dedans
    call initialisation_topologie
    
    ! Préparer les données nécessaires à la résolution de l'équation de Poisson
    call construire_donnees_resolution(taille_tableau_u, donnees_resolution)
    
    ! Itérer jusqu'à ce que l'erreur globale soit suffisamment faible
    do
        erreur_locale = iteration_resolution_jacobi(donnees_resolution)
        
        call echanger_valeurs_limites(donnees_resolution)
        
        erreur_globale = calcul_erreur_globale(erreur_locale)
        
        if (world_rank == root_process_world_rank) then
            print *, "Erreur globale : ", erreur_globale
        end if
        
        if (erreur_globale <= erreur_globale_desiree) exit
    end do
    
    ! Stocker le résultat du calcul dans un fichier avec MPI-IO (attention aux valeurs limites !)
    if (world_rank == root_process_world_rank) then
        print *, "Le calcul a convergé, écriture des résultats..."
    end if
    call stocker_resultats(donnees_resolution, nom_fichier)
    
    ! Libérer les données allouées dynamiquement
    call detruire_donnees_resolution(donnees_resolution)
    
    ! Finalization de la gestion MPI
    call finalisation_gestion_mpi
    
end program
