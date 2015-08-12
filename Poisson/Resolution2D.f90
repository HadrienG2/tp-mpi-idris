!**********************************************************************
!**********************************************************************
! Ce module contient un solveur adapté à des équations de Poisson
! bidimensionnelles avec conditions aux limites
!
! Les conflits entre conventions MPI et Fortran pour les coordonnées
! de tableaux et de topologies peuvent rendre le code assez dur à
! suivre, voir les schémas joints pour mieux visualiser tout ça.
!**********************************************************************
!**********************************************************************

module resolution_2D

    use equation_specifique
    use mpi
    use topologie_cartesienne
    implicit none
    
    ! Type de données MPI associé aux valeurs de u (type u_kind)
    integer :: u_mpi_type

    ! Etiquette utilisée pour les échanges de messages entre voisins, auquel on ajoute un numéro entre 0 et 3 pour
    ! différencier les quatre types d'échanges de données possibles (bordure x-inf, x-sup, y-inf, y-sup)
    integer, parameter :: tag_base_echange_voisins = 42
    
    ! Opération de réduction utilisée pour le calcul de l'erreur globale
    integer :: op_reduction_erreur

    ! Les types de données MPI suivants servent à l'échange de données de u "aux limites" entre voisins.
    ! Les données envoyées n'arrivent pas au même endroit que les données reçues, donc deux types sont nécessaires
    ! pour chaque bordure du domaine de u que possède un processus.
    ! 
    ! ATTENTION : Ici, on suit la convention de MPI et pas celle de Fortran : bordures horizontales, puis verticales
    type types_MPI_bordure_partagee
        integer :: envoi, reception
    end type
    type types_bordures_sur_dimension
        type(types_MPI_bordure_partagee) :: inferieur, superieur
    end type
        
    ! Le type suivant contient l'ensemble des données dont le solveur a besoin pour travailler
    type donnees_resolution_poisson
        integer, dimension(nb_dimensions) :: taille_tableau_u_global
        real(kind=u_kind), dimension(:,:), allocatable :: tableau_u_local, second_membre
        real(kind=u_kind), dimension(nb_dimensions) :: pas_simulation
        real(kind=u_kind) :: c0, c1, c2
        type(types_bordures_sur_dimension), dimension(nb_dimensions) :: types_bordures_domaine
    end type

contains

    ! Cette routine prépare les données locales à chaque processus dont le solveur d'équation de Poisson
    ! a besoin pour travailler : pas en coordonnées et coefficients de Jacobi associés, tableau de valeurs de u...
    subroutine construire_donnees_resolution(taille_tableau_u, &
                                             donnees_resolution)
        integer, dimension(nb_dimensions), intent(in) :: taille_tableau_u
        type(donnees_resolution_poisson), intent(out) :: donnees_resolution
        
        real(kind=u_kind), dimension(nb_dimensions) :: carre_pas_simulation
        
        ! Noter la taille du tableau global de u
        donnees_resolution%taille_tableau_u_global(:) = taille_tableau_u(:)
        
        ! Déterminer le pas en coordonnées et les coefficients de la méthode de Jacobi
        donnees_resolution%pas_simulation(:) = (coord_fin(:) - coord_debut(:))/(taille_tableau_u(:) - 1)
        carre_pas_simulation(:) = donnees_resolution%pas_simulation(:)**2
        donnees_resolution%c0 = (carre_pas_simulation(1) * carre_pas_simulation(2)) /           &
                                    (2 * (carre_pas_simulation(1) + carre_pas_simulation(2)))
        donnees_resolution%c1 = carre_pas_simulation(1)**(-1)
        donnees_resolution%c2 = carre_pas_simulation(2)**(-1)
        
        ! Préparer le tableau de valeurs de la fonction u local à chaque processus
        call construire_tableau_u_local(donnees_resolution)
                                        
        ! Calculer le second membre de l'équation
        call construire_second_membre(donnees_resolution)
                                        
        ! Définir des types de données MPI pour l'échange de valeurs limites avec les voisins du processus
        call creer_types_echange_voisin(donnees_resolution)
    end subroutine construire_donnees_resolution
    
    ! Finalisation des données locales créées pour les besoins de la résolution de l'équation
    subroutine detruire_donnees_resolution(donnees_resolution)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
        
        ! Libération des tableaux alloués dynamiquement
        deallocate(donnees_resolution%tableau_u_local, donnees_resolution%second_membre)
    end subroutine detruire_donnees_resolution
    
    ! Application d'une itération de la méthode de Jacobi à l'équation de Poisson, avec calcul simultané de l'erreur
    ! (ici choisie comme max des écarts relatifs absolus entre la valeur ancienne et nouvelle de u en chaque point)
    function iteration_resolution_jacobi(donnees_resolution) result(erreur)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
        real(kind=u_kind) :: erreur
        
        integer :: i, j, premiere_ligne, derniere_ligne, premiere_colonne, derniere_colonne
        
        real(kind=u_kind) :: ancien_u
        
        ! Récupération de la taille de la région active du tableau de valeurs de u
        premiere_ligne = lbound(donnees_resolution%tableau_u_local, 1) + nb_valeurs_limites
        derniere_ligne = ubound(donnees_resolution%tableau_u_local, 1) - nb_valeurs_limites
        premiere_colonne = lbound(donnees_resolution%tableau_u_local, 2) + nb_valeurs_limites
        derniere_colonne = ubound(donnees_resolution%tableau_u_local, 2) - nb_valeurs_limites
        
        ! Itération de la méthode de Jacobi
        ! (NOTE : si la définition de l'erreur est modifiée, penser à modifier aussi la méthode de réduction)
        erreur = 0.0
        !$omp parallel default(none) &
        !$omp shared(premiere_colonne, derniere_colonne, premiere_ligne, derniere_ligne, donnees_resolution, erreur) &
        !$omp private(ancien_u)
        !$omp do reduction(max:erreur)
        do j = premiere_colonne, derniere_colonne
            do i = premiere_ligne, derniere_ligne
                ancien_u = donnees_resolution%tableau_u_local(i, j)
                donnees_resolution%tableau_u_local(i, j) = donnees_resolution%c0*(  &
                    donnees_resolution%c1 * (                                       &
                        donnees_resolution%tableau_u_local(i+1, j) +                &
                        donnees_resolution%tableau_u_local(i-1, j)                  &
                    ) + donnees_resolution%c2 * (                                   &
                        donnees_resolution%tableau_u_local(i, j+1) +                &
                        donnees_resolution%tableau_u_local(i, j-1)                  &
                    ) - donnees_resolution%second_membre(i, j) )
                erreur = max(erreur, abs((donnees_resolution%tableau_u_local(i, j) - ancien_u) / &
                                         donnees_resolution%tableau_u_local(i, j)))
            end do
        end do
        !$omp end do
        !$omp end parallel
    end function iteration_resolution_jacobi
    
    ! Echanger les valeurs aux limites calculées pour l'équation de Poisson entre processus voisins
    subroutine echanger_valeurs_limites(donnees_resolution)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
    
        integer :: dim_actuelle, tag_base_dim_actuelle, premiere_ligne, premiere_colonne
        
        premiere_ligne = lbound(donnees_resolution%tableau_u_local, 1)
        premiere_colonne = lbound(donnees_resolution%tableau_u_local, 2)
        do dim_actuelle = 1, nb_dimensions
            tag_base_dim_actuelle = tag_base_echange_voisins + (2*(dim_actuelle-1))
            call MPI_SENDRECV(donnees_resolution%tableau_u_local(premiere_ligne, premiere_colonne),         &
                              1,                                                                            &
                              donnees_resolution%types_bordures_domaine(dim_actuelle)%inferieur%envoi,      &
                              rang_voisins(dim_actuelle)%inferieur,                                         &
                              tag_base_dim_actuelle,                                                        &
                              donnees_resolution%tableau_u_local(premiere_ligne, premiere_colonne),         &
                              1,                                                                            &
                              donnees_resolution%types_bordures_domaine(dim_actuelle)%inferieur%reception,  &
                              rang_voisins(dim_actuelle)%superieur,                                         &
                              tag_base_dim_actuelle,                                                        &
                              comm_topologie,                                                               &
                              MPI_STATUS_IGNORE,                                                            &
                              return_code)
            call MPI_SENDRECV(donnees_resolution%tableau_u_local(premiere_ligne, premiere_colonne),         &
                              1,                                                                            &
                              donnees_resolution%types_bordures_domaine(dim_actuelle)%superieur%envoi,      &
                              rang_voisins(dim_actuelle)%superieur,                                         &
                              tag_base_dim_actuelle+1,                                                      &
                              donnees_resolution%tableau_u_local(premiere_ligne, premiere_colonne),         &
                              1,                                                                            &
                              donnees_resolution%types_bordures_domaine(dim_actuelle)%superieur%reception,  &
                              rang_voisins(dim_actuelle)%inferieur,                                         &
                              tag_base_dim_actuelle+1,                                                      &
                              comm_topologie,                                                               &
                              MPI_STATUS_IGNORE,                                                            &
                              return_code)
        end do
    end subroutine echanger_valeurs_limites
    
    ! Calcul de l'erreur globale par réduction des erreurs locales
    function calcul_erreur_globale(erreur_locale) result(erreur_globale)
        real(kind=u_kind), intent(in) :: erreur_locale
    
        real(kind=u_kind) :: erreur_globale
        
        call MPI_ALLREDUCE(erreur_locale,       &
                           erreur_globale,      &
                           1,                   &
                           u_mpi_type,          &
                           op_reduction_erreur, &
                           comm_topologie)
    end function calcul_erreur_globale
    
    ! Sauvegarde des résultats locaux dans un fichier qui contiendra le tableau de valeurs complet de u, bords inclus
    subroutine stocker_resultats(donnees_resolution, nom_fichier)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
        character(len=*), intent(in) :: nom_fichier
        
        integer, dimension(nb_dimensions) :: position_topo_tableau, dimensions_topo_tableau
        integer, dimension(nb_dimensions) :: origine_sous_tableau, taille_sous_tableau, origine_absolue
        
        integer :: dim_actuelle
        integer :: file_descriptor, local_subarray_type, global_subarray_type
        
        integer(kind=MPI_OFFSET_KIND) :: zero_bytes = 0
        
        ! Changement de repère de la topologie MPI vers la topologie du tableau
        call conversion_topologie_tableau(position_topologie,       &
                                          dimensions_topologie,     &
                                          position_topo_tableau,    &
                                          dimensions_topo_tableau)
        
        ! Tous les processus sauvegardent leur région active (exprimée avec les conventions des sous-tableaux MPI)...
        origine_sous_tableau(:) = (/ nb_valeurs_limites,    &
                                     nb_valeurs_limites /)
        taille_sous_tableau(:) = (/ size(donnees_resolution%tableau_u_local, 1) - 2*nb_valeurs_limites,     &
                                    size(donnees_resolution%tableau_u_local, 2) - 2*nb_valeurs_limites /)
        
        ! ...et ceux qui sont au bord ou dans un coin du domaine sauvegardent aussi les valeurs aux limites de la solution
        do dim_actuelle = 1, nb_dimensions
            if (position_topo_tableau(dim_actuelle) == 0) then
                origine_sous_tableau(dim_actuelle) = origine_sous_tableau(dim_actuelle) - nb_valeurs_limites
                taille_sous_tableau(dim_actuelle) = taille_sous_tableau(dim_actuelle) + nb_valeurs_limites
            end if
            if (position_topo_tableau(dim_actuelle) == dimensions_topo_tableau(dim_actuelle)-1) then
                taille_sous_tableau(dim_actuelle) = taille_sous_tableau(dim_actuelle) + nb_valeurs_limites
            end if
        end do
        
        ! On déduit de l'origine du sous-tableau la position à laquelle on démarre dans le tableau global
        origine_absolue(:) = (lbound(donnees_resolution%tableau_u_local)-1+nb_valeurs_limites) &
                             + origine_sous_tableau(:)
        
        ! Définir un type MPI correspondant à la plage de lignes et colonnes choisie dans le tableau local
        call MPI_TYPE_CREATE_SUBARRAY(2,                                            &
                                      shape(donnees_resolution%tableau_u_local),    &
                                      taille_sous_tableau,                          &
                                      origine_sous_tableau,                         &
                                      MPI_ORDER_FORTRAN,                            &
                                      u_mpi_type,                                   &
                                      local_subarray_type,                          &
                                      return_code)
        call MPI_TYPE_COMMIT(local_subarray_type, return_code)
        
        ! Définir un type MPI localisant le tableau local dans le tableau global
        call MPI_TYPE_CREATE_SUBARRAY(2,                                            &
                                      donnees_resolution%taille_tableau_u_global,   &
                                      taille_sous_tableau,                          &
                                      origine_absolue,                              &
                                      MPI_ORDER_FORTRAN,                            &
                                      u_mpi_type,                                   &
                                      global_subarray_type,                         &
                                      return_code)
        call MPI_TYPE_COMMIT(global_subarray_type, return_code)
        
        ! Ouvrir et effacer le fichier dans lequel on va écrire les résultats
        call MPI_FILE_OPEN(comm_topologie,                      &
                           nom_fichier,                         &
                           MPI_MODE_WRONLY + MPI_MODE_CREATE,   &
                           MPI_INFO_NULL,                       &
                           file_descriptor,                     &
                           return_code)
        call MPI_FILE_SET_SIZE(file_descriptor, zero_bytes, return_code)
                           
        ! Définir une vue adaptée aux besoins de chaque processus
        call MPI_FILE_SET_VIEW(file_descriptor,         &
                               zero_bytes,              &
                               u_mpi_type,              &
                               global_subarray_type,    &
                               "native",                &
                               MPI_INFO_NULL,           &
                               return_code)
        
        ! Ecrire les données
        call MPI_FILE_WRITE_ALL(file_descriptor,                    &
                                donnees_resolution%tableau_u_local, &
                                1,                                  &
                                local_subarray_type,                &
                                MPI_STATUS_IGNORE,                  &
                                return_code)
        
        ! Fermer le fichier
        call MPI_FILE_CLOSE(file_descriptor, return_code)
    end subroutine stocker_resultats
    
    
    ! *** Les routines suivantes sont internes à ce module, et ne sont pas conçues pour être utilisées de l'extérieur ***
    
    ! Séparation du tableau "global" de valeurs de la fonction u en plusieurs sous-tableaux locaux à chaque processus,
    ! comprenant chacun les cases-limites et le contenu requis par le problème, et dont les indices sont choisis
    ! pour illustrer la position d'un processus dans la topologie.
    !
    ! Le processus du coin inférieur droit de la topologie reçoit ainsi le coin inférieur-droit du tableau global, etc.
    subroutine construire_tableau_u_local(donnees_resolution)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
        
        integer, dimension(nb_dimensions) :: position_topo_tableau, dimensions_topo_tableau
        integer, dimension(nb_dimensions) :: taille_active_tab_u, taille_active_tab_u_local, taille_tab_u_local
        integer, dimension(nb_dimensions) :: decalage_u_local
        
        integer :: i, premiere_ligne, derniere_ligne, premiere_colonne, derniere_colonne
        
        ! Changement de repère de la topologie MPI vers une topologie plus adaptée à la division d'un tableau Fortran
        call conversion_topologie_tableau(position_topologie,       &
                                          dimensions_topologie,     &
                                          position_topo_tableau,    &
                                          dimensions_topo_tableau)
        
        ! Détermination de la taille et la position dans u de chaque sous-tableau de u (cases limites non incluses)
        taille_active_tab_u(:) = donnees_resolution%taille_tableau_u_global(:) - nb_valeurs_limite_u(:)
        do i = 1, nb_dimensions
            ! Cas exact si le tableau global est parfaitement divisible par la dimension de la topologie
            taille_active_tab_u_local(i) = taille_active_tab_u(i) / dimensions_topo_tableau(i)
            decalage_u_local(i) = position_topo_tableau(i) * taille_active_tab_u_local(i)
            
            ! Dans le cas où la divisibilité n'est pas parfaite, les premiers processus reçoivent plus de points
            if (position_topo_tableau(i) < modulo(taille_active_tab_u(i), dimensions_topo_tableau(i))) then
                taille_active_tab_u_local(i) = taille_active_tab_u_local(i) + 1
                decalage_u_local(i) = decalage_u_local(i) + position_topo_tableau(i)
            else
                decalage_u_local(i) = decalage_u_local(i) + modulo(taille_active_tab_u(i), dimensions_topo_tableau(i))
            end if
        end do
        
        
        ! Déterminer la taille du tableau de u local de chaque processus, cases limites incluses
        taille_tab_u_local(:) = taille_active_tab_u_local(:) + nb_valeurs_limite_u(:)
        
        ! Allouer ce tableau, et supposer initialement u égale à sa valeur limite
        premiere_ligne = (1 - nb_valeurs_limites) + decalage_u_local(1)
        derniere_ligne = premiere_ligne + (taille_tab_u_local(1) - 1)
        premiere_colonne = (1 - nb_valeurs_limites) + decalage_u_local(2)
        derniere_colonne = premiere_colonne + (taille_tab_u_local(2) - 1)
        allocate(donnees_resolution%tableau_u_local(premiere_ligne:derniere_ligne,      &
                                                    premiere_colonne:derniere_colonne))
        donnees_resolution%tableau_u_local(:,:) = valeur_u_limite
    end subroutine construire_tableau_u_local
    
    ! Cette routine convertit des caractéristiques de topologie MPI 2D (premier indice représentant les lignes,
    ! de gauche à droite et de bas en haut) vers une topologie équivalente représentant mieux le découpage d'un
    ! tableau Fortran (premier indice représentant les colonnes, de haut en bas et de gauche à droite)
    pure subroutine conversion_topologie_tableau(position_topologie,        &
                                                 dimensions_topologie,      &
                                                 position_topo_tableau,     &
                                                 dimensions_topo_tableau)
        integer, dimension(nb_dimensions), intent(in) :: position_topologie, dimensions_topologie
        integer, dimension(nb_dimensions), intent(out) :: position_topo_tableau, dimensions_topo_tableau
        
        position_topo_tableau(:) = (/ dimensions_topologie(2) - 1 - position_topologie(2),  &
                                      position_topologie(1) /)
        dimensions_topo_tableau(:) = (/ dimensions_topologie(2),    &
                                        dimensions_topologie(1) /)
    end subroutine conversion_topologie_tableau
    
    ! Cette routine construit le second membre d'une équation de Poisson dont on a déjà le tableau de valeurs local
    subroutine construire_second_membre(donnees_resolution)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
        
        integer :: i, j, premiere_ligne, derniere_ligne, premiere_colonne, derniere_colonne
        
        ! Le second membre est calculé pour les indices "actifs" du tableau de valeurs de u, hors valeurs limites
        premiere_ligne = lbound(donnees_resolution%tableau_u_local, 1) + nb_valeurs_limites
        derniere_ligne = ubound(donnees_resolution%tableau_u_local, 1) - nb_valeurs_limites
        premiere_colonne = lbound(donnees_resolution%tableau_u_local, 2) + nb_valeurs_limites
        derniere_colonne = ubound(donnees_resolution%tableau_u_local, 2) - nb_valeurs_limites
        
        ! Allocation du tableau pour le stockage du second membre
        allocate(donnees_resolution%second_membre(premiere_ligne:derniere_ligne,        &
                                                  premiere_colonne:derniere_colonne))
                
        ! Remplissage de ce dernier avec les valeurs du second membre correspondant aux plages d'abscisse considérées
        forall (i = premiere_ligne:derniere_ligne, j = premiere_colonne:derniere_colonne)
            donnees_resolution%second_membre(i, j) = second_membre(                             &
                coord_debut(1) + donnees_resolution%pas_simulation(1)*(i-premiere_ligne),       &
                coord_debut(2) + donnees_resolution%pas_simulation(2)*(j-premiere_colonne) )
        end forall
    end subroutine construire_second_membre
    
    ! Définition des types dérivés MPI permettant d'échanger les valeurs limites de U entre voisins, sans échanger
    ! les valeurs de coin telles que (0, 0) qui sont inutiles au calcul.
    ! Les positions des bordures, en envoi et en réception, sont prises du point de vue du processus appelant
    subroutine creer_types_echange_voisin(donnees_resolution)
        type(donnees_resolution_poisson), intent(inout) :: donnees_resolution
        
        integer, dimension(nb_dimensions) :: taille_tableau_u_local
        integer, dimension(nb_dimensions) :: taille_bordure_axe_horizontal, taille_bordure_axe_vertical
        integer, dimension(nb_dimensions) :: position_bordure_envoi, position_bordure_reception
        
        integer :: return_code
        
        ! Récupération du type de données MPI associé à u_kind
        call MPI_TYPE_CREATE_F90_REAL(u_significant_digits,     &
                                      u_max_power_of_ten,       &
                                      u_mpi_type,               &
                                      return_code)
                                      
        ! Création d'une opération de réduction adaptée à ce type
        call MPI_OP_CREATE(reduction_erreur,    &
                           .true.,              &
                           op_reduction_erreur, &
                           return_code)
        
        ! Calculer la taille des bordures du tableau local sur les axes horizontaux et verticaux : leur longueur
        ! est égale à la taille active du tableau (nb de valeurs vraiment calculées, hors cases limites), et leur
        ! largeur au nombre de cases limites.
        taille_tableau_u_local(:) = shape(donnees_resolution%tableau_u_local)
        taille_bordure_axe_horizontal(:) = (/ taille_tableau_u_local(1) - nb_valeurs_limite_u(1),   &
                                              nb_valeurs_limites /)
        taille_bordure_axe_vertical(:) = (/ nb_valeurs_limites,                                     &
                                            taille_tableau_u_local(2) - nb_valeurs_limite_u(2) /)
        
        ! Création des types pour l'échange de données sur la bordure horizontale inférieure
        position_bordure_envoi(:) = (/ nb_valeurs_limites,      &
                                       nb_valeurs_limites /)
        position_bordure_reception(:) = (/ nb_valeurs_limites,                                  &
                                           taille_tableau_u_local(2) - nb_valeurs_limites /)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                            &
                                      taille_tableau_u_local,                                       &
                                      taille_bordure_axe_horizontal,                                &
                                      position_bordure_envoi,                                       &
                                      MPI_ORDER_FORTRAN,                                            &
                                      u_mpi_type,                                                   &
                                      donnees_resolution%types_bordures_domaine(1)%inferieur%envoi, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(1)%inferieur%envoi, return_code)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                                &
                                      taille_tableau_u_local,                                           &
                                      taille_bordure_axe_horizontal,                                    &
                                      position_bordure_reception,                                       &
                                      MPI_ORDER_FORTRAN,                                                &
                                      u_mpi_type,                                                       &
                                      donnees_resolution%types_bordures_domaine(1)%inferieur%reception, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(1)%inferieur%reception, return_code)
        
        ! Création des types pour l'échange de données sur la bordure horizontale supérieure
        position_bordure_envoi(:) = (/ nb_valeurs_limites,                                  &
                                       taille_tableau_u_local(2) - 2*nb_valeurs_limites /)
        position_bordure_reception(:) = (/ nb_valeurs_limites,  &
                                           0 /)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                            &
                                      taille_tableau_u_local,                                       &
                                      taille_bordure_axe_horizontal,                                &
                                      position_bordure_envoi,                                       &
                                      MPI_ORDER_FORTRAN,                                            &
                                      u_mpi_type,                                                   &
                                      donnees_resolution%types_bordures_domaine(1)%superieur%envoi, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(1)%superieur%envoi, return_code)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                                &
                                      taille_tableau_u_local,                                           &
                                      taille_bordure_axe_horizontal,                                    &
                                      position_bordure_reception,                                       &
                                      MPI_ORDER_FORTRAN,                                                &
                                      u_mpi_type,                                                       &
                                      donnees_resolution%types_bordures_domaine(1)%superieur%reception, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(1)%superieur%reception, return_code)
                                      
        ! Création des types pour l'échange de données sur la bordure verticale inférieure
        position_bordure_envoi(:) = (/ taille_tableau_u_local(1) - 2*nb_valeurs_limites,    &
                                       nb_valeurs_limites /)
        position_bordure_reception(:) = (/ 0,                       &
                                           nb_valeurs_limites /)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                            &
                                      taille_tableau_u_local,                                       &
                                      taille_bordure_axe_vertical,                                  &
                                      position_bordure_envoi,                                       &
                                      MPI_ORDER_FORTRAN,                                            &
                                      u_mpi_type,                                                   &
                                      donnees_resolution%types_bordures_domaine(2)%inferieur%envoi, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(2)%inferieur%envoi, return_code)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                                &
                                      taille_tableau_u_local,                                           &
                                      taille_bordure_axe_vertical,                                      &
                                      position_bordure_reception,                                       &
                                      MPI_ORDER_FORTRAN,                                                &
                                      u_mpi_type,                                                       &
                                      donnees_resolution%types_bordures_domaine(2)%inferieur%reception, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(2)%inferieur%reception, return_code)
                                      
        ! Création des types pour l'échange de données sur la bordure verticale inférieure
        position_bordure_envoi(:) = (/ nb_valeurs_limites,      &
                                       nb_valeurs_limites /)
        position_bordure_reception(:) = (/ taille_tableau_u_local(1) - nb_valeurs_limites, &
                                           nb_valeurs_limites /)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                            &
                                      taille_tableau_u_local,                                       &
                                      taille_bordure_axe_vertical,                                  &
                                      position_bordure_envoi,                                       &
                                      MPI_ORDER_FORTRAN,                                            &
                                      u_mpi_type,                                                   &
                                      donnees_resolution%types_bordures_domaine(2)%superieur%envoi, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(2)%superieur%envoi, return_code)
        call MPI_TYPE_CREATE_SUBARRAY(2,                                                                &
                                      taille_tableau_u_local,                                           &
                                      taille_bordure_axe_vertical,                                      &
                                      position_bordure_reception,                                       &
                                      MPI_ORDER_FORTRAN,                                                &
                                      u_mpi_type,                                                       &
                                      donnees_resolution%types_bordures_domaine(2)%superieur%reception, &
                                      return_code)
        call MPI_TYPE_COMMIT(donnees_resolution%types_bordures_domaine(2)%superieur%reception, return_code)
    end subroutine creer_types_echange_voisin
    
    ! Opération de réduction utilisée pour le calcul de l'erreur globale (ici, un maximum)
    subroutine reduction_erreur(invec, inoutvec, length, datatype)
        real(kind=u_kind), dimension(length), intent(in) :: invec
        real(kind=u_kind), dimension(length), intent(inout) :: inoutvec
        integer, intent(in) :: length, datatype
        
        if (datatype /= u_mpi_type) call abort
        inoutvec(:) = max(invec(:), inoutvec(:))
    end subroutine reduction_erreur
end module resolution_2D
