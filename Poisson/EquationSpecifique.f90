!**********************************************************************
!**********************************************************************
! Ce module encapsule toutes les problématiques "métiers" liées à
! l'équation de Poisson spécifique résolue, à savoir :
!       d^2(u)/d(x^2) + d^2(u)/d(y^2) = f(x, y) dans ]0,1[^2
!       u(x,y) = 0 aux limites de l'intervalle
!       f(x,y) = 2 * (x^2-x + y^2-y)
!
! La solution analytique est connue, c'est x^2-x + y^2-y. Cela donne
! un moyen de vérifier les résultats obtenus.
!**********************************************************************
!**********************************************************************

module equation_specifique

    implicit none

    ! On traite ici une équation à deux dimension. Si cette donnée est modifiée, tout le code ci-dessous
    ! doit être réécrit car il est très spécifique au cas deux dimensions.
    integer, parameter :: nb_dimensions = 2
    
    ! Ce paramètre définit le type de réel recommandé pour les valeurs de u
    integer, parameter :: u_significant_digits = 5, u_max_power_of_ten = 1
    integer, parameter :: u_kind = selected_real_kind(p=u_significant_digits, r=u_max_power_of_ten)
    
    ! Ce paramètre détermine les bornes de l'intervalle de travail dans l'espace réel
    real(kind=u_kind), dimension(nb_dimensions), parameter :: coord_debut = (/ 0.0, 0.0 /), &
                                                              coord_fin   = (/ 1.0, 1.0 /)
    
    ! On a besoin d'une valeur limite de u à chaque borne de l'intervalle
    integer, parameter :: nb_valeurs_limites = 1
    integer, dimension(nb_dimensions), parameter :: nb_valeurs_limite_u = (/ 2*nb_valeurs_limites, 2*nb_valeurs_limites /)
    real(kind=u_kind), parameter :: valeur_u_limite = 0.d0
    
    ! La fonction u considérée ici n'est pas périodique
    logical, dimension(nb_dimensions), parameter :: periodicite_u = (/ .FALSE., .FALSE. /)
    
contains

    ! Le second membre de l'équation est donné en fonction de coordonnées d'espace réel
    pure function second_membre(x, y) result(valeur)
        real(kind=u_kind), intent(in) :: x, y
        real(kind=u_kind) :: valeur
        
        valeur = 2.0*(x**2-x + y**2-y)
    end function second_membre

end module equation_specifique
