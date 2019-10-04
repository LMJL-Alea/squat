#' Dynamic Time Warping for Quaternion Time Series
#'
#' This function evaluates the Dynamic Time Warping (DTW) distance between two
#' quaternion time series (QTS).
#'
#' @param s1 A QTS stored as an object of class \code{\link[onion]{quaternion}}.
#' @param s2 A QTS stored as an object of class \code{\link[onion]{quaternion}}.
#' @param distance_only A boolean specifyung whether to only compute distance
#'   (no backtrack, faster). Default is \code{FALSE}.
#'
#' @return An object of class \code{\link[dtw]{dtw}} storing the dynamic time
#'   warping results.
#' @export
#'
#' @examples
#' # Pour générer 2 séries de quaternions :
#' s1.brut <- onion::rquat(15, rand = "norm")
#' s1 <- s1.brut / Mod(s1.brut)
#' s2.brut <- onion::rquat(20, rand = "norm")
#' s2 <- s2.brut / Mod(s2.brut)
#'
#' DTW(s1, s2)
DTW <- function(s1, s2, distance_only = FALSE) {
  M <- GetCostMatrix(as.matrix(s1), as.matrix(s2))
  dtw::dtw(M, distance.only = distance_only)
}

oldDTW <- function(t1, t2, s1, s2, distance = geodesic, size = 1) {
  # Contrainte 1 : Boundary
  # Calcul distance entre les 2 premiers éléments des 2 séries
  dist_map <- distance(s1[1] , s2[1])

  # Contrainte 2 : Step Size
  # Choix de la taille de la zone de mapping
  z <- as.matrix(expand.grid(0:size, 0:size)[-1, 2:1])
  colnames(z) <- NULL

  # Initialisation du premier élément de la fonction de warping
  warp <- rbind(c(1, 1))

  # Initialisation de l'algo DTW
  k <- 0

  # Algo de DTW

  while (warp[k + 1, 1] < length(t1) | warp[k + 1, 2] < length(t2)) {

    # Tant que les derniers points des 2 séries n'ont pas été mappés
    k <- k + 1
    pos <- z
    pos[, 1] <- pos[, 1] + warp[k, 1]
    pos[, 2] <- pos[, 2] + warp[k, 2]

    # Vérifier que la zone de recherche ne dépasse pas la matrice de distance
    cond1 <- pos[, 1] <= length(t1)
    cond2 <- pos[, 2] <= length(t2)
    z <- z[cond1 & cond2, , drop = FALSE]

    # Calcul des distances dans la zone de mapping
    dz <- rep(NA, nrow(z))
    for (i in 1:nrow(z)) {
      ind <- z[i, ] + warp[k, ]
      dz[i] <- distance(s1[ind[1]], s2[ind[2]])
    }

    # Sélection de la distance minimale dans la zone de mapping
    dist_map <- c(dist_map , min(dz))

    # Sélection du mapping correpondant à la distance minimale
    warp <- rbind(warp, z[which.min(dz), ] + warp[k, ])
  }

  list(
    distance = sum(dist_map),
    warping = tibble::tibble(
      label1 = warp[, 1],
      label2 = warp[, 2],
      distance = dist_map
    )
  )
}
