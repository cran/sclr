# Datasets included with sclr

#' Simulated one-titre antibody data
#'
#' A simulated dataset containing 6000 independent observations on antibody
#' titres and the corresponding infection status. The data was simulated to
#' resemble real influenza infection and haemagglutinin titre data.
#'
#' @format A data frame with 6000 observations and 4 variables: 
#' \describe{
#'   \item{logHI}{haemagglutinin-inhibiting (HI) titre. 
#'   True simulated titre on a log scale.} 
#'   \item{logHIcens}{HI censored (observed) titre.
#'   The titre value on a log scale that would be observed in a real dataset
#'   with a typical HI assay.}
#'   \item{logHImid}{Midpoint of the interval (on a log scale)
#'   to which observed HI values are censored.}
#'   \item{status}{influenza infection status. 1 -
#'   infected. 0 - not infected}
#'  }
#'
#' @section Model:
#'
#'   The model behind the simulation was
#'
#'   \deqn{\lambda * (1 - f(\beta_0 + \beta_1 * HI))}
#'
#'   Where
#'
#'   \itemize{ \item \eqn{f} - Inverse logit function \item \eqn{\lambda} = 0.25
#'   \item \eqn{\beta_0} = 7.5 \item \eqn{\beta_1} = 2 }
#'   
"sclr_one_titre_data"

#' Simulated two-titre antibody data
#'
#' A simulated dataset containing 6000 independent observations on antibody
#' titres and the corresponding infection status. The data was simulated to
#' resemble real influenza infection and haemagglutinin + neuraminidase titre
#' data.
#'
#' @format A data frame with 6000 observations and 6 variables: 
#' \describe{
#'   \item{logHI}{haemagglutinin-inhibiting (HI) titre. 
#'   True simulated titre on a log scale.} 
#'   \item{logHIcens}{HI censored (observed) titre.
#'   The titre value on a log scale that would be observed in a real dataset
#'   with a typical HI assay.}
#'   \item{logHImid}{Midpoint of the interval (on a log scale)
#'   to which observed HI values are censored.}
#'   \item{logNI}{neuraminidase-inhibiting titre. True simulated titre on a
#'   log scale.} 
#'   \item{logNIcens}{haemagglutinin-inhibiting censored (observed)
#'   titre. The titre value that would be observed in a real dataset with a
#'   typical NI assay.} 
#'   \item{status}{influenza infection status. 1 - infected.
#'   0 - not infected} 
#'  }
#'
#' @section Model:
#'
#'   The model behind the simulation was
#'
#'   \deqn{\lambda * (1 - f(\beta_0 + \beta_1 * HI + \beta_2 * NI))}
#'
#'   Where
#'
#'   \itemize{ \item \eqn{f} - Inverse logit function \item \eqn{\lambda} = 0.25
#'   \item \eqn{\beta_0} = -7.5 \item \eqn{\beta_1} = 2 \item \eqn{\beta_2} = 2
#'   }
#'   
"sclr_two_titre_data"