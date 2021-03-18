#' @name bayes.pick.winner
#' @title Bayesian pick the winner design
#' @description This function will run Shiny application for the Bayesian pick the winner design in randomized phase II trials.
#'              two-arm stage one and two and a two stage mulitple arm design.
#' @export

bayes.pick.winner <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "BayesianPickWinner")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `BayesianPickWinner`.", call. = FALSE)
  }

  shiny::runApp(paste(appDir,'/bayes_pick_winner_shiny_application.R',sep = ''), launch.browser = T,
                host = getOption( "127.0.0.1"))
}


