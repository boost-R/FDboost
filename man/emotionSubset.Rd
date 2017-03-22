\name{emotionSubset}
\alias{emotionSubset}
\docType{data}
\title{EEG and EMG recordings in a computerised gambling study}
\description{
  
   To analyse the functional relationship between electroencephalography (EEG) and facial electromyography (EMG), Gentsch et al. (2014) simultaneously recorded EEG and EMG signals from 24 participants while they were playing a computerised gambling task. The given subset contains aggregated observations of 23 participants. Curves were averaged over each subject and each of the 8 study settings, resulting in 23 times 8 curves.
  
}
\usage{data("emotionSubset")}
\format{
  A list with the following 10 variables.
  \describe{
    \item{\code{power}}{factor variable with levels \emph{high} and \emph{low}}
    \item{\code{game\_outcome}}{factor variable with levels \emph{gain} and \emph{loss}}
    \item{\code{control}}{factor variable with levels \emph{high} and \emph{low}}
    \item{\code{time\_long}}{time points in long format}
    \item{\code{subject}}{factor variable with 23 levels}
    \item{\code{curveid}}{id for identification of each curve}
    \item{\code{EEG}}{AsIs matrix; EEG signal in wide format}
    \item{\code{EMG}}{AsIs matrix; EMG signal in wide format}
    \item{\code{s}}{time points for the functional covariate}
    \item{\code{t}}{time points for the functional response}
  }
}
\details{
  The aim is to explain potentials in the EMG signal by study settings as well as
  the EEG signal (see Ruegamer et al., 2016).
}
\source{
  
  Gentsch,  K.,  Grandjean,  D.  and  Scherer,  K.  R.  (2014)  Coherence  explored  between  emotioncomponents:  Evidence from event-related potentials and facial electromyography.BiologicalPsychology, 98, 70 – 81.
  
  Ruegamer D, Brockhaus S, Gentsch K, Scherer K, Greven S (2016). Detecting Synchronisation in
EEG- and EMG-Signals via Boosted Functional Historical Models.
  
}
\keyword{datasets}