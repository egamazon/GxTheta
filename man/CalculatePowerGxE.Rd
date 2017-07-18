\name{CalculatePowerGxE}
\alias{CalculatePowerGxE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CalculatePowerGxE
}
\description{
This function simulates an environmental variable E that is correlated with global ancestry. The correlation is a function of sigma. It performs the AIT test and computes power at a given alpha (pval.cutoff). 
}
\usage{
CalculatePowerGxE(save.path, pval.cutoff, E.effect.size, genetic.effect.size, global.anc.effect, interaction.effect.size, n.idvs, sigma.max, sigma.step, n.iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{save.path}{
path to save results to
}
  \item{pval.cutoff}{
alpha for computing power
}
  \item{E.effect.size}{
effect size of environmental variable
}
  \item{genetic.effect.size}{
genetic effect size
}
  \item{global.anc.effect}{
effect size of global ancestry
}
  \item{interaction.effect.size}{
effect size of interaction
}
  \item{n.idvs}{
number of individuals
}
  \item{sigma.max}{
maximum standard deviation used for the normal distribution when creating environmental variables correlated with global ancestry
}
  \item{sigma.step}{
how much to increase sigma by starting from 0 to \code{sigma.max} when creating environmental variables correlated with global ancestry
}
  \item{n.iter}{
number of iterations for the simulation
}
}
\details{
%% ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
