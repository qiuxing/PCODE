\name{randomization}
\alias{randomization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generating block design randomizations used for clinical trials.}
\description{
  Generate block design randomizations used for clinical trials.
}
\usage{
randomization(ntotal, ngroups = 2, nblocks = 1, verbose = FALSE,
group.labels = LETTERS[1:ngroups], block.group.sizes = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntotal}{ The total number of patients to be recruited. }
  \item{ngroups}{ Number of control/treatment groups used in this
    study.  Default to 2 (control v. treatment).}
  \item{nblocks}{Number of blocks to be used.  Default to 1 (complete
    randomization design).}
  \item{verbose}{If \code{verbose=TRUE}, print the block size table.}
  \item{group.labels}{Optional, labels (names) of the treatment groups.
  Default to letters "A", "B", \dots, etc.}
  \item{block.group.sizes}{Optional. A table containing block sizes.  Rows:
    blocks; columns: treatment groups.  Elements: number of patients in
    each block/group.  Sum of this table must equal \code{ntotal}. If
    not given, the default value is a matrix with every entry equaling
    \code{ntotal / (ngroups * nblocks)}.  Or an error message if
    \code{ntotal} is not a multiple of \code{ngroups * nblocks}.}
}
\details{
  This function generates a randomization list of treatments, possibly
  with a blocked design.
}
\value{
  \item{X}{A list of treatment group IDs that can be used in a
    randomized clinical trial.}
  \item{block.group.sizes}{A table containing block sizes.  Rows:
    blocks; columns: treatment groups.  Elements: number of patients in
    each block/group.}
}
\references{Some text book on clinical trials.}
\author{Xing Qiu}
\seealso{blockrand {blockrand}}
\examples{
## Default behavior, control/treatment, completely randomization
## design. Total number of patients: 20.  In each group: 10.
randomization(20)

## Three treatment groups with labels "Placebo", "Low", and "High".  5
## blocks of size 12 (4 for each treatment).  Total number of patients:
## 60.  Verbose output.
randomization(60, 3, 5, TRUE, c("Placebo", "Low", "High"))

## Customized blocks.  Needed because ntotal=78 is not a very good number
## for block randomization.  78 will be roughly divided into 3 blocks of
## size 20, each containing 10 controls, and 10 treatments, plus 1 block
## with 18 patients, 9 controls and 9 treatments.

randomization(78, ngroups=2, nblocks=4, group.labels=c("Control",
"Treatment"), block.group.sizes=matrix(c(rep(10, 6), rep(9,2)), nrow=4, byrow=TRUE))
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{randomization}
\keyword{block}% __ONLY ONE__ keyword per line
