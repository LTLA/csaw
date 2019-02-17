\name{csawUsersGuide}
\alias{csawUsersGuide}
\alias{csaw}

\title{View csaw user's guide}
\description{Finds the location of the user's guide and opens it for viewing.}

\usage{
csawUsersGuide()
} 

\value{The user's guide is opened on the system's default browser.}

\details{
The csaw package is designed for de novo detection of differentially bound regions from ChIP-seq data. 
It provides methods for window-based counting, normalization, filtering and statistical analyses via edgeR.
The user guide for this package can be obtained by running this function.

The user's guide is obtained online from the \pkg{csawUsersGuide} package.
We do not build the guide as a vignette in \pkg{csaw} itself due to the time-consuming nature of the code when run on realistic case studies.
This does mean, though, that this function requires an internet connection to work properly.
}

\seealso{
\code{\link{browseURL}}
}

\author{Aaron Lun}

\examples{
if (interactive()) csawUsersGuide()
}

\keyword{documentation}
