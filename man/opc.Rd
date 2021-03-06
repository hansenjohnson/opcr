% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opc.R
\docType{data}
\name{opc}
\alias{opc}
\title{Processed OPC downcast}
\format{
A nested tibble with 148 rows and 8 variables:
\describe{
\item{scan}{the scan number, which increases with each data record}
\item{timer}{the timer count, which is sent every 0.5 seconds since the unit was powered on and resets after 4095 values}
\item{atten}{light attenuance}
\item{depth}{instrument depth in meters}
\item{flag}{quality control flag, with zero meaning good. See \code{opc_flag()} for the other definitions}
\item{secs}{the number of seconds elapsed since the instrument was powered on}
\item{time}{the datetime since the instrument was powered on}
\item{volume_filtered}{the volume of water that has passed through the OPC since the previous record, in cubic meters}
\item{esd}{a nested list of the particle sizes, in Equivalent Spherical Diameter (ESD; mm) detected during this data record}
}
}
\source{
hansen.johnson@dal.ca
}
\usage{
opc
}
\description{
A dataset containing a processed, trimmed OPC downcast (cast id \verb{2018_33}) that
was collected in the southern Gulf of St Lawrence in August, 2018.
}
\keyword{datasets}
