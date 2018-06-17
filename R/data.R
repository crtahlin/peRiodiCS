#' @title Viral etiology, seasonality and severity of hospitalized patients with severe acute respiratory infections in the Eastern Mediterranean Region, 2007–2014
#' 
#' @description 
#' Data about infections with different viruses across several years. 
#' 
#' For more information see Source and References section. 
#' 
#' @format A data frame with variables:
#' \describe{
#' \item{UniqueID}{record identification number}
#' \item{Enrolled}{Did the patient consent and enroll in the study?:	1=Yes, 0=No}
#' \item{Country}{Country of enrollment:	Egypt, Jordan, Oman, Qatar, Yemen}
#' \item{EpiYear}{Year of enrollment:	Integers (2007-2014)}
#' \item{EpiMonth}{Month of enrollment:	Integers (1-12)}
#' \item{EpiWeek}{Week of enrollment:	Integers (1-53)}
#' \item{Interval}{Number of days between onset of symptoms and hospitalization:	Integer}
#' \item{Stay}{Number of days between hospitalization and outcome:	Integer}
#' \item{Sex}{Sex: 1=Female, 0=Male}
#' \item{AgeGrp}{Age group:	1=<1 year, 2=1-4 years, 3=5-49 years, 4=50+ years}
#' \item{ChronicDis}{Does the patient have any pre-existing chronic disease?:	1=Yes, 0=No}
#' \item{OxTherapy}{Did the patient receive oxygen therapy during hospitalization?:	1=Yes, 0=No}
#' \item{Ventilated}{Was the patient ventilated during hospitalization?:	1=Yes, 0=No}
#' \item{ICU}{Was the patient admitted to an intensive care unit during hospitalization?:	1=Yes, 0=No}
#' \item{Outcome}{What was the patient"'"s final hospitalization outcome?:	1=Discharge, 2=Transfer, 3=Death}
#' \item{RSV}{Results for respiratory syncytial virus: 1=Positive, 0=Negative}
#' \item{AdV}{Results for adenovirus:	1=Positive, 0=Negative}
#' \item{hMPV}{Results for human metapneumovirus:	1=Positive, 0=Negative}
#' \item{hPIV1}{Results for human parainfluenzavirus type 1: 1=Positive, 0=Negative}
#' \item{hPIV2}{Results for human parainfluenzavirus type 2: 1=Positive, 0=Negative}
#' \item{hPIV3}{Results for human parainfluenzavirus type 3: 1=Positive, 0=Negative}
#' \item{Flu}{Results for influenza:	1=Positive, 0=Negative}
#' }
#' 
#' @source \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0180954}
#' 
#' @references 
#' Horton, Katherine C. AND Dueger, Erica L. AND Kandeel, Amr AND Abdallat, Mohamed AND El-Kholy, Amani AND Al-Awaidy, Salah AND Kohlani, Abdul Hakim AND Amer, Hanaa AND El-Khal, Abel Latif AND Said, Mayar AND House, Brent AND Pimentel, Guillermo AND Talaat, Maha (2017). Viral etiology, seasonality and severity of hospitalized patients with severe acute respiratory infections in the Eastern Mediterranean Region, 2007–2014. PLOS ONE, 12, 1-17.
#' 
"viral_east_mediteranean"