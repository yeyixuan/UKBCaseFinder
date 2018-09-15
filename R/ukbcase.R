#' ukbcase: find cases in UK Biobank
#'
#' @title ukbcase
#'
#' @author yixuan ye
#'
#' \code{ukbcase} is a tool for identifying patients in UK biobank given the definition (ICD or OPCS) of
#' disease phenotypes. Here we combine 3 data sources: hospital in-patient episode records (hesin),
#' death records, and cancer records to firstly identify all the patients ID. And by comparing the date
#' of all the records, we then identify the earlist onset date for each patient.
#'
#' For more details about UK Biobank data sources related to health outcomes, see:
#' \url{http://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=100091}
#'
#' @rdname ukbcase
#'
#' @aliases ukbcase
#'
#' @export
#'
#' @param icd10 String array. International Classification of Diseases, 10th revision. For example, the icd10 for
#' myocardial infarction is I21-I24,I25, written as c('I211','I212',...).
#' @param icd9 String array. International Classification of Diseases, 9th revision. For example, the icd9 for
#' myocardial infarction is 410.
#' @param oper4 String array. Operative procedures which could also be used to imply specific diseases.
#' For example, the oper4 K40-K46, K49, K50, K75 imply the onset of MI, written as c('K40','K41',...).
#' @param histology String array. Histology code which could be used to imply specifc cancer. Meaningful only when
#' cancerdata is available, and 'cancer==Ture'.
#' @param hesin main hospital record dataset.
#' @param hesin_diag10 hospital record subdataset, containing diag10 except main diagnosis.
#' @param hesin_oper4 hospital record subdataset, containing oper4 except main opsc records.
#' @param deathdata death dataset, eid, icd10 and diagnosis time should be recorded, and named as
#' ('eid','case_of_death','date_of_death'), long formate dataset is required.
#' @param cancerdata cancer dataset, eid, cancertype and time should be recorded, cancertype could be
#' recorded either by icd10,icd9 or histology of the cancer, the corresponding names of column should
#' be ('eid','date,'icd10','icd9','histology'), long formate dataset is required.
#' @param icd10main logical; if TRUE, the main icd10 is considered.
#' @param icd10sec logical; if TRUE, the secondary icd10 is considered.
#' @param icd9main logical; if TRUE, the icd9 is considered.
#' @param oper4main logical; if TRUE, the main oper4 is considered.
#' @param oper4sec logical; if TRUE, the secondary oper4 is considered.
#' @param death logical; if TRUE, the death dataset is considered.
#' @param cancer logical; if TRUE, the cancer dataset is considered.
#'
#' @examples
#'
#' ## simulate hesin, hesin_diag10, hesin_oper4 dataset
#' hesin <- data.frame(eid=1:100,record_id=1:100,
#'   epistart=sample(c("2018-06-14",'2018-01-01'),100,replace=TRUE),
#'   diag_icd10=sample(c('I211','I212'),100,replace = TRUE), diag_icd9=sample(c('410','411'),100,replace=TRUE),
#'   oper4=sample(c('K400','K401'),100,replace = TRUE))
#' hesin_diag10 <- hesin[sample(1:100,80,replace = TRUE),c('eid','record_id','epistart')]
#' hesin_diag10$diag_icd10 <- sample(c('I211','I212','I213','I214'),80,replace = TRUE)
#' hesin_oper4 <- hesin[sample(1:100,80,replace = TRUE),c('eid','record_id','epistart')]
#' hesin_oper4$oper4 <- sample(c('K400','K401','K402','K403'),80,replace = TRUE)
#'
#' ## select cases by definition "icd10=I211 or oper4=K400"
#' icd10 <- 'I211'
#' oper4 <- 'K400'
#' case <- ukbcase(hesin=hesin,hesin_diag10=hesin_diag10,hesin_oper4=hesin_oper4,icd10=icd10,oper4=oper4)
#' summary(case)
#'

ukbcase <- function(icd10=NULL,icd9=NULL,oper4=NULL,histology=NULL,
                    hesin=NULL,hesin_diag10=NULL,hesin_oper4=NULL,
                    deathdata=NULL,cancerdata=NULL,
                    icd10main=T,icd10sec=T,icd9main=T,oper4main=T,oper4sec=T,
                    death=T,cancer=T) {

  case <- NULL # initialization
  if(is.null(icd10) & is.null(icd9) & is.null(oper4) & is.null(histology)) {stop("check the definition of disease")}

  if(((icd10main | icd10sec | icd9main) & (!is.null(icd10))) == T){
    if(is.null(hesin)){stop('main hesin dataset is required')}
    case_diag <- ukbcase_diag(icd10,icd9,hesin,hesin_diag10,icd10main,icd10sec,icd9main)
    case <- case_diag
  }

  if(((oper4main | oper4sec) & (!is.null(oper4))) == T){
    if(is.null(hesin)){stop('main hesin dataset is required')}
    case_oper <- ukbcase_oper(oper4,hesin,hesin_oper4,oper4main,oper4sec)
    if(is.null(case)){case <- case_oper}
    else{case <- merge(case_diag,case_oper,all=T)}
  }

  if(death==T & (!is.null(deathdata))){
    deathcase <- deathdata[which(deathdata$cause_of_death %in% icd10),]
    deathcase$death_description <- 'death'
    deathcase$date_of_death <- as.Date(deathcase$date_of_death)
    names(deathcase) <- c('eid','epistart','diag_icd10','death_description')
    if(is.null(case)){case <- deathcase}
    else{case <- merge(case, deathcase, all=T)}
  }

  if(cancer==T & (!is.null(cancerdata))){
    if(is.null(icd10) & is.null(icd9) & is.null(histology)) {stop("check the definition of disease")}
    cancercase <- cancerdata[unique(c(which(cancerdata$icd10 %in% icd10),which(cancerdata$icd9 %in% icd9),
                                      which(cancerdata$histology %in% histology))),c('eid','date','icd10','icd9','histology')]
    cancercase$date <- as.Date(cancercase$date)
    names(cancercase) <- c('eid','epistart','diag_icd10','diag_icd9','cancer_histology')
    if(is.null(case)){case <- cancercase}
    else{case <- merge(case, cancercase, all=T)}
  }

  # get unique case with the earlist onset date
  case <- case[!is.na(case$epistart),]
  colid<-sapply(unique(case$eid),function(x){
    col<-which(case$eid==x & case$epistart==min(case$epistart[which(case$eid==x)]))[1]
    return(col)
  })
  case <- case[colid,]
  case <- case[!is.na(case$eid),]
  return(case)
}

###################### diag
ukbcase_diag <- function(icd10=NULL,icd9=NULL,hesin=NULL,hesin_diag10=NULL,
                         icd10main=T,icd10sec=T,icd9main=F) {
  hesin_diag <- NULL # initialization
  if(icd10main==F & icd10sec==F & icd9main==F) {stop('at least one icd code should be given when using function ukbcase_diag')}

  # main icd10 from hesin dataset
  if (icd10main==T & (!is.null(icd10))) {
    hesin_diag_1 <- hesin[which(hesin$diag_icd10 %in% icd10),c('eid','record_id','diag_icd10','epistart')]
    hesin_diag_1$epistart <- as.Date(hesin_diag_1$epistart)
    hesin_diag <- hesin_diag_1
  }

  # secondary icd10 from hesin_diag10 subdataset
  if(icd10sec==T & (!is.null(hesin_diag10))) {
    hesin_diag_2 <- hesin_diag10[which(hesin_diag10$diag_icd10 %in% icd10),c('eid','record_id','diag_icd10')]
    # find the date of each record from the main hesin dateset
    hesin_diag_2_date <- hesin[which(hesin$record_id %in% hesin_diag_2$record_id),]
    hesin_diag_2 <- merge(hesin_diag_2_date,hesin_diag_2, by='record_id',all=F)
    hesin_diag_2$epistart <- as.Date(hesin_diag_2$epistart)
    hesin_diag_2 <- hesin_diag_2[,c('eid.x','record_id','diag_icd10.y','epistart')]
    names(hesin_diag_2) <- c('eid','record_id','diag_icd10','epistart')

    if(is.null(hesin_diag)){
      hesin_diag <- hesin_diag_2
    }else{
      hesin_diag <- merge(hesin_diag_1,hesin_diag_2,all=T)
      names(hesin_diag) <- c('eid','record_id','diag_icd10','epistart')
    }
  }

  # main icd9 from hesin dataset
  if (icd9main==T & (!is.null(icd9))) {
    hesin_diag_3 <- hesin[which(hesin$diag_icd9 %in% icd9),c('eid','record_id','diag_icd9','admidate')]
    hesin_diag_3$admidate <- as.Date(hesin_diag_3$admidate)
    names(hesin_diag_3) <- c('eid','record_id','diag_icd9','epistart')
    
    if(is.null(hesin_diag)){
      hesin_diag <- hesin_diag_3
    }else{
      hesin_diag <- merge(hesin_diag,hesin_diag_3,all=T)
    }
  }

  # select the epi with earlist date
  hesin_diag <- hesin_diag[which(!is.na(hesin_diag$eid)),]
  colid<-sapply(unique(hesin_diag$eid),function(x){
    col<-which(hesin_diag$eid==x & hesin_diag$epistart==min(hesin_diag$epistart[which(hesin_diag$eid==x)]))[1]
    return(col)
  })
  hesin_diag <- hesin_diag[colid,]
  hesin_diag <- hesin_diag[which(!is.na(hesin_diag$eid)),]

  return(hesin_diag)
}

###################### oper
ukbcase_oper <- function(oper4=NULL,hesin=NULL,hesin_oper4=NULL,
                         oper4main=T,oper4sec=T) {
  hesin_oper <- NULL # initialization
  if(oper4main==F & oper4sec==F) {stop('at least one oper4 code should be given when using function ukbcase_oper')}
  # main oper4 from hesin dataset
  if (oper4main==T) {
    hesin_oper_1 <- hesin[which(hesin$oper4 %in% oper4),c('eid','record_id','oper4','epistart')]
    hesin_oper_1$epistart <- as.Date(hesin_oper_1$epistart)
    hesin_oper <- hesin_oper_1
  }

  # secondary oper4 from hesin_oper4 subdataset
  if(oper4sec==T & (!is.null(hesin_oper4))) {
    hesin_oper_2 <- hesin_oper4[which(hesin_oper4$oper4 %in% oper4),c('eid','record_id','oper4')]

    # find the date of each record from the main hesin dateset
    hesin_oper_2_date <- hesin[which(hesin$record_id %in% hesin_oper_2$record_id),]
    hesin_oper_2 <- merge(hesin_oper_2_date,hesin_oper_2, by='record_id',all=F)
    hesin_oper_2$epistart <- as.Date(hesin_oper_2$epistart)
    hesin_oper_2 <- hesin_oper_2[,c('eid.x','record_id','oper4.y','epistart')]
    names(hesin_oper_2) <- c('eid','record_id','oper4','epistart')

    if(is.null(hesin_oper)){
      hesin_oper <- hesin_oper_2
    }else{
      hesin_oper <- merge(hesin_oper,hesin_oper_2,all=T)
    }
  }

  # uniqueness, earlist onset
  hesin_oper <- hesin_oper[which(!is.na(hesin_oper$eid)),]
  colid<-sapply(unique(hesin_oper$eid),function(x){
    col<-which(hesin_oper$eid==x & hesin_oper$epistart==min(hesin_oper$epistart[which(hesin_oper$eid==x)]))[1]
    return(col)
  })
  hesin_oper <- hesin_oper[colid,]
  names(hesin_oper) <- c('eid','record_id','oper4','epistart')
  hesin_oper <- hesin_oper[which(!is.na(hesin_oper$eid)),]

  return(hesin_oper)
}
