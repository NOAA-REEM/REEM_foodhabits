#' Temperature scaling consumption function; Fish Bioenergetics model
#' This function calculates the temperature and weight specific consumption rate for an individual fish
#'
#' Rcode developed by: Kirstin Holsman
#' kirstin.holsman@noaa.gov
#' This function calculates the weight and tempature specific physiological paramaters given a set of parameters and data inputs.
#' The model below is based on the "Wisconsin" Fish Bioenergetics model (Hanson et al. 1997)
#'
#' Based on the Fish Bioenergetics Model 3.0, see update 4.0 by Deslauriers, D.,
#' Chipps, S.R., Breck, J.E., Rice, J.A. & Madenjian, C.P. 2017. Fish Bioenergetics 4.0:
#' An R- based modeling application. Fisheries, 42(11): 586–596.
#'
#' @param par is a list of parameters used for the model
#' @param data is a list with a daily observation of weight, temperature, prey energy density, and predator energy density
#' @param data$TempC  Temperature
#' @param Ceq is a par list object (par$Ceq) that specifies the type of consumption equation to use: 0 = specify function using data$fTCmodel, 1= exponential (type 1 from Fish Bioenergetics), 2 = type 2, 3 = type 3
#' @keywords Temperature, scaling, Respiration
#' @export fTC_fun
#' @examples
#' c_data<-ebs_data
#' c_data$fTCmodel<-function(TempC){-.5*TempC}
#' fTC_fun(par=plk_par,data=c_data)

fTC_fun<-function(par,data){

  TempC  <-  data$TempC
  Ceq    <-  par$Ceq

  if(Ceq==0){
    fTc  <-  data$fTCmodel(TempC)
  }else if  (Ceq==1){
    fTc  <-  exp(par$QC*TempC)
  }else if (Ceq==2){
    Tco  <-  par$Tco
    Tcm  <-  par$Tcm
    Yc   <-  log(par$QC)*(Tcm-Tco+2)
    Zc   <-  log(par$QC)*(Tcm-Tco)
    Vc   <-  (Tcm-TempC)/(Tcm-Tco)
    Xc   <-  ((Zc^2)*(1+((1+40/Yc)^0.5))^2)/400
    fTc  <-  (Vc^Xc)*exp(Xc*(1-Vc))
  }else if (Ceq==3){
    CK1  <-  par$CK1
    CK4  <-  par$CK4
    Tcl  <-  par$Tcl
    Tco  <-  par$Tco
    Tcm  <-  par$Tcm
    G2   <-  (1/(Tcl-Tcm))*log((0.98*(1-CK4))/(CK4*0.02))
    L2   <-  exp(G2*(Tcl-TempC))
    Kb   <-  (CK4*L2)/(1+CK4*(L2-1))
    G1   <-  (1/(Tco-par$QC))*log((0.98*(1-CK1))/(CK1*0.02))
    L1   <-  exp(G1*(TempC-par$QC))
    Ka   <-  (CK1*L1)/(1+CK1*(L1-1))
    fTc  <-  Ka*Kb

  }else{

    message("ERROR: Ceq does not match criteria [0,3]")

  }

  return( fTc = fTc )

}
