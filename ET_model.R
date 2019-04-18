require(Evapotranspiration)
setwd("E:/R Programing/ET")
climatedata_shanghai <- read.csv("climatedata_shanghai.csv", header = T )
constants_shanghai <- as.list(read.csv("constants_shanghai.csv", header = T))

modeldata <- ReadInputs(varnames = c("Temp","Rs","RH","uz"),
                   climatedata_shanghai, 
                   constants_shanghai, 
                   stopmissing=c(10,10,3),
                   timestep = "subdaily",
                   interp_missing_days = FALSE, 
                   interp_missing_entries = FALSE, 
                   interp_abnormal = FALSE, 
                   missing_method = NULL, 
                   abnormal_method = NULL)


result_Penman <- ET.Penman(modeldata,constants_shanghai,solar = "data", 
                    wind= "yes", windfunction_ver = "1948", 
                    alpha = 0.08, z0=0.001)

result_PenmanMonteith <- ET.PenmanMonteith(modeldata, constants_shanghai, wind="yes",solar= "data")

result_PriestleyTaylor <- ET.PriestleyTaylor(modeldata, constants_shanghai, solar="data", alpha=0.23)

result_PenPan <- ET.PenPan(modeldata, constants_shanghai, solar="data", alpha=0.23)

result_Radsum <- ET.Radsum(modeldata, constants_shanghai, solar="data")

##
accurate <- function(x,x.hat,k,output = TRUE)
{
        n <- length(x)
        SST <- sum((x - mean(x))^2)
        SSE <- sum((x - x.hat)^2)
        MSE <- SSE/(n - k)
        RMSE <- sqrt(MSE)
        MAPE <- 100*mean(abs((x - x.hat)/x))
        MPE <- 100*mean((x - x.hat)/x)
        MAE <- mean(abs(x - x.hat))
        ME <- mean(x - x.hat)
        R2 <- 1 - SSE/SST
        ADJ.R2 <- 1 - (n - 1)*(1 - R2)/(n - k)
        z <- embed(x,2)
        RW.R2 <- 1 - (n - 1)*SSE/(n*sum(z[,1] - z[,2] - mean(z[,1] - z[,2])))
        AIC <- n*log(SSE/n) + 2*k
        SBC <- n*log(SSE/n) + k*log(n)
        APC <- ((n + k)/(n*(n - k)))*SSE
        result <- c(SST,SSE,MSE,RMSE,MAPE,MPE,MAE,ME,R2,ADJ.R2,RW.R2,AIC,SBC,APC)
        names(result) <- c("SST","SSE","MSE","RMSE","MAPE","MPE","MAE","ME","R.squared",
                           "R.adj.squared","RW.R.squared","AIC","SBC","APC")
        if (output) {
                cat("MAPE", "RMSE","R.squared","AIC","\n")
                cat( MAPE, RMSE, R2, AIC,"\n")
        }
        accurate <- result
}

dat <- read.table("clipboard", header = T)
accurate(dat$ETa, dat$ETc_Radsum,k=0)
accurate(dat$ETa, dat$ETc_Penman,k=0)
accurate(dat$ETa, dat$ETc_FAOPM,k=0)
accurate(dat$ETa, dat$ETc_PT,k=0)



