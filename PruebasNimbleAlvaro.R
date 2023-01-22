library(nimble)   

load("E:/Mi unidad/Docencia/2022-2023/TFM/Luis/R/DatosPCA.rda")

Modelo_EjercicioN <- nimbleCode(
  {
    for(i in 1:293)
    {
      # y[i, 1:7] recoge el valor para cada uno de las 293 localizaciones y 7 CP
      # Asumiremos una distribucion normal multivariante
      # mu[i, ] es el vector de medias para cada localizacion i
      # El valor de mu[i, ] vendrá determinado por mus[, 1], mus[, 2] o mus[, 3] en funcion del cluster
      # Seguimos trabajando como en el caso univariante para la distribucion de z
      
      y[i, 1:7] ~ dmnorm(mus[1:7, z[i]],wi_tau[1:7, 1:7])
      # y[i, 1:7] ~ dmnorm(mus[1:7, z[i]], tau[1:7, 1:7, z[i]])
     
      # y[i, 1:7] ~ dmnorm(aux_mus[1:7,z[i]], aux_tau[1:7, 1:7])
      
      # y[i, 1:7] ~ dmnorm(mus_i[1:7], tau_i[1:7, 1:7])
      # mus_i[1:7] <- mus[1:7, z[i]]
      # tau_i[1:7, 1:7] <- tau[1:7, 1:7, z[i]]
      # mu[i] <- mus[1:7, z[i]]
      # z[i] <- 1
      wi_tau[1:7, 1:7] <- w[i]*tau[1:7, 1:7]
      # wi_tau[1:7, 1:7] <- w[i]*tau[1:7, 1:7, z[i]]
      w[i] ~ dgamma(2,2)
      omega[1:3] ~ ddirich(dirch_alpha[1:3])
      z[i] ~ dcat(omega[1:3])
    }
    for (h in 1:3) 
    {
      # La distribución Wishart es idonea para estimadores de varianza, como en este caso tau
      # Sus parametros recogen la matriz de varianzas y los grados de libertad
      # En un principio era suficiente con indicar 7 grados de libertad
      # mus[1:7, h] ~ dmnorm(mu_0[1:7],tau_h[1:7, 1:7])
      # tau_h[1:7, 1:7] <- tau[1:7, 1:7, h]
      # tau[1:7, 1:7, h] ~ dwish(wish_V[1:7,1:7], 8)
      mus[1:7, h] ~ dmnorm(mu_0[1:7],tau_0[1:7, 1:7])
      # tau[1:7, 1:7, h] ~ dwish(wish_V[1:7,1:7], 8)
    }
    
    tau[1:7, 1:7] ~ dwish(wish_V[1:7,1:7], 8)
    
  }
)

set.seed(1)
# tau_ini <- array(rep(0,7*7*3), c(7, 7, 3))  
# for (i in 1:7){
#   for (h in 1:3){
#     tau_ini[i,i,h]=1
#   }
# }
tau_ini <- diag(1,7)
z_ini=sample(c(1,2,3),293,replace=T)
Iniciales_EjercicioN <- function()
{
  list(tau = tau_ini,
       z = z_ini,
       mus = matrix(rep(0,7*3),ncol=3),
       omega = c(1/3,1/3,1/3),
       w = rep(1,293))
}

Parametros_EjercicioN <- c("omega", "mus", "z", "tau", "w")
nimbleOptions(showCompilerOutput = F)
Resul_EjercicioN <- nimbleMCMC(Modelo_EjercicioN, data = list(y = y), 
                               constants = list(mu_0=rep(0,7),
                                                wish_V=diag(1,7),
                                                dirch_alpha=c(1,1,1),
                                                tau_0=diag(0.01,7)),
                               inits = Iniciales_EjercicioN,
                               monitors = Parametros_EjercicioN, 
                               thin = 10, niter = 20000, nburnin = 4000, nchains = 1, 
                               summary = TRUE, WAIC = TRUE)

# Comprobamos la unimodalidad de las distribuciones a posteriori para comprobar convergencia (a falta del RHat)

# Por ejemplo
plot(density(Resul_EjercicioN$samples[,"mus[1, 1]"]))

# Para la z, variable categorica, alguna tabla de frecuencias para ver si abunda una categoria
table(Resul_EjercicioN$samples[,"z[1]"])
