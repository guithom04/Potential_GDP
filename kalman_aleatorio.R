# pacotes 
library(sidrar)
library(FKF)
library(KFAS)

# dados
# pib
gdp = sidrar::get_sidra(api = "/t/6613/n1/all/v/all/p/all/c11255/90707/d/v9319%202")
GDP = ts(gdp$Valor, start = c(1996,01), frequency = 4)
# taxa de desemprego
un = sidrar::get_sidra(api = "/t/6396/n1/all/v/4099/p/all/c2/6794/d/v4099%201")
UN = ts(un$Valor, start = c(2012,01), frequency = 4)


# estimar


library(KFAS)

# Define matrices
Z <- matrix(c(1, 0), nrow=1)
T <- matrix(c(1, 1, 0, 1), ncol=2)
R <- matrix(c(0, 1), ncol=1)
Q <- matrix(c(0, 0, 0, NA), ncol=2)  # We'll estimate the variance for the growth rate

# Set up the state space model
model <- SSModel(GDP ~ -1 + SSMcustom(Z=Z, T=T, Q=Q, R=R), H=NA)

# Estimate the model using maximum likelihood
fit <- fitSSM(model, method="BFGS")

# Check the results
print(fit)



??SSMcustom()





library(KFAS)

# Define matrices
Z <- matrix(c(1, 0), nrow=1)
T <- matrix(c(1, 1, 0, 1), ncol=2)
R <- matrix(c(0, 1), nrow=2)
Q <- matrix(c(0, 0, 0, NA), ncol=2)  # We'll estimate the variance for the growth rate

# Set up the state space model using SSMcustom
model <- SSModel(GDP ~ -1 + SSMcustom(Z=Z, T=T, R=R, Q=Q), H=NA)

# Estimate the model using maximum likelihood
fit <- fitSSM(model, method="BFGS")

# Check the results
print(fit)




# 


# Carregar o pacote necessário
library(dlm)

# Suponha que y seja seus dados de PIB e u seja seus dados de desemprego
y <- PIB
u <- UN

# Criar função de construção para o modelo dlm
buildModel <- function(par) {
  dlm(
    FF = matrix(c(1, 0, 0, 0,
                  0, par[1], 0, 1), 2, 4),
    GG = matrix(c(1, 0, -1, 0,
                  1, 1, -1, 0,
                  0, 0, 1, 0,
                  0, 0, 0, par[2]), 4, 4),
    V = diag(c(exp(par[3]), exp(par[4])), 2),
    W = diag(exp(par[5:8]), 4)
  )
}

# Parâmetros iniciais para otimização
initPar <- c(-0.2, 0.03, log(1e-2), log(1e-2), log(1e-2), log(1e-2), log(1e-2), log(1e-2))

# Estimativa de Máxima Verossimilhança
fit <- dlmMLE(y = cbind(y, u), parm = initPar, build = buildModel)

# Modelo estimado
dlmFit <- buildModel(fit$par)
dlmSmoothed <- dlmSmooth(cbind(y, u), mod = dlmFit)

# PIB potencial estimado
y_pot <- ts(dlmSmoothed$s[-1, 1], start = start(y))

# Ver resultados
plot(y_pot)




# kfas


ssm <- function(model, signalToNoise) {
  
  # get trend
  trend <- attr(model, "trend")
  
  # obtain system variances
  Z <- model$SSModel$Z
  T <- model$SSModel$T
  R <- model$SSModel$R
  Q <- model$SSModel$Q
  a1 <- model$SSModel$a1
  P1 <- model$SSModel$P1
  P1inf <- model$SSModel$P1inf
  H <- model$SSModel$H
  stateNames <- colnames(model$SSModel$T)
  
  # get rid of trend variance
  if (trend != "RW1") {
    varNames <- model$loc$variableRow[model$loc$sysMatrix == "Q"]
    name_delete <- stateNames[grepl("trend", stateNames)][1]
    index_delete <- which(R[name_delete, , ] == 1)
    R <- R[, -index_delete, ]
    Q <- Q[-index_delete, -index_delete, ]
    model$loc <- model$loc[!(model$loc$sysMatrix == "Q" & model$loc$variableRow == "trend"), ]
  }
  
  # signal to noise ratio specified
  if (!is.null(signalToNoise)) {
    model$loc <- model$loc[!(model$loc$sysMatrix == "Q" & grepl("trend", model$loc$variableRow)), , drop = FALSE]
    if (trend != "RW1") {
      varNames <- model$loc$variableRow[model$loc$sysMatrix == "Q"]
      index <- varNames != "trend"
      Q <- Q[index, index]
      R <- R[, index]
    }
  }
  
  # ----- state space model
  if (inherits(model, "TFPmodel")) {
    modelSS <- SSModel(cbind(logtfp, cubs) ~ -1 + SSMcustom(Z = Z, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf, state_names = stateNames),
                       H = H,
                       data = model$tsl
    )
  } else if (inherits(model, "NAWRUmodel")) {
    modelSS <- SSModel(cbind(ur, pcInd) ~ -1 + SSMcustom(Z = Z, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf, state_names = stateNames),
                       H = H,
                       data = model$tsl
    )
  } else if (inherits(model, "KuttnerModel")) {
    modelSS <- SSModel(cbind(loggdp, dinfl) ~ -1 +
                         +SSMcustom(Z = Z, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf, state_names = stateNames),
                       H = H,
                       data = model$tsl
    )
  }
  
  # return
  model$SSModel <- modelSS
  model
}
