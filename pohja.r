# CTRL+Shift+C kommentoi tai poistaa kommentit valituilta riveilt?
library(forecast)
# Asetetaan ty?hakemisto SET Working Directory
# Muokkaa tarvittaessa oman ty?hakemistosi polku
# RStudiossa t?m?n voi tehd? my?s valikosta Session->Set working directory
setwd('Z:/Documents')

# Luetaan s?hk?nkulutus- ja l?mp?tiladata, hyp?t??n headerrivin yli
eletemp = read.table(file = "sahko.csv",
                     sep = ";",
                     dec = ",",
                     skip = 1,
                     col.names = c('kWh','Celcius'))

# S?hk?nkulutus aikasarjaksi
ele = ts(eletemp$kWh[1:816], start = 1, frequency = 24)
elemodel=auto.arima(ele)
elemodel
# L?mp?tila kahdeksi aikasarjaksi: 816 ensimm?ist? havaintoa k?ytet??n mallin estimointiin
# ja 24 viimeist? havaintoa ennustamiseen.
temp = ts(eletemp$Celcius, start = 1, frequency = 24)
temp816 = ts(eletemp$Celcius[1:816], start = 1, frequency = 24)
# start parametrina vektori: 817. tunti = 35. p?iv?n ensimm?inen tunti
temp24 = ts(eletemp$Celcius[817:840], start = c(35,1), frequency = 24)
par(mfrow=c(2,1))
# Plotataan aikasarjat
ts.plot(ele,
     xlab = "aika/vrk",
     ylab = "kulutus/kWh")
                    

ts.plot(temp816,temp24,
        xlab = "aika/vrk",
        ylab = expression(~degree~C),
        col = c("black", "blue"))

# M??ritell??n 2x2 plottausruudukko.
par(mfrow=c(2,2))
# Jos haluat katsoa kuvia 2x2 matriisin sijaan yksitellen, niin
# par(mfrow=c(1,1))

# Plotataan autokorrelaatio-, osittaisautokorrelaatio- ja ristikorrelaatiofunktiot.
acf(ele, lag.max=168)
acf(ele, lag.max=168, type = "partial")
acf(temp, lag.max=168)
acf(temp, lag.max=, type = "partial")
# Piirret??n ristikorrelaatiofunktio omaan kuvaan
par(mfrow=c(1,1))
ccf(ele,temp, lag.max=168)
par(mfrow=c(2,2))
# Stationarisoidaan aikasarjat. M??rittele parametrit d,S,D
# Huomaa, ett? s?hk?nkulutuksen ja l?mp?tilan aikasarjojen differointien asteiden ei v?ltt?m?tt? tarvitse olla samoja.
d = 1
d =1  # Differoinnin kertaluku d
S1 = 24
S2 = 24     # Kausidifferoinnin jakso S

D = 1 # Kausidifferensoinnin kertaluku D
dele = ele
dtemp = temp
if (d2 > 0) {
  dele = diff(dele, lag = 1, differences = d)
  dtemp = diff(dtemp, lag = 1, differences = d)
}
if (D > 0) {
  dele = diff(dele, lag = S1, differences = D)
  dele = diff(dele, lag = 168, differences = D)
  dtemp = diff(dtemp, lag = S2, differences = D)
}

# Differoitujen aikasarjojen autokorrelaatio-, osittaisautokorrelaatio- ja ristikorrelaatiofunktiot.
acf(dele, lag.max=816)
acf(dele, lag.max=816, type = "partial")

acf(dtemp, lag.max=816)
acf(dtemp, lag.max=816, type = "partial")
ccf(dele, dtemp, lag.max=816)

# Estimoidaan malli ja lasketaan ennusteet ilman ulkoista muuttujaa.
p = 1
q = 1
P = 2
Q = 0
malli = arima(ele,
              order = c(p,d,q),
              seasonal = list(order = c(P, D, Q),period = S),
              method = "CSS")
enne = predict(malli, n.ahead = 24)
plot(enne)
# Estimoidaan malli l?mp?tilan kanssa. M??r?? l?mp?tilan mahdollinen viive L.
L = 1
tempestimointi = eletemp$Celcius[1:(816-L)]
tempennuste = eletemp$Celcius[(816-L+1):(816-L+24)]
eleestimointi = ts(eletemp$kWh[(1+L):816], start = 1, frequency = 24)
malli2 = arima(eleestimointi,
               order = c(p,d,q),
               seasonal = list(order = c(P, D, Q), period = S),
               xreg = tempestimointi,
               method = "CSS")
enne2 = predict(malli2,
                n.ahead = 24,
                newxreg = tempennuste)

# Esimerkki Portmanteau-testist?. Onko residuaaliaikasarjan alussa nollia?
Box.test(malli$residuals,
         lag = 20,
         type = "Ljung-Box",
         fitdf = p + q + P + Q)

# Palautetaan plottaus normaaliin 1x1 ruutuun
par(mfrow=c(1,1))

# Plotataan kuva sahkonkulutusaikasarjasta, mallin (1) sovitteesta,
# ennusteesta ja ennusteen 95 %:n eli 1,96*sigma-luottamusv?leist?.
ts.plot(ele,
        ele - malli$residuals,
        enne$pred,
        enne$pred + 1.96*enne$se,
        enne$pred - 1.96*enne$se,
        col = c("black", "red", "blue", "blue", "blue"),
        main  = "Sovite ja ennuste")

# Plotataan kuva pelk?st? ennusteesta.
ts.plot(enne$pred,
        enne$pred + 1.96*enne$se,
        enne$pred - 1.96*enne$se,
        col = c("black", "blue", "blue"),
        main = "Ennuste ja  95 %:n luottamusv?lit")

# Kirjoitetaan ennuste ja luottamusv?lit .csv-tiedostoon, jonka voi avata Excelill?.
output = cbind(enne$pred,
               enne$pred + 1.96*enne$se,
               enne$pred - 1.96*enne$se)
write.csv2(output, file = 'ennuste.csv')
