setwd('Z:/Downloads')
# Luetaan sähkönkulutus- ja lämpötiladata, hypätään headerrivin yli
eletemp = read.table(file = file.choose(),
                     sep = ";",
                     dec = ",",
                     skip = 1,
                     col.names = c('kWh','Celcius'))

# Sähkönkulutus ja lämpötila aikasarjoiksi
ele = ts(eletemp$kWh[1:816], start = 1, frequency = 24)
temp = ts(eletemp$Celcius, start = 1, frequency = 24)
temp816 = ts(eletemp$Celcius[1:816], start = 1, frequency = 24)
temp24 = ts(eletemp$Celcius[817:840], start = c(35,1), frequency = 24)

# Plotataan aikasarjat
ts.plot(ele, xlab = "aika/vrk", ylab = "kulutus/kWh")
ts.plot(temp816, temp24, xlab = "aika/vrk", ylab = expression(~degree~C), col = c("black", "blue"))

# Plotataan autokorrelaatio-, osittaisautokorrelaatio- ja ristikorrelaatiofunktiot
par(mfrow=c(2,2))
acf(ele, lag.max=816)
acf(ele, lag.max=168, type = "partial")
acf(temp, lag.max=814)
acf(temp, lag.max=168, type = "partial")
par(mfrow=c(1,1))
ccf(ele,temp, lag.max=816)


#Differoidut aikasarjat (viikottainen ja päivittäinen trendi):
dele = diff(diff(diff(ele, lag = 168), lag = 24), lag = 1)
dtemp = diff(diff(diff(temp, lag = 168), lag = 24), lag = 1)

#Tulostus:
par(mfrow=c(2,2))
acf(dele, lag.max=816, main="ACF: Differoitu sähkö")
acf(dele, lag.max=816, type = "partial", main="PACF: Differoitu sähkö")
acf(dtemp, lag.max=816, main="ACF: Differoitu lämpötila")
acf(dtemp, lag.max=816, type = "partial", main="PACF: Differoitu lämpötila")

#Ristikorrelaatio
par(mfrow=c(1,1))
ccf(dele, dtemp, lag.max=168, main="CCF: Differoidut sähkö ja lämpötila")

#Lähdetään rakentamaan mallia:
# Poistetaan viikkovaihtelu sähkönkulutuksesta manuaalisesti
dele168 = diff(ele, lag = 168, differences = 1)

#Valitaan kertaluvut:
p = 2
d = 0
q = 1
P = 1
D = 1
Q = 1

#Esitmoidaan malli ilman lämpötilaa (sarima)
malli = arima(dele168,
              order = c(p,d,q),
              seasonal = list(order = c(P, D, Q), period = 24),
              method = "CSS")
#Lämpötilan viive
L = 2

# Kohdistetaan sähkö ja lämpötila L tunnin viiveen mukaisesti
ele_mod = ts(eletemp$kWh[(1+L):816], frequency = 24)
temp_mod = ts(eletemp$Celcius[1:(816-L)], frequency = 24)

#Differoidaan molemmista aikasarjoista viikottainen trendi
dele168_mod = diff(ele_mod, lag = 168)
dtemp168_mod = diff(temp_mod, lag = 168)

# Estimoidaan malli 2, xreg-muuttujana differoitu lämpötila(Sarimax)
malli2 = arima(dele168_mod,
               order = c(p,d,q),
               seasonal = list(order = c(P, D, Q), period = 24),
               xreg = dtemp168_mod,
               method = "CSS")

#Ennustetta varten tarvitaan (t-L) hetken lämpötilat
temp_uusi = eletemp$Celcius[(817-L):(840-L)]
#niistä täytyy vähentää 168h takaiset arvot, jotta xreg pysyy differoituna!
temp_uusi_vanha = eletemp$Celcius[(817-L-168):(840-L-168)]
dtemp_uusi = temp_uusi - temp_uusi_vanha

# Lasketaan ennusteet mallille 2
enne2 = predict(malli2, n.ahead = 24, newxreg = dtemp_uusi)

#Tehdään ennuste:
# Koska syötimme ARIMA:lle differoitua dataa, predict antaa differoituja ennusteita.
# Palautetaan todellinen taso lisäämällä ennusteeseen viikon (168h) takainen toteutunut arvo.
ennuste_final = c(1:24)
clyla_final = c(1:24)
clala_final = c(1:24)

for (x in 1:24) {
  # Haetaan toteutunut kWh täsmälleen 168 tuntia (1 viikko) aiemmin
  historia_arvo = eletemp$kWh[816 - 168 + x]
  ennuste_final[x] = historia_arvo + enne2$pred[x]
  clyla_final[x]   = ennuste_final[x] + 1.96 * enne2$se[x]
  clala_final[x]   = ennuste_final[x] - 1.96 * enne2$se[x]
}

# Tehdään R:n aikasarjaobjekteja plottausta varten
ennuste_ts = ts(ennuste_final, start = c(35,1), frequency = 24)
clyla_ts   = ts(clyla_final, start = c(35,1), frequency = 24)
clala_ts   = ts(clala_final, start = c(35,1), frequency = 24)

#Diagnostiikka: Onko residuaaleissa enää säännönmukaisuutta?
#Haetaan residuaalit
acf(malli2$residuals, lag.max=816, main="Malli2: Residuaalien ACF")
Box.test(malli2$residuals,
         lag = 48,
         type = "Ljung-Box",
         fitdf = p + q + P + Q)

#Samat mallille ilman lämpötilan viivettä:
acf(malli$residuals, lag.max=816, main="Malli: Residuaalien ACF")
Box.test(malli$residuals,
         lag = 48,
         type = "Ljung-Box",
         fitdf = p + q + P + Q)

# Plotataan kuva pelkästä ennusteesta luottamusväleineen
ts.plot(ennuste_ts,
        clyla_ts,
        clala_ts,
        col = c("black", "blue", "red"),
        main = "Lämpötilamallin ennuste ja 95 % luottamusvälit",
        ylab = "kWh")

# Kirjoitetaan LOPULLISET (integroidut) ennusteet Exceliä varten csv:ksi!
output = cbind(ennuste_final, clyla_final, clala_final)


# Kohdistetaan sovite oikeaan kohtaan aikasarjaa (alkaa viipeiden L+168 jälkeen)
L = 2
alku_indeksi <- 1 + L + 168 
sovite_arvot <- ele[alku_indeksi:816] - malli2$residuals
sovite_ts <- ts(sovite_arvot, start = c(8, 3), frequency = 24)

par(mfrow=c(2,2))
# Plotataan kuva käyttäen palautettuja ts-objekteja (ennuste_ts, clyla_ts, clala_ts)
ts.plot(ele,
        sovite_ts,
        ennuste_ts,
        clyla_ts,
        clala_ts,
        col = c("black", "red", "green", "blue", "purple"),
        main  = "Sähkönkulutus: Sovite ja ennuste (SARIMAX)",
        ylab = "Kulutus (kWh)",
        xlab = "Aika (vrk)",
        lty = c(1, 1, 1, 2, 2),
        xlim = c(34, 36))

# Palautetaan malli1:n ennusteet todelliselle kWh-tasolle
ennuste_final1 = c(1:24)
clyla_final1 = c(1:24)
clala_final1 = c(1:24)

for (x in 1:24) {
  # Haetaan toteutunut kWh täsmälleen 168 tuntia (1 viikko) aiemmin
  historia_arvo = eletemp$kWh[816 - 168 + x]
  
  ennuste_final1[x] = historia_arvo + enne$pred[x]
  clyla_final1[x]   = ennuste_final1[x] + 1.96 * enne$se[x]
  clala_final1[x]   = ennuste_final1[x] - 1.96 * enne$se[x]
}

# Tehdään R:n aikasarjaobjekteja plottausta varten (malli1)
ennuste_ts1 = ts(ennuste_final1, start = c(35,1), frequency = 24)
clyla_ts1   = ts(clyla_final1, start = c(35,1), frequency = 24)
clala_ts1   = ts(clala_final1, start = c(35,1), frequency = 24)

# Kohdistetaan malli1:n sovite oikeaan kohtaan (alkaa viipeen 168 jälkeen)
alku_indeksi1 <- 1 + 168 
sovite_arvot1 <- ele[alku_indeksi1:816] - malli$residuals
sovite_ts1 <- ts(sovite_arvot1, start = c(8, 1), frequency = 24)

# Plotataan kuva mallille 1
ts.plot(ele,
        sovite_ts1,
        ennuste_ts1,
        clyla_ts1,
        clala_ts1,
        col = c("black", "red", "green", "blue", "purple"),
        main  = "Sähkönkulutus: Sovite ja ennuste (SARIMA ilman lämpötilaa)",
        ylab = "Kulutus (kWh)",
        xlab = "Aika (vrk)",
        lty = c(1, 1, 1, 2, 2),
        xlim = c(34, 36))

output = cbind(ennuste_final, clyla_final, clala_final)

write.csv2(output, file = "ennuste_output.csv", row.names = FALSE)
