library(pgirmess)

# abrir um dataframe com 3 colunas (Latitude, longitude e variável de interesse)
xy <- na.omit(xy)
coords <- xy[, c("lat","lon")] #Latitude e longitude nas colunas 1 e 2
corr_log <- correlog(coords, xy$extent, nbclass = 7) #resultados do I de Moran. Medida é a variável de interesse
plot(corr_log) #correlograma
