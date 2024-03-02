######################################################################
## Biotecnología
## Ejercicios de introducción a R

######################################################################

## Se pide:

## (a) Leer el marco o tabla de datos contenidos en el fichero
## femaleMiceWeights.csv y guardarlo en la variable miceWt.

miceWt <- read.csv("femaleMiceWeights.csv", header= TRUE, sep=";", dec=",")

## (b) Consultar los nombres de las columnas de la tabla/marco de
## datos anterior e indicar cual contiene el peso.

names(miceWt)

## (c) ¿Cuál es el peso de la decimosegunda muestra?

miceWt[12,"Bodyweight"]

## (d) Utiliza el operador $ para extraer la columna del peso y vuelve
## a consultar el de la decimosegunda muestra.

miceWt$Bodyweight[12]

## (e) ¿Cuántas muestras hay en la tabla miceWt?

nrow(miceWt)

## (f) Calcular el peso medio de los ratones que han seguido la dieta
## "high fat" (anotada como "hf" en la tabla).

mean(miceWt$Bodyweight[miceWt$Diet =="hf"])


## (g) Calcula la diferencia en el peso medio de los ratones según la
## dieta que hayan seguido.

abs(mean(miceWt$Bodyweight[miceWt$Diet =="hf"]) - 
      mean(miceWt$Bodyweight[miceWt$Diet =="chow"]))

ref <- abs(diff(tapply(miceWt$Bodyweight,mean)))

## (h) En el fichero femaleControlsPopulation.csv tenemos una muestra
## de ratones que no han seguido la dieta hf y queremos utilizarla
## para determinar si la diferencia anterior es fruto de la dieta
## seguida. Para ello seleccionar al azar dos poblaciones de 12
## ratones y calculamos la diferencia entre las medias de sus pesos.

female_control <- read.csv("femaleControlsPopulation.csv", sep="", header= TRUE, dec=",")

## Para la selección aleatoria de datos consultar la ayuda de la 
## función sample.

#Buscamos la diferencia entre grupos de muestras aleatorias de ambos ficheros
abs(mean(sample(female_control, 12, TRUE)) - mean(sample(female_control,12,TRUE)))


 ## Repetimos la operación 10000 veces (el vector difWtAzar contendrá
## dichos valores tras evaluar las siguientes expresiones) ¿Qué
## proporción de veces dicha diferencia supera a la del experimento?

difWAzar <- numeric(1e4) #arreglo numerico con 1000 iteraciones

for (i in 1:1e4){
  difWAzar[i] <-abs(mean(sample(female_control, 12, TRUE)) - mean(sample(female_control,12,TRUE)))
}

difWAzar

## El análisis de la significancia del resultado suele acompañarse
## de un histograma, marcando el valor de referencia y el percentil 95/99.
## Usualmente, por encima del 95 % se considera relativamente significativo
## y muy/realmente significativo si sobrepasa el 99 %

d <- density(difWAzar)
hist(difWAzar, frec =F)
lines(d, col="red")
abline(v = ref, col = "green") # aqui se encuentra la diferencia obtenida entre las dos dietas

p9 <- quantile(difWAzar, c(0.95,0.99))
abline(v = p9, col="blue") #dibuja una linea vertical por los cuartiles definidos


## Al estar entre el 95% y el 99% es relativamente significativo
## sin llegar a poder considerarse muy significativo.
