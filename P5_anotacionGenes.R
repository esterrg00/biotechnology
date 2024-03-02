######################################################################
## Biotecnología
## Anotación de genes

######################################################################
## Consideremos una situación más simple. Consideremos anotados los
## distintos codones que aparecen en uno de los marcos de lectura de
## una secuencia de ADN con los símbolos del alfabeto A = {M, P, X} (M
## := codón de inicio, P:= codón de parada, X := cualquier otro
## codón).

## Los siguientes ejercicios van encaminados a localizar los distintos
## ORF existentes en secuencias sobre dicho alfabeto. Para las
## distintas pruebas utilizaremos las secuencias x1, x2, x3 y x4 que
## puede encontrar en el fichero P5_anotacionGenes.RData.
######################################################################
## (1) Localizar los codones de inicio, aquellos etiquetados con M, de
## la secuencia x1.

## Solución: están en las posiciones 2, 7 y 13

inicioX1 = c()
for (i in 1:length(x1)){
    if(x1[i] == "M"){
        inicioX1 <- c(inicioX1,i)
    }
}

inicioX1

## SOLUCION CLASE:

## posMX1 <- which(x1=="M")



######################################################################
## (2) Localizar el primer codón de parada posterior al primer codón
## de inicio de x1. ¿Y el de x3?

## Solución: está en la posición 11, el de x3 está en la posición 12
## (pese a que hay un codón de parada en la posición 5, el primer
## codón de inicio de x3 está en la posición 6, por lo que hay que
## buscar uno posterior).

terminal1X1 <- c()
for (i in inicioX1[1]:length(x1)){
    if(x1[i] == "P"){
        terminal1X1 <- c(i)
        break
    }
}

terminal1X1

inicioX3 <- c()
for (i in 1:length(x3)){
    if(x3[i] == "M"){
        inicioX3 <- c(inicioX3,i)
    }
}
inicioX3

terminal1X3 <- c()
for (i in inicioX3[1]:length(x3)){
    if(x3[i] == "P"){
        terminal1X3 <- c(i)
        break
    }
}
terminal1X3



######################################################################
## (3) Localizar el primer codón de parada posterior a cada codón de
## inicio de x1.

## Están en las posiciones 11, 11 y 16

######### SOLUCIÓN CLASE

codonParada <- function(x, i){
    posMx <- which(x=="M") # posiciones de codones inicio
    MiX <- posMx[i] #codon de inicio en la posición i
    n <- length(x)
    for(j in (MiX+1):n){
        if(x[j] =="P"){
            break
        }
    }
    return(j)
}


posMX1 <- which(x1=="M")
res <- c()
for (k in 1:length(posMX1)){
    res <- c(res, codonParada(x1, k))
}

res


## Otra forma más efectiva

posM <- which(x1=="M")
paradas <- which(x1=="P")
# Nos quedamos con el primer elemento de los codones de parada cuya posicion sea 
# superior a la posicion que tiene dicho codon de inicio
posP <- mapply(function(posM){paradas[paradas > posM][1]}, posM)
posP


######################################################################
## (4) Resolver el ejercicio anterior para las secuencias x2, x3 y
## x4. Si no hay un codón de parada posterior a uno de inicio,
## devolver el tamaño de la secuencia.

## Nota: Puede ser de utilidad la función is.na

# Para x2
posM2 <- which(x2=="M")
paradas2 <- which(x2=="P")
posP2 <- mapply(function(pM2){paradas2[paradas2 > pM2][1]}, posM2)
posP2
x2

# Para x3
posM3 <- which(x3=="M")
paradas3 <- which(x3=="P")
posP3 <- mapply(function(pM3){paradas3[paradas3 > pM3][1]}, posM3)
posP3
x3

# Para x4
posM4 <- which(x4=="M")
paradas4 <- which(x4=="P")
posP4 <- mapply(function(pM4){paradas4[paradas4 > pM4][1]}, posM4)
posP4

# Para asegurarnos de si hay valores NA (no disponible), devolvamos el tamaño 
# de la posicion siguiente
posP4NA <- mapply(function(pP4){if (is.na(pP4)){
    length(x4)+1} else{pP4}
}, posP4)
posP4NA

    # OTRA OPCIÓN del de arriba: si vemos que ya hemos cambiado uno, cambiamos el resto
## no está bien; revisar ##
posP4NAv2 <- mapply(function(pP4){if (is.na(paradas4 > posM4[1])){
    length(x4)+1} else{paradas4[paradas4 > posM4[1]]}
}, posM4)
posP4NAv2


######################################################################
## (5) Definir una función localizaORF que, dada una secuencia sobre
## el alfabeto A devuelva un marco de datos con la información
## relativa a los distintos ORF existentes en dicha secuencia:
## * cInicial: posición del codón inicial
## * nCodones: número de codones que contiene (incluidos el codón
##   inicial y el de parada)

## > localizaORF(x1)
##   cInicial nCodones
## 1        2       10
## 2        7        5
## 3       13        4
## > localizaORF(x2)
##   cInicial nCodones
## 1        6        7
## 2       10        3
## 3       17        4
## > localizaORF(x3)
##   cInicial nCodones
## 1        6        7
## > localizaORF(x4)
##   cInicial nCodones
## 1       18        4
## 2       32       22
## 3       39       15
## 4       50        4


localizaORF <- function(s){
    
    inicio = c() # almacena las posiciones de los codones de inicio
    for (i in 1:length(s)){
        if(s[i] == "M"){
            inicio <- c(inicio,i)
        }
    }
    
    todas <- c() #almacena en la posición k todos los codones terminales que se 
               # encuentren tras el codón de inicio k
    for (j in 1:length(inicio)){
        paradaj <- c()
        if(j+1 < length(inicio)){
            for(m in s[j]:s[j+1]){
                if(s[m] == "M"){
                    paradaj <- c(paradaj,m)
                }
            }
        } else{
            for(m in s[j]:s[length(s)]){
                if(s[m] == "M"){
                    paradaj <- c(paradaj,m)
                }
            }
        }
        todas[[j]] <- paradaj
    }
    
    
}

## SOLUCIÓN CLASE
## falla algo; revisar

localizaORF <- function(x){
    l <- length(x)
    inicios <- which(x=="M") 
    paradas <- which(x=="P")
    # función anónima que la aplicamos sobre las posiciones de los codones de 
    # inicio; para codon de inicio, nos qudamos con su codon de parada y lo 
    # guardamos en fin; si la posicion no existe, lo guardaremos en la posicion 
    # l que es la longitud de la secuencia
    finales <- mapply(function(i){
        fin <- paradas[paradas > i][1]
        ifelse(is.na(fin), l, fin)},
        inicios)
    # contiene tantos registros como codones de inicio, donde mostraremos lo 
    # que se pide
    data.frame(cInicial= inicios,
               nCodones = finales -inicios-1)
}

localizaORF(x1)
localizaORF(x2)
localizaORF(x3)
localizaORF(x4)


######################################################################
## Generamos una secuencia al azar sobre el alfabeto A con una
## probabilidad 1/64 de ser M y 3/64 de ser P.

m1 <- sample(c("M", "P", "X"), 610046, replace = TRUE, 
             prob = c(1/64, 3/64, 60/64))
m1

######################################################################
## Mostramos gráficamente, con un histograma, la longitud de los
## distintos ORF de dicha secuencia.

azar <- localizaORF(m1)
library("lattice")

histogram(azar$nCodones, type = "percent",
          ylab = "Porcentaje",
          xlab = "Longitud",
          main = paste("ORF en una secuencia aleatoria de ",
                       length(m1), " codones
M(1/64, 3/64, 60/64)", sep = ""),
          xlim = c(-20, 250), ylim = c(-5, 45),
          panel = function (x, ...) {
              panel.histogram(x, ...)
              panel.abline(h = c(1,5), col = "red")
              panel.text(c(-11, -11), c(1.6, 5.6), cex = 0.75,
                         col = "red", labels = c("1%", "5%"))
})


######################################################################
## (6) Obtener un vector con la secuenciación del genoma completo del
## de H. Influenzae (ID NC_000907) de la base de datos Nucleotide del
## NCBI.

## > length(H.influenzae)
## 1830138
## > head(H.influenzae, 18)
## "T" "A" "T" "G" "G" "C" "A" "A" "T" "T" "A" "A" "A" "A" "T" "T" "G"
## "G"

library(seqinr)
H.influenzae <-
    getSequence(read.fasta("sequence.fasta"), forceDNAtoLower=FALSE)[[1]]


head(H.influenzae)

######################################################################
## (7) Definir una función codones que, dado un vector con una
## secuencia de nucleótidos, devuelva un vector con la primera lectura
## de codones en dicha cadena (los codones del primer marco de
## lectura directo).

## Ejemplo:
## > c1.H.influenzae <- codones(H.influenzae)
## > head(c1.H.influenzae)
## [1] "TAT" "GGC" "AAT" "TAA" "AAT" "TGG"


seq(1,7, by =3)

codones <- function(x,l = length(x)){
    ic <- seq(1,l-2, by=3)
    mapply(function(i){ c2s(x[i:(i-2)]) }, ic)
}

codones(H.influenzae)

head(H.influenzae)
codones(head(H.influenzae))

######################################################################
## (8) Definir una función etiqueta.codones que, dado un vector con
## una secuencia de codones, devuelva el correspondiente vector sobre
## el alfabeto A = {M, P, X} indicando si es un codón de inicio, ATG,
## uno de parada (TAA, TAG o TGA) o cualquier otro.

## > head(etiqueta.codones(c1.H.influenzae))
## "X" "X" "X" "P" "X" "X"


# recibe una lista x de CODONES
etiqueta.codones <- function(x){
    res <- rep ("X", length(x))
    res[x=="ATG"] <- "M"
    res[x=="TAA" | x=="TAG" | x=="TGA"] <- "P"
    res
}


etiqueta.codones(codones(head(H.influenzae,30)))


######################################################################
## La función calcularORF, dada una secuencia de ADN y un tamaño
## mínimo (n. de codones); devuelve un marco de datos con la
## información relativa a todos los ORF con tamaño superior al mínimo
## proporcionado. Las columnas del marco de datos serán:

## * mLectura: el marco de lectura (+1, +2, +3, -1, -2 -3),
## * nInicial: la posición, dentro de la secuenciación del genoma, del
##   primer nucleótido del codón "ATG" que marca el inicio del ORF y, 
## * nCodones: el número de codones que componen el ORF (sin incluir
##   el de parada).

## Además, el marco de datos está ordenado (de mayor a menor) por el
## número de codones y, en caso de empate, por la posición del codón
## de inicio.

## Comparar los resultados con los obtenidos con la herramienta
## ORFfinder de NCBI (https://www.ncbi.nlm.nih.gov/orffinder/)

calculaORF.Marco <- function (x, lo, lx, r, m) {
    cx <- codones(x, lx)
    ex <- etiqueta.codones(cx)
    lc <- length(cx)
    inicios <- which(ex == "M")
    paradas <- which(ex == "P")
    finales <- mapply(function (i) {
                          fin <- paradas[paradas > i][1]
                          ifelse(is.na(fin), lc, fin)},
                      inicios)
    orf <- (finales - inicios + 1) > lo
    data.frame(marco = rep(m, sum(orf)),
               nInicial = if (r >= 0) { inicios[orf] * 3 - 2 + r }
                          else {lx - inicios[orf] * 3 + 2 - r },
               nCodones = finales[orf] - inicios[orf])
}

calcularORF <- function (x, lo) {
    lm <- length(x)
    p1 <- calculaORF.Marco(x, lo, lm, 0, "+1")
    p2 <- calculaORF.Marco(x[-1], lo, lm-1, 1, "+2")
    p3 <- calculaORF.Marco(x[-c(1,2)], lo, lm-2, 2, "+3")

    C.x <- rev(comp(x, forceToLower = FALSE))

    n1 <- calculaORF.Marco(C.x, lo, lm, - 1, "-1")
    n2 <- calculaORF.Marco(C.x[-1], lo, lm-1, -1, "-2")
    n3 <- calculaORF.Marco(C.x[-c(1,2)], lo, lm-2, -1, "-3")

    tabla <- rbind(p1,p2,p3,n1,n2,n3)
    tabla[order(tabla$nCodones, tabla$nInicial, decreasing = TRUE),]
}

orf.H.influenzae <- calcularORF(H.influenzae, 150)
head(orf.H.influenzae[orf.H.influenzae$nInicial <= 5e4,], 30)

