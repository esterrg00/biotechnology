######################################################################
## Biotecnología
## Modelos ocultos de Markov
######################################################################
## Modelo sobre un conjunto de secuencias, no vacío, cuyos símbolos
## pertenecen a un determinado alfabeto A = {a1, ..., am} en las que
## se observan n tipos distintos de región (fragmentos a lo largo de
## los cuales la composición se mantiene estable).

## Cada variable E{j} representa el estado (o tipo de región, zona,
## ...) al que pertenece la observación (símbolo, nucleótido ...) que
## ocupa la posición j en una secuencia: val(E{j}) = {r1, ..., rn}

## Cada variable X{j} representa la observación (símbolo, nucleótido
## ...) que ocupa la posición j en una secuencia: val(X{j}) = {a1,
## ..., am}

## Para manejar un modelo oculto de Markov con n posibles estados
## (tipos de región) y m posibles observaciones (símbolos que componen
## las secuencias) dispondremos de una lista con tres elementos:

## * pi: un vector con las probabilidades iniciales de las distintas
##       regiones. El nombre de cada dato será el correspondiente
##       tipo de región.

## r1               rn
## P(E{1} = r1) ... P(E{1} = rn)

## * A: una matriz de tamaño nxn con las probabilidades de transición
##      entre los distintos tipos de región (los nombres de las filas
##      y columnas serán los correspondientes tipos de región).

##     r1                             rn
## r1  P(E{j+1} = r1 | E{j} = r1) ... P(E{j+1} = rn | E{j} = r1)
##     ...
## rn  P(E{j+1} = r1 | E{j} = rn) ... P(E{j+1} = rn | E{j} = rn)


## * B: una matriz de tamaño nxm con las probabilidades de los
##      distintos símbolos en cada una de las posibles regiones (los
##      nombres de las filas serán los correspondientes tipos de
##      región y los nombres de las columnas los símbolos del
##      alfabeto).

##     a1                           am
## r1  P(X{j} = a1 | E{j} = r1) ... P(X{j} = am | E{j} = r1)
##     ...
## rn  P(X{j} = a1 | E{j} = rn) ... P(X{j} = am | E{j} = rn)

######################################################################
## (1) Descargar el genoma del bacterio fago lambda, identificado con
## NC_001416 de la base de datos NCBI, y cargar la función
## local.composition definida en una práctica anterior (se encuentra
## en el fichero P6_HMM.RData).
######################################################################

library(seqinr)
fago <- getSequence(read.fasta("sequence.fasta" ,forceDNAtolower = FALSE)[[1]])
fago



######################################################################
## (2) Crear un gráfico para comparar el contenido global y local en
## GC del genoma del bacterio fago lambda usando una ventana de
## longitud 500 y un desplazamiento de 50.
######################################################################


bfl <- local.composition(fago, GC, 500, 60)
plot(bfl$positions, bfl$local.results, xlab = "", ylab = "GC",
     main = "NC_001416", type="l", col ="blue", lty =1, ylim= c(0,1))
abline(h = GC(fago), col = "red", lty = 1)


######################################################################
## (3) Definir la lista ricoCG_AT con la representación del siguiente
## modelo oculto de Markov.

## Estados: cgR y atR (representarán dos tipos de zonas que se suceden
## a lo largo del genoma; unas ricas en símbolos C y G, y otras ricas
## en símbolos A y T)

## Ambos estados son equiprobables al inicio (es decir, la primera
## zona de un genoma puede ser, indistintamente, de cualquiera de los
## dos tipos).

## Por otro lado, probabilidad de cambiar de estado es 0.0002
## (probabilidad de que una zona rica en Cs y Gs cambie a una zona
## rica en As y Ts; y viceversa).

## Por último, las probabilidades de los distintos nucleótidos en cada
## uno de los estados (puede observarse como la probabilidad de los
## símbolos C y G es mayor en las zonas ricas en Cs y Gs, idem para el
## otro tipo de zona):

## P(X{j} = A| E{j} = cgR) = 0.2462   P(X{j} = A| E{j} = atR) = 0.2700
## P(X{j} = C| E{j} = cgR) = 0.2476   P(X{j} = C| E{j} = atR) = 0.2084
## P(X{j} = G| E{j} = cgR) = 0.2985   P(X{j} = G| E{j} = atR) = 0.1980
## P(X{j} = T| E{j} = cgR) = 0.2077   P(X{j} = T| E{j} = atR) = 0.3236

# > ricoCG_AT
# $pi
# cgR atR 
# 0.5 0.5 
# $A
#        cgR    atR
# cgR 0.9998 0.0002
# atR 0.0002 0.9998
# $B
#          A      C      G      T
# cgR 0.2462 0.2476 0.2985 0.2077
# atR 0.2700 0.2084 0.1980 0.3236

reg <- c("cgR", "atR") #lista con el nombre de las regiones
simb <- c("A", "C", "G", "T")

piRico <- c( 0.5, 0.5)
names (piRico) <- reg

piRico
ARico <- matrix(c(0.9998, 0.0002, 0.0002, 0.9998), nrow=2,
                dimnames = list(reg,reg))
ARico

BRico <- matrix(c(0.2462, 0.2700, 0.2476, 0.2084, 0.2985, 0.1980, 0.2077, 0.3236)
    , nrow =2, dimnames = list(reg, simb))

BRico

ricoCG_AT <- list(pi = piRico, A = ARico, B = BRico); ricoCG_AT

ricoCG_AT



######################################################################
## (4) Definir una función probSecuenciaEstados que dado un modelo
## oculto de Markov (una lista según lo descrito) y una secesión de
## tipos de región, calcule la probabilidad de la misma.

## Utilizar esta función para calcular la probabilidad, según el
## modelo anterior, de que los tres primeros símbolos de un fragmento
## de ADN pertenezcan a una zona rica en Cs y Gs, y el cuarto no.


probSecuenciaEstados <- function(mom, sEst){
    l <- length(sEst)
    # prob de una region multiplicada por cada una de las probabilidades de 
    # transicion del resto de simbolos
    # funcion prod multiplica (igual que sum suma jaja xd)
    res <- mom[sEst[1]]*prod(mom$A[matrix(c(sEst[-l], sEst[-1]),ncol=2)])
    names(res) <- NULL # se quita el "nombre" que tiene asociada la probabilidad
    # para quedarnos sólo con el valor numérico que es lo que nos interesa
    res
}

secEstados <- c("cgR", "cgR", "cgR", "atR")
probSecuenciaEstados(ricoCG_AT,secEstados)

##ricoCG_AT["cgR"] * ricoCG_AT["cgR","cgR"]^2 *
 ##   ricoCG_AT["cgR", "atR"]


######################################################################
## El algoritmo de Viterbi se define como sigue:

## Entrada: un modelo oculto de Markov y una secuencia
##          de símbolos, s1, ..., sl, 
## Salida: La secuencia de tipos de región más probable para la
##         secuencia de símbolos dada.

## nu{1}[rj] = B[rj,s1] * pi[rj] para 1 <= j <= n
## Para k desde 2 a l:
##   nu{k}[rj] = B[rj,sk] * max{nu{k-1}[ri] * A[ri,rj] | 1 <= i <= n}
##   pr{k}[rj] = argmax{nu{k-1}[ri] * A[ri,rj] | 1 <= i <= n}
## q{l} = argmax{nu{l}[rj] | 1 <= j <= n }
## Para k desde l-1 hasta 1
##   q{k} = pr{k+1}[q{k+1}]
## Devolver la sucesión de tipos de región q{1}, ... q{l}

## La función Viterbi, dado un modelo oculto de Markov y una sucesión
## de símbolos S, calcula la sucesión de tipos de región más probable
## aplicando el algoritmo de viterbi.

viterbi <- function (mom, obs) {
    l <- length(obs) # longitud secuencia
    regiones <- names(mom$pi) 
    ## pr es donde calcula DÓNDE se ha alcanzado el máximo en posicion k
    ## es decir, en qué región se ha alcanzado el máximo de probable
    pr <- matrix(nrow = l-1, ncol = length(regiones),
                 dimnames = list(2:l, regiones))
    nu <- mom$B[,obs[1]] * mom$pi #calculo de nu iniciales
    for (k in 2:l) { # para el resto
        aux <- nu * mom$A
        nu <- mom$B[,obs[k]] * apply(aux, 2, max)
        pr[as.character(k),] <- apply(aux, 2, which.max) # se asocia a cada nu 
                                            #la region mas probable
    }
    q <- c()
    q[l] <- regiones[which.max(nu)]
    for (k in (l-1):1) {
        q[k] <- regiones[pr[as.character(k+1),q[k+1]]]
    }
    q
}

## viterbi(ricoCG_AT, fago[1:500]) # esto daria problemas porque no estamos
## trabajando con la version logarítmica y no seremos capaces de solucionar el problema

######################################################################
## (5) Utilizar la función viterbi para calcular el tipo de región más
## probable para los cuatro primeros símbolos del genoma del
## Enterobacteria phage lambda. Y para los 40 primeros. Y para los 200
## primeros. Y los 400 primeros.

## Resultado: todos de una zona rica en Cs y Gs. El aumento en la
## cantidad de observaciones modifica la aproximación a una zona rica
## en As y Ts.

q4 <- viterbi(ricoCG_AT,fago[1:4]); q4
q40 <- viterbi(ricoCG_AT,fago[1:40]); q40
q200 <- viterbi(ricoCG_AT,fago[1:200]); q200
q400 <- viterbi(ricoCG_AT,fago[1:400]); q400

######################################################################
## (6) Pretendemos utilizar la función viterbi para calcular la
## sucesión de estados más probable para el genoma del Enterobacteria
## phage lambda. ¿Qué problema nos encontraremos para realizar el
## cálculo pedido?

## Longitud de las secuencias: conforme aumenta el tamaño de la secuencia
## las probabilidades van siendo más epqueñas, y llega un punto en el que 
## la precisión del programa hace que el resultado sea  0, cuando esto no 
## debe de ser posible, ya que siempre debemos de considerar la probabilidad
## que obtengamos aunque sea muy pequeña


## Versión con log-probabilidades sería la solución


######################################################################
## (7) Modificar adecuadamente la función viterbi para que nos
## permita realizar el cálculo solicitado en el ejercicio anterior.

## CAMBIOS: 
## Cálculo de logaritmos
## Las multiplicaciones pasan a sumas

viterbiLog <- function (mom, obs) {

    l <- length(obs) # longitud secuencia
    regiones <- names(mom$pi) 
    pr <- matrix(nrow = l-1, ncol = length(regiones),
                 dimnames = list(2:l, regiones))
    B <- log(mom$B) # matriz de prob de los simbolos
    A <- log (mom$A) # matriz de prob de los estados
    pi <- log(mom$pi)
    nu <- B[,obs[1]] + pi #calculo de nu iniciales
    for (k in 2:l) { # para el resto
        aux <- nu * A
        nu <- B[,obs[k]] + apply(aux, 2, max)
        pr[as.character(k),] <- apply(aux, 2, which.max) # se asocia a cada nu 
        #la region mas probable
    }
    q <- c()
    q[l] <- regiones[which.max(nu)]
    for (k in (l-1):1) {
        q[k] <- regiones[pr[as.character(k+1),q[k+1]]]
    }
    q
}



######################################################################
## (8) Superponer al gráfico anterior el tipo de región identificado.

bfl <- local.composition(fago, GC, 500,50)

plot(bfl$positions, bfl$local.results, xlab = "", ylab = "GC",
     main= "NC_001416", type="l",col= "blue", lty = 1, ylim=c(0,1))

atR <-  which(regiones=="atR")
cgR <- which(regiones =="cgR")

valGCatR <- GC(fago[atR])
points(atR, rep(valGCatR, length(atR)), col ="orange", pch=".")
text(0, valGCarR)


######################################################################
