######################################################################
## Biotecnología
## Análisis de variabilidad a lo largo de una secuencia
######################################################################
## Repaso de R
######################################################################
## * Carga de bibliotecas adicionales
## library(nombreBiblioteca) o a través del menú Packages
## * Variables y asignación 
## variable <- valor
## * Definición de funciones
## nombreFunción <- function(datos, ..., arg1 = valor1, ...)
## VECTORES: Todos los elementos de la misma naturaleza
## * Selección:
##   - Por los índices de las posiciones: vector[vectorÍndices]
##   - Por exclusión de las posiciones: vector[-vectorÍndices]
## MATRICES: reorganización de los datos de un vector
## * Selección:
## - Por los índices de las filas y columnas:
## matriz[, índicesColumnas] => columnas completas
## matriz[índicesFilas, ] => filas completas
## matriz[índicesFilas, índicesColumnas] => submatriz
## LISTAS: Cada elemento puede tener una naturaleza diferente.
## * Selección:
## - Por el índice del elemento: lista[[índice]]
## - Por el nombre del elemento (no siempre disponible): lista$nombre
## GRÁFICOS: n.color, n.línea y n.punto son un entero indicando el
## color, el tipo de línea y el tipo de punto, respectivamente (por
## defecto 1 en todos los casos: color negro, línea continua y círculo
## respectivamente).
## * En general: main (título), xlab e ylab (etiquetas de los ejes, x
##   e y respectivamente) y xlim e ylim (dimensiones de los ejes, x e
##   y respectivamente).
## * Puntos:
##   plot(coord.x, coord.y, col = n.color, pch = n.punto)
## * Líneas (concatenación de segmentos entre los puntos dados por las
##   coordenadas): 
##   plot(coord.x, coord.y, type = "l", col = n.color, lty = n.línea)
## * Histograma (de densidad):
##   hist(datos, freq = FALSE)
## * Diagramas de caja y bigote
##   boxplot(datos1, datos2, ..., col = vector.colores,
##           names = vector.nombres)
## Añadidos a un gráfico:
## * Líneas verticales a distancias vector.x del origen
##   abline(v = vector.x, col = vector.color, lty = vector.línea)
## * Líneas horizontales a distancias vector.y del origen
##   abline(h = vector.y, col = vector.color, lty = vector.línea)
## * Líneas
##   lines(x.coord, y.coord, col = n.color, lty = n.línea)
## * Puntos
##   points(x.coord, y.coord, col = n.color, pch = n.punto)
## * Leyenda
##   legend(x.pos, y.pos, vector.nombres, col = vector.colores, 
##          lty = vector.líneas)
######################################################################
## (1) Definir una función, prop.alfabeto, que dados una secuencia de
## símbolos, el alfabeto o parte del mismo, alph, y opcionalmente el
## tamaño de dicha secuencia; calcule la proporción que, de cada uno
## de los símbolos de alf, hay en la secuencia. 


prop.alfabeto <- function (sec, alph, tam = NULL){
    if (is.null(tam) == TRUE){
        tam <- length(sec)
    }
    m <- matrix(nrow = length(alph), ncol = 1)
    colnames(m) <- c("proportion")
    rownames(m) <- alph
    
    for (i in 1:length(alph)){
        c <- 0 ## contador de veces
        for (j in 1:tam){
            if (alph[i] == sec[j]){
                c <- c+1 
            }
        }
        m[i,1] <- c / tam
    }
    return (m)
    
}



#alph <- c(unique(sec1))
#mejor utilizar factor para indicar aquellos aas que nos interesan
#alph <- factor(sec1, levels = c("A","T","G", "C"))

prop.alfabeto(sec1, c("A","T","G", "C"))


######################################################################
## (2) Calcular la proporción que, de cada nucleótido, hay en la
## secuencia sec1, y la que hay en los fragmentos de longitud 30
## inicial y sucesivos con un desplazamiento de 5 símbolos cada vez
## (hasta tres fragmentos).

DNA <- s2c("ACGT")
prop.alfabeto(sec1, DNA)

prop.alfabeto(sec1[1:30], DNA, 30)
prop.alfabeto(sec1[6:35], DNA, 30)
prop.alfabeto(sec1[11:40], DNA, 30)

#####################################################################
## La siguiente función calcula (para una secuencia, seq) la sucesión
## de valores que se obtiene al aplicar una función, func, a sucesivos
## fragmentos de seq de longitud window.length con un desplazamiento
## offset.

local.composition <-
    function (seq, func, window.length, offset, ...) {
        global <- func(seq, ...)
        low.extrems <- seq(1, length(seq) - window.length + 1,
                         by = offset) #hace las divisiones de los fragmentos
        local <-
            mapply(function (i) {
                func(seq[i:(i + window.length - 1)], ...)},
                low.extrems)
        list(local.results = local, global.result = global,
             positions = low.extrems)
    }

######################################################################
## (3) Calcular, combinando las dos funciones anteriores, la
## variabilidad en la proporción de nucleótidos a lo largo de sec1
## (fragmentos de tamaño 30 y valor de desplazamiento 5); y la
## variabilidad en su contenido GC (la cantidad de pares
## Guanina-Citosina presentes en una molécula de ADN).
## Nota: Consultar la ayuda de la función GC de la biblioteca seqinr.

local.composition(sec1, prop.alfabeto, 30, 5, DNA, 30)

local.composition(sec1, GC, 30, 5)


######################################################################
## Haemophilus influenzae (NC_000907) es una bacteria gram-negativa
## que fue erróneamente considerada la causa de la gripe hasta que en
## 1933 se descubrió la causa viral de la misma.

## El genoma de H. influenzae fue el primer genoma en secuenciarse de
## un organismo de vida libre por el grupo de Craig Venter en 1995.
######################################################################

######################################################################
## (4) Descargar su genoma de la base de datos de NCBI (como un
## fichero fasta) y, combinando las funciones read.fasta y
## getSequence de la biblioteca seqinr, obtener el vector genomaHI
## con la sucesión de nucleótidos que lo componen.

genomaHI <- getSequence(read.fasta("NC_000907.fasta",
                                   forceDNAtolower = F))[[1]]
genomaHI

## > length(genomaHI)
## 1830138



######################################################################
## (5) Calcular, para dicho genoma, el contenido GC global y su
## variabilidad a lo largo del mismo utilizando una ventana de
## longitud 20000 (aprox. un 1% del tamaño de la secuencia) y un
## desplazamiento de 2000 (un 10% del fragmento en cada
## ocasión). Representar gráficamente los resultados obtenidos.

dataGC <- local.composition(genomaHI, GC, 20000, 2000); dataGC
plot(dataGC$positions, dataGC$local.results, xlab = "", ylab= "GC",
     main= "NC_000907", type ="l", col="blue", lty=1)
abline(h = dataGC$global.result, col ="red", lty =1)

######################################################################

######################################################################
## El contenido en GC de un genoma, la proporción de los nucleótidos g
## y c, es una característica muy específica de cada organismo en
## concreto.

## La transferencia horizontal de genes hace referencia a la
## transmisión de genes entre organismos vivos (comúnmente de
## diferentes especies) sin que uno sea el progenitor del otro.

## Debido a que el contenido en GC es una de las características más
## específicas de cada especie, la identificación de zonas en el
## genoma cuyo contenido en GC diverge significativamente del
## contenido GC global puede indicar la existencia de un evento de
## transferencia horizontal de genes.

## ¿Se observa en el gráfico anterior una desviación significativa?
## ¿Puede esto indicar transferencia horizontal de genes?
######################################################################

######################################################################
## (6) Identificar los puntos en los que se observa la mayor
## desviación.


abline(v=1550000, col = "green")
abline(v=1600000, col = "green")



######################################################################

######################################################################
## (7) Comprobar, visualizando los datos y con un test de
## Shapiro-Wilk, si el contenido GC local en los distintos tramos
## sigue una distribución normal.

#posiciones donde esta el pico nose de donde lo saka xd
frag <- dataGC$local.results[766:806]
resto <- dataGC$local.results[-(766:806)]


shapiroFrag <- shapiro.test(frag)
shapiroResto <- shapiro.test(resto)

(1- shapiroFrag$p.value)*100
(1 - shapiroResto$p.value) *100

# Aunque para el fragmento el resultado no muy significativo( por
# poco, 99.8%) para eñ resto si, por lo que podemos afirmar que no sigue
# una distribucion normal
######################################################################

######################################################################
## (8) Determinar, utilizando un test adecuado, si el contenido GC
## local en el fragmento tiene una distribución con una media
## significativamente distinta (resp. mayor) que en el resto.


# usas el wilcox.test (igual que en 7)

## En ambos casos el resultado es muy significativo

######################################################################

######################################################################
## (9) Visualizar, utilizando diagramas de caja y bigote, el resultado
## obtenido en el apartado anterior. 

## - Un asterisco p-valor < 0.05
## - Dos asteriscos p-valor < 0.01
## - Tres asteriscos p-valor < 0.001

boxplot(frag, resto, main ="NC_000907", col = c("green", "blue"),
        ylim = c(0.33,0.49), names = c("Fragmento", "Resto"))
text(1, 0.48, "***", cex = 2)

# Un asterisco pvalor < 0.05
# Dos asteriscos pvalor < 0.01
# Tres asteriscos pvalor < 0.001

######################################################################

######################################################################
## (10) Realizar una análisis parecido para el genoma de
## Methanocaldococcus jannaschii (NC_000909). Es una archaea
## metanógena termofílica que habita ventanas hidrotermales; crece
## usando como fuente de energía dióxido de carbono e hidrógeno y
## produciendo metano como producto secundario de su metabolismo. El
## grupo de Craig Venter fue el primero en secuenciar su genoma en
## 1996; constituyó el primer genoma de archaea en secuenciarse
## completamente. La secuenciación de su genoma produjo evidencias
## claves para la existencia de los tres dominios de la vida (Archaea,
## Bacteria y Eukarya).





######################################################################
