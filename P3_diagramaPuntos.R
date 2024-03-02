######################################################################
## Biotecnología
## Alineación de secuencias
######################################################################
## La construcción de un diagrama de puntos es el método más simple
## para visualizar las posibles coincidencias que pueden darse entre
## dos secuencias sobre un mismo alfabeto.

## El método consiste en dibujar un punto por cada par de puntos de
## cada una de las secuencias que coincida.

## El fichero P3_diagramaPuntos.pdf contiene una imagen de los
## gráficos esperados para los dos distintos ejemplos considerados. El
## fichero P3_secuencias.RData las distintas secuencias utilizadas
## para generar los mismos, y los correspondientes gráficos.

## El objetivo del ejercicio es definir una función diagramaPuntos, lo
## más genérica posible, que permita realizar los gráficos ejemplo1,
## ..., ejemplo5 a partir de las secuencias correspondientes (y los
## distintos parámetros que se quieran personalizar).
## ¡Ojo! El último gráfico tarda bastante en visualizarse.

######################################################################
## Ejemplos: 
######################################################################

## ejemplo1 <-
##     diagramaPuntos(c("U","N","A","C","A","B","R","A","L","O","C","A"),
##                    c("A","B","R","A","C","A","D","A","B","R","A"),
##                    ejes = TRUE, tam = 7)

## ejemplo2 <-
##     diagramaPuntos(sec1, sec2,
##                    "Primera secuencia", "Segunda Secuencia",
##                    título = "Coincidencias", ejes = TRUE, tam = 6)

## De la base de datos de proteínas, Uniprot, utilizando los ID
## O18381 y Q66SS1 obtenemos dm y hs, a partir de los cuales
## generamos el gráfico: Gen Eyeless de la Drosophila Melanogaster
## versus Gen Box 6 del Homo Sapiens.

## ejemplo3 <-
##     diagramaPuntos(dm, hs,
##                    "Drosophila Melanogaster", "Homo Sapiens",
##                    título = "Gen eyeless DM vs Box 6 HS",
##                    ejes = FALSE, tam = 3)

######################################################################
## (1) Prueba a descargar los correspondientes ficheros .fasta
## (ficheros de texto con un formato estándar con la secuenciación de
## una molécula: proteína, mRNA, DNA, ...). Utiliza la función
## read.fasta de la biblioteca seqinr para obtener una lista con las
## distintas secuencias contenidas en dichos ficheros. Cada elemento
## de dicha lista es una estructura SeqFasta, entre cuyos elementos se
## encuentra un vector con la correspondiente secuenciación; además de
## variada información extraída de la correspondiente anotación que
## acompaña a cada secuencia. Utiliza la función getSequence para
## obtener una lista de vectores, los caracteres que componen cada una
## de las secuencias contenidas en el fichero.

## Nota: Al contener cada fichero .fasta la secuenciación de una única
## molécula, la lista obtenida al utilizar tanto la función read.fasta
## como la función getSequence contendrá un único elemento. Utiliza la
## selección en listas para obtener un vector con la correspondiente
## sucesión de caracteres.
######################################################################
library(seqinr)

secuencias_dm <- read.fasta("dm.fasta",forceDNAtolower = FALSE)
dm <- getSequence(secuencias_dm[[1]])
length(unique(dm))

secuencias_hs <- read.fasta("hs.fasta", forceDNAtolower = FALSE)
hs <- getSequence(secuencias_hs[[1]])



######################################################################

######################################################################
## Más ejemplos: 
######################################################################

## En este caso, de la base de datos de proteínas del NCBI, utilizando
## los ID NP_000511 y BAC38018 obtenemos hH y mH, a partir de los
## cuales realizamos la comparativa gráfica de la composición de la
## hexosaminidase subunit alpha (HEXA) vs hexosaminidase A.

## ejemplo4 <-
##     diagramaPuntos(hH, mH, "Homo Sapiens", "Mus musculus",
##                    título = "Coincidencias en la proteína hexosaminidase A",
##                    ejes = F, tam = 3)

## Por último, descargando la composición del correspondiente gen para
## ejemplo anterior.

## ejemplo5 <-
##     diagramaPuntos(hNH, mNH, "Homo Sapiens", "Mus musculus",
##                    título = "Coincidencias en el mRNA
## (hexosaminidase A)", ejes = F, tam = .3)

######################################################################
## (2) Calcular las posiciones de las letras coincidentes para las
## secuencias "UNACABRALOCA" y "ABRACADABRA".


s1 <- s2c("UNACABRALOCA")
s2 <- s2c("ABRACADABRA")

# creamos unos vectores en los que almacenaremos, respectivamente, las posiciones donde coinciden
# las letras de s1 en s2 y viceversa
r1 <- c()
r2 <- c()

for (i in 1:length(s1)){
  for(j in 1:length(s2)){
    if(s1[i] == s2[j]){
      r1<- c(r1,i)
      r2 <- c(r2,j)
    }
  }
}

r1
r2

tabla <- data.frame(X=r1, Y = r2); tabla


## Suponiendo que la primera es la que se utilizaría para el eje X una
## solución sería:

##      X    Y  
##  1   3    1  
##  2   3    4  
##  3   3    6  
##  4   3    8  
##  5   3   11  
##  6   4    5  
##  7   5    1  
##  8   5    4  
##  9   5    6  
## 10   5    8  
## 11   5   11  
## 12   6    2  
## 13   6    9  
## 14   7    3  
## 15   7   10  
## 16   8    1  
## 17   8    4  
## 18   8    6  
## 19   8    8  
## 20   8   11  
## 21  11    5  
## 22  12    1  
## 23  12    4  
## 24  12    6  
## 25  12    8  
## 26  12   11  



######################################################################
## (3) Crear un diagrama de puntos que refleje las coincidencias
## anteriores (ver ejemplo1).


plot(data.frame(x=r1,y = r2), xlim=c(0,length(s1)),ylim=c(0,length(s2)),
     xlab = "", ylab= "",
     pch = ".", cex= 5, col = "blue", xaxt = "n", yaxt="n",
     main= "Matriz de puntos")
axis(side= 1, at = 1:length(s1), labels=s1)
axis(side= 2, at = 1:length(s2), labels=s2)

######################################################################
## (4) Definir la función diagramaPuntos que, dadas dos secuencias
## sobre un mismo alfabeto (y el resto de argumentos que se consideren
## necesario), permita realizar el diagrama de puntos correspondiente.

library("lattice")

diagramaPuntosClase <- function(s1, s2, ny ="Y", nx = "X", titulo = 
                             "Matriz de Puntos", ejes= F, tam = 1.2){
  r1 <- c()
  r2 <- c()
  
  for (i in 1:length(s1)){
    for(j in 1:length(s2)){
      if(s1[i] == s2[j]){
        r1<- c(r1,i)
        r2 <- c(r2,j)
      }
    }
  }
  
  plot(data.frame(x=r1,y = r2), xlim=c(0,length(s1)),ylim=c(0,length(s2)),
       xlab = "", ylab= "",
       pch = ".", cex= 5, col = "blue", xaxt = "n", yaxt="n",
       main= "Matriz de puntos")
  axis(side= 1, at = 1:length(s1), labels=s1)
  axis(side= 2, at = 1:length(s2), labels=s2)
  
  
  
}
  

###### de mi huerta de aqui pa bajo #######
diagramaPuntos <- function (s1, s2, ...){
  r1 <- c()
  r2 <- c()
  
  for (i in 1:length(s1)){
    for(j in 1:length(s2)){
      if(s1[i] == s2[j]){
        r1<- c(r1,i)
        r2 <- c(r2,j)
      }
    }
  }
  
  plot(data.frame(x=r1,y = r2), xlim=c(1,length(s1)),ylim=c(1,length(s2)),
       xlabel = s1, ylabel=s2, ...)
}


######################################################################


