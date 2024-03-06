
## IMPORTACIÓN DE LIBRERÍAS ##

library(seqinr)
library(Biostrings)


######################################################################
## Apartado 1: 2 ptos
######################################################################
## Definir una función NeedlemanWunschTODAS que devuelva todas las
## alineaciones globales óptimas que podrían obtenerse al aplicar el
## algoritmo de Needleman-Wunsch.

## Nota: Puede utilizarse la función anterior como referencia.

## Con los datos del ejemplo anterior, las siguientes serían todas las
## alineaciones globales óptimas que deberían obtenerse.

## AAAC
## AG-C

## AAAC
## -AGC

## AAAC
## A-GC

## Nota: Aunque no es el caso con las secuencias y la función de
## puntuación propuestas como ejemplo, debe tenerse en cuenta que en
## general no todas las alineaciones globales óptimas tienen que
## coincidir en tamaño.

## Solución:

## NO ESTÁ TERMINADO. He intentado plantear cómo podría ser una posible solución
# pero no sabía cómo continuar o cómo enfocarlo de otra forma. 

# He creado una función auxiliar buscaOps que devuelve los índices, en caso de  
# haber repeticiones, o el índice, si no las hay, correpondiente con la/s 
# puntuación/es máximas. Esta función servirá de forma interna en la función
# que se pide en este apartado

buscaOps <- function(ops, maxOps){
  # Recorremos las tres posibles opciones de puntuación del vector ops
  for (i in 1: length(ops)){
    rep <- duplicated(ops) # vector de valores booleanos que devuelve TRUE si el
    # valor de esa posición ya ha aparecido antes
    k <- 1 # iterador
    c <- 0 # contador de veces que aparece TRUE
    while (k <= length(rep)){
      if (rep[k] == TRUE){
        c <- c+1
      }
      k <- k+1
    }
    
    # Si sólo hay un valor repetido, devolveremos un vector con dos índices
    if (c == 1){
      # Comprobamos si maxOps en la posición i de rep es TRUE y devolvemos las
      # posiciones que le corresponden a este valor repetido
      if (maxOps == ops[i] && rep[i]== TRUE){
        if(i < length(ops)){
          if(ops[i] == ops[i+1]){
            a <- i; b <- i+1
          } else{
            a <- i; b <- i-1
          }
        } else{
          if(ops[i] == ops[i-1]){
            a <- i; b <- i-1
          } else{
            a <- i; b <- i-2
          }
        }
        return (c(a,b))
      }
    } else if (c > 1) {
      # Si aparece más de una vez, todos los valores posibles son iguales, por
      # tanto, devolvemos todas las posiciones como posibles caminos para formar
      # alineaciones.
      return(c(1,2,3))
    }
  }
  # Si no hay repetidos, devolvemos la posición que le corresponde al máximo
  return(c(maxOps))
  
} 

NeedlemanWunshTODAS <- function (ref, comp, mS, gap) {
  alf <- colnames(mS)
  if (all(unique(ref) %in% alf) &
      all(unique(comp) %in% alf)) {
    lr <- length(ref) ; lc <- length(comp)
    mVal <- matrix(ncol = lr + 1, nrow = lc + 1)
    # Creamos una matriz de rastreo para cada posible opción
    mRas1 <- matrix(ncol = lr + 1, nrow = lc + 1) # caso diagonal (match)
    mRas2 <- matrix(ncol = lr + 1, nrow = lc + 1) # caso arriba (insertion)
    mRas3 <- matrix(ncol = lr + 1, nrow = lc + 1) # caso derecha (deletion)
    # Inicialización de todas las matrices
    mVal[1, 1] <- 0 ; mRas1[1, 1] <- mRas2 <- mRas3 <- 0
    mVal[1, -1] <- cumsum(rep(gap, lr))
    mRas1[1, -1] <- mRas2[1, -1] <- mRas3[1, -1] <- 2
    mVal[-1, 1] <- cumsum(rep(gap, lc))
    mRas1[-1, 1] <- mRas2[-1, 1] <- mRas3[-1, 1] <- 3
    for (i in 2:(lc+1)) {
      for (j in 2:(lr+1)) {
        ops <- c(mVal[i - 1, j - 1] + mS[comp[i - 1], ref[j - 1]],
                 mVal[i, j - 1] + gap, mVal[i - 1, j] + gap)
        mVal[i, j] <- max(ops)
        
        # Buscamos si hay duplicados; 
        bRep <- buscaRep(ops, max(ops))
        if (!is.null(bRep)){
          # tam_mRas almacenará los índices de opciones para las posibles
          # alineaciones
          tam_mRas <- bRep
          # Por cada repetición, tendremos una matriz de rastreo diferente, por
          # lo que le asignaremos a cada una su correspondiente índice en caso
          # de que se encuentre entre los repetidos
          for (k in 1: length(tam_mRas)){
            if(bRep[k] == 1){
              mRas1[i,j] <- 1
            } else if (bRep[k] == 2){
              mRas2[i,j] <- 2
            } else{
              mRas3[i,j] <- 3
            }
          }
        } else{
          mRas1[i,j] <- mRas2[i,j] <- mRas3[i,j] <- NULL
          mRas <- c(mRas1[i,j],mRas2[i,j],mRas3[i,j])
          for (n in 1:length(mRas)){
            if(which.max(ops) == n){
              mRas[n] <- which.max(ops)
            }
          }
        }
      }}
    alin <- matrix(nrow = 2) ; i <- lc + 1 ; j <- lr + 1
    opsi <- c(-1, 0, -1) ; opsj <- c(-1, -1, 0)
    
    # Comprobamos si mRas contiene un vector de valores o un único valor
    if (length(mRas[i,j]) > 1){
      todasAlin <- c() # lista que almacenará cada una de las alineaciones
      # Para cada valor que contenga el vector de mRas[i,j], creamos una 
      # alineación diferente
      for (k in 1:length(mRas[i,j])){
        sel <- mRas[i, j][k] ; nRef <- c(NA, ref); nComp <- c(NA, comp)
        while (sel != 0) {
          ops <- matrix(c(nRef[j], nComp[i], nRef[j],
                          "-", "-", nComp[i]),
                        byrow = T, ncol = 2)
          alin <- cbind(ops[sel,], alin)
          i <- i + opsi[sel] ; j <- j + opsj[sel]
          sel <- mRas[i, j][k]
        }
        totalAlin[[k]] <- alin[,-ncol(alin)]
      }
      
      ## de aquí para abajo, no se ha modificado nada del que se proporciona ##
      
      sel <- mRas[i, j] ; nRef <- c(NA, ref); nComp <- c(NA, comp)
      while (sel != 0) {
        ops <- matrix(c(nRef[j], nComp[i], nRef[j],
                        "-", "-", nComp[i]),
                      byrow = T, ncol = 2)
        alin <- cbind(ops[sel,], alin)
        i <- i + opsi[sel] ; j <- j + opsj[sel]
        sel <- mRas[i, j]
      }
      list(valor = mVal[lc + 1, lr + 1],
           alineacion = alin[,-ncol(alin)])
    } else {
      stop("Alfabeto(s) incorrectos")        
    }}
  
}

##NeedlemanWunshTODAS(s2c("AAAC"), s2c("AGC"), punt, -2)

## Ejemplo de uso y posible salida obtenida:
## > NeedlemanWunshTODAS(s2c("AAAC"), s2c("AGC"), punt, -2)
## $valor
## [1] -1
## $alineaciones
## $alineaciones[[1]]
##      [,1] [,2] [,3] [,4]
## [1,] "A"  "A"  "A"  "C" 
## [2,] "A"  "G"  "-"  "C" 
## $alineaciones[[2]]
##      [,1] [,2] [,3] [,4]
## [1,] "A"  "A"  "A"  "C" 
## [2,] "A"  "-"  "G"  "C" 
## $alineaciones[[3]]
##      [,1] [,2] [,3] [,4]
## [1,] "A"  "A"  "A"  "C" 
## [2,] "-"  "A"  "G"  "C" 

######################################################################
## Apartado 2: 0'8 ptos
######################################################################
## Obtener, a partir de los códigos de acceso de la base de datos NCBI
## proporcionados, la secuenciación de distintas macromoléculas; y
## definir una función de puntuación. Para ello:
######################################################################
## (2.1) Obtener un fichero en formato FASTA con las secuencias
## codificantes del genoma del organismo asignado y definir, a partir
## de su contenido, la lista de vectores de aminoácidos
## correspondientes a cada una de dichas secuencias.


## Solución:

sec_genome <- read.fasta("NZ_CP058907_genome.fasta" ,forceDNAtolower = FALSE)
## Cantidad de secuencias descargadas
length(sec_genome)

NZ_CP058907 <- getSequence(sec_genome[[1]])
## Tamaño de la única secuencia que contiene el genoma
lnz <- length(NZ_CP058907); lnz

## Vector de aminoácidos de la secuencia del genoma
symb <- c(NZ_CP058907[1:lnz])
## Primeros 25 aminoácidos
symb[1:25]




## Ejemplo de datos obtenidos (la información dependerá del organismo
## asignado a cada grupo)
## Cantidad de secuencias descargadas
## > lcs <- length(cs_organism); lcs
## 3906
## Tamaño de la primera de ellas
## > length(cs_organism[[1]])
## 507
## Tamaño de la última
## > length(cs_organism[[lcs]])
## 47
## Primeros 25 aminoácidos de la proteína correspondiente a la
## secuencia 1500.
## > cs_organism[[1500]][1:25]
##  [1] "M" "P" "S" "G" "E" "P" "S" "T" "A" "G" "H" "F" "E" "H" "L" "P" "R" "G" "S"
## [20] "F" "G" "R" "I" "L" "S"

######################################################################
## (2.2) Obtener un fichero en formato FASTA con la secuencia de la
## proteína, producida por el organismo anterior, asignada y definir,
## a partir de su contenido, el vector de aminoácidos correspondiente.

## Solución:

sec_protein <- read.fasta("WP_011157096_protein.fasta" ,forceDNAtolower = FALSE)
## Cantidad de secuencias descargadas
length(sec_protein)

WP_011157096 <- getSequence(sec_protein[[1]])
## Tamaño de la única secuencia que contiene la proteína
lwp <- length(WP_011157096); lwp

## Vector de aminoácidos de la secuencia de la proteína
symb_protein <- c(WP_011157096[1:lwp])
## Primeros diez aminoácidos
symb_protein[1:10]



## Ejemplo de datos obtenidos (la información dependerá de la proteína
## asignada a cada grupo).
## Tamaño de la proteína
## > length(protein_organism)
## 128
## Primeros 10 aminoácidos
## > protein_organism[1:10]
## "M" "P" "K" "S" "F" "Y" "D" "A" "V" "G"

######################################################################
## (2.3) Definir una función, Punt_B62, que dados dos aminoácidos
## cualesquiera devuelva el valor asociado a los mismos según la matriz
## de sustitución BLOSUM62.

#install.packages("Biostrings")
#library("Biostrings") ## La matriz BLOSUM62 forma parte de los datos
## disponibles en la biblioteca Biostrings

data(BLOSUM62)

## Solución:

Punt_B62 <- function(a1, a2){
  
  b62 <- BLOSUM62
  
  # La opción de dos huecos superpuestos no debe considerse porque estamos
  # estudiando la alineación entre un par de secuencias y nunca se va a dar 
  # este caso.
  
  if (a1 == "-"){
    a1 <- "*"
  } else if (a2 == "-"){
    a2 <- "*"
  }
  
  return (b62[a1, a2])
}

## Ejemplos de uso:
## > Punt_B62("V", "I") 
## 3
## > Punt_B62("R", "A")
## -1
## > Punt_B62("-", "W")
## -4

Punt_B62("V", "I") 
Punt_B62("R", "A")
Punt_B62("-", "W")

######################################################################
## Apartado 3: 2 ptos
######################################################################
## Localizar entre las secuencias codificantes del organismo asignado
## cuál se corresponde con la proteína descargada. Para ello:
######################################################################
## (3.1) Definir una función, NeedlemanWunschFUNC, que dadas dos
## secuencias de aminoácidos y una función de puntuación devuelva la
## valoración de una alineación global óptima, obtenida al aplicar el
## algoritmo Neddleman-Wunsch, entre las dos secuencias.

## Solución:

NeedlemanWunschFUNC <- function(s1, s2, func){
  ## Obtenemos la matriz de puntuacion para poder introducirla como mS en
  ## la función NeedlemanWunsch
  alph <- unique(c(s1, s2))
  
  # Creación de la matriz vacía que representará la matriz de sustitución
  mS <- matrix(data=NA, ncol = length(alph), nrow = length(alph)) 
  colnames(mS) <- rownames(mS) <- alph
  
  for (i in 1:nrow(mS)){
    for (j in 1:ncol(mS)){
      mS[i,j] <- func(alph[i],alph[j])
    }
  }
  gap <- func("-",s1[[2]])  # puntuación del hueco para un símbolo cualquiera
  
  return(NeedlemanWunsh(s1, s2, mS, gap)[[1]] )
  
}


## Ejemplos de uso:
## > NeedlemanWunschFUNC(c("E", "G", "H", "I", "L"), 
##                       c("E", "L", "G", "H", "I"),
##                       function (x, y) { if (x == y) {2} else {-2}})
## 4
## > NeedlemanWunschFUNC(protein_organism, cs_organism[[1]], Punt_B62)
## -1063

NeedlemanWunschFUNC(c("E", "G", "H", "I", "L"), 
                    c("E", "L", "G", "H", "I"),
                    function (x, y) { if (x == y) {2} else {-2}})

NeedlemanWunschFUNC(protein_organism, cs_organism[[1]], Punt_B62)

######################################################################
## (3.2) Definir una función, localiza, que dadas una secuencia de
## aminoácidos de referencia y una lista de secuencias de aminoácidos
## localice en la lista aquella con la que, al comparar la de
## referencia con ella utilizando la función anterior (con la función
## de puntuación Punt_B62), se obtenga el mejor resultado. La función
## debe devolver un diagrama de puntos que permita visualizar el grado
## de correspondencia entre las mismas. Pueden incluirse argumentos
## adicionales si se estiman necesarios.

## Nota: El proceso puede tardar bastantes minutos dependiendo del
## organismo asignado. Se recomienda comprobar la corrección de la
## solución con algunos datos de prueba preparados con tal finalidad
## (que deben entregarse).

## Solución:


localiza <-  function(ref, lsec, título = NULL, nx = NULL){
  punt <- c() # almacenará las puntuaciones con cada secuencia de la lista
  for (i in 1:length(lsec)){ # recorremos la lista de secuencias
    punt.i <- NeedlemanWunschFUNC(ref, lsec[[i]], Punt_B62)
    punt <- c(punt, punt.i) # actualizamos el vector
  }
  k <- which.max(punt) # posición de la puntuación más alta
  comp <- lsec[[k]] # guardamos la secuencia que está en la posición k
  
  # Diagrama de puntos
  
  # r1 y r2: vectores en los que se almacenan, respectivamente, las posiciones 
  # donde coinciden las letras de ref en comp y viceversa
  r1 <- c()
  r2 <- c()
  lr <- length(ref); lc <- length(comp)
  for (i in 1:lr){
    for(j in 1:lc){
      if(ref[i] == comp[j]){
        r1<- c(r1,i)
        r2 <- c(r2,j)
      }
    }
  }
  
  if (is.null(título) == T){
    título <- "Diagrama de Puntos"
  } else if (is.null(nx) == T){
    nx <- "Secuencia de referencia"
  }
  
  # Para que muestre el número de la secuencia elegida utilizamos sprintf
  yn <- sprintf("Secuencia comparada %d", k)
  
  plot(data.frame(x=r1,y = r2), xlim=c(0,lr),ylim=c(0,lc),
       xlab = nx, ylab= yn,
       pch = ".", cex= 5, col = "blue", xaxt = "n", yaxt="n",
       main = título)
  axis(side= 1, at = 1:length(ref), labels=ref)
  axis(side= 2, at = 1:length(comp), labels=comp)
  
}


## Ejemplo de uso con dos argumentos adicionales, título para el
## diagrama de puntos y título para el eje X (usando la biblioteca de
## gráficos lattice). Obsérvese que el índice de la secuencia
## localizada en la lista aparece como parte del título del eje Y.
## > localiza(protein_organism, cs_organism,
##            título = "Mycobacterium tuberculosis H37Rv",
##            nx = "Hemoglobin GlbO")
## Resultado esperado: diagramaPuntos.png

localiza(protein_organism, cs_organism, 
         título = "Mycobacterium tuberculosis H37Rv",
         nx = "Hemoglobin GlbO")


######################################################################
## Apartado 4: 3'7 ptos
######################################################################
## Analizar el parecido entre la proteína asignada y un conjunto de 10
## proteínas similares con un grado de coincidencia entre sí de más de
## un 50%. Los códigos de acceso de cada una de las proteínas del
## conjunto son de la base de datos UniProt.
######################################################################
## (4.1) Unificar en un único fichero con formato FASTA, a partir de
## los códigos de acceso del conjunto de 10 proteínas asignado, la
## secuenciación de las mismas y a partir de éste obtener la lista de
## vectores de aminoácidos correspondientes a cada una de dichas
## secuencias.

## Solución:

proteinas10 <- read.fasta("10proteinas.fasta",forceDNAtolower = FALSE)
p10_sec <- getSequence(proteinas10)
# Primeras 6 secuencias del conjunto de 10 proteínas
head(p10_sec); length(p10_sec)


######################################################################
## (4.2) Definir una función, SmithWaterman, que dadas dos secuencias
## de aminoácidos y una función de puntuación devuelva la valoración
## de una alineación local óptima, obtenida al aplicar el algoritmo
## Smith-Waterman, entre las dos secuencias; y que, opcionalmente,
## muestre dicha alineación en pantalla.

SmithWaterman <- function(ref, comp, func, alin = F){
  
  alfRef <- c(unique(ref))
  alfComp <- c(unique(comp))
  alf <- c(unique(c(alfRef,alfComp)))
  lr <- length(ref) ; lc <- length(comp)
  mVal <- matrix(ncol = lr + 1, nrow = lc + 1)
  mRas <- matrix(ncol = lr + 1, nrow = lc + 1)
  
  # Tomamos un símbolo cualquiera para conocer el valor de penalización que 
  # tiene con cualquier símbolo del alfabeto (gap)
  gap <- func("-",ref[[2]]) 
  
  # La primera columna y la primera fila de mVal se completan con 0
  # Utilizamos el 4 como representación del asterisco en mRas
  mVal[, 1] <- mVal[1, ] <- 0 ; mRas[, 1] <- mRas[1,] <- 4
  
  
  # i recorre los símbolos de la secuencia que se compara (filas)
  # j recorre los símbolos de la secuencia de referencia (columnas)
  for (i in 2:(lc +1)){ 
    for (j in 2:(lr +1)){ 
      ops <- c(mVal[i - 1, j - 1] + func(ref[[j-1]], comp[[i-1]]),
               mVal[i, j - 1] + gap, mVal[i - 1, j] + gap, 0)
      mVal[i, j] <- max(ops)
      mRas[i, j] <- which.max(ops) # asignamos la posición de la puntuación que
      # tiene la puntuación mayor de las opciones
    }
  }
  
  al <- matrix(nrow = 2) 
  maxVal <- 0
  # Buscamos las posiciones donde está la puntuación más alta de mVal
  for (i in 2:(lc+1)){
    for (j in 2:(lr+1)){
      if(maxVal < mVal[i,j]){
        maxVal <- mVal[i,j]
        maxVali <- i
        maxValj <- j
      }
    }
  }
  
  i <- maxVali; j <- maxValj # indices de la posición con la puntuación máxima
  opsi <- c(-1, 0, -1) ; opsj <- c(-1, -1, 0)
  sel <- mRas[i, j]; nRef <- c(NA, ref); nComp <- c(NA, comp)
  
  while (sel < 4){
    ops <- matrix(c(nRef[j], nComp[i], nRef[j],
                    "-", "-", nComp[i],
                    nRef[j], nComp[i]),
                  byrow = T, ncol = 2)
    al <- cbind(ops[sel,], al)
    i <- i + opsi[sel] ; j <- j + opsj[sel]
    lastSymb <- c(i,j)
    sel <- mRas[i, j]
  }
  
  # En caso de que queramos que devuelva la alineación
  if(alin == T){
    alinRef <- c(al[1, -ncol(al)])
    alinComp <- c(al[2,-ncol(al)])
    
    # Con seq conseguimos los índices para segmentar la alineación final obtenida
    low.extrems <- seq(1, length(alinRef), 60) 
    
    # Obtenemos los fragmentos para la secuencia de referencia
    linesRef <- c()
    for (i in 1:length(low.extrems)){
      if (i == length(low.extrems)){
        sRef <- c2s(alinRef[low.extrems[[i]]: length(alinRef)])
      } else{
        sRef <- c2s(alinRef[low.extrems[[i]]: low.extrems[[i+1]]-1])
      }
      linesRef[[i]] <- sRef # en la posición i del vector incluimos cada uno de  
      # los fragmentos correspondientes
    }
    
    # Lo mismo para la secuencia de comparación
    linesComp <- c()
    for (j in 1:length(low.extrems)){
      if (j == length(low.extrems)){
        sComp <- c2s(alinComp[low.extrems[[j]]: length(alinComp)])
      } else{
        sComp <- c2s(alinComp[low.extrems[[j]]: low.extrems[[j+1]]-1])
      }
      linesComp[[j]] <- sComp
    }
    
    print(sprintf("Query %d", lastSymb[[2]]))
    print(sprintf("Sbjct %d", lastSymb[[1]]))
    print("   ")
    
    
    # Cálculo de las posiciones a las que llega cada uno de los fragmentos de la 
    # alineación local obtenida.
    
    # Los contadores comienzan a contar por la posición anterior a la 
    # correspondiente con el primer símbolo que tiene la alineación local con 
    # respecto a la secuencia original. De esta forma, el primer símbolo que
    # encuentre en la alineación local coincidirá realmente con su posición
    # en la secuencia original.
    
    # Para los fragmentos de la secuencia de referencia
    n <- c() # vector que almacena las posiciones que queremos
    a <- lastSymb[[2]] -1 # contador
    for (m in 1:length(linesRef)){
      s1 <- s2c(linesRef[[m]])
      for (p in 1:length(s1)){
        if (s1[[p]] != "-"){ # contamos todos los símbolos menos los huecos 
          # ya que no pertenecen a la secuencia original
          a <- a+1 # actualización del contador
        }
      }
      n[[m]] <- a
    }
    
    # Lo mismo para la secuencia de comparación
    q <- c()
    b <- lastSymb[[1]] - 1 # contador
    for(m in 1: length(linesComp)){
      s2 <- s2c(linesComp[[m]])
      for (p in 1:length(s2)){
        if (s2[[p]] != "-"){
          b <- b+1
        }
        q[[m]] <- b
      }
      
    }
    
    
    # Imprimimos cada par de secuencias juntas junto con sus posiciones numéricas
    
    # Como la longitud del vector linesRef es la misma que la de las posiciones,
    # (porque saldrá un índice por cada línea) podemos utilizar el mismo índice
    # para recorrer el vector que contiene las posiciones. Para la línea
    # correspondiente con la secuencia de referencia, el vector n, y para la 
    # secuencia de comparación, el vector q.
    
    for (i in 1: length(linesRef)){ 
      print(sprintf(paste(linesRef[[i]],"   %d"), n[[i]]))
      print(sprintf(paste(linesComp[[i]],"   %d
                           "), q[[i]]))
    }
  }
  
  return(maxVal)
}


SmithWaterman(protein_organism, protein_set[[1]],
              Punt_B62, alin = T)
## Ejemplo de salida esperada: 
## > SmithWaterman(protein_organism, protein_set[[1]],
##                 Punt_B62, alin = T)
## Query 2
## Sbjct 12

## PKSFYDAVGGAKTFDAIVSRFYAQVAEDEVLRRVYPEDDLAGAEERLRMFLEQYWGGPRT   61 
## P-SFYEQVGGHDTFRRLVDAFYRGVAADPVLKPMYPEEDLGPAAERLTLFLEQYWGGPGT   70 

## YSEQRGHPRLRMRHAPFRISLIERDAWLRCMHTAVASIDSETLDDEHRRELLDYLEMAAH   121 
## YSAERGHPRLRMRHMPFRVNPDARDRWLTHMRAAVDELGLPPLQEE---TLWSYLERAAF   127 

## SLVNSPF   128 
## AMVNT-F   133 

## [1] 336

## Descripción de la salida: Por defecto, se muestran 60 caracteres
## por línea. Los números iniciales indican la posición en cada de las
## secuencias originales del símbolo inicial de cada uno de los
## fragmentos con los que se construye la alineación local óptima
## obtenida. Y los números que aparecen al final de cada fila, para
## cada tramo de 60 de la alineación, indican la posición que ocuparía
## en cada una de las secuencias originales el último símbolo que
## aparece en el correspondiente tramo representado del fragmento con
## el que se construye la alineación.
## Esto es, de la secuencia de referencia (el vector protein_organism)
## el fragmento que forma parte de la alineación local óptima obtenida
## sería el formado por los aminoácidos que van desde la posición 2 a
## la posición 128 de dicha secuencia. Puesto que son más de 60
## símbolo, aparece en tres tramos: el primero sería el
## correspondiente a las posiciones que van desde la posición 2 a la
## 61, el segundo desde la posición 62 hasta la 121 y el tercero desde
## la 122 a la 128. Puede observarse que no se ha incluido ningún
## hueco en este prim er fragmento.
## En cuanto a la secuencia con la que se compara (el
## primer vector de la lista protein_set). El fragmento que forma
## parte de la alineación sería el formado por los aminoácidos que van
## desde la posición 12 hasta la 133 (nótese que es de menor tamaño
## que el fragmento de la secuencia de referencia). Su primer tramo
## estaría compuesto por los aminoácidos entre las posiciones 12 y 70
## (nótese que sólo son 59 aminoácidos, ya que se ha incluido un hueco
## al inicio), el segundo tramo desde la posición 71 a la 127 y el
## tercero desde la 128 hasta la 133.
## Por último, además de mostrarse en la pantalla la alineación
## (dividida en distintos tramos tal y como se ha descrito para
## facilitar su lectura), se obtiene como resultado la valoración de
## dicha alineación.

######################################################################
## (4.3) Calcula las valoraciones de alineaciones locales óptimas
## entre la proteína de referencia (considerar como tal la obtenida en
## el apartado 2.2) y cada una de las proteínas del conjunto de 10
## asignado.

## Solución:

local_alin <- c()
for (i in 1: length(p10_sec)){
  local_alin <- c(local_alin, SmithWaterman(WP_011157096, p10_sec[[i]], 
                                            Punt_B62))
}
# Primeras 6 valoraciones obtenidas con la alineación local
head(local_alin)



######################################################################
## (4.4) Genera 500 secuencias aleatorias de N aminoácidos, siendo N
## en cada ocasión un valor aleatorio dentro del rango de tamaños de
## las proteínas del conjunto de 10 asignado. Calcula la valoración de
## una alineación local óptima entre cada una de las secuencias
## aleatorias generadas y la proteína de referencia. Determina si las
## valoraciones del apartado anterior pueden considerarse
## estadísticamente significativas en base a estos resultados.

## Solución:

# El vector tam almacenará los tamaños de cada proteína del conjunto de 10 
# asignado.
tam <- c()
for (i in 1:length(p10_sec)){
  tam <- c(tam, length(p10_sec[[i]]))
}

# El vector simb almacenará todos los símbolos de las proteínas del 
# conjunto de 10 asignado.
simb <- c()
for (j in 1:length(p10_sec)){
  for (k in 1: length(p10_sec[[j]])){
    simb <- c(simb, p10_sec[[j]][[k]])
  }
}
# Obtenemos el alfabeto
simb <- c(unique(simb))

# Creamos otro vector Nvalues que contenga el rango de tamaños que incluye a los
# almacenados en el vector tam. Utilizamos las funciones min y max.
Nvalues <- c(min(tam): max(tam)); Nvalues

# Generación de las secuencias aleatorias. Utilizamos sample tanto para obtener
# la sucesión de símbolos aleatoria como para seleccionar el tamaño aleatorio de
# la misma.
s <- 1
secAleat <- c() # lista que almacenará cada secuencia aleatoria generada
while (s < 501){
  secAleat[[s]] <- sample(simb, sample(Nvalues,1), TRUE)
  s <- s+1
}

# Primeras 6 secuencias aleatorias generadas
head(secAleat); length(secAleat)


# Generación de la puntuación de alineación local entre cada secuencia aleatoria
# y la proteína wP_011157096

local_alinAleat <- c() # vector que almacenará todas las puntuaciones obtenidas
for (i in 1:length(secAleat)){
  local_alinAleat <- c(local_alinAleat, SmithWaterman(WP_011157096, 
                                                      secAleat[[i]], Punt_B62))
}

# Planteamos un contraste de hipótesis para determinar si siguen, o no, las
# secuencias generadas aleatoriamente una distribución normal, con un nivel
# de confianza del 95% (alpha = 1 - 0.95)

# Definimos como hipótesis nula H0 que estas puntuaciones siguen una 
# distribución normal y utilizamos el p-valor que obtengamos con shapiro.test 
# para aceptar o rechazarla. Si supera el umbral, la aceptaremos, si no, se 
# rechazará y se dará por válida la hipótesis H1 (no siguen una distribución
# normal)

alpha = 1 - 0.95

# Análisis normalidad de las puntuaciones de las secuencias generadas
# aleatoriamente

H0 <- shapiro.test(local_alinAleat)$p.value > alpha; H0
if(H0 == TRUE){
  print("Se acepta H0: siguen una distribución normal")
} else{
  print("Se rechaza H0: no siguen una distribución normal")
}

hist(local_alinAleat)

# No supera el umbral, por tanto, deberíamos de rechazar la hipótesis nula, pero
# tal y como se ve en el histograma, sí presenta una representación similar
# a la campana de Gauss. Por tanto, rechazamos este método de estudio y 
# optamos por otro.


# En concreto, he decidido comparar las puntuaciones obtenidas con las 
# secuencias generadas aleatoriamente, y estudiar el porcentaje de puntuaciones
# que son superiores e inferiores a las puntuaciones obtenidas con las
# secuencias correspondientes al conjunto de 10 proteínas asignado.

# De esta forma, procedemos a calcular el porcentaje de valores que no supera
# a las valoraciones obtenidas con el conjunto de 10 proteínas y el porcentaje
# que es inferior.

inf <- 0 # contador de las veces que la puntuación es inferior a la 
# puntuación más pequeña obtenida con las secuencias de 10 proteínas
for (k in 1: length(local_alinAleat)){
  if(local_alinAleat[k] < min(local_alin)){
    inf <- inf +1
  }
}


# Como lo que buscamos es, determinar si las valoraciones obtenidas con las
# secuencias del conjunto de 10 proteínas asignado son estadísticamente
# significativas, establecemos un umbral con que el podamos afirmar que esto 
# sucede. En este caso, consideraremos un nivel de confianza del 98%. 
# Así, pues:

# Pasamos a porcentaje inf
pInf <- (inf/length(local_alinAleat)) *100; pInf

pInf > 95

# Como el 100% de las puntuaciones obtenidas por alineación local con
# secuencias generadas de forma aleatoria supera el umbral del 95% que se 
# había fijado, podemos afirmar que los resultados obtenidos con las secuencias
# pertenecientes al conjunto de 10 proteínas son significativas estadísticamente, 
# ya que no resultan en su conjunto, valores similares obtenidos con secuencias 
# que sean producto del azar.


######################################################################


