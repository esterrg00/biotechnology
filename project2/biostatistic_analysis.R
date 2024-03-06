## Uso perfiles para caracterizar las secuencias de una alineación
## múltiple.
######################################################################
## Algunas funciones de interés
######################################################################
## Funciones para la lectura de datos a partir del contenido de un
## archivo: scan y read.table

## Funciones para la manipulación de textos con R: grep, nchar, paste,
## sprintf, substr, strsplit, regex, gregexpr

## Funciones para la manipulación de textos (biblioteca seqinr):
## count, s2c, c2s

## Funciones para la aplicación de funciones: lapply, sapply, wapply,
## mapply, apply, tapply, aggregate, by
######################################################################
## Definir un conjunto de funciones que permitan:
## (1) (3 puntos) Obtener los datos de una alineación a partir del contenido de
##     un fichero de texto. Debe incluirse una descripción del tipo de
##     datos utilizado para contener toda información sobre la
##     alineación (alfabeto, alineación, secuencia, ...).
## Formato del fichero de texto:
## -------------------------------------------------------------------
## Alfabeto: {L[,L]*}
## Alineación:
##     [{L|-}]+
##     [{L|-}]+
##     [{L|-}]+
##     [{L|-}]+
##     [{L|-}]+
## Secuencia: [L]+
## -------------------------------------------------------------------
##     * La primera línea del fichero contiene el texto "Alfabeto: "
##       y, a continuación, entre "{" y "}" una sucesión de letras (al
##       menos una, pero la cantidad puede variar de un fichero a
##       otro) separadas por comas.
##     * La segunda línea del fichero únicamente contiene el texto
##       "Alineación:"
##     * Entre las líneas 3 y 7 están las secuencias que componen
##       la alineación: sucesiones de letras del alfabeto, sin
##       espacios entre ellas, que pueden contener el símbolo "-". Una
##       secuencia por línea, todas de la misma longitud (aunque esta
##       puede variar de un fichero a otro).
##     * La línea 8 contiene el texto "Secuencia: " y, a continuación,
##       una sucesión de letras del alfabeto, sin separación entre
##       ellas, cuya longitud puede variar de un fichero a otro.


## IMPORTACIÓN DE LIBRERÍAS ##
library(seqinr)

# Con esta función se pretende obtener una lista con la que se pueda acceder
# a cada una de las partes que conforman el fichero: alfabeto, alineación y 
# secuencia.

obtenerDatos <- function(fichero){
  #En la lista datos almacenaremos el alfabeto, la alineación y la secuencia
  datos <- list(NA, NA, NA)
  names(datos) <- c("alf", "alin", "sec")
  
  # Obtención del alfabeto
  fAlf <- scan(fichero, what=character(), nlines = 1, sep = ",") # línea del 
  # alfabeto sin "pulir"
  alf <- c() # vector que almacenará el alfabeto
  alf[1] <- substr(fAlf[1], 12, 12) # nos quedamos sólo con el primer símbolo
  # nos quedamos sólo con el elemento en la posición 12, porque la primera
  # línea del fichero siempre comenzará de la forma "Alfabeto: {"
  
  # Obtención del resto de símbolos
  for (i in 2:length(fAlf)){
    alf[i] <- substr(fAlf[i], 2,2) # eliminamos el espacio y nos quedamos sólo 
    # con el símbolo que nos interesa
  }
  # Guardamos en la lista
  datos[["alf"]]<- alf
  
  
  # Obtención de la alineación
  # Realizamos la lectura del fichero a partir de la línea 3 (nos saltamos las
  # dos primeras) y leemos hasta la línea 7 (7-3= 5)
  fAlin <- scan("datos03.txt", what=character(), skip = 2, nlines = 5)
  
  # En alin almacenamos en cada posición un vector con cada uno de
  # los símbolos que componen cada una de las líneas de la alineación 
  # del fichero.
  alin <- c()
  for (j in 1:length(fAlin)){
    alin[[j]]<- s2c(fAlin[j]) 
  }
  
  # Guardamos en la lista
  datos[["alin"]] <- alin
  
  # Obtención de la secuencia
  # Realizamos la lectura del fichero de la línea 8, que es la que contiene 
  # la secuencia.
  fSec <- scan("datos03.txt", what=character(), skip = 7, nlines = 1)
  # Como obtenemos dos elementos, porque por defecto considera como separador un
  # espacio, nos quedamos sólo con el segundo elemento de ese vector, que será
  # la secuencia. Además, también la pasamos a un vector con s2c()
  sec <- s2c(fSec[2])
  # Guardamos en la lista
  datos[["sec"]] <- sec
  
  return(datos)
  
}

# De esta forma, obtenemos una lista con cada uno de los componentes que nos
# proporciona el fichero de texto.

datos <- obtenerDatos("datos03.txt")

## (2) (3 puntos) Construir un HMM (la primera aproximación vista en teoría y
##     aplicando suavizado, únicamente, en la estimación de la matriz
##     de observaciones) a partir de la información obtenida sobre una
##     alineación y por tanto, según el tipo de datos descrito en el
##     apartado anterior. Debe incluirse una descripción del tipo de
##     datos utilizado para contener toda la información sobre el HMM.

# Los elementos que componen un HMM son un conjunto de estados, est; 
# una matriz de transición, mt; una lista que contiene los vectores 
# iniciales para cada región e, pi; un alfabeto, alf, y una matriz
# de observaciones, mo.

# Así pues, planteamos la funcion hmm() que dada una alineación y el alfabeto
# correspondiente a la misma, devuelva un modelo oculto de Markov el cual 
# está compuesto por los siguientes elementos:
# - Conjunto de estados, est; 
# - Matriz de transición, mt;
# - Lista que contiene las probabilidades iniciales para cada estado e, pi;
# - Alfabeto de la alineación, alf;
# - Matriz de observaciones, mo

# Por tanto, a partir de las variables de entrada, la función devolverá
# la matriz de transición, mt; la matriz de observaciones, mo; y la lista
# que contenga las probabilidades inciciales, pi.

alf <- datos$alf # alfabeto
alin <- datos$alin # alineación de secuencias

# Definimos una función estados() que dada una alineacion, devuelva una lista
# que contenga el conjunto de estados asociado (conj) y el estado que tiene cada
# columna (col) 
estados <- function(al){
  est <- c() # almacenará en la posición de cada línea de la alineacion
  # un vector con cada una de las regiones a la que pertenece cada uno
  # de los símbolos de esa línea
  
  for (a in 1:length(al)){ # recorremos cada sencuencia de la alineación
    lest <- c() # creamos un vector para cada secuencia que almacene los estados 
    # a los que pertenece
    for (i in 1:length(al[[a]])){ # recorremos cada símbolo de cada secuencia
      if (al[[a]][[i]] == "-"){ # si es hueco
        lest[i] <- "I" # lo asignamos como inserción
      } else { # en cualquier otro caso
        lest[i] <- "M" # lo asignamos como coincidencia
      }
    }
    est[[a]] <- lest # guardamos en cada posición correspondiente para cada
    # secuencia, la lista con los estados para cada uno de
    # sus símbolos
  }
  
  estAlin <- c() # vector que tendrá los estados de la alineacion en su conjunto
  # Pasamos a matriz, porque es más fácil de recorrer cada columna para definir
  # el estado mayoritario (I o M) en cada columna, y poder conseguir el conjunto
  # de estados que definen nuestro modelo
  
  mt <- matrix(nrow=5, ncol=length(est[[1]]))
  for (i in 1:length(est)){
    for (j in 1:length(est[[i]])){
      mt[i,j] <- est[[i]][[j]]
    }
  }
  
  estCol <- c() # almacena vectores con los estados que hay en cada columna de 
  # la alineación
  for (b in 1:ncol(mt)){
    estCol[[b]] <- c(mt[, b])
  }
  
  estMasRep <- c() # vector que almacena, para cada columna el estado que más
  # se repite
  for (p in 1:length(estCol)){
    cM <- 0 # contador M
    cI <- 0 # contador I
    # Recorremos cada "estado" en cada columna y realizamos el conteo de "M"
    # e "I" para poder conocer cuál es el "estado" mayoritario
    for (k in 1:length(estCol[[p]])){ 
      if (estCol[[p]][k] == "M"){
        cM <- cM+1
      } else{
        cI <- cI+1
      }
    }
    if (cM < cI){
      estMasRep[p] <- "I"
    } else{
      estMasRep[p] <- "M"
    }
  }
  
  estFin <- c() # vector que almacena el conjunto de estados de la alineación
  contadorM <- 1 # contador para regiones M
  contadorI <- 0 # contador para regiones I
  # Recorremos la lista de los estados más repetidos en cada columna
  for (i in 1:length(estMasRep)){ 
    # Si el iterador i no corresponde con el último elemento
    if (i < length(estMasRep)){ 
      # Si el elemento de en posición i es región "I" 
      if(estMasRep[i] == "I"){
        # Comprobamos si la región siguiente es también "I"
        if(estMasRep[i] == estMasRep[i+1]){ # si coinciden
          # Renombramos en ambas posiciones el estado con el contador I
          estFin[i] <- estFin[i+1] <- paste("I", contadorI, sep="")
        } else{ # si no coincide la region en posicion i con su siguiente (i+1)
          # Renombramos sólo en esa posición el estado con el contador M
          estFin[i] <- paste("I", contadorM-1, sep="")
          # El valor que le asignamos será el mismo valor que esté almacenando 
          # el contador de regiones M - 1, ya que nombramos las regiones I con 
          # el índice de la región M a la que precede. Debemos de restar una 
          # posición porque el contador M se incrementa cada vez que se asigna.
          contadorI <- contadorM-1 # actualización del contador I por la último
          # cifra asignada, para poder ir teniendo un
          # seguimiento de las regiones que ya han sido cifradas
        }
      } else{ #Para el caso en el que sea región "M"
        # Renombramos a esa región con su índice adecuado con ayuda del contadorM
        estFin[i] <- paste("M", contadorM, sep="") 
        contadorM <- contadorM+1 # actualización del contador
      }
    } else{ # Si estamos en la última posición, no habrá región siguiente, y
      # por tanto procedemos de la siguiente forma
      if(estMasRep[i] == "I"){ # si es región "I"
        # Renombramos ese estado con la posición de la región a la que precede
        estFin[i] <- paste("I", i-1,sep="") 
      } else{ # si es "M"
        # Renombramos a esa región con su índice adecuado con el contadorM
        estFin[i] <- paste("M", contadorM,sep="")
      }
    }
  }
  
  return(list(conj = unique(estFin), col = estFin))
}

# Definimos una función buscaAnterior() que dado un vector noNulos, si 
# en la alineacion, el elemento en la posición pos es "-" , tenemos que 
# buscar el estado anterior no nulo a esa posición. 
# Con esta funcion podemos acceder a la posicion anterior donde hay un estado.

buscaAnterior<- function(vector, pos){
  pVector <- 0
  # Localizamos la posición en la que se encuentra pos en el vector para poder
  # buscar desde ahí hacia atrás el estado precesor a él en noNulo, que
  # corresponderá con el estado en que realmente se corresponde su transición
  for (i in 1: length(vector)){ # recorremos el vector de entrada
    if(vector[i] == pos){ # si coincide la posición
      pVector <- i # la guardamos en la variable pVector
    }
  }
  # Recorremos hacia atrás el vector desde la posicion que acabamos de guardar
  for (j in (pVector-1):1){ 
    anterior <- vector[j] # guardamos el valor que corresponde al elemento 
    # justo anterior a pVector
    break
  }
  return (anterior) # lo devolvemos
}


# Función que permite crear un modelo oculto de Markov a partir de una
# alineación y un alfabeto.
hmm <- function(alineacion, alfabeto){
  # Definimos algunas variables con datos que nos serán útiles
  est <- estados(alineacion)$conj # conjunto de estados de la alineación
  todosEst <- estados(alineacion)$col # estado al que pertenece cada columna 
  # de la alineación
  le <- length(est) # tamaño del conjunto de estados
  lAlin <- length(alineacion) # tamaño de la alineación (número de secuencias
  # que la forman)
  lAlinSec <- length(alineacion[[1]]) # tamaño de una línea cualquiera de la
  # alineación (todas tendrán el mismo tamaño, por eso no importa cuál coger)
  
  # Trasladamos la alineación de secuencias a una matriz para poder trabajar
  # más facilmente
  mcol <- matrix(nrow=lAlin, ncol=lAlinSec)
  for (i in 1:lAlin){
    for (j in 1:lAlinSec){
      mcol[i,j] <- alin[[i]][[j]]
    }
  }
  
  ############### Probabilidades de inicio ###############
  # Queremos conocer P(E1 = e), siendo e cada uno de los estados de est.
  
  piVector <- c() # vector en el que almacenaremos para cada secuencia en qué
  # estado se encuentra su primer símbolo, así como su posición.
  
  # Recorremos cada posición de la matriz, para poder en cuenta en qué estado ha 
  # comenzado cada una de las secuencias de la alineación
  for (i in 1:nrow(mcol)){  
    for (j in 1:ncol(mcol)){
      # Guardamos en la posición i de piVector el estado correspondiente a la 
      # posición en donde se ha encontrado el primer símbolo diferente a "-",
      # así como esta posición.
      if(mcol[i,j] != "-"){ 
        piVector[[i]] <- c(todosEst[j], j)
        break
      }
    }
  }
  
  # Ahora vemos el número de repetición de cada uno de estos valores, y dividimos
  # entre el total de secuencias para poder hallar la probabilidad de incio.
  # Estas probabilidades las dejaremos almacenadas en el vector pi, que tendrá
  # el tamaño del conjunto de estados definido, est. En caso de que no se haya
  # dado algún estado para el inicio de una secuencia, su probabilidad inicial
  # será nula.
  
  pi <- c() # probabilidades iniciales
  for (e in 1:le){ # para cada estado en el conjunto de estados
    pi[e] <- 0 # le asignamos, en principio, una frecuencia de aparición nula
    for (p in 1:length(piVector)){ # para cada estado observado en piVector
      if (est[e] == piVector[[p]][1]){ # si aparece el estado e en piVector
        pi[e] <- pi[e] + 1 # incrementamos en una unidad su frecuencia
      }
    }
  }
  # Dividimos cada una de los valores entre el total de observaciones (que 
  # coincide con el número de secuencias de la alineación, las filas de la 
  # matriz mcol)
  for (q in 1:length(pi)){
    pi[q] <- pi[q]/nrow(mcol) # probabilidad
  }
  
  # Así obtenemos un vector que en cada una de sus posiciones i contiene la
  # probabilidad incial de aparición del estado en la posición i del vector est.
  names(pi) <- est # asignamos a cada posición el estado correspondiente
  
  
  
  ############### Matriz de transicion de estados ###############
  
  # Queremos conocer P(E[i+1]=q | E[i]=e), donde e,q pertenecen al conjunto de
  # estados de la alineación.
  
  # Tenemos que recorrer cada secuencia, y, a partir de conocer la región en
  # la que se encuentra su primer símbolo y su posición (obtenido con piVector),
  # debemos de estudiar desde ahí en adelante la región a la que pertenece el
  # símbolo siguiente distinto a "-" de esta posición.
  
  # Es decir, si definimos r como cada una de las posiciones que hay en
  # piVector: podemos saber que la secuencia r de la alineación comienza  
  # en el estado piVector[[r]][1]
  
  lpiV <- length(piVector) # tamaño de piVector
  transicionesVector <- c() # lista que contendrá para cada secuencia, un vector
  # con los estados siguientes a cada uno de sus estados
  
  for(r in 1:lpiV){ # para cada elemento de piVector
    primerEst <- piVector[[r]][2] # guardamos en primerEst la posición donde  
    # está el primer símbolo en la secuencia r
    tPos <- c() # vector donde están las transiciones para cada secuencia
    for (b in primerEst:lAlinSec){ # para cada letra desde pos hasta el final-1
      # (-1 para que podamos acceder a la última posición)
      if(b < lAlinSec){ # verificamos que no estamos en la última posición
        if (mcol[r, b+1] != "-"){ # comprobamos que en esa coordenada no hay "-"
          tPos[b] <- todosEst[b+1] # guardamos el estado al que pertenece el 
          # siguiente símbolo
        }
      } else{ # en caso de estar en la última posición
        tPos[b] <- NA 
      }
    }
    transicionesVector[[r]] <- tPos # guardamos en la posición r de la lista el
    # vector resultante (tPos)
  }
  
  
  # Creación de la matriz de transición con ceros
  mt <- matrix(0, nrow=le, ncol =le)
  
  # A continuación procedemos a hacer el conteo de cada transición que se ha 
  # dado en la alineación.
  for (t in 1:length(transicionesVector)){ # recorremos cada secuencia   
    for (m in 1:length(transicionesVector[[t]])-1){ 
      # Guardamos en noNulos las posiciones de las regiones que no son NA
      noNulos <- which(!is.na(transicionesVector[[t]]))
      
      # Verificamos que la posición m está en noNulos, es decir, hay un 
      # siguiente. Además, verificamos que en esa posición de la secuencia 
      # hay un símbolo y no un hueco.
      if(m%in%noNulos){
        # Guardamos cada uno de esos estados en variables para facilitar el
        # manejo de estos valores
        estSig <- transicionesVector[[t]][m] # estado siguiente
        
        # Si en la alineación, el estado de referencia de la secuencia t en la 
        # posición m es un hueco, buscamos el estado anterior no nulo con
        # la función definida buscaAnterior.
        if (mcol[t, m] == "-"){ 
          anterior <- buscaAnterior(noNulos, m)
          estRef <- transicionesVector[[t]][anterior]
        } else{
          estRef <- todosEst[m] # estado de referencia
        }
        # Debemos conocer la posicion que ocupan estos estados en el vector est
        # para poder actualizar la coordenada [estRef, estSig] de la matriz
        # de transición.
        i <- which(est == estRef)
        j <- which(est == estSig)
        # Una vez localizada la coordenada, actualizamos su valor
        mt[i,j] <- mt[i,j] + 1
      }
    }
  }
  
  # Procedemos a realizar la proporción de cada una de las transiciones. Para
  # ello, tenemos que dividir entre el número de símbolos que hay en cada estRef.
  
  nSimbEstados <- c() # lista que almacena para cada estado, el número de
  # símbolos pertenecientes al alfabeto (será el denominador
  # para cada estRef, es decir, las filas de la matriz mt)
  
  # Es decir, tendrá tantos elementos como estados haya en nuestra alineación.
  
  for(region in 1:le){ # para cada estado contar el número de símbolos
    apariciones <- c() # almacena las apariciones por cada estado
    for(m in 1:length(todosEst)){
      if(todosEst[m] == est[region]){ # si coinciden los estados en los que 
        # estamos mirando
        # Guardamos en apariciones la suma del valor que contenga, en caso de
        # que tuviera alguno debido a ya haber aparecido ese estado antes, con
        # el tamaño del vector que contiene las posiciones de aquellos elementos
        # de la columna m de la matriz mcol (alineación de las secuencias) que
        # sean distintas al hueco
        apariciones <- sum(apariciones, length(which(mcol[,m] != "-"))) 
      }
    }
    # Guardamos este valor en la posición correspondiente a su estado
    nSimbEstados[[region]] <- apariciones
  }
  
  # Recorremos cada coordenada(fila,columna) de la matriz mt, y sólo a
  # aquellas coordenadas en las que la fila (estRef (Ei)) coincide con k, las 
  # dividimos por el valor resultante obtenido en nSimbEstados en la posición k.
  for (k in 1:length(nSimbEstados)){
    for(fila in 1:nrow(mt)){
      for(columna in 1:ncol(mt)){
        if(fila == k){
          mt[fila,columna] <- mt[fila,columna]/nSimbEstados[[k]]
        }
      }
    }
  }
  
  # Asignamos nombres a las filas y columnas para una mejor lectura
  dimnames(mt) <- list(est,est)
  
  
  ############### Matriz de observaciones ###############
  
  # Queremos conocer P(S[i]=a | E[i]=e), donde e pertenece al conjunto de
  # estados de la alineación, y a pertenece al alfabeto de la alineación
  
  # La idea que tenemos que buscar es ir recorriendo cada una de los símbolos
  # de las columnas e ir guardando el número de veces que aparecen cada una con
  # respecto a un determinado estado (columnas).
  
  # Primero, creamos un vector para cada estado que almacene todos los símbolos 
  # que aparecen en él (buscar aquellos que sean distintos de "-").
  
  simbolosEst <- c() # lista que almacenará en cada posición (correspondiente 
  # con cada uno de los estados de la alineación), los 
  # símbolos que aparecen en el mismo estado
  for (e in 1:le){ # para cada estado
    sEst <- c() # creamos un vector que almacene los símbolos de ese estado
    for(j in 1: ncol(mcol)){ #recorremos cada una de las columnas de la alineación
      if (est[e] == todosEst[j]){ # si el estado coincide con el que le corresponde
        # a la columna que estamos consultando
        leidos <- which(mcol[,j] != "-") # nos quedamos con los distintos de "-"
        # obtenemos las "filas" de la matriz
        for(q in 1:length(leidos)){ # para cada uno de ellos
          sEst <- c(sEst, mcol[leidos[q],j]) #lo añadimos al vector de ese estado
        }
      }
    }
    simbolosEst[[e]] <- sEst # asignamos el vector sEst a la posición del 
    # estado que hemos consultado
  }
  
  # Explicación de simbolosEst: por ejemplo, simbolosEst[[2]] contiene los 
  # símbolos que hay en el estado correspondiente con est[2]. En este caso,
  # correspondería con el estado M1.
  
  # Creamos la matriz de observaciones, mo
  # En las filas tenemos cada uno de los estados de la alineación, y en las
  # columnas encontramos cada uno de los símbolos del alfabeto asociado a la
  # alineación.
  mo <- matrix(0, nrow=le, ncol =length(alfabeto))
  
  # Ahora debemos de recorrer cada posible combinación, y realizar el conteo
  # de ocurrencia de aparición de los símbolos
  
  for (l in 1:length(simbolosEst)){ # recorremos la lista simbolosEst
    for(t in 1:length(simbolosEst[[l]])){ # recorremos cada elemento de cada 
      # posición de simbolosEst
      for(i in 1:nrow(mo)){ # recorremos las filas (estados) de la matriz
        if(est[l] == est[i]){ # si coincide ese estado con el que estamos 
          # consultando en simbolosEst
          for(j in 1:ncol(mo)){ # recorremos las columnas (símbolos) de la matriz
            # Si coindice ese símbolo con el que aparece en el vector de la 
            # posición t de la lista simbolosEst
            if(simbolosEst[[l]][t] == alfabeto[j]){ 
              mo[i,j] <- mo[i,j] + 1 # incrementamos el valor que tenía en esa
              # coordenada (pareja de Estado, Simbolo)
            }
          }        
        }
      }
    }
  }
  
  
  # Como tenemos que aplicar suavizado, incrementaremos cada coordenada en una 
  # unidad, ya que estaremos introduciendo un caso ficticio por cada símbolo
  # del alfabeto asociado a la alineación.
  # De esta forma, estamos considerando, aunque sea con valores muy pequeños, 
  # la posibilidad de que aparezca cualquier símbolo del alfabeto en 
  # cualquier estado.
  
  mo <- mo + 1
  
  # Ahora debemos de dividir cada ocurrencia contada por el total de símbolos  
  # que aparecen en cada estado. Para ello, utilizaremos el tamaño de cada
  # vector que contiene la lista simbolosEst. Como hay que aplicar suavizado,
  # a este valor tendremos que sumarle el total de casos ficticios que hemos
  # introducido, es decir, el tamaño del alfabeto.
  
  for (l in 1:length(simbolosEst)){ # recorremos la lista simbolosEst
    for(i in 1:nrow(mo)){
      if(est[l] == est[i]){
        mo[i,] <- mo[i,] / (length(simbolosEst[[l]])+length(alfabeto))
      }
    }
  }
  
  # Asignamos nombres a los índices de la matriz para una mejor lectura
  dimnames(mo) <- list(est,alfabeto)
  
  # Devolvemos en una lista cada uno de los conjuntos de probabilidades
  # calculadas
  return(list(pi = pi, mt = mt, mo= mo))
  
}


## (3) (3 puntos) Calcular la logProbabilidad, utilizando adecuadamente la
##     versión normalizada del algoritmo de avance, de la última
##     secuencia del archivo según el modelo anterior. Hay que tener
##     en cuenta que no hay límite para la longitud de las secuencias
##     de la alineación y de última secuencia ¿es significativo el
##     resultado obtenido?.

# Con el algoritmo de avance normalizado, lo que estamos calculando son
# probabilidades condicionadas, es decir, P(E[t]= e | s[1]...s[t]), siendo e 
# cualquier estado del conjunto de estados, y s cualquier símbolo perteneciente 
# al alfabeto de la alineación que estamos estudiando.
# Debemos saber, que para calcular la logProbabilidad de la secuencia a partir
# de los datos que se vayan obteniendo con este algoritmo equivalen a la 
# suma de los gamma calculados.

## NO FUNCIONA ## No sé qué puede estar pasando. He trabajado con las listas
# de igual forma que lo he hecho en el resto de ejercicios y he estado 
# haciendo pruebas por separado y sí se incluyen los elementos en 
# alphaPrima[[k]] y alphaNorm[[k]], pero al meterlos en los bucles, me sale el
# siguiente fallo: Error in `*tmp*`[[k]] : subíndice fuera de  los límites.

avanceNorm <- function(modelo, secuenciaObservada){
  tam <- length(secuenciaObservada) # longitud de la secuencia
  conjEstados <- names(modelo$pi) # estados del modelo
  mTrans <- modelo$mt # matriz de transición
  mObs <- modelo$mo # matriz de observaciones
  pi <- modelo$pi # probabilidades iniciales
  
  alphaPrima<- c() # lista que almacenará los valores de alpha
  # Realizamos el primer paso a parte
  alphaPrima[[1]] <- pi + mObs[,secuenciaObservada[1]]
  alphaNorm <- c() # lista que almacenará los valores de alpha normalizados
  gamma <- c() # vector que almacenará los inversos de la suma de los alphaPrima
  gamma[1] <- 1 / sum(alphaPrima[[1]])
  alphaNorm[[1]] <- alphaPrima[[1]] * gamma[1]
  
  for(k in 2:tam){
    alphaPrima[[k]] <- c()
    alphaNorm[[k]] <- c()
    for(j in 1: length(conjEstados)){
      # Para cada estado en la posición j que consultemos, calculamos el valor
      # de la probabilidad de que es
      sumatorio <- sum(mTrans[,j]*alphaNorm[[k-1]])
      alphaPrima[[k]][j] <- mObs[j,secuenciaObservada[k]] * sumatorio
    }    
    gamma[k] <- 1/ sum(alphaPrima[[k]])
    alphaNorm[[k]] <- alphaPrima[[k]] * gamma[k] 
  }
  
  # Devolvemos la suma de los logaritmos de las variables gamma para poder
  # conocer la probabilidad de la secuencia
  return(sum(log(gamma)))
  
}

# A pesar de no poder visualizar los resultados, a continuación se ha 
# realizado el procedimiento llevado a cabo para poder considerar
# significativo o no el resultado obtenido

# Para estudiar el nivel de significación de la secuencia observada,
# procederemos a generar una gran cantidad de secuencias aleatorias de 
# igual tamaño. En concreto, generaremos 500. De esta forma, podremos 
# estudiar si existe una gran diferencia entre la logProbabilidad 
# obtenida con la secuencia observada, en comparación con la de las 
# secuencias aleatorias.

obs <- datos$sec # observación que vamos a estudiar
mom <- hmm (alin, alf) # modelo oculto de Markov que utilizaremos

# Para evaluar esto, se realizará un sesgo en el cual comprobaremos si
# la logP de las secuencias aleatorias, por ejemplo, el 95%, difiere con
# la logP de la secuencia observada.

# Si obtenemos que la diferencia es superior al 95%, podremos considerar
# que realmente, el resultado obtenido con el algoritmo de avance normalizado
# sí es significativo; mientras que, si no supera este umbral, podemos estimar
# que el resultado no es lo suficientemente significativo como para tenerlo
# en consideración en nuestro estudio.

s <- 1 # iterador
aleatorias <- c() # lista que almacenará cada secuencia aleatoria generada
while (s < 501){
  # Añadimos en la posición correspondiente s de la lista aleatorias cada una
  # de las secuencias generadas aleatoriamente.
  aleatorias[[s]] <- sample(alf, length(obs), TRUE)
  s <- s+1
}

head(aleatorias) # primeras 6 secuencias generadas aleatoriamente

# A continuación procedemos a aplicar el algoritmo de avance normalizado para
# cada una de estas secuencias, empleando el mismo Modelo Oculto de Markov.

logPAleatorias <- c() # vector que almacenará todas las logP obtenidas para
# cada secuencia generada aleatoriamente
for (i in 1:length(aleatorias)){
  logPAleatorias[[i]] <- avanceNorm(mom, aleatorias[[i]])
}

# Guardamos en una variable la logP de que ocurra la secuencia observada.
logPObservada <- avanceNorm(mom, obs)


# Procedemos a realizar un conteo de aquellas secuencias aleatorias cuya 
# logP es inferior a la logP de la secuencia observada
contador <- 0 # contador que representa el número de veces de coincidencia
for (p in 1:length(logPAleatorias)){
  if (logPAleatorias[[p]] < logPObservada) {
    contador <- contador + 1
  }
}

# Pasamos los resultados a porcentaje
porcentajeContador <- (contador/length(aleatorias))*100

# Comprobamos si supera el umbral
porcentajeContador > 95 

# Si lo supera: el resultado obtenido con el algoritmo de Avance normalizado
# para la secuencia observada podemos considerarlo significativo, ya que,
# generando secuencias al azar de igual tamaño, la logProb de ocurrencia es
# menor que la logP de la secuencia observada, y por tanto, habría que
# tener en cuenta en nuestro estudio esta secuencia, puesto que no es "normal"
# que se obtenga una probabilidad de ocurrencia tan elevada para secuencias de 
# ese tamaño.

# En caso de no superar este umbral, NO podemos considerar significativo el
# resultado obtenido.
