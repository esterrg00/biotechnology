######################################################################
## Biotecnología
## Ejercicios de introducción a R

######################################################################
## Ejercicio 1. Disponemos de datos sobre la fecundidad (número de
## huevos puestos por día durante los primeros 14 días de vida) de 25
## hembras de tres líneas génicas diferentes de la mosca de la fruta
## Drosophila melanogaster. Estas líneas están etiquetadas como
## Resistente y Susceptible ya que se alimentaron selectivamente
## buscando resistencia y susceptibilidad a DDT, respectivamente. Se
## usa una línea no seleccionada como control. Estos datos están
## almacenados en el fichero drosophila.txt

## Se pide:

## (1) Leer los datos anteriores y almacenarlos en un vector de datos.

## Nota: Para obtener un vector con todo los datos contenidos en
## un fichero de texto: scan (ver Ayuda).

fecundidad <- scan(file="drosophila.txt")


## (2) Genera un vector con la línea génica de cada muestra (los datos
## están ordenados de tal forma que, los 25 primeros se corresponden a
## la línea Resistente, los siguientes 25 a la línea Susceptible y los
## restantes 25 a la línea Control.


linea <- rep(c("Resistente", "Susceptible", "Control"), each=25)

             
## (3) Calcula la media y la desviación típica de la fecundidad de
## todas las hembras. 

mean(fecundidad)
sd(fecundidad)


## (4) Determina el grupo (línea génica) al que pertenece la hembra
## más fecunda. 


linea[which.max(fecundidad)]

## (5) Determina el grupo al que pertenece la hembra menos fecunda. 


linea[which.min(fecundidad)]

## (6) Representa en un histograma la distribución de la fecundidad en
## el conjunto total de las muestras.

hist(fecundidad)

## (7) Almacena la fecundidad de cada uno de los grupos en vectores
## diferentes 

resistente <- fecundidad[1:25]
susceptible <- fecundidad[26:50]
control <- fecundidad[51:75]

##Otra forma si valores no están ordenados

resistentes <- fecundidad[linea=="Resistente"]
susceptibles <- fecundidad[linea=="Susceptible"]
controls <- fecundidad[linea=="Control"]


## (8) Calcula la media y la desviación típica de la fecundidad para
## cada grupo 

mean(resistente)
sd(resistente)

mean(susceptible)
sd(susceptible)

mean(control)
sd(control)

## (9) Representa en un gráfico de barras la información obtenida en
## el apartado anterior.
medias <- c(mean(resistente),mean(susceptible),mean(control))
varianzas <- c(sd(resistente), sd(susceptible), sd(control))
vals <- barplot(medias, main="Fecundidad en la mosca Drosophila (media de huevos)",
                names=c("Resistencia", "Susceptible", "Control"), xlab="Linea génica", 
                ylab="Media del n. de huevos", col=2:4,ylim=c(0,50))


arrows(vals, medias- varianzas, vals, medias+varianzas, angle = 90, code =3)

## (10) Representa en un gráfico de caja y bigotes la información
## relativa a los cuartiles de las distribución de la fecundidad en
## cada línea.

boxplot(fecundidad[linea=="Resistente"], fecundidad[linea=="Susceptible"], 
        fecundidad[linea=="Control"], names= c("Resistente","Susceptible","Control"),
        xlab ="Línea génica")


######################################################################
## Ejercicio 2: Obtener el genoma completo de: Haemophilus influenzae
## Rd KW20.

## Para ello acceder a la página web de GenBank
##                 http://www.ncbi.nlm.nih.gov/genbank/
## y elegir la opción Nucleotide, y en la pestaña de la derecha Send
## elegimos: complete record, file y formato FASTA.

## (1) Guardar el genoma obtenido en un fichero de texto:
## hinfluenzae.fasta

## (2) Utilizando las funciones read.fasta y getSequence de la biblioteca 
## seqinr obtener un vector con la sucesión de nucleótidos del genoma del
## Haemophilus influenzae Rd KW20.

## read.fasta (biblioteca seqinr): Lista con las secuencias contenidas
## en un fichero con formato fasta.
## getSequence (biblioteca seqinr): Lista con los vectores de
## caracteres contenidos en una lista de secuencias.

library(seqinr)

secuencias <- read.fasta("hinfluenzae.fasta")

## (3) Calcular las frecuencias absolutas y relativas de los
## nucleótidos: A, C, G y T en el vector hinfluenzae

hinfluenzae <- getSequence(secuencias)[[1]]
length(hinfluenzae)

unique(hinfluenzae) #para decir los nucleotidos q componen la secuencia

table(hinfluenzae) #frecuencias absolutas

dna <- c("a", "c", "g","t")

factor(hinfluenzae, dna)

frecAbs <- table(factor(hinfluenzae, dna))
frecAbs
  

frecRel <- frecAbs/length(hinfluenzae)

#frecRel <- frecAbs / sum(frecAbs) #para tener en cuenta sólo las de a c t g


## (4) Definir una función proporciones que, dada una secuencia sobre
## un alfabeto (por defecto, el alfabeto {A, C, G T}), devuelva las
## frecuencias relativas de los símbolos del correspondiente alfabeto.

proporciones <- function(vector, alfabeto){
  frecAbs <- table(factor(vector, alfabeto))
  frecAbs/sum(frecAbs)
}

prueba <- sample(rep(c("A", "B"), each = 10))
proporciones(prueba, c("A", "B"))
##   A   B 
## 0.5 0.5 

## Comprueba su funcionamiento con los datos del ejercicio anterior.

proporciones(hinfluenzae, dna)

## (5) Obtener la lista de sub-secuencias codificantes del genoma de
## Haemophilus influenzae Rd KW20 (eligiendo la opción Conding
## Sequences y utilizando, de nuevo las funciones read.fasta y 
## getSequence).

secuencias_CS <- read.fasta("hinfluenzae_CS.fasta")
hinfluenzae_CS <- getSequence(secuencias_CS)

## (6) Calcular las frecuencias absolutas y relativas de los
## nucleótidos para el conjunto de secuencias codificantes.

#Frecuencias relativas
frecs <- sapply(hinfluenzae_CS, proporciones, dna) 
#devuelve para cada fila un nucleotido, y cada columna es una subsecuencia;
#como no queremos eso, hacemos la traspuesta, y solo mostramos unas pocas, en este caso
# las frecuencias con respecto a las tres primeras subsecuencias
t(frecs[,1:3])

#frecuencias absolutas
frecuenciasAbs <- function(vector, alfabeto){
  table(factor(vector, alfabeto))
}

frecsAbs <- sapply(hinfluenzae_CS,frecuenciasAbs,dna)
t(frecsAbs[,1:3])


######################################################################
