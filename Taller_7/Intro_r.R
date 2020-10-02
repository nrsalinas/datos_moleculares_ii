#############################
# Introducción sucinta a R
#############################

#################
# Preliminares 
#################
# Esto es un comentario, no será interpretado por R
# Comentarios están precedidos por una `#`
# `help`: ayuda sobre funciones

help(help)
help(sum)


###########################
# Tipos de datos básicos
###########################
#  `logical`, `character` y `float`

help(typeof)
typeof(2)
typeof("Hola") # "Hola" es equivalente a 'Hola'
typeof(TRUE)

############################################
# Operadores y matemática básica */+- ++ 
###########################################

2 + 2
3 ** 6
9 ** 0.5 # atajo de raíz cuadrada

##############
# Expresiones
##############
# Comandos sintacticamente correctos que producen un resultado

2 + 2
sum(2,3,5)
"e" == "r"

######################
#	Variables
######################
# Declaración de variables y palabras clave
# El resultado de cualquier expresión puede ser guardado en una variable

c = 1
a = 1 + 2
s = sum(2,3,4,5)
?reserved

############################
# 	Estructuras de datos
############################

### Vectores
# Conjunto de datos. Son accesados a través del operador de acceso `[]`

av = c(1, 2, TRUE, "Hola")
av[3] 	# TRUE
av = c(av, 4)	# Añade `4` al final del vector `av` 

### Listas
# Mini bases de datos. Conjunto de datos organizados en pares "llave" - "valor". Las llaves deben ser únicas, los valores no necesariamente. 

mi.lista = list(nombre="Fred", esposa="Wilma", no.vastagos=3, edad.vastagos=c(4,7,9))
mi.lista[[1]]
mi.lista$nombre


### DataFrames
# Tablas con multiples funcionalidades

ta = read.csv("seed_plants.csv") # Funciona para csv estricto
?read.table
head(ta)
rownames(ta)
colnames(ta)
ta[1,3]
mean(ta$Occurrence)

#########	Ejercicio	#########
# Leer el archivo seed_plant_sp.csv. Primero examinarlo en un editor de texto y luego leer la ayuda de la función `read.table`

##################
#	Funciones
#################
# variable que ejecuta una serie de comandos, recibe inputs y produce un output

cir.area = function( diam ){	# declaración de la función
	rad = diam / 2				#--
	ar = (rad ** 2) * 3.14159	# Cuerpo de la función
	return (ar)					#--
}

#########	Ejercicio	#########
# Completar la siguiente función para el cálculo del volúmen de un cono
# parámetros de entrada: `diam` y `alt`

vol.cone = function( ){
	# Completar código
}

####################
# Control de flujo
####################
# Tomar decisiones durante la ejecución del código

if (0 == 1){
	print('yes')
} else {
	print("no")
}

if (0 == 1){
	print('yes')
} else if(0 == 2) {
	print("maybe")
} else {
	print('nooo')
}

for (i in 1:10) {
	print(i)
}


##############
# Paquetes
##############
# Conjunto de datos y funciones distribuidos por terceros

install.packages("ape")
library(ape)

# Info del programa: https://cran.r-project.org/web/packages/ape/index.html
# Allí se encontrará código, página web oficial del proyect y manual

mis.accs = c("AF404817", "U48604", "GU176622", "U48612", "GU176623", "AF519552", "AF519556", "GU176625", "GU176626", "U48613", "GU176627", "AY520786", "AY520804", "U48609", "GU176629", "U48582", "U48611", "U48599", "U48581", "U48601", "U48600", "U48597", "U48608", "GU176630", "U48607", "U48606", "U48605", "GU176631", "AF393440", "GU176632", "GU176633", "GU176634", "GU176635", "GU176636", "GU176637", "AF382741", "AF106811", "AF133752")

mis.seqs = read.GenBank(mis.accs)
names(mis.seqs)
attr(mis.seqs, "species")


new.names = paste(attr(mis.seqs, "species"), names(mis.seqs), sep="_")
names(mis.seqs) = new.names
write.FASTA(mis.seqs, "test.fasta")

##########		Ejercicio
# Descargar el archivo `ficus.csv`, leerlo como una DataFrame y descargar las sequencias del gen ETS


