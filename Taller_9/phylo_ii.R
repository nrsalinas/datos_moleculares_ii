###############################################
# Información filogenética en R, segunda parte
###############################################

library(ape)
library(phangorn)

#########################
# Inferencia filogenética 
# Métodos de distancia
#########################

matkal = read.FASTA("matk_al.fst")

# matriz de distancias genéticas
matkdist = dist.ml(matkal)

# Arbol de Neighboor Joining
njtree = NJ(matkdist)

# Graficando el árbol
plot(njtree)

# No está enraizado
is.rooted(njtree)

# Si solo se tiene un grupo externo, se puede enraizar utilizardo dicho nombre
njtree = root(njtree, outgroup="Meliosma herbertii", resolve.root=TRUE)

# Pero en este caso hay varios grupos externos, así que es mejor seleccionar el nodo en la gráfica
plot(njtree)
njtree = root(njtree, interactive=TRUE, resolve.root=TRUE)
# se pincha la base de las astéridas


########	Ejercicio
# Inferir el árbol de NJ para el alineamiento de rbcL

#########################
# Inferencia filogenética 
# Parsimonia
#########################

# Conversion to phyDat object of phangorn
matkphydat = phyDat(matkal)
partree = optim.parsimony(njtree, matkphydat)

# adicionalmente se puede utilizar el algoritmo Ratchet
# dicha función tiene muchos más parámetros, pero todos opcionales
partree = pratchet(matkphydat)

########	Ejercicio
# Inferir el árbol de parsimonia para el alineamiento de rbcL

####################################
# Manipulación de árboles
####################################

# Lectura y escritura de árboles en formato Newick
write.tree(njtree, "arbol_nj.newick")
tr = read.tree("arbol_nj.newick")

#Inspección de las terminales
tr$tip.label


# Graficar los árboles
plot(tr)

plot(tr, no.margin=TRUE)

plot(tr, no.margin=TRUE, underscore=TRUE)

plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001)

plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001, cex=0.7)

tr = ladderize(tr)

####### Ejercicio
# Inserte un título y un subtítulo a la gráfica
# con la función `title`


#####################################
# Comparar árboles filogenéticos
#####################################
# Suponiendo que t1 y t2 son los árboles a comparar
# y tienen idénticas terminales

t1 = rtree(10) # arbol aleatorio de 10 terminales
t2 = rtree(10) # arbol aleatorio de 10 terminales

astab = cbind(t1$tip.label, t1$tip.label)
cophyloplot(t1, t2, assoc = astab, length.line = 5, space = 40, gap = 2)

# Guardar gráficas de los árboles en archivos

pdf("arbol.pdf")
plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001, cex=0.7)
dev.off()

png("arbol.png", 600, 600)
plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001, cex=0.7)
dev.off()


#################   Ejercicio
# Graficar la comparación topológica entre los árboles de NJ de matK y rbcL y guardarla como una gráfica

