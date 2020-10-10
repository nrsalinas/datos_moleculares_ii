#################################
# Información filogenética en R
#################################

install.packages("ape")
library(ape)

# Info del programa: https://cran.r-project.org/web/packages/ape/index.html
# Allí se encontrará código, página web oficial del proyect y manual

mis.accs = c("AF404817", "U48604", "GU176622", "U48612", "GU176623", "AF519552", "AF519556", "GU176625", "GU176626", "U48613", "GU176627", "AY520786", "AY520804", "U48609", "GU176629", "U48582", "U48611", "U48599", "U48581", "U48601", "U48600", "U48597", "U48608", "GU176630", "U48607", "U48606", "U48605", "GU176631", "AF393440", "GU176632", "GU176633", "GU176634", "GU176635", "GU176636", "GU176637", "AF382741", "AF106811", "AF133752")

mis.seqs = read.GenBank(mis.accs)
names(mis.seqs) # Los nombres de las secuencias son los códigos de acceso 
attr(mis.seqs, "species") # Pero los nombres taxonómicos están en el atributo "species"


new.names = paste(attr(mis.seqs, "species"), names(mis.seqs), sep="_")
names(mis.seqs) = new.names
write.FASTA(mis.seqs, "test.fasta")

##########		Ejercicio
# Descargar el archivo `ficus.csv`, leerlo como una DataFrame y descargar las sequencias del gen ETS

#######################
# Alineamientos
#######################

mialg = clustal(mis.seqs)

# También se pueden visualizar los alineamientos
x = mialg[, 1:50]
image(x)

x = mialg[, length(mialg[1,]) - (50: 0)]
image(x)

# Para editarlos (cortar inicios/finales del alineamiento
short = mialg[, 1:(length(mialg[1,]) - 46)]
x = short[, length(short[1,]) - (50: 0)]
image(x)

# Finalmente se pueden guardar en formato fasta
write.FASTA(short, "alin.fasta")

##########		Ejercicio
# Descargar todas las secuencias de la tabla ficus.csv y alinear cada secuencia por separado. Guardar los alineamientos en un archivo fasta.

##############################################
# Seleccion de modelos para inferencia por ML
##############################################

library(phangorn) # Si no está preinstalada, instalar

itsdat = phyDat(itsal) # Alineamiento de ape tiene que ser convertido a objeto phyDat de phangorn

itsmt = modelTest(itsdat)

# Ordenar deacuerdo a AIC para seleccionar el mejor modelo
itsmt[order(itsmt$AICc),]

#########	Ejercicio
# Seleccionar el mejor modelo para cada partición del dataset de Ficus

# Concatenación
concat = cbind(itsal, etsal, fcbal, fill.with.gaps=TRUE, check.names=TRUE)

# Guardar alineamiento a formato phylip
write.phyDat(midat, "align.phy", "phylip")

#################	Ejercicio 
# Concatenar los alineamientos y guardar la matrix para un análisis en RAxML


#########################
# Arboles filogeneticos
#########################

# matriz de distancias genéticas
itsdist = dist.ml(itsdat)

# Arbol de Neighboor Joining
njtree = NJ(itsdist)

# Reenraizando el árbol con el grupo externo
njtree = root(njtree, outgroup="Ficus_ingens", resolve.root=TRUE)

########	Ejercicio
# Inferir arboles de NJ para todos los genes del dataset de Ficus

# Lectura y escritura de árboles en formato Newick
write.tree(njtree, "arbol_nj.newick")
tr = read.tree("arbol_nj.newick")

#Inspección de las terminales
tr$tip.label


# Graficar los árboles
plot(tr)

plot(tr, no.margin=TRUE)

plot(tr, no.margin=TRUE, underscore=TRUE)

tr = drop.tip(tr, "Ficus_ottoniifolia_2")

plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001)

plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001, cex=0.7)

tr = ladderize(tr)

# Comparar árboles filogenéticos
# suponiendo que t1 y t2 son los árboles a comparar
# y tienen idénticas terminales

at = cbind(t1$tip.label, t1$tip.label)
cophyloplot(t1, t2, assoc = ass, length.line = 5, space = 40, gap = 2)

# Guardar gráficas de los árboles en archivos
pdf("arbol.pdf")
plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001, cex=0.7)
dev.off()

png("arbol.png", 600, 600)
plot(tr, no.margin=TRUE, underscore=TRUE, label.offset=0.001, cex=0.7)
dev.off()




