####################################################
#
#	Diversificación de linajes en R
#
####################################################

library(geiger)
library(phytools)

# Simulacion de proceso de Yule: no extincion
t0 = pbtree(0.3, 0, n=50)
plot(t0)
axisPhylo()


# Simular un arbol bird-death, lambda = 0.3, mu = 0.1, 50 terminales  
t0 = pbtree(0.3, 0.1, n=50))
plot(t0)
axisPhylo()


# Interaccion variables de linajes presentes y unidades de tiempo. 
# El criterio de interrupcion de la simulacion debe senalar alguno de aquellas variables
t0 = pbtree(0.3, 0.1, t=100))
plot(t0)
axisPhylo()


# ¿Como es el arbol sin linajes extintos
t0p = drop.extinct(t0)
plot(t0p)



# LTT plots - linages a traves del tiempo
ltt.plot(t0)

# Usualmente se transforma el eje y a logaritmo natural
ltt.plot(t0, log="y")

# se puede anadir una linea recta como patron esperado
# las coordenadas deben ser verificadas siempre
lines(c(-10,0), c(1,50), lty=2)


# Tambien se pueden graficar los LTT plots de varios arboles al tiempo
t1 = pbtree(0.3, 0, n=50)
mltt.plot(t0, t1, log="y")

# generacion de distribuciones nulas
treesYule = pbtree(0.3, 0, n=50, nsim=100)
obj = ltt(treesYule, plot=FALSE)
plot(obj, col="grey")
ltt(t0, add=TRUE, col="blue")

# Comportamiento modelos y datos empiricos
# Arbol filogenetico de percas (Etheostoma y Percina)
darters = read.tree("darters.newick")
plot(darters, show.tip.label=FALSE, type='fan')

# Comparar LTT con patron esperado para proceso Yule
ltt.plot(darters, log="y")
lines(c(-25,0), c(1,200), lty=2, col='blue')

# Grafica LTT observado vs distribucion nula Yule
obsLambda = phytools:::qb(darters) # calculo rapido de especiacion observada
trees = pbtree(b=obsLambda, n=Ntip(darters), t=max(nodeHeights(darters)), nsim=100, method='direct')
plot(obj, col="grey")
ltt(darters, add=TRUE, col="blue")


# Ajuste de modelos
yule.model = fit.yule(darters)

# Incorporando el submuestreo en el ajuste del modelo
samp = 201/216
yule.model = fit.yule(darters, rho=samp)
bd.model = fit.bd(darters, rho=samp)

# Seleccionando el modelo que mejor se ajusta a los datos
AIC(yule.model,bd.model)

# Detectando cambios en la diversificación de la filogenia

md = medusa(darters, criterion="aic", model = "yule")
plot(md, show.tip.label=FALSE)

# en el objeto `md` está la información de los cambios en
# las tasas de diversificación. La función `nodelabels`
# puede usarse para ubicar los cambios

nodelabels("Cambio", 309)

#######################################################
#
#					Ejercicio
#
# Explorar la dinámica de diversificación del 
# género Anolis. Para ello se utilizarán la filogenia
# calibrada del grupo presentada por Poe et al. (2017,
# accesible en https://doi.org/10.1093/sysbio/syx029).
# Los resultados deben ser presentados en el miniartículo
# de la sesión.
#
# El miniartículo de esta semana debe tener una longitud
# máxima de una página tamaño carta, excluyendo bibliografía 
# y gráficas de apoyo
#
##########################################################


#cargar datos de Anolis
anoles = read.nexus("anolis_tree.nexus")
plot(anoles)

# Es necesario eliminar outgroups
outs = c("Bplumifrons","Pmarmoratus","Pscapulatus","Ugallardoi")
anoles = drop.tip(anoles, outs)
