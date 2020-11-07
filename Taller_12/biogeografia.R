# Cargar los paquetes necesarios para el analisis
library(GenSA)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

# Archivo que contiene el arbol datado en formato Newick
trfn = "Psychotria_5.2.newick" 

# Distribuciones geograficas de los terminales
geogfn = "Psychotria_geog.data"

# Variable con el numero maximo de areas en las que puede distribuirse un taxa
max_range_size = 4

# modelo DEC (predefinido)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA" # FALSE or "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDEC = bears_optim_run(BioGeoBEARS_run_object)


# modelo DEC+J
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA" # FALSE or "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDECJ = bears_optim_run(BioGeoBEARS_run_object)


# modelo DIVA-like
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA" # FALSE or "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDIVA = bears_optim_run(BioGeoBEARS_run_object)

# modelo DIVA-like+J
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA" # FALSE or "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
dstart = resDIVA$outputs@params_table["d","est"]
estart = resDIVA$outputs@params_table["e","est"]
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDIVAJ = bears_optim_run(BioGeoBEARS_run_object)

# modelo BAYAREA-like
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA" # FALSE or "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resBAYAREA = bears_optim_run(BioGeoBEARS_run_object)

# modelo BAYAREA-like+J
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA" # FALSE or "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
dstart = resBAYAREA$outputs@params_table["d","est"]
estart = resBAYAREA$outputs@params_table["e","est"]
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resBAYAREAJ = bears_optim_run(BioGeoBEARS_run_object)

###########################################################################
# Comparacion de modelos: AIC y valores de verosimilitud
# La comparacion con AIC solo se puede en modelos que tienen un subset 
# de parametros, para los demas casos se debe revisar el valor de la
# verosimilitud  
###########################################################################
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECJ)
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
stats = AICstats_2models(LnL_1, LnL_2, 3, 2)
stats
stats[,c("LnLalt", "LnLnull", "AIC1", "AIC2")]

LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVAJ)
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVA)
stats = AICstats_2models(LnL_1, LnL_2, 3, 2)
stats[,c("LnLalt", "LnLnull", "AIC1", "AIC2")]

LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREAJ)
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREA)
stats = AICstats_2models(LnL_1, LnL_2, 3, 2)
stats[,c("LnLalt", "LnLnull", "AIC1", "AIC2")]

########################################
# Grafica de areas ancestrales
########################################
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
tr = read.tree(trfn)

plot_BioGeoBEARS_results(resDECJ, "DEC+J", addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(resDECJ, "DEC+J", addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


#####################################################################
#
#							Ejercicio
#
# El siguiente ejercicio analítico será el tema para el reporte de
# la sesión. Como en los reportes anteriores, éste debe tener una
# longitud de una página, a excepción de gráficas y bibliografía.
#
# Realizar la selección de modelo biogeográfico para el dataset de
# Aguirre-Santoro et al. (2020. Botanical Journal of the Linnean
# Society 192: 609–624. https://doi.org/10.1093/botlinnean/boz087).
# El análisis se debe realizar sobre un subset de los terminales:
#
# Aechmea_dactylina	
# Aechmea_gracilis	
# Aechmea_napoensis	
# Aechmea_tayoensis	
# Billbergia_viridiflora	
# Bromelia_flemingii	
# Bromelia_karatas	
# Bromelia_trianae	
# Hohenbergia_rosea	
# Orthophytum_gurkenii	
# Ronnbergia_veitchii	
# Wittmackia_limae	
# Wittmackia_urbaniana	
#
# Para ello debe reducir las terminales del árbol filogenético
# (wittmackia.newick) usando la función `keep.tip`. También debe 
# reducir la información de las distribuciones geográficas
# (wittmackia_areas.txt). En esta última debe asegurarse que las 
# áreas no utilizadas en la codificación de la matriz reducida deben 
# ser eliminadas. 
#
# Una vez realizada la selección del modelo, se realizarán los 
# siguientes experimentos:
#
# 1. Efecto de la distribución del grupo externo. En este dataset, 
# las tres especies de Bromelia son el grupo externo. Las tres
# especies tienen rangos de distribución distintos. De esta manera,
# realizará tres análisis, utilizando solo una especie de Bromelia
# en cada uno y observará su influencia en la reconstrucción de los
# rangos ancestrales en el grupo interno.
#
# 2. Efecto del número máximo de áreas conjuntas. La estructura de
# la matriz de evolución geográfica depende del número máximo de
# áreas que puede ocupar un linaje en un periódo determinado. Para
# evaluar el efecto de este parámetro en la estimación biogeográfica,
# se realizará análisis con diferentes valores máximos de áreas
# conjuntas (1, 2, 3 y 4).
#
######################################################################





