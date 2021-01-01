## get_data ##
# get and process example data used in opcr
library(opcr)
opc = opc_process('../gosl_habitat/data/raw/2019b/opc/OPC027.D00')
save(opc,file='data/opc.rda')
