## get_data ##
# get and process example data used in opcr
library(opcr)
opc = opc_process('../gosl_habitat/data/raw/2018/opc/OPC033.D00')
save(opc,file='data/opc.rda')
