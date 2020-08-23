#!/usr/bin/Rscript 

traits=read.csv ("tpg2plantgenome2015080073-sup-0003.csv")

tr="tuber_eye_depth";
write.csv (traits [, c("Name", tr)], paste0("gwaspoly-trait-", tr, ".csv"), quote=F, row.names=F)

tr="total_yield"
write.csv (traits [, c("Name", tr)], paste0("gwaspoly-trait-", tr, ".csv"), quote=F, row.names=F)

tr="tuber_shape"
write.csv (traits [, c("Name", tr)], paste0("gwaspoly-trait-", tr, ".csv"), quote=F, row.names=F)

tr="tuber_length"
write.csv (traits [, c("Name", tr)], paste0("gwaspoly-trait-", tr, ".csv"), quote=F, row.names=F)
	
