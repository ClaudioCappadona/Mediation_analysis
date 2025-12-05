data_1$mode_delivery <- ifelse(data_1$mode_delivery %in% c("3", "4"), "3", as.character(data_1$mode_delivery))
data_1$mode_delivery <- factor(data_1$mode_delivery)

data_1$HDP_bin <- as.factor(ifelse(data_1$HDP == 0, 0, 1))
data_1$parity_m_bin <- as.factor(ifelse(data_1$parity_m == 0, 0, 1))
