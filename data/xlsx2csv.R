
data_file = "ivsc_3cmtct_shed3_param.xlsx"
file_name = "ivsc_3cmtct_shed3_param.csv"
d = xlsx::read.xlsx(data_file, 1)
write.csv(d, file= file_name, row.names=FALSE, na="NA")
