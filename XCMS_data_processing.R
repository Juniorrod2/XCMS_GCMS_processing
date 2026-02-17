library(xcms)
library(Spectra)
library(MsExperiment)

register(SnowParam(4))

#Loading files

mzFiles <- list.files("C:/Users/Edmil/Desktop/Thermo Raw files",pattern = "mzXML")


anno_df <- data.frame(sample_name=sub(x=basename(mzFiles),pattern = ".mzXML",replacement = ""),
                      group=c("A","A","A","B","B"))
#reading spectra

raw_data <- readMsExperiment(spectraFiles = paste("C:/Users/Edmil/Desktop/Thermo Raw files/",mzFiles,sep = ""),
                             sampleData = anno_df)

#Extracting and ploting the chromatogram

chrom <- chromatogram(raw_data)

plot(chrom)

spectra(data)

chr_raw <- chromatogram(raw_data, mz = c(65,72), rt = c(20.2*60,28*60))
plot(chr_raw)


raw_data |>
  filterRt(rt = c(27.2*60,28*60)) |>
  filterMz(mz = c(65,72)) |>
  plot(type = "XIC")


#Peak Detection


mfp <- MatchedFilterParam(
  fwhm = 30,          # largura média do pico cromatográfico (segundos)
  max = 5E^6,         # intensidade máxima esperada
  snthresh = 10,      # threshold sinal/ruído (GC-MS geralmente ≥ 5)
)

data <- findChromPeaks(raw_data, param = mfp)

chromPeaks(data)|>View()

