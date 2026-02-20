library(xcms)
library(Spectra)
library(MsExperiment)

register(SnowParam(8))

#Loading files

files_location = "GC spectra/"
mzFiles <- list.files(files_location,pattern = "mzXML")


spectra_metadata <- data.frame(sample_name=sub(x=basename(mzFiles),pattern = ".mzXML",replacement = ""),
                      group=rep(c("CD","UC"),each=10))
#reading spectra

raw_data <- readMsExperiment(spectraFiles = paste(files_location,mzFiles,sep = ""),
                             sampleData = spectra_metadata)

sampleData(raw_data)

#Extracting and ploting the chromatogram

full_chromatogram <- chromatogram(raw_data,aggregationFun = "max")

group_color_pattern <- rep(c("red","blue"),each=10)

plot(full_chromatogram,col=group_color_pattern)

# Cut region without signal

raw_data <- filterRt(raw_data,c(700,2200))

full_chromatogram <- chromatogram(raw_data,aggregationFun = "max")

plot(full_chromatogram,col=group_color_pattern)


total_ion_count <- spectra(raw_data) |>
  tic() |>
  split(f = fromFile(raw_data))

boxplot(total_ion_count, col = group_color_pattern,
        ylab = "intensity", main = "Total ion current")



#Peak Detection

#Extraction of a single ion peak (EIC) for inspection
EIC_RT <- c(1645,1655)
EIC_mz <- c(80,85)

EIC_data <- chromatogram(raw_data, rt=EIC_RT,mz=EIC_mz)

plot(EIC_data,col=group_color_pattern)

filterRt(raw_data,rt = EIC_RT) |>
filterMz(mz = EIC_mz) |>
plot(type = "XIC")


#Testing the peak detection of the selected EIC

MF_parameters <- MatchedFilterParam(
  fwhm = 20,          # largura média do pico cromatográfico (segundos)
  max = 5,         # intensidade máxima esperada
  snthresh = 5,      # threshold sinal/ruído (GC-MS geralmente ≥ 5)
)

detected_peaks <- findChromPeaks(EIC_data, param = MF_parameters)

chromPeaks(detected_peaks)|>View()

plot(detected_peaks,col=group_color_pattern)

#Peak Detection algorithm

detected_peaks <- findChromPeaks(raw_data, param = MF_parameters)

chromPeaks(detected_peaks)|>View()

# Peak alignment

aligned_peaks <- adjustRtime(detected_peaks,param = ObiwarpParam(binSize = 0.6))

full_chromatogram <- chromatogram(aligned_peaks,aggregationFun = "max")

plot(full_chromatogram,col=group_color_pattern)

# Checking alignment on test peak

EIC_RT <- c(1645,1655)
EIC_mz <- c(80,85)

EIC_data <- chromatogram(aligned_peaks,rt=EIC_RT,mz=EIC_mz)

plot(EIC_data,col=group_color_pattern)

#Peak correspondence

pdp <- PeakDensityParam(sampleGroups = sampleData(aligned_peaks)$group,
                        minFraction = 0.4, bw = 30)

grouped_data <- groupChromPeaks(aligned_peaks, param = pdp)


#Gap filling

gap_filled_values <- fillChromPeaks(grouped_data, param = ChromPeakAreaParam())

#Final feature matrix

Feature_matrix <- featureValues(gap_filled_values, value = "into")

features_metadata <- featureDefinitions(gap_filled_values)
