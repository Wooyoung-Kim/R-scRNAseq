##############################
### 1. 환경 설정 및 라이브러리 로드
##############################

# 라이브러리 설치 경로 지정 (필요에 따라 수정)
.libPaths("/usr/lib/R/library/4.3")

# 데이터 전처리, 시각화, single-cell 분석에 필요한 패키지 로드
library(tidyverse)       # 데이터 조작 및 시각화 (ggplot2, dplyr 등 포함)
library(Seurat)          # Single-cell RNA-seq 분석 패키지
library(DoubletFinder)   # 더블렛(doublet) 검출 도구 (McGinnis et al., 2019)
library(SeuratDisk)      # Seurat 객체 형식 변환 (h5Seurat, h5ad 등)
library(scCustomize)     # Seurat 전용 추가 시각화 기능
library(qs)              # R 객체의 빠른 직렬화/저장을 위한 패키지

##############################
### 2. 10X 데이터 읽기 및 Seurat 객체 생성
##############################

# 각 샘플의 10X Genomics 데이터 로드 (SoupX 전처리된 데이터)
# Chow 다이어트: Sedentary (SD) 및 Training (TR) 그룹
Chow_sd_1 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554964_SoupX/")
Chow_sd_2 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554965_SoupX/")
Chow_sd_3 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554966_SoupX/")

Chow_TR_1 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554967_SoupX/")
Chow_TR_2 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554968_SoupX/")
Chow_TR_3 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554969_SoupX/")
Chow_TR_4 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554970_SoupX/")

# HFD (High-Fat Diet) 다이어트: Sedentary (SD) 및 Training (TR) 그룹
HFD_sd_1 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554971_SoupX/")
HFD_sd_2 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554972_SoupX/")
HFD_sd_3 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554973_SoupX/")

HFD_TR_1 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554974_SoupX/")
HFD_TR_2 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554975_SoupX/")
HFD_TR_3 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554976_SoupX/")
HFD_TR_4 <- Read10X(data.dir = "~/Exercise_Full/scWAT/GSM5554977_SoupX/")

# 각 데이터를 Seurat 객체로 생성 (최소 3개 세포 및 200개 유전자 기준 필터링)
Chow_sd_1 <- CreateSeuratObject(counts = Chow_sd_1, min.cells = 3, min.features = 200, project = "scWAT_Chow_SD_1")
Chow_sd_2 <- CreateSeuratObject(counts = Chow_sd_2, min.cells = 3, min.features = 200, project = "scWAT_Chow_SD_2")
Chow_sd_3 <- CreateSeuratObject(counts = Chow_sd_3, min.cells = 3, min.features = 200, project = "scWAT_Chow_SD_3")

Chow_TR_1 <- CreateSeuratObject(counts = Chow_TR_1, min.cells = 3, min.features = 200, project = "scWAT_Chow_TR_1")
Chow_TR_2 <- CreateSeuratObject(counts = Chow_TR_2, min.cells = 3, min.features = 200, project = "scWAT_Chow_TR_2")
Chow_TR_3 <- CreateSeuratObject(counts = Chow_TR_3, min.cells = 3, min.features = 200, project = "scWAT_Chow_TR_3")
Chow_TR_4 <- CreateSeuratObject(counts = Chow_TR_4, min.cells = 3, min.features = 200, project = "scWAT_Chow_TR_4")

HFD_sd_1 <- CreateSeuratObject(counts = HFD_sd_1, min.cells = 3, min.features = 200, project = "scWAT_HFD_SD_1")
HFD_sd_2 <- CreateSeuratObject(counts = HFD_sd_2, min.cells = 3, min.features = 200, project = "scWAT_HFD_SD_2")
HFD_sd_3 <- CreateSeuratObject(counts = HFD_sd_3, min.cells = 3, min.features = 200, project = "scWAT_HFD_SD_3")

HFD_TR_1 <- CreateSeuratObject(counts = HFD_TR_1, min.cells = 3, min.features = 200, project = "scWAT_HFD_TR_1")
HFD_TR_2 <- CreateSeuratObject(counts = HFD_TR_2, min.cells = 3, min.features = 200, project = "scWAT_HFD_TR_2")
HFD_TR_3 <- CreateSeuratObject(counts = HFD_TR_3, min.cells = 3, min.features = 200, project = "scWAT_HFD_TR_3")
HFD_TR_4 <- CreateSeuratObject(counts = HFD_TR_4, min.cells = 3, min.features = 200, project = "scWAT_HFD_TR_4")

##############################
### 3. 함수 정의: Doublet 검출 및 데이터 처리
##############################

# Function: doublet_finder_analysis
# 설명: PCA 결과를 바탕으로 최적의 차원 수를 결정한 후,
#        UMAP, 이웃 탐색, 클러스터링을 수행하고 DoubletFinder를 이용하여
#        더블렛을 검출합니다.
# 참고: McGinnis et al., 2019, Cell Systems. (PMID: 31068753)
doublet_finder_analysis <- function(seurat_obj, multiplet_rate = NULL) {
  # PCA 표준편차를 이용하여 분산 기여도 계산
  stdv <- seurat_obj[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  
  # 두 기준(누적 분산 90% 초과 및 분산 기여도의 급격한 감소)으로 최적 차원 수 결정
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(which((percent_stdv[1:(length(percent_stdv) - 1)] - percent_stdv[2:length(percent_stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)
  
  # UMAP, 이웃 찾기 및 클러스터링 수행
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:min_pc)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:min_pc)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
  
  # 최적의 pK 값 찾기: 파라미터 sweep을 통해 후보 pK 값 평가
  sweep_list <- paramSweep(seurat_obj, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)
  optimal_pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>% 
    dplyr::select(pK) %>% 
    as.numeric(as.character(.))
  
  # Homotypic doublet 비율 추정
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  
  # Multiplet rate가 제공되지 않은 경우 10x Genomics 기준값을 이용하여 추정
  if (is.null(multiplet_rate)) {
    message('multiplet_rate not provided... estimating from dataset cell count')
    
    multiplet_rates_10x <- data.frame(
      'Multiplet_rate' = c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
    )
    
    multiplet_rate <- multiplet_rates_10x %>% 
      dplyr::filter(Recovered_cells < nrow(seurat_obj@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% 
      dplyr::select(Multiplet_rate) %>% 
      as.numeric(as.character(.))
    
    message(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # 예상 더블렛 수 계산 (Homotypic 보정 적용)
  nExp_poi <- round(multiplet_rate * nrow(seurat_obj@meta.data))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  
  # DoubletFinder 실행
  seurat_obj <- doubletFinder(seu = seurat_obj, 
                              PCs = 1:min_pc, 
                              pK = optimal_pk, 
                              nExp = nExp_poi_adj)
  
  # DoubletFinder 결과 컬럼명을 "doublet_finder"로 변경
  colnames(seurat_obj@meta.data)[grepl('DF.classifications.*', colnames(seurat_obj@meta.data))] <- "doublet_finder"
  
  return(seurat_obj)
}

# Function: processing
# 설명: 각 Seurat 객체에 대해 기본 전처리(메타토큰 계산, 정규화, 변수 유전자 검출, 스케일링, PCA 수행)
#        후, 더블렛 검출을 실행하여 Singlet만 남깁니다.
processing <- function(seurat_obj) {
  # 미토콘드리아 유전자 비율 계산 (패턴: "^mt-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # (필요시) 추가 필터링: 유전자 수, UMI 수, 미토콘드리아 비율 기준 적용
  # seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
  
  # 데이터 정규화 및 변수 유전자 탐색
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 모든 유전자에 대해 스케일링 진행
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  # PCA 수행 (변수 유전자 사용)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  # 더블렛 검출 수행 후, Singlet만 선택
  seurat_obj <- doublet_finder_analysis(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = doublet_finder == "Singlet")
  
  return(seurat_obj)
}

##############################
### 4. 각 샘플 전처리 수행
##############################

# Chow 다이어트 Sedentary 그룹
Chow_sd_1 <- processing(Chow_sd_1)
Chow_sd_2 <- processing(Chow_sd_2)
Chow_sd_3 <- processing(Chow_sd_3)

# Chow 다이어트 Training 그룹
Chow_TR_1 <- processing(Chow_TR_1)
Chow_TR_2 <- processing(Chow_TR_2)
Chow_TR_3 <- processing(Chow_TR_3)
Chow_TR_4 <- processing(Chow_TR_4)

# HFD 다이어트 Sedentary 그룹
HFD_sd_1 <- processing(HFD_sd_1)
HFD_sd_2 <- processing(HFD_sd_2)
HFD_sd_3 <- processing(HFD_sd_3)

# HFD 다이어트 Training 그룹
HFD_TR_1 <- processing(HFD_TR_1)
HFD_TR_2 <- processing(HFD_TR_2)
HFD_TR_3 <- processing(HFD_TR_3)
HFD_TR_4 <- processing(HFD_TR_4)

##############################
### 5. 샘플 통합 및 메타데이터 추가
##############################

# 동일 조건 내에서 샘플 병합 (add.cell.ids를 통해 셀 ID 구분)
Chow_SD <- merge(x = Chow_sd_1, y = c(Chow_sd_2, Chow_sd_3), 
                 add.cell.ids = c("SkM_Chow_SD_1", "SkM_Chow_SD_2", "SkM_Chow_SD_3"))
Chow_TR <- merge(x = Chow_TR_1, y = c(Chow_TR_2, Chow_TR_3, Chow_TR_4), 
                 add.cell.ids = c("SkM_Chow_TR_1", "SkM_Chow_TR_2", "SkM_Chow_TR_3", "SkM_Chow_TR_4"))
HFD_SD <- merge(x = HFD_sd_1, y = c(HFD_sd_2, HFD_sd_3), 
                add.cell.ids = c("SkM_HFD_SD_1", "SkM_HFD_SD_2", "SkM_HFD_SD_3"))
HFD_TR <- merge(x = HFD_TR_1, y = c(HFD_TR_2, HFD_TR_3, HFD_TR_4), 
                add.cell.ids = c("SkM_HFD_TR_1", "SkM_HFD_TR_2", "SkM_HFD_TR_3", "SkM_HFD_TR_4"))

# 각 그룹에 조건, 조직, 다이어트, 운동 여부 메타데이터 추가
Chow_SD@meta.data$Condition <- "Chow_SD"
Chow_TR@meta.data$Condition <- "Chow_TR"
HFD_SD@meta.data$Condition <- "HFD_SD"
HFD_TR@meta.data$Condition <- "HFD_TR"

Chow_SD@meta.data$tissue <- "SkM"
Chow_TR@meta.data$tissue <- "SkM"
HFD_SD@meta.data$tissue <- "SkM"
HFD_TR@meta.data$tissue <- "SkM"

Chow_SD@meta.data$diet <- "Chow"
Chow_TR@meta.data$diet <- "Chow"
HFD_SD@meta.data$diet <- "HFD"
HFD_TR@meta.data$diet <- "HFD"

Chow_SD@meta.data$exercise <- "Sedentary"
Chow_TR@meta.data$exercise <- "Training"
HFD_SD@meta.data$exercise <- "Sedentary"
HFD_TR@meta.data$exercise <- "Training"

# 전체 데이터셋(scWAT)으로 병합 후 레이어 통합 (JoinLayers)
scWAT <- merge(x = Chow_SD, y = c(Chow_TR, HFD_SD, HFD_TR))
scWAT <- JoinLayers(scWAT)

# 샘플 구분을 위해 기본 식별자 지정
Idents(scWAT) <- "orig.ident"

# 그래픽 장치 초기화 (필요시)
dev.off()

##############################
### 6. 추가 필터링 및 최종 저장
##############################

# QC 기준: 유전자 수, UMI 수, 미토콘드리아 비율 필터 적용
scWAT <- subset(scWAT, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5 & nCount_RNA > 200)

# 특정 유전자 (Mki67) 발현이 0인 세포만 선택 (세포 증식 마커 제거)
Mki67_expression <- FetchData(scWAT, vars = "Mki67")
scWAT <- subset(scWAT, subset = Mki67 == 0)

# 최종 Seurat 객체 저장 (RDS 형식)
saveRDS(scWAT, "~/Exercise_Full/RDS/SkM.RDS")

##############################
### 참고 및 대체 처리 코드 (주석 처리됨)
##############################
# 아래 코드는 필요에 따라 추가 전처리나 파일 포맷 변환을 위해 사용될 수 있습니다.
# scWAT <- JoinLayers(scWAT)
# scWAT <- NormalizeData(scWAT)
# scWAT <- FindVariableFeatures(scWAT, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(scWAT)
# scWAT <- ScaleData(scWAT, features = all.genes)
# scWAT <- RunPCA(scWAT, features = VariableFeatures(object = scWAT))
# a <- nrow(scWAT@meta.data)
# scWAT <- doublet_finder_analysis(scWAT)
# scWAT <- subset(scWAT, subset = doublet_finder == "Singlet")
# nrow(scWAT@meta.data)
# a
# scWAT <- JoinLayers(scWAT)
# scWAT <- Convert_Assay(seurat_object = scWAT, convert_to = "V3")
# SaveH5Seurat(scWAT, overwrite = TRUE, verbose = TRUE, filename ="~/Exercise_Full/scWAT.h5Seurat")
# Convert("~/Exercise_Full/scWAT.h5Seurat", dest = "h5ad")
