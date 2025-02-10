##############################################################
# 1. 라이브러리 및 경로 설정
##############################################################
# 라이브러리 경로 설정 (R 버전 4.3 라이브러리 경로 지정)
.libPaths("/usr/lib/R/library/4.3")

# 필수 라이브러리 로드
library(dplyr)
library(Seurat)
library(tidyverse)
library(scibet)

##############################################################
# 2. Seurat 객체 불러오기
##############################################################
# 저장된 RDS 파일에서 각 조직의 Seurat 객체 불러오기
vWAT  <- readRDS("~/Exercise_Full/RDS/vWAT.RDS")
scWAT <- readRDS("~/Exercise_Full/RDS/scWAT.RDS")
SkM   <- readRDS("~/Exercise_Full/RDS/SkM.RDS")

##############################################################
# 3. Seurat 객체 전처리 함수 정의 (sce)
##############################################################
# sce() 함수는 데이터를 정규화, 변수 유전자 선정, 스케일링, PCA, 
# 최근접 이웃 찾기, 클러스터링, UMAP 임베딩을 수행합니다.
sce <- function(x) {
  sce <- NormalizeData(x, verbose = FALSE)                # 데이터 정규화
  sce <- FindVariableFeatures(sce, verbose = FALSE)         # 변수 유전자 선택
  sce <- ScaleData(sce, verbose = FALSE)                    # 데이터 스케일링
  sce <- RunPCA(sce, verbose = FALSE)                       # PCA 수행
  sce <- FindNeighbors(sce, dims = 1:30, verbose = FALSE)   # 최근접 이웃 찾기 (1~30 PC 사용)
  sce <- FindClusters(sce, resolution = 0.5, verbose = FALSE)  # 클러스터링 (해상도 0.5)
  sce <- RunUMAP(sce, dims = 1:30, verbose = FALSE)         # UMAP 임베딩
  return(sce)
}

# RNA assay가 기본으로 사용되도록 layers 결합 (JoinLayers 함수 사용)
vWAT  <- JoinLayers(vWAT, assay = "RNA")
scWAT <- JoinLayers(scWAT, assay = "RNA")
SkM   <- JoinLayers(SkM, assay = "RNA")

# 각 객체 전처리 수행
vWAT  <- sce(vWAT)
scWAT <- sce(scWAT)
SkM   <- sce(SkM)

##############################################################
# 4. Scibet를 이용한 세포 타입 어노테이션
##############################################################
# --- 4-1. vWAT 객체 ---
# RNA 데이터 및 메타데이터 추출 후 결합
vWAT_input <- vWAT@assays$RNA$data %>% t() %>% as.data.frame()
vWAT_meta  <- vWAT@meta.data %>% select(seurat_clusters)
colnames(vWAT_meta) <- "label"
vWAT_input <- cbind(vWAT_input, vWAT_meta)

# --- 4-2. scWAT 객체 ---
scWAT_input <- scWAT@assays$RNA$data %>% t() %>% as.data.frame()
scWAT_meta  <- scWAT@meta.data %>% select(seurat_clusters)
colnames(scWAT_meta) <- "label"
scWAT_input <- cbind(scWAT_input, scWAT_meta)

# --- 4-3. SkM 객체 ---
SkM_input <- SkM@assays$RNA$data %>% t() %>% as.data.frame()
SkM_meta  <- SkM@meta.data %>% select(seurat_clusters)
colnames(SkM_meta) <- "label"
SkM_input <- cbind(SkM_input, SkM_meta)

# --- 모델 불러오기 및 프로세싱 ---
# Fat 조직 모델과 Muscle 모델을 CSV 파일에서 불러옴
model_fat    <- readr::read_csv("~/Organ_cross_2/Scibet_reference/GSE109774_Fat_Smart-seq2_scibet_core.csv")
model_Muscle <- readr::read_csv("~/Organ_cross_2/Scibet_reference/GSE109774_Limb_Muscle_10x_scibet_core.csv")

# 모델 전처리 (pro.core 함수 사용)
model_fat    <- pro.core(model_fat)
model_Muscle <- pro.core(model_Muscle)

# --- vWAT에 scibet 어노테이션 적용 ---
ori_label_WAT <- vWAT_input$label
vWAT_input    <- vWAT_input[,-ncol(vWAT_input)]  # 마지막 열(라벨) 제거
prd_vWAT      <- LoadModel(model_fat)
label_vWAT    <- prd_vWAT(vWAT_input)
vWAT@meta.data$Scibet <- label_vWAT

# --- scWAT에 scibet 어노테이션 적용 ---
ori_label_WAT <- scWAT_input$label
scWAT_input   <- scWAT_input[,-ncol(scWAT_input)]
prd_scWAT     <- LoadModel(model_fat)
label_scWAT   <- prd_scWAT(scWAT_input)
scWAT@meta.data$Scibet <- label_scWAT

# --- SkM에 scibet 어노테이션 적용 ---
ori_label_skm <- SkM_input$label
SkM_input     <- SkM_input[,-ncol(SkM_input)]
prd_SkM       <- LoadModel(model_Muscle)
label_SkM     <- prd_SkM(SkM_input)
SkM@meta.data$Scibet <- label_SkM

##############################################################
# 5. scType를 이용한 세포 타입 어노테이션 (Adipose/Muscle)
##############################################################
# scType 어노테이션에 필요한 라이브러리 및 함수 로드
lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = TRUE)

# gene set 준비 함수 및 어노테이션 함수 GitHub에서 소스 코드로 로드
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

#### 5-1. scType 어노테이션: scWAT (Adipose tissue) ####
# DB 파일 및 Tissue 설정
db_    <- "~/Exercise_Full/RDS/ScTypeDB_Adipose.xlsx"
Tissue <- "Adipose"  # 예: Immune system, Pancreas, Liver 등

# gene set 준비
gs_list <- gene_sets_prepare(db_, Tissue)

# Seurat 객체 버전 확인 (v4/v5에 따라 count 데이터 접근 방법이 다름)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(scWAT[["RNA"]])))
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# 스케일링된 RNA 데이터 추출
scRNAseqData_scaled <- if (seurat_package_v5) {
  as.matrix(scWAT[["RNA"]]$counts)
} else {
  as.matrix(scWAT[["RNA"]]@counts)
}

# sctype_score 함수 실행: positive 및 negative gene sets 사용
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, 
  scaled = TRUE, 
  gs = gs_list$gs_positive, 
  gs2 = gs_list$gs_negative
)

# 각 클러스터별 점수 집계 및 상위 cell type 도출
cL_resutls <- do.call("rbind", lapply(unique(scWAT@meta.data$seurat_clusters), function(cl){
  es.max.cl <- sort(
    rowSums(es.max[, rownames(scWAT@meta.data[scWAT@meta.data$seurat_clusters == cl, ])]),
    decreasing = TRUE
  )
  head(data.frame(
    cluster = cl, 
    type    = names(es.max.cl), 
    scores  = es.max.cl, 
    ncells  = sum(scWAT@meta.data$seurat_clusters == cl)
  ), 10)
}))

# 각 클러스터별 최고 점수 cell type 선택
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# 낮은 점수의 클러스터는 "Unknown"으로 지정 (임계값: ncells/4)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[, 1:3])

# scType 어노테이션을 meta.data에 저장
scWAT@meta.data$scType <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  scWAT@meta.data$scType[scWAT@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# UMAP 시각화: scibet와 scType 어노테이션 비교
DimPlot(scWAT, reduction = "umap", group.by = "Scibet")
DimPlot(scWAT, reduction = "umap", group.by = "scType")

#### 5-2. scType 어노테이션: vWAT (Adipose tissue) ####
# Seurat 객체 버전 확인
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(vWAT[["RNA"]])))
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# 스케일링된 RNA 데이터 추출
scRNAseqData_scaled <- if (seurat_package_v5) {
  as.matrix(vWAT[["RNA"]]$counts)
} else {
  as.matrix(vWAT[["RNA"]]@counts)
}

# sctype_score 함수 실행
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, 
  scaled = TRUE, 
  gs = gs_list$gs_positive, 
  gs2 = gs_list$gs_negative
)

# 클러스터별 점수 집계
cL_resutls <- do.call("rbind", lapply(unique(vWAT@meta.data$seurat_clusters), function(cl){
  es.max.cl <- sort(
    rowSums(es.max[, rownames(vWAT@meta.data[vWAT@meta.data$seurat_clusters == cl, ])]),
    decreasing = TRUE
  )
  head(data.frame(
    cluster = cl, 
    type    = names(es.max.cl), 
    scores  = es.max.cl, 
    ncells  = sum(vWAT@meta.data$seurat_clusters == cl)
  ), 10)
}))

# 각 클러스터별 최고 점수 cell type 선택
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# 낮은 점수 클러스터 "Unknown" 지정
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[, 1:3])

# scType 어노테이션 저장
vWAT@meta.data$scType <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  vWAT@meta.data$scType[vWAT@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# UMAP 시각화: scibet와 scType 어노테이션 비교
DimPlot(vWAT, reduction = "umap", group.by = "Scibet")
DimPlot(vWAT, reduction = "umap", group.by = "scType")

#### 5-3. scType 어노테이션: SkM (Muscle tissue) ####
# DB 파일 및 Tissue 설정 (Muscle)
db_    <- "~/Exercise_Full/RDS/ScTypeDB_Muscle.xlsx"
Tissue <- "Muscle"

# gene set 준비
gs_list <- gene_sets_prepare(db_, Tissue)

# Seurat 객체 버전 확인
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(SkM[["RNA"]])))
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# 스케일링된 RNA 데이터 추출
scRNAseqData_scaled <- if (seurat_package_v5) {
  as.matrix(SkM[["RNA"]]$counts)
} else {
  as.matrix(SkM[["RNA"]]@counts)
}

# sctype_score 함수 실행
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, 
  scaled = TRUE, 
  gs = gs_list$gs_positive, 
  gs2 = gs_list$gs_negative
)

# 클러스터별 점수 집계
cL_resutls <- do.call("rbind", lapply(unique(SkM@meta.data$seurat_clusters), function(cl){
  es.max.cl <- sort(
    rowSums(es.max[, rownames(SkM@meta.data[SkM@meta.data$seurat_clusters == cl, ])]),
    decreasing = TRUE
  )
  head(data.frame(
    cluster = cl, 
    type    = names(es.max.cl), 
    scores  = es.max.cl, 
    ncells  = sum(SkM@meta.data$seurat_clusters == cl)
  ), 10)
}))

# 각 클러스터별 최고 점수 cell type 선택
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# 낮은 점수 클러스터 "Unknown" 지정
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[, 1:3])

# scType 어노테이션 저장
SkM@meta.data$scType <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  SkM@meta.data$scType[SkM@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# UMAP 시각화: scibet와 scType 어노테이션 비교
DimPlot(SkM, reduction = "umap", group.by = "Scibet")
DimPlot(SkM, reduction = "umap", group.by = "scType")

##############################################################
# 6. scMayoMap를 이용한 추가 세포 타입 어노테이션
##############################################################
# scMayoMap 관련 라이브러리 로드
library(scMayoMap)
library(MAST)

#### 6-1. vWAT: scMayoMap 어노테이션 ####
# vWAT의 마커 유전자 탐색 (MAST 방법 사용)
seurat.markers_vWAT <- FindAllMarkers(vWAT, method = 'MAST')

# scMayoMap 객체 생성 (데이터베이스: scMayoMapDatabase, tissue: 'adipose tissue')
scMayoMap.obj <- scMayoMap(
  data = seurat.markers_vWAT, 
  database = scMayoMapDatabase, 
  tissue = 'adipose tissue'
)
plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)

# 수동으로 클러스터 아이덴티티 재설정
res_1 <- c("Pericyte",
           "ADSc",
           "ADSc",
           "ADSc",
           "T cells",
           "NK cells",
           "Platelet",
           "Adipocyte",
           "Mast cells",
           "Dendritic cells",
           "Epithelial cells",
           "Epithelial cells",
           "T cells",
           "B cells",
           "B cells",
           "Endothelial cells",
           "memory T cells",
           "Endothelial cells",
           "Dendritic cells",
           "memory T cells",
           "Neuron",
           "Neuron",
           "Dendritic cells",
           "Mast cells",
           "Basophils")

# 클러스터 레벨 확인 및 정렬 (현재 factor levels를 숫자형으로 정렬 후 재설정)
current_levels <- levels(vWAT)
sorted_levels <- sort(as.numeric(current_levels))
levels(vWAT) <- as.character(sorted_levels)

# 클러스터별 이름 매핑 후 아이덴티티 재설정
names(res_1) <- levels(vWAT)
vWAT <- RenameIdents(vWAT, res_1)
vWAT@meta.data$scMayomap <- vWAT@active.ident

# UMAP 시각화 (scMayomap과 scType 비교)
DimPlot(vWAT, reduction = "umap", label = TRUE, group.by = "scMayomap")
DimPlot(vWAT, reduction = "umap", label = TRUE, group.by = "scType")

#### 6-2. scWAT: scMayoMap 어노테이션 ####
# scWAT의 마커 유전자 탐색
seurat.markers_scWAT <- FindAllMarkers(scWAT, method = 'MAST')

# scMayoMap 객체 생성 (tissue: 'adipose tissue')
scMayoMap.obj <- scMayoMap(
  data = seurat.markers_scWAT, 
  database = scMayoMapDatabase, 
  tissue = 'adipose tissue'
)
plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)

# 수동으로 클러스터 아이덴티티 재설정
res_1 <- c("memory T cells",
           "memory T cells",
           "B cells",
           "Epithelial cells",
           "B cells",
           "Mast cells",
           "Adipocyte",
           "memory T cells",
           "memory T cells",
           "NK cells",
           "B cells",
           "Macrophages",
           "Dendritic cells",
           "memory T cells",
           "B cells",
           "Dendritic cells",
           "Neuron",
           "B cells",
           "Dendritic cells",
           "Dendritic cells")

# 클러스터 레벨 정렬 및 재설정
current_levels <- levels(scWAT)
sorted_levels <- sort(as.numeric(current_levels))
levels(scWAT) <- as.character(sorted_levels)
names(res_1) <- levels(scWAT)
scWAT <- RenameIdents(scWAT, res_1)
scWAT@meta.data$scMayomap <- scWAT@active.ident

# UMAP 시각화 (scMayomap과 scType 비교)
DimPlot(scWAT, reduction = "umap", label = TRUE, group.by = "scMayomap")
DimPlot(scWAT, reduction = "umap", label = TRUE, group.by = "scType")

#### 6-3. SkM: scMayoMap 어노테이션 ####
# SkM의 마커 유전자 탐색
seurat.markers_SkM <- FindAllMarkers(SkM, method = 'MAST')

# scMayoMap 객체 생성 (tissue: 'muscle')
scMayoMap.obj <- scMayoMap(
  data = seurat.markers_SkM, 
  database = scMayoMapDatabase, 
  tissue = 'muscle'
)
plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)

# 수동으로 클러스터 아이덴티티 재설정
res_1 <- c("Fibroblast",
           "Fibroblast",
           "Fibroblast",
           "B cells",
           "Fibroblast",
           "Fibroblast",
           "Schwann cells",
           "Endothelial cells",
           "Fibroblast",
           "Smooth muscle cells",
           "Myocytes",
           "Neutrophils",
           "Schwann cells",
           "T cells",
           "T cells",
           "Myocytes",
           "Schwann cells",
           "Fibroblast",
           "Schwann cells",
           "MSCs",
           "Neutrophils",
           "Macrophages",
           "Macrophages",
           "Monocytes",
           "B cells",
           "Myocytes",
           "Satellite cells",
           "Fibroblast",
           "Macrophages",
           "Schwann cells")

# 클러스터 레벨 정렬 및 재설정
current_levels <- levels(SkM)
sorted_levels <- sort(as.numeric(current_levels))
levels(SkM) <- as.character(sorted_levels)
names(res_1) <- levels(SkM)
SkM <- RenameIdents(SkM, res_1)
SkM@meta.data$scMayomap <- SkM@active.ident

# UMAP 시각화 (scMayomap과 scType 비교)
DimPlot(SkM, reduction = "umap", label = TRUE, group.by = "scMayomap")
DimPlot(SkM, reduction = "umap", label = TRUE, group.by = "scType")

##############################################################
# 7. 최종 Seurat 객체 저장
##############################################################
saveRDS(vWAT,  "~/Exercise_Full/RDS/vWAT.RDS")
saveRDS(scWAT, "~/Exercise_Full/RDS/scWAT.RDS")
saveRDS(SkM,   "~/Exercise_Full/RDS/SkM.RDS")
