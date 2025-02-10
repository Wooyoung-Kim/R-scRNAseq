# -----------------------------------------------------------
# 환경 설정 및 라이브러리 로드
# -----------------------------------------------------------
# R 라이브러리 경로 지정
.libPaths("/usr/lib/R/library/4.3")

# 필수 패키지 로드
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)

# future 패키지의 전역 변수 메모리 제한 증가 (대규모 데이터 분석 대비)
options(future.globals.maxSize = 1e10)

# -----------------------------------------------------------
# 1. 데이터 불러오기 및 병합
# -----------------------------------------------------------
# 각 조직(조직별 RDS 파일)에서 Seurat 객체 불러오기
vWAT <- readRDS("~/Exercise_Full/RDS/vWAT.RDS")
scWAT <- readRDS("~/Exercise_Full/RDS/scWAT.RDS")
SkM   <- readRDS("~/Exercise_Full/RDS/SkM.RDS")

# 서로 다른 조직 데이터를 하나의 Seurat 객체(tissue)로 병합
tissue <- merge(x = vWAT, y = c(scWAT, SkM))

# 사용한 개별 객체 메모리 정리
rm(vWAT, scWAT, SkM)

# -----------------------------------------------------------
# 2. RNA 어세이 준비 및 조건별 분할
# -----------------------------------------------------------
# (사용자 정의 함수로 추정) RNA 레이어들을 하나로 결합
tissue <- JoinLayers(tissue, assay = "RNA")

# RNA 어세이를 'Condition' 메타데이터에 따라 분할
tissue[["RNA"]] <- split(tissue[["RNA"]], f = tissue$Condition)

# -----------------------------------------------------------
# 3. 전처리 및 비통합 분석 (Pre-integration)
# -----------------------------------------------------------
# 데이터 정규화, 고변이 유전자 선택, 스케일링 수행
tissue <- NormalizeData(tissue)
tissue <- FindVariableFeatures(tissue)
tissue <- ScaleData(tissue)

# PCA를 통해 차원 축소 수행
tissue <- RunPCA(tissue)

# PCA 결과를 기반으로 이웃 그래프 생성
tissue <- FindNeighbors(tissue, dims = 1:30, reduction = "pca")

# 해상도 0.5로 클러스터링 (비통합 클러스터)
tissue <- FindClusters(tissue, resolution = 0.5, cluster.name = "unintegrated_clusters")

# PCA 기반 UMAP 계산 (비통합 데이터 시각화)
tissue <- RunUMAP(tissue, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# -----------------------------------------------------------
# 4. 비통합 데이터 UMAP 시각화
# -----------------------------------------------------------
DimPlot(tissue, reduction = "umap.unintegrated", group.by = "Scibet")
DimPlot(tissue, reduction = "umap.unintegrated", group.by = "scType")
DimPlot(tissue, reduction = "umap.unintegrated", group.by = "scMayomap")
DimPlot(tissue, reduction = "umap.unintegrated", group.by = "Condition")

# -----------------------------------------------------------
# 5. SCTransform 및 Harmony를 이용한 통합 분석
# -----------------------------------------------------------
# SCTransform을 통해 분산 안정화 및 정규화 수행
tissue <- SCTransform(tissue)

# SCTransform 결과로 다시 PCA 수행
tissue <- RunPCA(tissue)

# PCA 결과를 바탕으로 UMAP 계산 (SCT 기반)
tissue <- RunUMAP(tissue, dims = 1:30)

# Harmony 통합 (SCT 정규화된 데이터를 기반으로)
tissue <- IntegrateLayers(object = tissue, 
                          method = HarmonyIntegration, 
                          normalization.method = "SCT", 
                          verbose = FALSE, 
                          new.reduction = "harmony")

# 통합된 (Harmony) 데이터를 기반으로 이웃 그래프 생성
tissue <- FindNeighbors(tissue, reduction = "harmony", dims = 1:30)

# 해상도 0.5로 클러스터링 (통합 클러스터)
tissue <- FindClusters(tissue, resolution = 0.5, cluster.name = "integrated_clusters")

# Harmony 결과를 이용하여 UMAP 계산 (통합 데이터 시각화)
tissue <- RunUMAP(tissue, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

# 통합 데이터 UMAP: Condition 및 Scibet에 따른 그룹 시각화
DimPlot(tissue, reduction = "umap", group.by = c("Condition", "Scibet"))

# 현재 Seurat 객체를 파일로 저장
qsave(tissue, "~/Exercise_Full/RDS/Tissue.qs")

# -----------------------------------------------------------
# 6. Reference Transfer: 레퍼런스 데이터를 통한 주석 전파
# -----------------------------------------------------------
# 레퍼런스 Seurat 객체(단일세포 아틀라스) 불러오기 및 업데이트
ref <- qread("~/Organ_cross_2/WAT_Muscle/GSE183288_Single_cell_atlas.qs")
ref <- UpdateSeuratObject(ref)

# tissue 객체에 대해 RNA 레이어 다시 결합 및 전처리
tissue <- JoinLayers(tissue, assay = "RNA")
tissue <- SetIdent(tissue, value = "integrated_clusters")
tissue <- NormalizeData(tissue)
DefaultAssay(tissue) <- "RNA"

# 레퍼런스와 query(tissue) 간의 전이 앵커(anchor) 찾기 (PCA 기준)
tissue.anchors <- FindTransferAnchors(reference = ref, query = tissue, dims = 1:30,
                                      reference.reduction = "pca", query.assay = "RNA")

# 레퍼런스의 cell_state_label을 이용하여 query에 레이블 전이
predictions <- TransferData(anchorset = tissue.anchors, 
                            refdata = ref$cell_state_label, 
                            dims = 1:30)
# 예측된 cell state 레이블을 metadata에 추가
tissue <- AddMetaData(tissue, metadata = predictions$predicted.id)

# 레퍼런스의 cell_type_label을 이용하여 query에 세포 유형 전이
predictions_2 <- TransferData(anchorset = tissue.anchors, 
                              refdata = ref$cell_type_label, 
                              dims = 1:30)
# cell_type 예측 결과를 별도 변수(celltype)로 지정
predictions_2$celltype <- predictions_2$predicted.id
tissue <- AddMetaData(tissue, metadata = predictions_2$celltype, col.name = "celltype")

# 통합 데이터(UAMP: Harmony)에서 예측된 cell state 레이블 시각화 (라벨 표시)
DimPlot(tissue, reduction = "umap.harmony", group.by = "predicted.id", label = TRUE)

# 업데이트된 tissue 객체 저장
qsave(tissue, "~/Exercise_Full/RDS/Tissue.qs")

# 필요 시 객체 재불러오기
tissue <- qread("~/Exercise_Full/RDS/Tissue.qs")

# -----------------------------------------------------------
# 7. Metadata 가공 및 세포 유형 분류
# -----------------------------------------------------------
# 분석에 필요한 메타데이터 열만 선택
tissue@meta.data <- tissue@meta.data %>% 
  select(orig.ident, nCount_RNA, nFeature_RNA, Condition, tissue, diet, exercise, 
         seurat_clusters, Scibet, scType, scMayomap, unintegrated_clusters, 
         integrated_clusters, predicted.id, celltype)

# 선택된 celltype의 고유 값 확인
print(unique(tissue@meta.data$celltype))

# 비면역 세포(non-immune)로 분류할 세포 유형 목록 정의
non_immune <- c("ASC", "Fibroblast", "EC", "Smooth_Muscle", "Glial", 
                "Epi", "FAP", "Epididymis", "Satellite", "FAP_Sca1-", 
                "Tenocyte", "Muscle_Fiber")

# celltype에 따라 각 세포를 "non-immune" 또는 "immune"으로 분류하여 새 메타데이터 열 추가
tissue@meta.data$Type <- ifelse(tissue@meta.data$celltype %in% non_immune, 
                                "non-immune", "immune")
