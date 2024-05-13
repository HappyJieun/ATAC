# 후성 유전체 분석을 위한 ATAC-sequencing 파이프라인 구축

> Period: 2023.09 ~ 2023.12
> 
> Subject: ATAC-sequencing



##### 23.12.20 last edit
---

## 0. Environment

+ Language : R

+ Editor : RStudio
---
## 1. Introduction

**Background**

닭의 배아에서 성별 결정 과정에 중요한 역할을 하는 비대칭성 골난성 발달에 대해 연구하고자 함.

---
## 2. Data Set

**Dataset Info.**

**(1) Samples**

Li J, Sun C, Zheng J, Li J, Yi G, Yang N. Time-Course Transcriptional and Chromatin Accessibility Profiling Reveals Genes Associated With Asymmetrical Gonadal Development in Chicken Embryos. Front Cell Dev Biol. 2022 Mar 8;10:832132. doi: 10.3389/fcell.2022.832132. PMID: 35345851; PMCID: PMC8957256.

**(2) Reference genome**

https://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/

---
## 3. Summary

**(1) Data Preprocessing**

- Trimmomatic, Alignment, Sorting, Duplicated remove, Merge 등 진행
- LHX9 유전자 발현양 비교

  
  ![image](https://github.com/HappyJieun/ATAC-seq/assets/166107244/20a7c2aa-6af2-4ce6-b554-c62b804e569a)


<br/>

**(2) Review**

- 4가지 성장 단계에서 LHX9 유전자 발현 차이 확인
![image](https://github.com/HappyJieun/ATAC-seq/assets/166107244/0a0f63fd-d43a-4eb3-b418-931eadc35a77)

- 유전자 발현과 염색체 접근성 사이의 상관관계 확인
- 성별 결정과 비대칭성 골난성 발달의 연관성 확인
