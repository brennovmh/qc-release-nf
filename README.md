# QC Release NF

Pipeline em Nextflow para controle de qualidade de dados brutos de sequenciamento antes da análise secundária.

## Visão geral

Este pipeline foi desenvolvido para apoiar a triagem inicial de corridas de NGS, automatizando a geração de métricas de qualidade e relatórios visuais.

A versão atual inclui:

- FastQC
- fastp
- MultiQC
- classificação automática das amostras em PASS, REVIEW ou FAIL
- relatório visual em HTML
- relatório visual em PDF

## Objetivo

O objetivo do pipeline é auxiliar a decisão de liberação de amostras para análise secundária, reunindo em um único fluxo informações técnicas de qualidade e um resumo interpretável da corrida.

## Estrutura do pipeline

O pipeline executa as seguintes etapas:

1. FastQC para avaliação inicial dos FASTQs
2. fastp para métricas adicionais de qualidade
3. QC_GATE para classificação automática das amostras
4. QC_VISUAL_REPORT para geração do relatório HTML da corrida
5. QC_REPORT_PDF para geração do PDF da corrida
6. MultiQC para consolidação técnica dos resultados

## Requisitos

Antes de executar o pipeline, é necessário ter instalado:

- Nextflow
- Docker
- Git

## Estrutura esperada do projeto

```bash
qc_release_nf/
├── main.nf
├── nextflow.config
├── README.md
├── samplesheet.example.csv
├── docker/
│   └── weasyprint/
│       └── Dockerfile
