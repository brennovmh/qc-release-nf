# QC_CHECK

Pipeline em Nextflow para controle de qualidade de dados de sequenciamento (NGS) antes da análise secundária.

---

## Visão geral

Este pipeline foi desenvolvido para automatizar a triagem inicial de corridas de sequenciamento, integrando métricas de qualidade e gerando relatórios interpretáveis para apoio à decisão de liberação de amostras.

A versão atual inclui:

- Avaliação de qualidade com FastQC
- Métricas adicionais com fastp
- Consolidação técnica com MultiQC
- Classificação automática das amostras (PASS, REVIEW, FAIL)
- Relatório visual da corrida (HTML e PDF)

---

## Objetivo

O pipeline tem como objetivo fornecer um fluxo reprodutível e padronizado para avaliação da qualidade de dados brutos, permitindo identificar rapidamente amostras aptas ou não para análise downstream.

---

## Arquitetura do pipeline

O pipeline segue as etapas abaixo:

1. FastQC  
   Avaliação inicial da qualidade dos arquivos FASTQ  

2. fastp  
   Geração de métricas adicionais, como taxa de Q30 e conteúdo de adaptadores  

3. QC_GATE  
   Classificação automática das amostras com base em critérios definidos  

4. QC_VISUAL_REPORT  
   Geração de relatório visual consolidado em HTML  

5. QC_REPORT_PDF  
   Conversão do relatório HTML em PDF  

6. MultiQC  
   Consolidação dos resultados técnicos em um relatório único  

---

## Requisitos

Certifique-se de ter os seguintes requisitos instalados:

- Nextflow >= 23  
- Docker  
- Git  

---

## Estrutura do projeto

```bash
qc_release_nf/
├── main.nf
├── nextflow.config
├── README.md
├── samplesheet.example.csv
├── docker/
│   └── weasyprint/
│       └── Dockerfile
```

---

## Formato do arquivo de entrada

O pipeline utiliza um arquivo CSV com o seguinte formato:

```csv
sample,fastq_1,fastq_2
sample1,/path/sample1_R1.fastq.gz,/path/sample1_R2.fastq.gz
sample2,/path/sample2_R1.fastq.gz,/path/sample2_R2.fastq.gz
```

---

## Instalação

### Clonar o repositório

```bash
git clone git@github.com:brennovmh/qc-release-nf.git
cd qc-release-nf
```

---

## Preparação do ambiente

### Construir o container para geração de PDF

```bash
docker build -t qc-weasyprint:1.0 docker/weasyprint
```

---

## Execução do pipeline

### Execução básica

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

### Execução com relatórios do Nextflow

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  -with-report results/execution_report.html \
  -with-timeline results/execution_timeline.html \
  -with-trace results/execution_trace.txt
```

### Reexecução usando cache

```bash
nextflow run main.nf -resume
```

---

## Saídas geradas

Após a execução, os principais resultados estarão em:

- `results/fastqc/`  
- `results/fastp/`  
- `results/summary/qc_summary.tsv`  
- `results/report/run_qc_report.html`  
- `results/report/run_qc_report.pdf`  
- `results/multiqc/multiqc_report.html`  

---

## Interpretação dos resultados

As amostras são classificadas automaticamente em:

**PASS**  
Amostra apta para análise secundária  

**REVIEW**  
Amostra que requer avaliação manual  

**FAIL**  
Amostra não apta para análise secundária  

---

## Critérios de avaliação

A classificação é baseada em métricas como:

- número total de reads  
- taxa de bases com qualidade Q30  
- percentual de adaptadores  
- módulos FAIL do FastQC  
- módulos WARN do FastQC  

Os limiares podem ser ajustados conforme o tipo de ensaio.

---

## Relatórios

### Relatório principal

Arquivo: `results/report/run_qc_report.html`

Contém:

- resumo da corrida  
- distribuição de PASS, REVIEW e FAIL  
- tabela consolidada por amostra  
- interpretação automática  

### Relatório em PDF

Arquivo: `results/report/run_qc_report.pdf`

Versão estática do relatório HTML, adequada para compartilhamento e documentação.

### MultiQC

Arquivo: `results/multiqc/multiqc_report.html`

Relatório técnico detalhado das métricas geradas pelas ferramentas.

---

## Boas práticas

- Sempre revisar amostras classificadas como REVIEW  
- Verificar consistência com MultiQC  
- Ajustar parâmetros conforme tipo de biblioteca  
- Evitar decisões automatizadas sem validação técnica  

---

## Problemas comuns

**MultiQC não encontrado**  
Garantir uso do profile docker  

**Erro na geração de PDF**  
Verificar se a imagem qc-weasyprint foi construída corretamente  

---

## Limitações

- Pipeline focado em FASTQ paired-end  
- Critérios de QC ainda genéricos  
- Não substitui avaliação técnica especializada  

---

## Próximos desenvolvimentos

- Relatórios individuais por amostra  
- Integração com LIMS  
- Ajuste dinâmico de thresholds  
- Interface gráfica  
- Integração com pipelines downstream  

---

## Contribuição

Contribuições são bem-vindas. Sugestões, melhorias e correções podem ser enviadas via pull request.

---

## Licença

Este projeto está licenciado sob os termos da licença MIT

---

## Autor

Brenno Martins  
Bioinformata
