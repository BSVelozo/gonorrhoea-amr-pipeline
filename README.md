Gonorrhoea AMR Typer

Pipeline Nextflow para tipagem molecular e detecção de resistência antimicrobiana em Neisseria gonorrhoeae a partir de dados Illumina paired‑end (FASTQ).

> Leve, rápido e sem montagem de genoma.
Utiliza mapeamento por BWA + SAMtools para identificar alelos de MLST, NG‑MAST e NG‑STAR.

Visão geral
Este pipeline foi desenvolvido para substituir abordagens pesadas (montagem de novo com SPAdes ou BLAST em contigs) por uma estratégia assembly‑free baseada em mapeamento direto das reads contra bancos de alelos curados.

 O que ele faz
- MLST – determina o Sequence Type (ST) a partir dos 7 genes housekeeping: `abcZ`, `adk`, `aroE`, `fumC`, `gdh`, `pdhC`, `pgm`. Suporta partial matching (não exige todos os 7 genes presentes).
- NG‑MAST – identifica alelos de `porB` e `tbpB` e atribui o ST correspondente.
- NG‑STAR – perfila alelos de resistência nos genes: `penA`, `mtrR`, `porB`, `ponA`, `gyrA`, `parC`, `23S`, e retorna o NG‑STAR ST.
- Controle de qualidade (opcional) – FastQC antes e depois da limpeza.
- Limpeza de reads (opcional) – CutAdapt em 3 níveis de rigor: leve, médio e rigoroso.
- Verificação de espécie (opcional) – mapeamento contra genoma de referência de N. gonorrhoeae.
- Saída consolidada – tabela binária `Yes`/`No`/`NA` por gene, mais os STs por tipagem.

 Por que ele é mais rápido e leve?
- Mapeamento com BWA‑MEM (rápido e preciso).
- Processamento paralelo por amostra com Nextflow.
- Bancos de alelos combinados em um único FASTA por tipagem (menos arquivos, menos I/O).
- Sem montagem de genoma – tudo direto das reads brutas.
- Cache do Nextflow – reexecuções reaproveitam etapas concluídas.

Requisitos

Ferramenta 
Versão mínima
Função 
Nextflow 
24.10.3
Orquestração
BWA 
0.7.17
Mapeamento 
SAMtools 
1.21
Manipulação de BAM
FastQC 
0.12.1
Controle de Qualidade (Opcional)
CutAdapt
4.6
Limpeza de Reads (Opcional)
Python 
3.10
Scripts auxiliares
bc
qualquer
Cálculos em shell
Conda 
qualquer
Ambiente (Opcional, mas recomendado)


Dependências Python: `pandas`, `biopython`, `numpy` (instaladas via Conda).

Estrutura do projeto

gonorrhoea-amr-pipeline/
├── main.nf                         Script principal do Nextflow
├── nextflow.config                 Configurações do executor
├── environment.yml                 Ambiente Conda
├── scripts/
│   ├── prepare_db.sh               Prepara bancos pyngoST
│   ├── prepare_pyngoST_db.py       Constrói FASTA combinados e perfis CSV
│   ├── extract_alleles.py          Extrai alelos e determina ST (partial matching)
│   ├── generate_binary_table.py    Gera tabela Yes/No/NA por amostra
│   ├── generate_coverage_report.py  Relatório de cobertura por gene
│   └── build_final_report.py       Consolida tabela final de todas as amostras
├── db/
│   ├── pyngoST/                    (fornecido pelo usuário) arquivos .fas e .tab
│   ├── mlst/                       (gerado) FASTA + profile MLST
│   ├── ng-mast/                    (gerado) FASTA + profile NG‑MAST
│   ├── ng-star/                    (gerado) FASTA + profile NG‑STAR
│   └── reference/                  (opcional) genoma de referência
├── data/                           (fornecido pelo usuário) FASTQs _{1,2}.fastq
└── results/                        (criado automaticamente) saída do pipeline
```

Instalação
1. Clonar o repositório

git clone https://github.com/seu-usuario/gonorrhoea-amr-pipeline.git
cd gonorrhoea-amr-pipeline

2. Criar o ambiente Conda
conda env create -f environment.yml
conda activate gonorrhoea-amr

3. Dar permissão de execução aos scripts

chmod +x scripts/.sh scripts/.py

 4. Corrigir quebras de linha (se você baixou de Windows)

sed -i 's/\r$//' scripts/.py scripts/.sh main.nf

Preparação dos bancos de dados
O pipeline utiliza o banco pré‑construído do pyngoST (alelos e perfis para MLST, NG‑MAST e NG‑STAR). Coloque todos os arquivos `.fas` e `.tab` em `db/pyngoST/`:

db/pyngoST/
├── abcZ.fas        ├── porB.fas         ├── MLST_profiles.tab
├── adk.fas         ├── TBPB.fas         ├── NGMAST_profiles.tab
├── aroE.fas        ├── penA.fas         └── NGSTAR_profiles.tab
├── fumC.fas        ├── 23S.fas
├── gdh.fas         ├── gyrA.fas
├── pdhC.fas        ├── parC.fas
└── pgm.fas         ├── ponA.fas
                    └── mtrR.fas

Executar a preparação

./scripts/prepare_db.sh

O script irá:

- Combinar todos os alelos de cada tipagem em um único FASTA (`mlst_combined.fa`, `ngmast_combined.fa`, `ngstar_combined.fa`).
- Converter os perfis `.tab` para `.csv` com cabeçalhos apropriados.
- Colocar os arquivos prontos em `db/mlst/`, `db/ng-mast/` e `db/ng-star/`.

 (Opcional) Baixar genoma de referência
Para habilitar a verificação de espécie:

mkdir -p db/reference
wget -O db/reference/neisseria_gonorrhoeae.fasta.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/845/GCF_000006845.1_ASM684v1/GCF_000006845.1_ASM684v1_genomic.fna.gz"
gunzip db/reference/neisseria_gonorrhoeae.fasta.gz

Execução do pipeline
 Comando básico (sem QC nem trimagem):
nextflow run main.nf --reads 'data/_{1,2}.fastq' --results_dir resultados

Com limpeza de reads e QC:
nextflow run main.nf \
    --reads 'data/_{1,2}.fastq' \
    --run_trimming true --trim_level medium \
    --run_fastqc true \
    --results_dir resultados

 Com parâmetros personalizados:
nextflow run main.nf \
    --reads 'data/_{1,2}.fastq' \
    --results_dir meus_resultados \
    --min_coverage 80 --min_depth 15

 Retomar execução interrompida
nextflow run main.nf --reads 'data/_{1,2}.fastq' --results_dir resultados -resume

Parâmetros
Entrada e saída
| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--reads` | Padrão dos arquivos FASTQ paired‑end (com wildcard) | `'data/_{1,2}.fastq'` |
| `--results_dir` | Diretório de saída | `'results'` |

Critérios de detecção
| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--min_coverage` | Cobertura mínima do gene para considerar presente (%) | `70` |
| `--min_depth` | Profundidade média mínima | `10` |

Bancos de dados
| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--mlst_fasta` | FASTA combinado MLST | `'db/mlst/mlst_combined.fa'` |
| `--mlst_profile` | Perfil MLST (CSV) | `'db/mlst/profiles.csv'` |
| `--ngmast_fasta` | FASTA combinado NG‑MAST | `'db/ng-mast/ngmast_combined.fa'` |
| `--ngmast_profile` | Perfil NG‑MAST | `'db/ng-mast/profiles.csv'` |
| `--ngstar_fasta` | FASTA combinado NG‑STAR | `'db/ng-star/ngstar_combined.fa'` |
| `--ngstar_profile` | Perfil NG‑STAR | `'db/ng-star/profiles.csv'` |
| `--ref_genome` | Genoma de referência (para verificação de espécie) | `'db/reference/neisseria_gonorrhoeae.fasta'` |

Qualidade e limpeza (opcionais)
| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--run_trimming` | Ativar limpeza com CutAdapt | `false` |
| `--trim_level` | Nível de rigor: `light`, `medium` ou `strict` | `'medium'` |
| `--run_fastqc` | Executar FastQC antes e depois da limpeza | `false` |

Parâmetros de trimagem por nível
| Nível | Qualidade (`-q`) | Comprimento mínimo (`-m`) | Uso recomendado |
|-------|-----------------|--------------------------|-----------------|
| `light` | 20 | 30 | Reads boas, remoção leve |
| `medium` | 25 | 50 | Recomendado |
| `strict` | 30 | 70 | Reads ruins ou com contaminação |

Saída gerada
Após a execução, a estrutura será:
resultados/
├── SRR10861750/                        uma pasta por amostra
│   ├── qc/
│   │   ├── raw/                        FastQC antes da limpeza
│   │   └── trimmed/                    FastQC depois da limpeza
│   ├── trimmed/                        FASTQs limpos (se trimming ativo)
│   ├── species_check/
│   │   ├── SRR10861750.species.flagstat
│   │   └── SRR10861750.species.stats.tsv   % mapeado + espécie
│   ├── mlst/
│   │   ├── bam/
│   │   │   ├── SRR10861750.bam
│   │   │   ├── SRR10861750.bam.bai
│   │   │   ├── SRR10861750.sam
│   │   │   ├── SRR10861750.flagstat
│   │   │   └── SRR10861750.depth.txt
│   │   ├── SRR10861750.mlst.alleles.tsv    alelos + cobertura + profundidade
│   │   └── SRR10861750.mlst.summary.tsv    ST + matches=X/Y
│   ├── ngmast/
│   │   └── ... (mesma estrutura)
│   ├── ngstar/
│   │   └── ... (mesma estrutura)
│   ├── SRR10861750.binary.tsv          tabela Yes/No/NA
│   └── SRR10861750.combined_report.tsv  cobertura consolidada
├── final_binary_table.tsv              ← TABELA FINAL (todas as amostras)
└── final_per_sample/
    ├── SRR10861750.binary.tsv
    └── `final_binary_table.tsv`

Tabela consolidada, uma linha por amostra, colunas por gene com `Yes`/`No`/`NA`, mais os STs:
| Sample | 23s | abcz | adk | gyra | mtrr | porb | tbpb | mlst_ST | ngmast_ST | ngstar_ST |
|--------|-----|------|-----|------|------|------|------|---------|-----------|-----------|
| SRR10861750 | Yes | Yes | Yes | Yes | Yes | No | Yes | 1901 | 32 | 45 |

Significado dos valores:
- `Yes` → alelo detectado com cobertura ≥ `min_coverage` e profundidade ≥ `min_depth`.
- `No` → gene presente no banco, mas nenhum alelo passou os limiares.
- `NA` → gene não se aplica (ex: `clonal_complex`, `cc`).

Colunas de ST:
- Valor numérico (ex: `1901`) → ST determinado com base no perfil.
- `Unknown` → nenhum perfil bateu; consulte o `.summary.tsv` correspondente e a coluna `MatchInfo` para ver quantos genes coincidiram.

 `.mlst.alleles.tsv` (exemplo)
| Gene | Allele | Coverage(%) | MeanDepth |
|------|--------|-------------|-----------|
| abcz | abcZ_126 | 98.5 | 42.3 |
| adk  | adk_481  | 100.0 | 55.1 |
| fumc | N/A      | 0.0   | 0.00 |

- Coverage(%) → percentual do gene coberto por pelo menos 1 read.
- MeanDepth → profundidade média (nº médio de reads por posição).

Fluxo de trabalho
FASTQ bruto
    │
    ├── (opcional) FastQC raw
    │
    ├── (opcional) CutAdapt (light/medium/strict)
    │       └── (opcional) FastQC trimmed
    │
    ├── (opcional) Verificação de espécie (BWA contra genoma ref)
    │
    ├── MAP_MLST   ──► BAM + BAI + SAM + depth + flagstat
    ├── MAP_NGMAST ──► BAM + BAI + SAM + depth + flagstat
    ├── MAP_NGSTAR ──► BAM + BAI + SAM + depth + flagstat
    │
    ├── EXTRACT_MLST   ──► .mlst.alleles.tsv + .mlst.summary.tsv
    ├── EXTRACT_NGMAST ──► .ngmast.alleles.tsv + .ngmast.summary.tsv
    ├── EXTRACT_NGSTAR ──► .ngstar.alleles.tsv + .ngstar.summary.tsv
    │
    └── FINAL_REPORT ──► final_binary_table.tsv + final_per_sample/.binary.tsv

Citação
Se utilizar este pipeline, cite:
- Nextflow: Di Tommaso et al. (2017) Nat. Biotechnol.
- BWA: Li & Durbin (2009) Bioinformatics
- SAMtools: Danecek et al. (2021) GigaScience
- pyngoST (banco de dados): repositório oficial em [https://github.com/sanger-bentley-group/pyngoST](https://github.com/sanger-bentley-group/pyngoST)
- NG‑STAR: [https://ngstar.canada.ca/](https://ngstar.canada.ca/)
- PubMLST: [https://pubmlst.org/neisseria/](https://pubmlst.org/neisseria/)

—
Desenvolvido para vigilância epidemiológica e pesquisa em resistência antimicrobiana de Neisseria gonorrhoeae


