<!--
================================================================================
POINT-IN-TIME HISTORICAL RECORD
Date: 2026-08-30
Historical branch: (feat)-siepe-ready-project
Historical HEAD: 0013c6d24d90070ee297ba6af7d12d125ebd1732
Authoritative current state: NO
Preserved for provenance.
Current decisions must be read from documentation/SCIENTIFIC_STATE.md and the
current documentation/decisions/ ADRs.
The body below is the original document, preserved without semantic rewriting.
================================================================================
-->

# Auditoria técnica somente‑leitura — Fase A do projeto pAgo

**Branch:** `(feat)-siepe-ready-project`
**HEAD:** `0013c6d24d90070ee297ba6af7d12d125ebd1732`
**Data da auditoria:** 2026‑08‑30
**Escopo:** reconstrução factual do estado final da Fase A. Nenhuma alteração foi feita no repositório.

---

## Convenção de rastreabilidade

Cada afirmação relevante recebe um selo:

| Selo | Significado |
|---|---|
| `[CÓDIGO]` | **CONFIRMADO NO CÓDIGO** — verificado lendo o fonte no HEAD atual. |
| `[VERSIONADO]` | **CONFIRMADO EM ARTEFATO VERSIONADO** — verificado num arquivo rastreado pelo Git (fixture, teste, doc, `.gitignore`, `.gitattributes`). |
| `[EXEC-LOCAL]` | **REPORTADO DA EXECUÇÃO LOCAL** — lido de manifestos/CSVs em `data/…` que existem no disco de trabalho mas **não estão versionados** (produzidos pela execução real no PyCharm). |
| `[INFERÊNCIA]` | Dedução do auditor a partir das evidências acima. |

Distinção central desta auditoria: **todo o código e os testes da Fase A estão versionados; nenhum artefato numérico da execução real (52.473 registros) está versionado.** Os números da execução completa vêm de manifestos locais não rastreados.

---

# 1. Estado do Git

## 1.1 Identificação

| Item | Valor | Selo |
|---|---|---|
| Branch atual | `(feat)-siepe-ready-project` | `[CÓDIGO]` |
| HEAD | `0013c6d24d90070ee297ba6af7d12d125ebd1732` | `[CÓDIGO]` |
| Autor do HEAD | Daniel Harley `<danielharleygoes@outlook.com>` | `[CÓDIGO]` |
| Data do HEAD | 2026‑08‑30 22:26:53 ‑0300 (= 2026‑08‑31T01:26:53Z) | `[CÓDIGO]` |
| `origin/(feat)-siepe-ready-project` | `0013c6d…` (idêntico ao local — branch publicada e sincronizada) | `[CÓDIGO]` |
| `master` / `origin/master` | `9d1e55a5e3819b9bbc943c464263ba4846b375e2` (iguais) | `[CÓDIGO]` |
| `merge-base master HEAD` | `9d1e55a…` (master é ancestral direto; branch 0 atrás, 23 à frente) | `[CÓDIGO]` |

## 1.2 Baseline da Fase A

A Fase A parte de **`0b4ed5a` — `(docs) record the NCBI retrieval performance plan and implementation`** (2026‑08‑15 21:01 ‑0300). `[INFERÊNCIA]` a partir de:

- Há um intervalo temporal limpo de **15 dias** entre `0b4ed5a` (2026‑08‑15) e o primeiro commit da Fase A `175a12a` (2026‑08‑30 19:47). `[CÓDIGO]`
- Todos os 13 commits de 2026‑08‑30 tocam **apenas** arquivos novos de `src/pago_pipeline/*` (5 estágios), `notebooks/10_dataset_audit.ipynb`, `tests/*`, `README.md`, `.gitignore` — nenhum toca módulos pré‑Fase‑A. `[CÓDIGO]`
- 9 dos 13 commits trazem o sufixo explícito `(Phase A)` na mensagem. `[CÓDIGO]`

`master..HEAD` contém 23 commits; os 10 primeiros (`22a967e` … `0b4ed5a`) são o trabalho de performance de retrieval NCBI **anterior** à Fase A. `[CÓDIGO]`

## 1.3 Commits da Fase A em ordem cronológica

| # | Hash | Data (‑03) | Assunto | Objetivo real (confirmado no código) |
|---|---|---|---|---|
| 1 | `175a12a` | 08‑30 19:47 | `(feat) add technical-only prefilter stage (Phase A)` | Cria `pago_technical_prefilter.py` + `_snapshot.py` + 2 testes. Partição técnica‑apenas de um metadata CSV. `[CÓDIGO]` |
| 2 | `8a200ff` | 08‑30 19:52 | `(feat) add provenance-honest derived protein FASTA stage (Phase A)` | Cria `derived_protein_fasta.py` + `_snapshot.py` + 2 testes. FASTA derivado de uma seleção de UIDs, com `artifact_type` próprio e revalidação por hash. `[CÓDIGO]` |
| 3 | `c66e467` | 08‑30 19:55 | `(feat) let SWeeP accept a derived protein FASTA snapshot (Phase A)` | Modifica `ncbi_fasta_snapshot.py` (parâmetro `allowed_artifact_types`) e `sweep_genes_snapshot.py` (aceita `derived_protein_fasta_snapshot`, compara `source_record_ids_sha256`). Novo teste `test_fasta_snapshot_artifact_types.py`. `[CÓDIGO]` |
| 4 | `599a111` | 08‑30 20:00 | `(feat) add query reference recall stage (Phase A)` | Cria `query_reference_recall.py` + `_snapshot.py` + fixture CSV (11 linhas) + 2 testes. Também adiciona as negações `!tests/fixtures/**` ao `.gitignore`. Versão inicial: matching só por accession. `[CÓDIGO]` |
| 5 | `f7e91ec` | 08‑30 20:03 | `(feat) add NCBI ESearch preflight stage (Phase A)` | Cria `ncbi_esearch_preflight.py` + `_snapshot.py` + 2 testes. ESearch com History + amostra + guarda `max_uid_count`. `[CÓDIGO]` |
| 6 | `b54c698` | 08‑30 20:06 | `(feat) add notebook 10 dataset audit + README … (Phase A)` | Cria `notebooks/10_dataset_audit.ipynb` (orquestração) e a seção "Annotation-enriched candidate set" no `README.md`. `[CÓDIGO]` |
| 7 | `8629d20` | 08‑30 20:06 | `(fix) use WebEnv (not webenv) for the ESearch preflight sample fetch` | Corrige 1 linha em `ncbi_esearch_preflight.py`: `Entrez.efetch(..., WebEnv=…)` (maiúsculo). Sem o fix, a busca de amostra falhava. `[CÓDIGO]` |
| 8 | `d21f25a` | 08‑30 20:23 | `(feat) expand query recall reference set; NaN for empty strata (Phase A)` | Introduz `RECALL_STRATA`, status `EVALUABLE`/`NOT_EVALUABLE`, recall `None` (não `0.0`) para estrato vazio. Fixture cresce. Ajustes em notebook 10 e testes. `[CÓDIGO]` |
| 9 | `0fd283a` | 08‑30 20:32 | `(feat) add PsPIWI-RE reference; stratify PIWI_RE by ago_family (Phase A)` | Estrato `PIWI_RE` passa a ser selecionado por `ago_family == "PIWI_RE"` (não por `clade`). Adiciona `WP_014597637.1`. `[CÓDIGO]` |
| 10 | `21793bc` | 08‑30 20:50 | `(feat) expand PIWI-RE reference set to 7 (1 experimental + 6 computational)` | +6 linhas PIWI‑RE `CURATED_COMPUTATIONAL`; cria `query_recall_reference_set_curation_notes.md`; ajusta README e teste. `[CÓDIGO]` |
| 11 | `98e66c6` | 08‑30 20:53 | `(test) drop PIWI_RE from allowed clade values; assert PIWI-RE rows are clade=UNRESOLVED` | Só `tests/test_query_reference_recall.py`: remove `"PIWI_RE"` do conjunto de `clade` permitido; assert de que toda linha `ago_family==PIWI_RE` tem `clade==UNRESOLVED`. `[CÓDIGO]` |
| 12 | `34c0ca6` | 08‑30 22:04 | `(fix) recall detail as string; drop 'retain' from CELL 10 exclusions; document RsAgo false miss` | Loader do recall lê `detail` com dtypes string; CELL 10 do notebook ignora `retain` ao listar exclusões; seção "Known limitation" nas notas de curadoria. `[CÓDIGO]` |
| 13 | `0013c6d` | 08‑30 22:26 | `(feat) recognize reference pAgos recovered under a different accession by identical sequence` | Matcher de identidade de sequência (`SEQUENCE_SHA256`), `sequence_sha256`/`sequence_length` nas 21 referências, `matching_strategy_sha256`, `snapshot_format_version 1.1`, dois readings de recall, 17 testes de lógica + 4 de snapshot, notebook 10 CELL 9/12, README, notas. `[CÓDIGO]` |

## 1.4 Arquivos adicionados / modificados na Fase A (baseline → HEAD)

`git diff --name-status 0b4ed5a..HEAD` `[CÓDIGO]`:

**Adicionados (A) — 20 arquivos** (1 notebook + 8 módulos `src/pago_pipeline` + 2 fixtures + 9 arquivos de teste):

```
notebooks/10_dataset_audit.ipynb

src/pago_pipeline/ncbi_esearch_preflight.py
src/pago_pipeline/ncbi_esearch_preflight_snapshot.py
src/pago_pipeline/query_reference_recall.py
src/pago_pipeline/query_reference_recall_snapshot.py
src/pago_pipeline/pago_technical_prefilter.py
src/pago_pipeline/pago_technical_prefilter_snapshot.py
src/pago_pipeline/derived_protein_fasta.py
src/pago_pipeline/derived_protein_fasta_snapshot.py

tests/fixtures/query_recall_reference_set.csv
tests/fixtures/query_recall_reference_set_curation_notes.md

tests/test_ncbi_esearch_preflight.py
tests/test_ncbi_esearch_preflight_snapshot.py
tests/test_query_reference_recall.py
tests/test_query_reference_recall_snapshot.py
tests/test_pago_technical_prefilter.py
tests/test_pago_technical_prefilter_snapshot.py
tests/test_derived_protein_fasta.py
tests/test_derived_protein_fasta_snapshot.py
tests/test_fasta_snapshot_artifact_types.py
```

**Modificados (M) — 4 arquivos:**

```
.gitignore                                  (+ negações de fixtures/resources; + 'tmp'; texto SWeeP)
README.md                                   (+ seção "Annotation-enriched candidate set")
src/pago_pipeline/ncbi_fasta_snapshot.py    (+ allowed_artifact_types, sem mudança de comportamento default)
src/pago_pipeline/sweep_genes_snapshot.py   (+ source_fasta_artifact_types; grava/compara source_fasta_record_ids_sha256)
```

**Removidos (D):** nenhum. `[CÓDIGO]`

**Totais do diff (`git diff --stat 0b4ed5a..HEAD`):** 24 arquivos (20 A + 4 M), +6.256 / −7 linhas. As 7 remoções estão em `ncbi_fasta_snapshot.py` (−1, refatoração de assinatura) e `sweep_genes_snapshot.py` (−6, refatoração de `_source_fasta_snapshot_identity_matches`) — nenhuma lógica pré‑existente foi deletada, só reindentada/estendida. `[CÓDIGO]`

**Não versionado, porém presente no working tree** (produzido pela execução real): `data/01-raw/{esearch_preflight,protein_uid_snapshots,protein_xml_snapshots}__annotation_enriched_candidate_set/` (untracked), `data/02-intermediate/{protein_metadata_csv,query_reference_recall}__annotation_enriched_candidate_set/`, `data/02-intermediate/derived_protein_fasta__annotation_enriched_proteome/`, `data/03-features/pago_technical_prefilter/` (ignorados por `.gitignore`). `[EXEC-LOCAL]`

---

# 2. Arquitetura final da Fase A

## 2.1 Pipeline real (ordem de execução no notebook 10)

```
CELL 3/4  configuração + paths
   │
CELL 5    ncbi_esearch_preflight      ──► preflight_report.json (+ sample_protein_uids.txt)
   │                                       guarda: Count, translated_query, WebEnv/QueryKey, QC de amostra
   │                                       bloqueia se Count > max_uid_count
   ▼
CELL 6    resolve_ncbi_protein_uid_snapshot            (módulo PRÉ-Fase-A, reusado)
   │        ESearch(history) → EFetch rettype=uilist   ──► protein_uids.txt  (sha256)
   ▼
CELL 7    resolve_ncbi_protein_xml_snapshot            (módulo PRÉ-Fase-A, reusado)
   │        EFetch rettype=gp&retmode=xml em lotes de 100, concorrência 4, resume  ──► protein_records.xml
   ▼
CELL 8    resolve_ncbi_protein_metadata_snapshot       (módulo PRÉ-Fase-A, reusado)
   │        XML → CSV achatado (148 colunas)  ──► protein_metadata.csv (+ qc_report.json)
   ├─────────────► CELL 9  query_reference_recall   (Fase A)  ──► summary.csv + detail.csv   [somente relatório; não filtra]
   ▼
CELL 10   pago_technical_prefilter     (Fase A)  ──► retained_protein_uids.txt + retained_records.csv
   │                                                 + excluded_technical.csv + prefilter_counts.csv
   ▼
CELL 11   derived_protein_fasta        (Fase A)  ──► protein_sequences.fasta (+ selection_report.json)
   │        dataset_kind = annotation_enriched_proteome
   ▼
CELL 12/13 audit_summary + exposição de variáveis para os notebooks a jusante
```

A ordem pedida no enunciado (`query → preflight → UID → XML → metadata → recall → prefilter → derived FASTA`) está **correta**, com uma precisão: **`query_reference_recall` e `pago_technical_prefilter` são ramos paralelos que ambos consomem o metadata snapshot**; o recall **não** está na cadeia que produz o FASTA. `[CÓDIGO]` (notebook 10 CELL 9 lê `METADATA_SNAPSHOT_ROOT_DIRECTORY`; CELL 11 encadeia de `TECHNICAL_PREFILTER_SNAPSHOT_ROOT_DIRECTORY`, não do recall).

## 2.2 Estágio por estágio

### E0 — Configuração (notebook 10, CELL 3)

- `DATASET_NAME = 'annotation_enriched_candidate_set'` `[VERSIONADO]`
- `SEARCH_QUERY = '(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])'` `[VERSIONADO]`
- `MAX_UID_COUNT = 250_000` (o default do módulo é `150_000`; o notebook eleva para 250k) `[VERSIONADO]` / `[CÓDIGO]`

### E1 — ESearch preflight

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Descrever o custo de um retrieval completo **antes** de baixá‑lo: quantos UIDs, como o NCBI traduziu a query, e uma QC de amostra. Bloqueia retrievals acidentalmente enormes. | `[CÓDIGO]` (docstring `run_ncbi_esearch_preflight`) |
| Entrada | `search_query`, `ncbi_email`, `ncbi_api_key?`, `max_uid_count`, `sample_size` | `[CÓDIGO]` |
| Processamento | 1× `Entrez.esearch(db="protein", term=q, retmax=0, usehistory="y")` → `Count`, `QueryTranslation`, `WebEnv`, `QueryKey`. Se `Count>0` e há History: `Entrez.efetch(rettype="uilist", retmax=sample_size, WebEnv=…, query_key=…)` → lista de UIDs de amostra; depois `Entrez.efetch(id=…, rettype="gb", retmode="xml")` → QC (contagem de registros, presença de sequência, UID extraível). Retries com backoff exponencial `5·2^k` s, 5 tentativas. A falha da amostra é best‑effort (registrada em `sample_fetch_error`, não fatal). | `[CÓDIGO]` |
| Saída | `preflight_report.json` (dump do dataclass `EsearchPreflightResult`), `sample_protein_uids.txt` | `[CÓDIGO]` |
| `artifact_type` | `ncbi_esearch_preflight` | `[CÓDIGO]` |
| `snapshot_format_version` | `1.0` | `[CÓDIGO]` |
| Snapshot root (notebook) | `data/01-raw/esearch_preflight__annotation_enriched_candidate_set` | `[VERSIONADO]` |
| Hashes/proveniência no manifesto | `output_files.{preflight_report_file,sample_uid_file}.sha256`; `search_query`; `translated_query`; `result_count`; `max_uid_count`; `exceeds_max_uid_count`; `history_web_env`; `history_query_key`; `retrieved_at_utc`; contadores de amostra; `python_version`; `biopython_version` | `[CÓDIGO]` |
| Módulos | `src/pago_pipeline/ncbi_esearch_preflight.py` (lógica), `…_snapshot.py` (persistência) | `[CÓDIGO]` |
| Célula | notebook 10, CELL 5 | `[VERSIONADO]` |
| Testes | `tests/test_ncbi_esearch_preflight.py` (4), `tests/test_ncbi_esearch_preflight_snapshot.py` (2) | `[VERSIONADO]` |
| Particularidade | `save_*` **sempre materializa** o relatório (mesmo se `exceeds_max_uid_count`); é `resolve_*` que **levanta `RuntimeError`** se exceder e `allow_exceeds_max_uid_count=False`. O snapshot `latest/` é considerado válido só se `manifest["search_query"] == search_query` pedido. | `[CÓDIGO]` |

### E2 — UID acquisition (módulo PRÉ‑Fase‑A, reusado sem alteração)

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Materializar a lista canônica e ordenada de `protein_uid` que define o dataset, com hash reproduzível. | `[INFERÊNCIA]` |
| Entrada | `search_query`, paginação (`page_size`), controles de falha, `ncbi_email/api_key` | `[VERSIONADO]` (CELL 6) |
| Processamento | `esearch_history_efetch_uilist`: 1 ESearch (History) + N EFetch `rettype=uilist` paginados; dedup + sort | `[EXEC-LOCAL]` (`uid_retrieval_strategy` no manifesto) |
| Saída | `protein_uids.txt` | `[EXEC-LOCAL]` |
| `snapshot_format_version` | `1.1` | `[EXEC-LOCAL]` |
| Snapshot root | `data/01-raw/protein_uid_snapshots__annotation_enriched_candidate_set` | `[VERSIONADO]` (CELL 4) |
| Hashes/proveniência | `protein_uids_sha256`; `ncbi_reported_result_count`; `normalized_protein_uid_count`; `raw_protein_uid_count`; `history_web_env/query_key`; telemetria de latência e `failure_counts` por estágio | `[EXEC-LOCAL]` |
| Célula | CELL 6 | `[VERSIONADO]` |
| Testes | os do módulo NCBI pré‑existente (fora do escopo Fase A) | `[INFERÊNCIA]` |

### E3 — XML acquisition (módulo PRÉ‑Fase‑A, reusado)

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Baixar o registro GenPept/XML completo de cada proteína, em lotes, com workspace resumível e telemetria. | `[INFERÊNCIA]` |
| Entrada | `source_uid_snapshot_root_directory`, `xml_batch_size`, `max_concurrent_requests`, `enable_batch_resume`, `purge_batch_workspace_on_success` | `[VERSIONADO]` (CELL 7) |
| Processamento | EFetch `rettype=gp`, `retmode=xml`, lotes de 100 UIDs, concorrência 4, rate‑limit no início da requisição, circuit breaker, parser XML em streaming | `[EXEC-LOCAL]` (`request_policy` no manifesto) |
| Saída | `protein_records.xml` consolidado + `protein_uids.txt` | `[EXEC-LOCAL]` |
| `artifact_type` / versão | `ncbi_protein_xml_snapshot` / `1.1` | `[EXEC-LOCAL]` |
| Snapshot root | `data/01-raw/protein_xml_snapshots__annotation_enriched_candidate_set` | `[VERSIONADO]` |
| Hashes/proveniência | `xml_file_sha256`; `batches[].xml_payload_sha256` (525 lotes); `consolidated_record_count`; `source_uid_sha256`; `source_uid_snapshot_manifest_sha256`; telemetria com `retry_count`, `failure_counts` (http_5xx, truncated_response, …) | `[EXEC-LOCAL]` |
| Célula | CELL 7 | `[VERSIONADO]` |

### E4 — Metadata (módulo PRÉ‑Fase‑A, reusado)

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Achatar o XML aninhado num CSV tabular de 1 linha por `protein_uid`, com QC. | `[INFERÊNCIA]` |
| Entrada | `source_xml_snapshot_root_directory` | `[VERSIONADO]` (CELL 8) |
| Processamento | XML → DataFrame (148 colunas: `gbseq__*`, `taxonomy__*`, `reference__*`, `feature__*`); QC de 5 checagens | `[EXEC-LOCAL]` |
| Saída | `protein_metadata.csv` + `qc_report.json` | `[EXEC-LOCAL]` |
| `artifact_type` / versão | `ncbi_protein_metadata_snapshot` / `1.0` | `[EXEC-LOCAL]` |
| Snapshot root | `data/02-intermediate/protein_metadata_csv__annotation_enriched_candidate_set` | `[VERSIONADO]` |
| Hashes/proveniência | `csv_file_sha256`; `column_count`; `columns[]`; `observed_feature_keys/qualifiers`; herda `search_query`/`translated_query`; `source_xml_snapshot_relative_path` | `[EXEC-LOCAL]` |
| Célula | CELL 8 | `[VERSIONADO]` |
| QC checagens | `protein_uid_has_no_duplicates`, `protein_uid_has_no_empty_values`, `row_count_matches_metadata_manifest`, `row_count_matches_source_xml`, `schema_matches_metadata_manifest` — todas `true` na execução real | `[EXEC-LOCAL]` |

### E5 — Query reference recall (Fase A) — ver Seção 8 para detalhamento

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Medir quantas pAgos/PIWI‑RE **conhecidas** o texto da query recupera, estratificado, para justificar o nome "annotation‑enriched" e não "universo". **Somente relatório**: não filtra nem altera o dataset. | `[CÓDIGO]` (docstring) |
| Entrada | `source_metadata_snapshot_root_directory`, `query_recall_reference_set_csv_path` (fixture versionada) | `[CÓDIGO]` |
| Processamento | Lê `protein_metadata.csv` (só `protein_uid`, `gbseq__accession_version`, `gbseq__sequence`); casa cada uma das 21 referências por hierarquia `EXACT_ACCESSION_VERSION → SAME_BASE_ACCESSION → SEQUENCE_SHA256 → NONE`; calcula 2 recalls por estrato | `[CÓDIGO]` |
| Saída | `reference_recall_summary.csv`, `reference_recall_detail.csv` | `[CÓDIGO]` |
| `artifact_type` / versão | `query_reference_recall` / `1.1` | `[CÓDIGO]` |
| Snapshot root | `data/02-intermediate/query_reference_recall__annotation_enriched_candidate_set` | `[VERSIONADO]` |
| Hashes/proveniência | `matching_strategy` + `matching_strategy_sha256`; `query_recall_reference_set_csv_sha256`; `source_metadata_csv_sha256`; `source_metadata_manifest_sha256`; `stratum_exact_recall`, `stratum_equivalent_recall`, `stratum_recall_status`; `output_files.*.sha256` | `[CÓDIGO]` |
| Célula | CELL 9 | `[VERSIONADO]` |
| Testes | `tests/test_query_reference_recall.py` (17), `tests/test_query_reference_recall_snapshot.py` (4) | `[VERSIONADO]` |

### E6 — Technical prefilter (Fase A) — ver Seção 9

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Remover **apenas** registros tecnicamente inutilizáveis; nunca por biologia, texto de anotação ou comprimento. | `[CÓDIGO]` |
| Entrada | metadata snapshot (`protein_uid`, `gbseq__sequence`, `gbseq__length` opcional) | `[CÓDIGO]` |
| Processamento | Primeira‑regra‑vence: `drop_unprocessable_record` (uid vazio) → `drop_technical_duplicate` (uid repetido) → `drop_missing_sequence` → `drop_invalid_sequence_characters` → `retain` (+ `length_warning` fora de [200, 2000]) | `[CÓDIGO]` |
| Saída | `retained_protein_uids.txt`, `retained_records.csv`, `excluded_technical.csv`, `prefilter_counts.csv` | `[CÓDIGO]` |
| `artifact_type` / versão | `pago_technical_prefilter` / `1.0` | `[CÓDIGO]` |
| Snapshot root (notebook) | `data/03-features/pago_technical_prefilter` (compartilhado, não sufixado por dataset) | `[VERSIONADO]` |
| Hashes/proveniência | `technical_prefilter_policy` + `technical_prefilter_policy_sha256`; `counts_by_decision`; `source_metadata_csv_sha256`; `source_metadata_manifest_sha256`; `source_metadata_row_count`; `output_files.*.sha256` | `[CÓDIGO]` |
| Célula | CELL 10 | `[VERSIONADO]` |
| Testes | `tests/test_pago_technical_prefilter.py` (7), `tests/test_pago_technical_prefilter_snapshot.py` (4) | `[VERSIONADO]` |

### E7 — Derived FASTA (Fase A) — ver Seção 10

| Aspecto | Detalhe | Selo |
|---|---|---|
| Por que existe | Materializar o proteoma retido como FASTA rastreável, com `artifact_type` próprio (não `ncbi_protein_fasta_snapshot`), para o SWeeP a jusante. | `[CÓDIGO]` |
| Entrada | metadata snapshot + selection snapshot (`retained_protein_uids.txt` do prefilter) + `record_selection_rule` + `record_selection_config_sha256` + `dataset_kind` | `[CÓDIGO]` |
| Processamento | Resolve os UIDs contra o metadata (ordem `as_selected`), exporta CSV temporário → `export_metadata_csv_to_fasta` | `[CÓDIGO]` |
| Saída | `protein_sequences.fasta`, `selection_report.json` | `[CÓDIGO]` |
| `artifact_type` / versão | `derived_protein_fasta_snapshot` / `1.0` | `[CÓDIGO]` |
| Snapshot root (notebook) | `data/02-intermediate/derived_protein_fasta__annotation_enriched_proteome` | `[VERSIONADO]` |
| Hashes/proveniência | `derived_from_artifact_type`, `derived_from_manifest_sha256`, `record_selection_rule`, `record_selection_config_sha256`, `record_order`, `source_record_ids_sha256`, `fasta_file_sha256`, `fasta_record_count`, `selection_report_file_sha256`, `source_metadata_*`, `source_selection_*` | `[CÓDIGO]` |
| Célula | CELL 11 | `[VERSIONADO]` |
| Testes | `tests/test_derived_protein_fasta.py` (7), `tests/test_derived_protein_fasta_snapshot.py` (4) | `[VERSIONADO]` |

### Adaptação transversal — SWeeP

`sweep_genes_snapshot.py` ganhou `SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES = ("ncbi_protein_fasta_snapshot", "derived_protein_fasta_snapshot")` e passa esse conjunto a `load_fasta_snapshot_by_directory`; o manifesto do SWeeP passa a gravar `source_fasta_artifact_type` e `source_fasta_record_ids_sha256`, e `_source_fasta_snapshot_identity_matches` compara também esse fingerprint quando ambos os lados o possuem. `ncbi_fasta_snapshot.py` ganhou o parâmetro `allowed_artifact_types` (default inalterado → sem mudança de comportamento para o caminho existente). `[CÓDIGO]`

---

# 3. Inventário arquivo por arquivo

## 3.1 Tabela

| Path | Tipo | Responsabilidade | Principais símbolos | Entradas | Saídas | Dependências | Testes |
|---|---|---|---|---|---|---|---|
| `src/pago_pipeline/ncbi_esearch_preflight.py` | código (lógica pura) | ESearch+History+amostra QC | `EsearchPreflightResult` (dataclass frozen), `run_ncbi_esearch_preflight`, `parse_esearch_history_response`, `parse_uilist_text`, `summarize_sample_xml`, `build_esearch_preflight_result`, `_run_with_retries` | `search_query`, `ncbi_email`, `ncbi_api_key`, `max_uid_count`, `sample_size` | `EsearchPreflightResult` | `Bio.Entrez`, `ncbi_api._configured_ncbi_entrez_urlopen`, `ncbi_xml_stream.extract_protein_uid_from_gbseq_element` | `test_ncbi_esearch_preflight.py` |
| `src/pago_pipeline/ncbi_esearch_preflight_snapshot.py` | código (I/O + manifesto) | Persistir/validar/reusar o preflight | `save_*`, `resolve_*`, `load_*_by_directory`, `load_latest_*`, `latest_*_is_available`, `list_saved_*`, `get_most_recent_*`, `_raise_if_exceeds_max_uid_count`, `_build_manifest_payload` | idem lógica + `snapshot_mode`, `snapshot_root_directory`, `allow_exceeds_max_uid_count` | dir imutável + `latest/` (`manifest.json`, `preflight_report.json`, `sample_protein_uids.txt`) | `ncbi_snapshot` (SnapshotMode, helpers), `storage` | `test_ncbi_esearch_preflight_snapshot.py` |
| `src/pago_pipeline/query_reference_recall.py` | código (lógica pura) | Recall estratificado do painel de referência | `QueryReferenceRecallResult`, `ReferenceMatchMethod` (Enum), `compute_query_reference_recall`, `normalize_protein_sequence`, `protein_sequence_sha256`, `build_matching_strategy_payload/_sha256`, `_recall`, `RECALL_STRATA`, `SEQUENCE_NORMALIZATION` | `reference_dataframe`, `retrieved_metadata_dataframe`, nomes de coluna | `QueryReferenceRecallResult` (summary df, detail df, dicts de recall) | `pandas`, `hashlib`, `json`, `re` | `test_query_reference_recall.py` |
| `src/pago_pipeline/query_reference_recall_snapshot.py` | código (I/O + manifesto) | Persistir/validar/reusar o recall | 7‑API padrão + `_build_manifest_payload`, `_validate_loaded_*` (rejeita `≠1.1` **e** `matching_strategy_sha256` divergente) | `snapshot_mode`, `snapshot_root_directory`, `source_metadata_snapshot_root_directory`, `query_recall_reference_set_csv_path` | dir imutável + `latest/` (`manifest.json`, `reference_recall_summary.csv`, `reference_recall_detail.csv`) | `ncbi_metadata_snapshot`, `ncbi_snapshot`, `query_reference_recall`, `storage` | `test_query_reference_recall_snapshot.py` |
| `src/pago_pipeline/pago_technical_prefilter.py` | código (lógica pura) | Partição técnica‑apenas | `PagoTechnicalPrefilterDecision` (Enum), `PagoTechnicalPrefilterResult`, `build_pago_technical_prefilter_partition`, `build_technical_prefilter_counts_dataframe`, `build_technical_prefilter_policy_payload/_sha256`, `DEFAULT_ALLOWED_RESIDUES`, `DEFAULT_TOLERATED_*`, `DEFAULT_LENGTH_WARNING_{MIN,MAX}` | `metadata_dataframe`, política (resíduos, colunas, banda de comprimento) | `PagoTechnicalPrefilterResult` (retained df, excluded df, contagens, uids retidos) | `pandas`, `hashlib`, `json` | `test_pago_technical_prefilter.py` |
| `src/pago_pipeline/pago_technical_prefilter_snapshot.py` | código (I/O + manifesto) | Persistir/validar/reusar o prefilter | 7‑API padrão + `_source_metadata_snapshot_identity_matches` | `snapshot_mode`, `snapshot_root_directory`, `source_metadata_snapshot_root_directory` | dir imutável + `latest/` (4 CSV/TXT + `manifest.json`) | `ncbi_metadata_snapshot`, `ncbi_snapshot`, `pago_technical_prefilter`, `storage` | `test_pago_technical_prefilter_snapshot.py` |
| `src/pago_pipeline/derived_protein_fasta.py` | código (lógica pura) | Resolver seleção de UIDs → FASTA rastreável | `DerivedFastaRecordOrder` (Enum), `DerivedFastaSelectionResult`, `build_derived_fasta_selection`, `parse_protein_uids_from_fasta_deflines`, `compute_record_ids_sha256_for_order`, `_compute_record_ids_sha256` | `metadata_dataframe`, `selected_protein_uids`, `record_order`, `drop_missing_uids` | `DerivedFastaSelectionResult` | `pandas`, `hashlib` | `test_derived_protein_fasta.py` |
| `src/pago_pipeline/derived_protein_fasta_snapshot.py` | código (I/O + manifesto) | Persistir/validar/reusar o FASTA derivado | 7‑API padrão + `_load_selection_snapshot`, `_source_identity_matches`, `_validate_loaded_*` (revalida deflines, contagem, `source_record_ids_sha256` por ordem) | `snapshot_mode`, roots de metadata e de seleção, `selection_artifact_type`, `selection_uid_list_file_name`, `record_selection_rule`, `record_selection_config_sha256`, `dataset_kind`, `sequence_line_width`, `record_order` | dir imutável + `latest/` (`protein_sequences.fasta`, `selection_report.json`, `manifest.json`) | `derived_protein_fasta`, `metadata_to_fasta.export_metadata_csv_to_fasta`, `ncbi_metadata_snapshot`, `ncbi_snapshot`, `storage` | `test_derived_protein_fasta_snapshot.py` |
| `src/pago_pipeline/ncbi_fasta_snapshot.py` | código (modificado) | +`allowed_artifact_types` (default = `("ncbi_protein_fasta_snapshot",)`) em `_validate_loaded_fasta_snapshot_payload`, `load_fasta_snapshot_by_directory`, `load_latest_fasta_snapshot`, `latest_fasta_snapshot_is_available` | `DEFAULT_FASTA_SNAPSHOT_ARTIFACT_TYPES` | — | — | — | `test_fasta_snapshot_artifact_types.py` |
| `src/pago_pipeline/sweep_genes_snapshot.py` | código (modificado) | Aceitar `derived_protein_fasta_snapshot` como fonte; gravar/comparar `source_fasta_record_ids_sha256` | `SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES`; param `source_fasta_artifact_types` em `save_*`, `latest_*_is_available`, `resolve_*`; `_source_fasta_snapshot_identity_matches` estendido | — | — | `ncbi_fasta_snapshot` | `test_fasta_snapshot_artifact_types.py` (parcial) |
| `notebooks/10_dataset_audit.ipynb` | notebook (orquestração) | Executar os 8 estágios em ordem, imprimir auditoria | 14 células (0 markdown + 13 código) | `.env` (`NCBI_EMAIL`, `NCBI_API_KEY`) | snapshots em `data/…`; `audit_summary` | todos os `*_snapshot` acima + `pandas`, `dotenv` | — (validação via `python -m unittest`, não via nbval) |
| `tests/fixtures/query_recall_reference_set.csv` | fixture (dados curados, versionada) | Painel de 21 referências pAgo/PIWI‑RE | 12 colunas; 21 linhas | — | — | — | consumida por ambos os testes de recall |
| `tests/fixtures/query_recall_reference_set_curation_notes.md` | documentação (versionada) | Justificar cada linha, método dos hashes, caso RsAgo, IDs históricos | — | — | — | — | — |
| `tests/test_*` (9 arquivos) | testes `unittest` | ver Seção 12 | — | — | — | `unittest`, `pandas`, `unittest.mock` | — |
| `README.md` | documentação (modificada) | +seção "Annotation-enriched candidate set", +item 11 na ordem dos notebooks | — | — | — | — | — |
| `.gitignore` | config (modificada) | +`!tests/fixtures/**`, +`!src/pago_pipeline/resources/**`, +`tmp` | — | — | — | — | — |

## 3.2 Explicação dos arquivos mais importantes

### `query_reference_recall.py` (370 linhas)

O núcleo científico da Fase A. Função pública `compute_query_reference_recall(*, reference_dataframe, retrieved_metadata_dataframe, retrieved_accession_column="gbseq__accession_version", retrieved_protein_uid_column="protein_uid", retrieved_sequence_column="gbseq__sequence", progress_callback=None)`. `[CÓDIGO]`

Fluxo interno `[CÓDIGO]`:
1. `_validate_reference_columns` exige `accession, ago_family, clade, reference_label_source, reference_label_evidence` (senão `RuntimeError`).
2. Constrói do lado recuperado: `retrieved_versioned_set` (conjunto de accession.version), `retrieved_by_bare` (accession sem versão → menor `(accession, uid)` lexicográfico), `versioned_to_uid`, e — se a coluna de sequência existe — `sequence_sha256_to_hits` (hash → lista **ordenada** de `(accession, uid)`).
3. Para cada referência: tenta na ordem `EXACT_ACCESSION_VERSION` → `SAME_BASE_ACCESSION` → `SEQUENCE_SHA256` → `NONE`. Registra `match_method`, `matched_accession`, `matched_protein_uid`, `sequence_match_count`.
4. `recovered = match_method != NONE`; `recovered_exact_accession = match_method ∈ {EXACT_ACCESSION_VERSION, SAME_BASE_ACCESSION}`.
5. Para cada estrato de `RECALL_STRATA` (`overall`, `LONG_A`, `LONG_B`, `SHORT`, `PIWI_RE`): `exact_recovered_count = Σ recovered_exact_accession`, `equivalent_recovered_count = Σ recovered`; `_recall(n, k)` devolve `(None, "NOT_EVALUABLE")` se `n==0`, senão `(k/n, "EVALUABLE")`.

`RECALL_STRATA` `[CÓDIGO]`: `(("overall", "overall_reference_recall", None, None), ("LONG_A", "long_a_reference_recall", "clade", "LONG_A"), ("LONG_B", …, "clade", "LONG_B"), ("SHORT", …, "clade", "SHORT"), ("PIWI_RE", "piwi_re_reference_recall", "ago_family", "PIWI_RE"))`. Note que o estrato PIWI‑RE seleciona por `ago_family`, não `clade`.

`build_matching_strategy_payload()` retorna `{strategy_kind:"query_reference_recall_matching", strategy_version:"2.0", match_hierarchy:[EXACT_ACCESSION_VERSION, SAME_BASE_ACCESSION, SEQUENCE_SHA256], sequence_normalization:"strip_all_whitespace_then_uppercase", sequence_hash:"sha256", multi_sequence_hit_representative:"min_accession_lexicographic", recall_readings:[exact_accession_recall, retrieval_equivalent_recall]}`; `build_matching_strategy_sha256()` = `sha256(json.dumps(payload, sort_keys=True, separators=(",",":")))` = `3460b048fc6de363ddf9282c2943a44c284e51d2c48092ca14426600b2871a08`. `[CÓDIGO]` (valor pinado em `tests/test_query_reference_recall.py::_EXPECTED_MATCHING_STRATEGY_SHA256`).

### `pago_technical_prefilter.py` (329 linhas)

`build_pago_technical_prefilter_partition` aplica as regras em **ordem reversa de prioridade** sobre uma `pd.Series` inicializada com `"retain"`, de modo que regras de maior prioridade sobrescrevem: primeiro `invalid_sequence`, depois `missing_sequence`, depois `duplicate`, por fim `missing_protein_uid` (a de maior prioridade). `[CÓDIGO]`

`_normalize_sequence_text(value)` = `"".join(str(value).split()).upper()` — remove todo whitespace e faz upper. O conjunto de caracteres aceitos é a união de `DEFAULT_ALLOWED_RESIDUES="ACDEFGHIKLMNPQRSTVWY"` + `DEFAULT_TOLERATED_AMBIGUOUS_RESIDUES="BZJXUO"` + `DEFAULT_TOLERATED_SEQUENCE_SYMBOLS="*-"`. `[CÓDIGO]`

`length_warning` = `numeric_length.isna() | (< 200) | (> 2000)`. **Nunca exclui.** `[CÓDIGO]`

Invariante interna: `len(retained) + len(excluded) == input_record_count`, senão `RuntimeError`. `[CÓDIGO]`

`build_technical_prefilter_policy_sha256()` = `032fbc727b68ceb97cdda00ca5764db3ddeebc0868058d0c73883c23771fa523` `[EXEC-LOCAL]` (registrado no manifesto real; o código produz o mesmo por serialização determinística).

### `derived_protein_fasta.py` (210 linhas)

`build_derived_fasta_selection`:
- normaliza `protein_uid` (strip); **exige unicidade** no metadata (`RuntimeError "must be unique"` se duplicado);
- remove UIDs vazios da lista pedida; detecta duplicatas na lista pedida (`RuntimeError "duplicates"` salvo `drop_missing_uids`);
- UIDs ausentes no metadata → `RuntimeError "absent from the metadata"` salvo `drop_missing_uids`;
- `record_order="sorted_by_uid"` reordena; `"as_selected"` preserva a ordem pedida;
- sequência vazia em qualquer registro selecionado → `RuntimeError "empty amino-acid sequence"`;
- `source_record_ids_sha256 = sha256(("\n".join(uids) + "\n").encode())` (string vazia → sha256 de `""`).

`parse_protein_uids_from_fasta_deflines` extrai o token `protein_uid=` de cada defline `>protein_uid=…|accession=…|length=…|organism=… <definição>`. `[CÓDIGO]` (formato da defline confirmado em `metadata_to_fasta` e no FASTA real).

### `notebooks/10_dataset_audit.ipynb`

Orquestração pura: cada célula chama um `resolve_*` e imprime. Todas as células usam `SnapshotMode.reuse_latest_or_create`. Reexecutar o notebook com snapshots válidos **não faz rede**. `[VERSIONADO]` (ver Seção 11).

---

# 4. Inventário função por função (módulos centrais)

## 4.1 `ncbi_esearch_preflight.py`

### `run_ncbi_esearch_preflight`
- **Assinatura:** `run_ncbi_esearch_preflight(*, search_query: str, ncbi_email: str, ncbi_api_key: str|None=None, max_uid_count: int=150_000, sample_size: int=200, max_retry_attempts: int=5, retry_backoff_seconds: float=5.0, ssl_ca_file: str|None=None, ssl_ca_directory: str|None=None) -> EsearchPreflightResult` `[CÓDIGO]`
- **Responsabilidade:** 1 ESearch(History) + (opcional) 1 fetch de UIDs de amostra + 1 fetch XML de amostra + QC.
- **Retorno:** `EsearchPreflightResult` (frozen dataclass, 20 campos).
- **Efeitos colaterais:** define `Entrez.email` e `Entrez.api_key` (estado global do Biopython); abre conexões HTTP ao NCBI; dorme durante backoff.
- **Validações:** `search_query` não vazio; `ncbi_email` obrigatório; `max_uid_count>0`; `sample_size>=0`.
- **Erros:** `ValueError` (args inválidos); `RuntimeError` (`"NCBI … failed after N attempts"`) se o ESearch falhar em todas as tentativas. A falha da **amostra** não propaga — é capturada em `sample_fetch_error`.
- **Chamada em:** `ncbi_esearch_preflight_snapshot.save_ncbi_esearch_preflight_snapshot`.
- **Testes:** `test_run_preflight_uses_history_and_sample` (mock de `Entrez`), indiretamente `test_ncbi_esearch_preflight_snapshot`.

### `parse_esearch_history_response`
- **Assinatura:** `(*, response_payload_bytes: bytes) -> dict` `[CÓDIGO]`
- **Retorno:** `{result_count:int, translated_query, history_web_env, history_query_key}` via `Entrez.read`.
- **Erros:** propaga erros de parsing do Biopython.
- **Testes:** exercitada por `test_run_preflight_uses_history_and_sample`.

### `parse_uilist_text`
- **Assinatura:** `(*, payload_text: str) -> list[str]` — split por linha, strip, remove vazias. `[CÓDIGO]`
- **Testes:** `test_parse_uilist_text` (`"100\n101\n\n 102 \n" → ["100","101","102"]`).

### `summarize_sample_xml`
- **Assinatura:** `(*, xml_payload_bytes: bytes) -> dict[str,int]` `[CÓDIGO]`
- **Retorno:** `{sample_record_count, sample_records_with_sequence, sample_records_missing_sequence, sample_records_with_extractable_uid}`.
- **Efeitos:** nenhum (parse `ElementTree` em memória).
- **Erros:** a extração de UID por registro é envolta em `try/except Exception: pass` (contador best‑effort).
- **Testes:** `test_summarize_sample_xml`.

### `build_esearch_preflight_result`
- **Assinatura:** `(*, search_query, retrieved_at_utc, result_count, translated_query, history_web_env, history_query_key, max_uid_count, sample_requested_count, sample_uid_list=None, sample_summary=None, sample_fetch_error=None) -> EsearchPreflightResult` `[CÓDIGO]`
- **Lógica‑chave:** `exceeds_max_uid_count = int(result_count) > int(max_uid_count)`; captura `sys.version` e `Bio.__version__`.
- **Testes:** `test_build_result_flags_exceeds_max_uid_count`.

### `_run_with_retries` (privada)
- Backoff exponencial `retry_backoff_seconds · 2^attempt`; após esgotar, `RuntimeError`. `[CÓDIGO]`

## 4.2 `ncbi_esearch_preflight_snapshot.py`

### `resolve_ncbi_esearch_preflight_snapshot`
- **Assinatura:** `(*, snapshot_mode, snapshot_root_directory, search_query=None, ncbi_email=None, ncbi_api_key=None, max_uid_count=150_000, sample_size=200, allow_exceeds_max_uid_count=False, ssl_ca_file=None, ssl_ca_directory=None, update_latest_directory=True) -> dict` `[CÓDIGO]`
- **Responsabilidade:** dispatch de 6 passos por `SnapshotMode` — `reuse_latest` (erro se ausente), `reuse_latest_or_create` (reusa se `latest_*_is_available`), senão valida modo e chama `save_*`.
- **Efeitos:** cria diretório de snapshot; pode escrever `latest/`.
- **Erros:** `FileNotFoundError` (reuse_latest sem snapshot); `ValueError` (modo inválido, ou `ncbi_email`/`search_query` ausente quando precisa rodar); **`RuntimeError`** via `_raise_if_exceeds_max_uid_count` se `exceeds_max_uid_count` e não `allow_exceeds_max_uid_count` — **em todos os caminhos** (reuse e create).
- **Chamada em:** notebook 10 CELL 5.
- **Testes:** `test_resolve_materializes_report_and_reuses`, `test_resolve_raises_when_result_count_exceeds_max_uid_count`.

### `latest_ncbi_esearch_preflight_snapshot_is_available`
- **Assinatura:** `(*, snapshot_root_directory, search_query=None) -> bool` `[CÓDIGO]`
- **Lógica:** existe `latest/manifest.json`; valida `artifact_type`, `snapshot_format_version`, e cada `output_files[*].sha256` contra o arquivo em disco; se `search_query` passado, exige `manifest["search_query"] == search_query`.
- **Erros:** captura `(FileNotFoundError, RuntimeError, OSError, ValueError)` → `False`.

### `save_ncbi_esearch_preflight_snapshot`
- Cria `snapshots/<ts>__q_<hash12>/`; grava `preflight_report.json` (via `write_json_atomic`), `sample_protein_uids.txt` (atômico, `\n`); monta manifesto; `write_json_atomic(manifest)`; marca `immutable_snapshot_complete=True`; se `update_latest_directory`, `_replace_latest_directory`. Em qualquer `Exception` antes de completo → `shutil.rmtree` do diretório parcial + re‑raise. `[CÓDIGO]`

### `_build_manifest_payload`
- Monta o dict com `artifact_type`, `snapshot_format_version`, `snapshot_created_at_utc`, `immutable_snapshot_*`, campos do preflight, `output_files` (`file_name`, `path`, `sha256`). `[CÓDIGO]`

### `_raise_if_exceeds_max_uid_count`
- `RuntimeError` com mensagem `"… above the configured max_uid_count …"` se `manifest["exceeds_max_uid_count"]` e não `allow_exceeds_max_uid_count`. `[CÓDIGO]`

## 4.3 `query_reference_recall.py`

### `compute_query_reference_recall`
- **Assinatura:** ver 3.2.
- **Parâmetros:** `reference_dataframe` (21 refs curadas), `retrieved_metadata_dataframe` (metadata achatado), nomes de coluna, `progress_callback`.
- **Retorno:** `QueryReferenceRecallResult(summary_dataframe, detail_dataframe, stratum_exact_recall: dict[str,float|None], stratum_equivalent_recall: dict[str,float|None], stratum_recall_status: dict[str,str], reference_count: int, exact_recovered_count: int, equivalent_recovered_count: int)`.
- **Efeitos colaterais:** nenhum (puro; opcionalmente chama `progress_callback`).
- **Validações:** colunas obrigatórias na referência (`RuntimeError`); coluna de accession no recuperado (`RuntimeError`).
- **Erros:** `RuntimeError` (colunas faltando).
- **Chamada em:** `query_reference_recall_snapshot.save_query_reference_recall_snapshot`.
- **Testes que exercitam o comportamento:** `test_exact_accession_version_match`, `test_same_base_accession_different_version`, `test_different_accession_identical_sequence_is_equivalent_only`, `test_different_sequence_is_not_recovered`, `test_exact_accession_takes_priority_over_sequence_hash`, `test_same_base_accession_takes_priority_over_sequence_hash`, `test_multiple_accessions_same_hash_resolved_deterministically`, `test_stratified_recall_and_both_readings`, `test_empty_stratum_is_not_evaluable_not_zero`, `test_missing_reference_columns_raise`, `test_rsago_is_recovered_by_sequence_identity_under_alias_accession`.

### `normalize_protein_sequence(value) -> str`
- `"".join(str(value).split()).upper()`; `None`/`NaN` → `""`. `[CÓDIGO]`
- **Testes:** `test_sequence_normalization_strips_whitespace_and_uppercases`.

### `protein_sequence_sha256(value) -> str`
- `sha256(normalize_protein_sequence(value).encode("utf-8")).hexdigest()`. `[CÓDIGO]`

### `build_matching_strategy_payload() / build_matching_strategy_sha256()`
- Payload literal (ver 3.2); sha256 determinístico. **Efeito de contrato:** o sha entra no manifesto e a mudança dele invalida snapshots de recall (só de recall). `[CÓDIGO]`
- **Testes:** `test_matching_strategy_sha256_is_pinned`.

### `_recall(reference_count, recovered_count) -> (float|None, str)`
- `(None, "NOT_EVALUABLE")` se `reference_count==0`; senão `(k/n, "EVALUABLE")`. `[CÓDIGO]`
- **Testes:** `test_empty_stratum_is_not_evaluable_not_zero`.

### `ReferenceMatchMethod(str, Enum)`
- `EXACT_ACCESSION_VERSION`, `SAME_BASE_ACCESSION`, `SEQUENCE_SHA256`, `NONE`. `_STRICT_ACCESSION_METHODS = {EXACT_ACCESSION_VERSION, SAME_BASE_ACCESSION}`. `[CÓDIGO]`

## 4.4 `query_reference_recall_snapshot.py`

### `resolve_query_reference_recall_snapshot`
- **Assinatura:** `(*, snapshot_mode, snapshot_root_directory, source_metadata_snapshot_root_directory, query_recall_reference_set_csv_path, update_latest_directory=True) -> dict` `[CÓDIGO]`
- **Retorno:** `{snapshot_directory, manifest_file_path, manifest, summary_file_path, detail_file_path, summary: DataFrame, detail: DataFrame}`.
- **Erros:** `FileNotFoundError` (reuse_latest sem snapshot; ou nenhum metadata snapshot reusável); `ValueError` (modo inválido).
- **Chamada em:** notebook 10 CELL 9.
- **Testes:** `test_resolve_computes_recall_against_committed_reference_set`, `test_reuse_and_invalidation_on_reference_set_change`, `test_snapshot_written_by_earlier_methodology_is_not_reused`, `test_recovery_by_sequence_identity_under_alias_accession`.

### `save_query_reference_recall_snapshot`
- Lê `protein_metadata.csv` com `usecols` restrito a `{protein_uid, gbseq__accession_version, gbseq__sequence}`, `dtype=str`; lê a fixture (`dtype=str`, `fillna("")`); chama `compute_query_reference_recall`; grava `reference_recall_summary.csv`, `reference_recall_detail.csv`; monta manifesto com `matching_strategy` + `matching_strategy_sha256` + `stratum_exact_recall` + `stratum_equivalent_recall` (`None`→JSON null) + `stratum_recall_status` + `query_recall_reference_set_csv_sha256` + `source_metadata_*`; `_replace_latest_directory`. Rollback por `rmtree` em exceção. `[CÓDIGO]`

### `_validate_loaded_query_reference_recall_payload`
- Rejeita `artifact_type != "query_reference_recall"`; **rejeita `snapshot_format_version != "1.1"`**; **rejeita `matching_strategy_sha256 != build_matching_strategy_sha256()`**; valida `output_files[*].sha256`. `[CÓDIGO]`
- **Testes:** `test_snapshot_written_by_earlier_methodology_is_not_reused` (patch do sha para `"0"*64` → indisponível).

### `latest_query_reference_recall_snapshot_is_available`
- Valida o payload; se `query_recall_reference_set_csv_path` passado, compara `query_recall_reference_set_csv_sha256`; se `source_metadata_snapshot_root_directory` passado, compara `source_metadata_manifest_sha256` contra o metadata `latest`. `[CÓDIGO]`
- **Testes:** `test_reuse_and_invalidation_on_reference_set_change` (append de 1 linha ao CSV → indisponível).

## 4.5 `pago_technical_prefilter.py`

### `build_pago_technical_prefilter_partition`
- **Assinatura:** `(*, metadata_dataframe, allowed_residues=…, tolerated_ambiguous_residues=…, tolerated_sequence_symbols=…, length_warning_min=200, length_warning_max=2000, protein_uid_column="protein_uid", sequence_column="gbseq__sequence", length_column="gbseq__length", progress_callback=None) -> PagoTechnicalPrefilterResult` `[CÓDIGO]`
- **Retorno:** `PagoTechnicalPrefilterResult(retained_records, excluded_records, counts_by_decision: dict[str,int], retained_protein_uids: tuple[str,...], input_record_count: int)`.
- **Efeitos:** nenhum (copia o df).
- **Validações:** `protein_uid` e `gbseq__sequence` obrigatórios (`RuntimeError "required technical prefilter input columns"`); `length_warning_min>=0` e `max>=min` (`ValueError`); invariante total/disjunto (`RuntimeError`).
- **Chamada em:** `pago_technical_prefilter_snapshot.save_*`.
- **Testes:** `test_retains_records_regardless_of_annotation_text_or_length`, `test_excludes_only_technical_problems`, `test_tolerated_ambiguous_and_symbol_residues_are_kept`, `test_partition_is_total_and_disjoint`, `test_missing_required_columns_raises`, `test_counts_dataframe_sums_to_input`.

### `build_technical_prefilter_policy_payload / _sha256`
- Payload com `policy_kind`, `policy_version:"1.0"`, conjuntos de caracteres **ordenados**, banda de comprimento, nomes de coluna, `exclusion_reasons` (4), `notes`. sha256 determinístico (`sort_keys`, `separators`). `[CÓDIGO]`
- **Testes:** `test_policy_sha256_is_stable_and_parameter_sensitive`.

### `build_technical_prefilter_counts_dataframe`
- 1 linha por `PagoTechnicalPrefilterDecision` (5): `decision, count, fraction, description`. `[CÓDIGO]`

## 4.6 `pago_technical_prefilter_snapshot.py`

### `resolve_pago_technical_prefilter_snapshot`
- **Assinatura:** `(*, snapshot_mode, snapshot_root_directory, source_metadata_snapshot_root_directory, update_latest_directory=True) -> dict` `[CÓDIGO]`
- **Retorno:** `{…, "retained_records": DataFrame, "prefilter_counts": DataFrame, paths…}`.
- **Erros:** `FileNotFoundError` (sem metadata snapshot); `ValueError` (modo inválido).
- **Chamada em:** notebook 10 CELL 10, e por `derived_protein_fasta` (como selection snapshot).
- **Testes:** `test_resolve_persists_outputs_and_manifest`, `test_reuse_latest_or_create_reuses_frozen_snapshot`, `test_latest_is_unavailable_when_metadata_snapshot_changes`, `test_validation_rejects_artifact_type_and_hash_mismatch`.

### `latest_pago_technical_prefilter_snapshot_is_available`
- Valida payload; exige `technical_prefilter_policy_sha256 == build_technical_prefilter_policy_sha256()`; se root de metadata passado, `_source_metadata_snapshot_identity_matches` (compara `source_metadata_manifest_sha256`). `[CÓDIGO]`

### `save_pago_technical_prefilter_snapshot`
- Lê `protein_metadata.csv` (`low_memory=False`, `dtype={protein_uid:"string"}`); particiona; grava `retained_protein_uids.txt` (atômico), `retained_records.csv`, `excluded_technical.csv`, `prefilter_counts.csv`; manifesto com política + `counts_by_decision` + proveniência de metadata + `source_xml_snapshot_relative_path`; `_replace_latest_directory`. Rollback por `rmtree`. `[CÓDIGO]`

## 4.7 `derived_protein_fasta.py`

### `build_derived_fasta_selection`
- **Assinatura:** ver 3.2 / Seção 4 acima.
- **Retorno:** `DerivedFastaSelectionResult(selected_metadata, resolved_protein_uids, requested_uid_count, resolved_uid_count, missing_uids, duplicate_requested_uids, empty_sequence_uids, record_order, source_record_ids_sha256)`.
- **Erros:** `RuntimeError` (colunas faltando; metadata com uid duplicado; duplicata na lista pedida; uid ausente; sequência vazia); `ValueError` (`record_order` inválido).
- **Testes:** `test_selection_preserves_requested_order`, `test_sorted_by_uid_order`, `test_record_ids_sha256_matches_helper`, `test_missing_uid_raises_unless_dropped`, `test_duplicate_request_raises`, `test_duplicate_metadata_protein_uid_raises`, `test_empty_sequence_raises`.

### `parse_protein_uids_from_fasta_deflines / compute_record_ids_sha256_for_order`
- Extração de `protein_uid=` por defline; recomputa o sha por ordem (`as_selected` mantém, `sorted_by_uid` ordena). Usadas pela revalidação no snapshot. `[CÓDIGO]`

## 4.8 `derived_protein_fasta_snapshot.py`

### `resolve_derived_protein_fasta_snapshot`
- **Assinatura:** `(*, snapshot_mode, snapshot_root_directory, source_metadata_snapshot_root_directory, source_selection_snapshot_root_directory, selection_artifact_type, selection_uid_list_file_name, record_selection_rule, record_selection_config_sha256, dataset_kind, sequence_line_width=60, record_order="as_selected", drop_missing_uids=False, update_latest_directory=True) -> dict` `[CÓDIGO]`
- **Retorno:** `{snapshot_directory, manifest_file_path, manifest, fasta_file_path}`.
- **Erros:** `FileNotFoundError` (sem metadata snapshot; manifest/uid‑list de seleção ausente); `RuntimeError` (`artifact_type` de seleção divergente; qualquer `RuntimeError` de `build_derived_fasta_selection`); `ValueError` (modo inválido).
- **Chamada em:** notebook 10 CELL 11.
- **Testes:** `test_resolve_produces_honest_provenance_manifest`, `test_validation_detects_reordered_or_tampered_fasta`, `test_latest_unavailable_when_selection_changes`, `test_reuse_latest_or_create_reuses_without_resaving`.

### `_validate_loaded_derived_protein_fasta_payload`
- Rejeita `artifact_type != "derived_protein_fasta_snapshot"`; rejeita versão `!= "1.0"`; verifica `fasta_file_sha256`; **reparse das deflines** → checa `fasta_record_count` e recomputa `source_record_ids_sha256` por `record_order` (`RuntimeError "record-id set … changed"`); verifica `selection_report_file_sha256`. `[CÓDIGO]`
- **Testes:** `test_validation_detects_reordered_or_tampered_fasta` (dois caminhos: hash mismatch e, com hash corrigido, record‑id set).

### `_source_identity_matches`
- Compara `source_metadata_manifest_sha256` e `source_selection_manifest_sha256`. `[CÓDIGO]` (bug de sintaxe corrigido nesta branch — `if not …():` seguido de `return False` / `return True`.)

## 4.9 Modificações em módulos pré‑existentes

### `ncbi_fasta_snapshot.py`
- `DEFAULT_FASTA_SNAPSHOT_ARTIFACT_TYPES = ("ncbi_protein_fasta_snapshot",)`. `[CÓDIGO]`
- Parâmetro `allowed_artifact_types` adicionado a `_validate_loaded_fasta_snapshot_payload`, `load_fasta_snapshot_by_directory`, `load_latest_fasta_snapshot`, `latest_fasta_snapshot_is_available`, **todos com default = a tupla acima** → o caminho NCBI existente não muda de comportamento.
- **Testes:** `test_fasta_snapshot_artifact_types.test_default_still_rejects_non_ncbi_fasta_snapshot`.

### `sweep_genes_snapshot.py`
- `SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES = ("ncbi_protein_fasta_snapshot", "derived_protein_fasta_snapshot")`. `[CÓDIGO]`
- Parâmetro `source_fasta_artifact_types` propagado por `save_sweep_genes_snapshot`, `latest_sweep_genes_snapshot_is_available`, `resolve_sweep_genes_snapshot`.
- Manifesto do SWeeP passa a gravar `source_fasta_artifact_type` e `source_fasta_record_ids_sha256` (do manifesto FASTA fonte).
- `_source_fasta_snapshot_identity_matches` refatorada: além dos dois hashes de arquivo, **quando ambos os lados têm** `source_fasta_record_ids_sha256`, exige igualdade (snapshots NCBI antigos não têm o campo → tolerado).
- **Testes:** `test_fasta_snapshot_artifact_types.test_sweep_supported_types_accept_derived_fasta_but_not_arbitrary`.

## 4.10 Integrações com módulos NCBI/snapshot pré‑existentes

Os módulos Fase A **reusam sem alterar**: `ncbi_snapshot` (`SnapshotMode`, `_coerce_snapshot_mode`, `_replace_latest_directory`, `build_snapshot_directory_name`, `get_most_recent_snapshot_directory`, `list_saved_snapshot_directories`, `_build_query_hash`), `ncbi_metadata_snapshot` (`load_latest_metadata_snapshot`, `load_metadata_snapshot_by_directory`), `ncbi_api` (`_configured_ncbi_entrez_urlopen`), `ncbi_xml_stream` (`extract_protein_uid_from_gbseq_element`), `metadata_to_fasta` (`export_metadata_csv_to_fasta`), `storage` (`read_json_file`, `read_text_lines_from_file`, `sha256_of_file`, `sha256_of_lines`, `write_json_atomic`). `[CÓDIGO]`

- `build_snapshot_directory_name(retrieved_at_utc, search_query)` → `"<ts sanitizado>__q_<sha256(search_query)[:12]>"`. `[CÓDIGO]`
- Cada `save_*` da Fase A que não é ligado a uma query NCBI usa um **literal** como "query": prefilter → `"pago_technical_prefilter"` (→ `q_bbdf63d86832`), recall → `"query_reference_recall"` (→ `q_8994498899ae`), derived FASTA → `f"derived_protein_fasta::{dataset_kind}"` (→ `q_c99c5e155233`). `[CÓDIGO]` / `[EXEC-LOCAL]` (hashes conferem com os diretórios reais).
- `_replace_latest_directory` **move o `latest/` antigo para o lado** (não apaga in‑place); em falha, restaura. `[CÓDIGO]`

---

# 5. Contratos de dados

Convenção: **campo obrigatório** = presente e não‑nulo por construção do produtor.

## 5.1 `preflight_report.json` — artefato `ncbi_esearch_preflight`

Um objeto por execução de preflight. Dump de `EsearchPreflightResult`. `[CÓDIGO]`

| Campo | Tipo | Significado | Obrigatório |
|---|---|---|---|
| `search_query` | str | query textual enviada ao NCBI | sim |
| `translated_query` | str\|null | como o NCBI expandiu a query (`QueryTranslation`) | não (null se ausente) |
| `result_count` | int | `Count` do ESearch — nº de UIDs que um retrieval completo traria | sim |
| `history_web_env` / `history_query_key` | str\|null | handles da NCBI History API | não |
| `retrieved_at_utc` | str ISO‑8601 Z | instante da consulta | sim |
| `max_uid_count` | int | teto configurado | sim |
| `exceeds_max_uid_count` | bool | `result_count > max_uid_count` | sim |
| `sample_requested_count` | int | `min(sample_size, result_count)` | sim |
| `sample_uid_list` | list[str] | UIDs de amostra (também em `sample_protein_uids.txt`) | sim (pode ser `[]`) |
| `sample_record_count` | int | registros no XML de amostra | sim |
| `sample_records_with_sequence` / `_missing_sequence` | int | QC de presença de `GBSeq_sequence` | sim |
| `sample_records_with_extractable_uid` | int | QC de extração de `protein_uid` | sim |
| `sample_fetch_error` | str\|null | mensagem se a amostra falhou (não fatal) | não |
| `python_version` / `biopython_version` | str | ambiente | sim |

**Invariantes:** `sample_records_with_sequence + sample_records_missing_sequence == sample_record_count`; `exceeds_max_uid_count == (result_count > max_uid_count)`. `[CÓDIGO]`
**Upstream:** nenhum (é a cabeça do pipeline). **Downstream:** o notebook usa `translated_query` e `result_count` no `audit_summary`; **não** alimenta diretamente a aquisição de UID (que refaz o próprio ESearch).

## 5.2 `manifest.json` do preflight

Superset do report + `artifact_type`, `snapshot_format_version:"1.0"`, `snapshot_created_at_utc`, `immutable_snapshot_directory_name`, `immutable_snapshot_relative_path`, `manifest_file_name`, `output_files.{preflight_report_file,sample_uid_file}.{file_name,path,sha256}`. `[CÓDIGO]`

## 5.3 `protein_uids.txt` — artefato UID snapshot (pré‑Fase‑A)

Uma linha = um `protein_uid` (GI numérico, string de dígitos). Ordenado, deduplicado. `[EXEC-LOCAL]`
**Manifesto (`snapshot_format_version:"1.1"`):** `protein_uids_sha256` (via `sha256_of_lines`, sem dedupe/sort adicional na verificação), `ncbi_reported_result_count`, `raw_protein_uid_count`, `normalized_protein_uid_count`, `deduplicate_uids`, `sort_uids`, `uid_retrieval_strategy`, `page_size`, `history_web_env/query_key`, telemetria. `[EXEC-LOCAL]`
**Invariante observada:** `raw == normalized == ncbi_reported == 52473`. `[EXEC-LOCAL]`
**Downstream:** XML snapshot consome via `source_uid_sha256` + `source_uid_snapshot_manifest_sha256`.

## 5.4 `protein_records.xml` — artefato XML snapshot (pré‑Fase‑A)

GenPept em XML (`<GBSet><GBSeq>…`), 1 `<GBSeq>` por proteína. `[EXEC-LOCAL]`
**Manifesto (`1.1`, `artifact_type:"ncbi_protein_xml_snapshot"`):** `xml_file_sha256`, `consolidated_record_count`, `batch_count`, `batch_size`, `batches[]` (cada um: `batch_index`, `batch_start_index`, `batch_end_index`, `protein_uid_count`, `reused_from_workspace`, `xml_payload_sha256`), `rettype:"gp"`, `retmode:"xml"`, `request_policy` (deadline 300s, circuit breaker: threshold 3 / cooldown 60s, `max_request_starts_per_second:8.0`, backoff 0.1→30s ×2, `max_concurrent_requests:4`, `reuse_http_connection:false`), `retrieval_telemetry` (por estágio: `request_count`, `retry_count`, `failure_counts:{deadline,http_429,http_5xx,other,rate_limit,response_validation,timeout,truncated_response}`, `reused_batch_count`), `source_uid_sha256`, `source_uid_snapshot_manifest_sha256`. `[EXEC-LOCAL]`
**Invariante:** `consolidated_record_count == requested_protein_uid_count == source_uid_count`. `[EXEC-LOCAL]`

## 5.5 `protein_metadata.csv` — artefato metadata snapshot (pré‑Fase‑A)

**1 linha = 1 proteína** (`protein_uid` único). 148 colunas: `protein_uid`, `gbseq__accession_version`, `gbseq__definition`, `gbseq__length`, `gbseq__organism`, `gbseq__sequence`, `gbseq__source_db`, `gbseq__primary_accession`, `taxonomy__raw`, `taxonomy__01..10`, `reference__*`, e ~100 `feature__*` (cds, gene, protein, region, site, source, …). `[EXEC-LOCAL]`
**Colunas usadas a jusante na Fase A:** `protein_uid`, `gbseq__accession_version`, `gbseq__sequence`, `gbseq__length`. `[CÓDIGO]`
**Manifesto (`1.0`):** `csv_file_sha256`, `column_count:148`, `columns[]`, `observed_feature_keys`, `observed_feature_qualifiers`, `max_taxonomy_depth:10`, herda `search_query`/`translated_query`, `source_xml_snapshot_relative_path`. `[EXEC-LOCAL]`
**`qc_report.json` (`ncbi_protein_metadata_csv_qc`):** 5 checagens booleanas (ver E4). Todas `true` na execução real. `[EXEC-LOCAL]`

## 5.6 `query_recall_reference_set.csv` — fixture versionada

**1 linha = 1 proteína de referência conhecida** (pAgo ou PIWI‑RE). 12 colunas, 21 linhas. `[VERSIONADO]`

| Coluna | Tipo | Significado | Obrigatório |
|---|---|---|---|
| `accession` | str (accession.version) | identificador NCBI da proteína de referência | sim |
| `protein_short_name` | str | rótulo humano (TtAgo, RsAgo, PsPIWI‑RE, …) | sim |
| `organism` | str | organismo de origem | sim |
| `ago_family` | `PAGO` \| `PIWI_RE` | identidade biológica de família | sim |
| `clade` | `LONG_A` \| `LONG_B` \| `SHORT` \| `UNRESOLVED` | clado MID‑PIWI (PIWI‑RE → sempre `UNRESOLVED`) | sim |
| `sequence_sha256` | str (64 hex) | `sha256(strip_whitespace+upper(sequência))` — habilita o matching por identidade | sim (todas as 21 preenchidas) |
| `sequence_length` | int (como str) | nº de resíduos após normalização (informativo) | sim |
| `uniprot_accession` | str | cross‑ref UniProt (pode ser vazio) | não |
| `reference_label_source` | str | citação/base do rótulo | sim |
| `reference_label_evidence` | `EXPERIMENTAL` \| `LITERATURE_PHYLOGENETIC` \| `CURATED_COMPUTATIONAL` \| `DATABASE_ANNOTATION` | força da evidência do rótulo | sim |
| `verification_status` | `verified` \| `provisional` | status de conferência do mapeamento accession↔estudo | sim |
| `notes` | str | justificativa livre | sim |

**Invariantes (asseguradas por `test_query_reference_recall.py`):** `[VERSIONADO]`
- `accession` casa `\.\d+$` (todas têm versão explícita);
- `clade.upper() ⊆ {LONG_A, LONG_B, SHORT, UNRESOLVED, NA, ""}` — **`PIWI_RE` não é valor de `clade`**;
- toda linha `ago_family == PIWI_RE` tem `clade == UNRESOLVED`;
- `LONG_A`, `LONG_B`, `SHORT` cada um com ≥1 linha;
- ≥21 linhas; `sequence_sha256` das 21 casa `^[0-9a-f]{64}$`, **todas distintas**; `sequence_length` inteiro positivo;
- contagens de estrato: `LONG_A=8`, `LONG_B=2`, `SHORT=4`, `PIWI_RE=7`;
- `reference_label_evidence ⊆` os 4 níveis; `verification_status ⊆ {verified, provisional}`;
- linha `ABP72561.1`: `clade==LONG_B`, `sequence_sha256=="cbdb6bb6…c5fb36"`, `sequence_length=="777"`.

**SHA‑256 da fixture (HEAD):** `1da47d3f36db8f5a328f446804262e4d93ab017a49d7122602703ad3dea7a70b`. `[VERSIONADO]`

**Upstream:** curada à mão (ver Seção 16). **Downstream:** `query_reference_recall` a compara com o metadata CSV.

## 5.7 `reference_recall_summary.csv` — artefato recall

**1 linha = 1 estrato.** 5 linhas (`overall`, `LONG_A`, `LONG_B`, `SHORT`, `PIWI_RE`). `[CÓDIGO]`

| Coluna | Tipo | Significado |
|---|---|---|
| `stratum` | str | nome do estrato |
| `metric_name` | str | `overall_reference_recall`, `long_a_reference_recall`, … |
| `reference_count` | int | nº de referências no estrato |
| `exact_recovered_count` | int | Σ `recovered_exact_accession` |
| `equivalent_recovered_count` | int | Σ `recovered` |
| `exact_accession_recall` | float\|NaN | `exact/ref` (NaN se `ref==0`) |
| `retrieval_equivalent_recall` | float\|NaN | `equiv/ref` |
| `recall_status` | `EVALUABLE` \| `NOT_EVALUABLE` | |

## 5.8 `reference_recall_detail.csv` — artefato recall

**1 linha = 1 proteína de referência.** 21 linhas. `[CÓDIGO]`

| Coluna | Tipo | Significado |
|---|---|---|
| `accession` | str | accession da referência (da fixture) |
| `clade` | str (upper) | clado da referência |
| `ago_family` | str (upper) | família da referência |
| `sequence_sha256` | str | hash da fixture (lowercased) |
| `recovered` | bool | casou por qualquer método ≠ NONE |
| `recovered_exact_accession` | bool | casou por EXACT ou SAME_BASE |
| `match_method` | str | `EXACT_ACCESSION_VERSION` \| `SAME_BASE_ACCESSION` \| `SEQUENCE_SHA256` \| `NONE` |
| `matched_accession` | str | accession recuperado que casou ("" se NONE) |
| `matched_protein_uid` | str | `protein_uid` correspondente |
| `sequence_match_count` | int | nº de registros recuperados com o mesmo hash (0 salvo em SEQUENCE_SHA256) |
| `reference_label_source` | str | herdado da fixture |
| `reference_label_evidence` | str | herdado da fixture |

**Manifesto do recall (`1.1`, `artifact_type:"query_reference_recall"`):** `matching_strategy` (objeto), `matching_strategy_sha256`, `reference_count`, `exact_recovered_count`, `equivalent_recovered_count`, `stratum_exact_recall` (dict, `None`→null), `stratum_equivalent_recall` (dict), `stratum_recall_status` (dict), `query_recall_reference_set_csv_path`, `query_recall_reference_set_csv_sha256`, `source_metadata_csv_sha256`, `source_metadata_manifest_sha256`, `source_metadata_snapshot_directory_name/relative_path`, `search_query`, `translated_query`, `output_files.*`. `[CÓDIGO]`

## 5.9 `retained_protein_uids.txt` — artefato prefilter

Uma linha = um `protein_uid` retido, **na ordem do metadata** (não reordenado). `[CÓDIGO]`
**Papel downstream:** é a `selection_uid_list_file_name` que o `derived_protein_fasta` consome; seu `sha256` (`7bea84d9…`) reaparece como `derived.source_record_ids_sha256`. `[EXEC-LOCAL]`

## 5.10 `retained_records.csv` / `excluded_technical.csv` — artefato prefilter

**1 linha = 1 registro do metadata**, com 3 colunas extras adicionadas: `technical_prefilter_decision`, `technical_prefilter_reason`, `length_warning` (bool). `retained` = todas com `decision == "retain"`; `excluded` = as demais. `[CÓDIGO]` (na execução real: 52.473 colunas 151 = 148 + 3; `excluded_technical.csv` vazio). `[EXEC-LOCAL]`

## 5.11 `prefilter_counts.csv` — artefato prefilter

**1 linha por decisão possível (5).** Colunas `decision, count, fraction, description`. `Σ count == input_record_count`. `[CÓDIGO]`

## 5.12 `protein_sequences.fasta` — artefato derived FASTA

**1 registro = 1 proteína.** Defline: `>protein_uid=<uid>|accession=<acc>|length=<n>|organism=<org_com_underscores> <definição livre>`; sequência em linhas de 60 col, minúscula. `[CÓDIGO]` / `[EXEC-LOCAL]`
**`selection_report.json`:** `dataset_kind`, `record_selection_rule`, `record_order`, `requested_uid_count`, `resolved_uid_count`, `missing_uids[]`, `duplicate_requested_uids[]`, `empty_sequence_uids[]`, `source_record_ids_sha256`, `fasta_record_count`, `skipped_missing_sequence_count`. `[CÓDIGO]`
**Manifesto (`1.0`, `artifact_type:"derived_protein_fasta_snapshot"`):** ver E7 / Seção 6. `[CÓDIGO]`

## 5.13 Estrutura comum de todos os `manifest.json` da Fase A

`artifact_type`, `snapshot_format_version`, `snapshot_created_at_utc`, `manifest_file_name`, `immutable_snapshot_directory_name`, `immutable_snapshot_relative_path`, `output_files` (mapa `chave → {file_name, path, sha256}`) — mais campos específicos do estágio e a proveniência do(s) snapshot(s) pai por `*_sha256`. `[CÓDIGO]`

---

# 6. Sistema de snapshots e proveniência

## 6.1 Estrutura de diretórios

```
<snapshot_root>/
  snapshots/
    <YYYY-MM-DDThh-mm-ssZ>__q_<sha256(query)[:12]>/   ← IMUTÁVEL: nunca reescrito
      manifest.json
      <arquivos de saída do estágio>
  latest/                                             ← cópia de conveniência, substituível
    manifest.json
    <mesmos arquivos>
```

`[CÓDIGO]` (`save_*` faz `immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)` — falha se já existe; nunca sobrescreve um snapshot imutável).

## 6.2 Ciclo de `save_*`

1. `mkdir(exist_ok=False)` do diretório imutável (timestamp + hash de query no nome).
2. Escreve cada saída **atomicamente** (`NamedTemporaryFile` no mesmo diretório → `replace`).
3. Monta o manifesto com `sha256_of_file` de cada saída + proveniência dos pais.
4. `write_json_atomic(manifest)`.
5. `immutable_snapshot_complete = True`.
6. Se `update_latest_directory`: `_replace_latest_directory` (stage num tempdir → move `latest/` antigo para o lado → move o stage para `latest/`; restaura em falha).
7. `except Exception` **antes** do passo 5: `shutil.rmtree` do diretório parcial + re‑raise (não deixa snapshot meia‑boca).

`[CÓDIGO]` (padrão idêntico nos 4 `_snapshot.py` da Fase A).

## 6.3 `manifest.json` e SHA‑256

Cada manifesto grava:
- o `sha256` de **cada arquivo de saída próprio** (`output_files`);
- o `sha256` do **manifesto de cada snapshot pai** (`source_<x>_manifest_sha256`);
- o `sha256` do **CSV/arquivo‑chave de cada pai** (`source_metadata_csv_sha256`, etc.);
- fingerprints derivados (`source_record_ids_sha256`, `record_selection_config_sha256`, `matching_strategy_sha256`, `technical_prefilter_policy_sha256`).

## 6.4 Compatibilidade / `reuse_latest_or_create`

`latest_<x>_snapshot_is_available(...)` retorna `True` **somente se todas** as condições valem `[CÓDIGO]`:

| Estágio | Condições de reuso |
|---|---|
| preflight | `latest/manifest.json` existe; `artifact_type` ok; `snapshot_format_version == "1.0"`; cada `output_files[*].sha256` bate com o disco; (se dado) `manifest["search_query"] == search_query`. |
| recall | idem + `snapshot_format_version == "1.1"` + `matching_strategy_sha256 == build_matching_strategy_sha256()` + (se dado) `query_recall_reference_set_csv_sha256` bate + (se dado) `source_metadata_manifest_sha256` bate com o metadata `latest`. |
| prefilter | idem base + `technical_prefilter_policy_sha256 == build_technical_prefilter_policy_sha256()` + (se dado) `source_metadata_manifest_sha256` bate. |
| derived FASTA | idem base (versão `"1.0"`) + `fasta_file_sha256` bate + reparse de deflines: `fasta_record_count` e `source_record_ids_sha256` (por `record_order`) batem + (se dado) `source_metadata_manifest_sha256` **e** `source_selection_manifest_sha256` batem. |

`resolve_*` em `reuse_latest_or_create`: se disponível → carrega `latest/` (rápido, sem rede/recompute); senão → `save_*` (cria snapshot novo).
`resolve_*` em `reuse_latest`: erro `FileNotFoundError` se indisponível.
`resolve_*` em `create_new`: sempre `save_*`.

## 6.5 Cadeia explícita de proveniência (execução real)

`[EXEC-LOCAL]` — todos os hashes abaixo lidos dos manifestos em `data/…` (não versionados):

```
SEARCH_QUERY = "(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])"
   │  sha256(query)[:12] = c9c3315a67a9   (usado por preflight, uid, xml, metadata)
   ▼
[E1] preflight  2026-08-30T23-57-17Z__q_c9c3315a67a9
     result_count = 52473 ; exceeds_max_uid_count = false (max 250000)
     translated_query = "(PIWI[All Fields] OR Argonaute[All Fields]) AND
                         ((\"Bacteria\"[Organism] OR \"Bacteria Latreille et al. 1825\"[Organism]) OR \"Archaea\"[Organism])"
   ▼  (o preflight NÃO passa hash adiante; a aquisição de UID refaz o ESearch)
[E2] uid        2026-08-30T23-59-40Z__q_c9c3315a67a9
     normalized_protein_uid_count = 52473
     protein_uids_sha256 = 949bcaabc6cae0a60238bcbd4ae88608dbaf6db1fa3557db5cd0ec96df5cb253
   ▼  source_uid_sha256 = 949bcaab… ; source_uid_snapshot_manifest_sha256 = feda9bb7…
[E3] xml        2026-08-31T00-00-11Z__q_c9c3315a67a9
     consolidated_record_count = 52473 ; batch_count = 525
     xml_file_sha256 = 777390433da45c16c10bf0958dca429f03132e517dd10b13311a192268c96169
   ▼  source_xml_snapshot_relative_path = snapshots/2026-08-31T00-00-11Z__q_c9c3315a67a9
[E4] metadata   2026-08-31T00-21-18Z__q_c9c3315a67a9
     row_count = 52473 ; column_count = 148
     csv_file_sha256      = b39096f339a4079750b0880d49ae7a9221d07a16943bf476b1e27786a6541510
     manifest_sha256      = d6ab54d551432082db5b25502b1aacd1495f6a4eeb2de7892cf6735f56c80994
   ├──────────────────────────────────────────────────────────────────────┐
   ▼                                                                       ▼
[E5] recall  2026-08-31T01-30-37Z__q_8994498899ae                     [E6] prefilter  2026-08-31T00-23-09Z__q_bbdf63d86832
     source_metadata_csv_sha256      = b39096f3…  ✔                        source_metadata_csv_sha256   = b39096f3…  ✔
     source_metadata_manifest_sha256 = d6ab54d5…  ✔                        source_metadata_manifest_sha256 = d6ab54d5…  ✔
     query_recall_reference_set_csv_sha256 = 1da47d3f…  (== fixture HEAD)  input_record_count = retained = 52473 ; excluded = 0
     matching_strategy_sha256 = 3460b048…                                  technical_prefilter_policy_sha256 = 032fbc72…
     reference_count 21 ; exact 20 ; equivalent 21                         retained_protein_uids.txt sha256 = 7bea84d983c58bff09b9dc6c539377982a939c2ce9545b3a35e4bb96acde8ef0
                                                                             │
                                                                             ▼
                                                       [E7] derived FASTA  2026-08-31T00-23-28Z__q_c99c5e155233
                                                            derived_from_artifact_type = pago_technical_prefilter
                                                            derived_from_manifest_sha256   = dc16489b4a09d31f62955e4ba184f73d8d6ff778067d1e9d5472293d6215514d
                                                            source_selection_manifest_sha256 = dc16489b…  ✔ (== derived_from)
                                                            record_selection_config_sha256 = 032fbc72…   ✔ (== prefilter policy)
                                                            source_record_ids_sha256       = 7bea84d9…   ✔ (== retained_protein_uids.txt)
                                                            source_metadata_csv_sha256     = b39096f3…   ✔
                                                            fasta_record_count = 52473
                                                            fasta_file_sha256 = f7bac76b982bfe71d621416df493dc77ecda6fd3ef29b9b314cf21838954ce43
```

**Todos os elos verificados batem.** A auditoria recomputou na hora `sha256_of_file` de `retained_records.csv`, `retained_protein_uids.txt` e `protein_sequences.fasta` em disco — idênticos aos manifestos. `[EXEC-LOCAL]`

## 6.6 O que impede reusar um snapshot incompatível — evidência real

Há **dois** snapshots de recall em disco `[EXEC-LOCAL]`:

| Dir | `snapshot_format_version` | `matching_strategy_sha256` | `query_recall_reference_set_csv_sha256` | `recovered_count` / readings |
|---|---|---|---|---|
| `2026-08-31T00-22-58Z__q_8994498899ae` (antigo) | `1.0` | ausente (`None`) | `e239679ed2f0c949f2589a3b91f74d9da09952faf24b92e6734e367b714d46bf` (CSV pré‑`0013c6d`) | `recovered_count: 20`; `stratum_recall.long_b = 0.5` |
| `2026-08-31T01-30-37Z__q_8994498899ae` (**latest**) | `1.1` | `3460b048…` | `1da47d3f…` (CSV do HEAD) | `exact 20 / equivalent 21`; `stratum_exact_recall.long_b = 0.5`, `stratum_equivalent_recall.long_b = 1.0` |

`[INFERÊNCIA]` A sequência real foi: (a) execução completa do pipeline ~00:00–00:23Z de 31/08 produziu o snapshot antigo de recall com a metodologia só‑accession e o CSV de 11+/sem `sequence_sha256`; (b) o commit `0013c6d` (01:26Z) trocou a metodologia (`1.0→1.1`, `matching_strategy_sha256`) **e** o CSV (novas colunas → novo sha `1da47d3f`); (c) ao reexecutar o notebook 10, `latest_query_reference_recall_snapshot_is_available` retornou `False` por **três** razões independentes (versão `≠1.1`, `matching_strategy_sha256` ausente, `csv_sha256` divergente), forçando `save_*` a produzir o snapshot `01-30-37Z`. Os demais estágios (`preflight/uid/xml/metadata/prefilter/derived`) **não** foram reexecutados porque seus módulos não mudaram em `0013c6d` e seus snapshots continuaram válidos — daí terem timestamp de 00:xxZ.

## 6.7 Papel de `matching_strategy_sha256` (introduzido em `0013c6d`)

- Gravado no manifesto do recall e **verificado em `_validate_loaded_*`**: um snapshot cujo `matching_strategy_sha256` difere do valor corrente do código é tratado como indisponível → recomputa. `[CÓDIGO]`
- Não é consumido por **nenhum** estágio a jusante (o recall é folha do grafo). Logo, mudar a metodologia de matching **invalida apenas o snapshot de recall** e nada mais. `[CÓDIGO]` / `[INFERÊNCIA]`
- Combinado com `snapshot_format_version == "1.1"` (checado em separado), garante que um snapshot escrito por metodologia antiga nunca seja lido como se fosse novo. `[CÓDIGO]`

---

# 7. Aquisição NCBI

## 7.1 Parâmetros efetivos no notebook 10 (CELL 3)

`[VERSIONADO]`

| Parâmetro | Valor | Estágio |
|---|---|---|
| `SEARCH_QUERY` | `(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])` | todos |
| `MAX_UID_COUNT` | `250_000` | preflight (guarda) |
| `PREFLIGHT_SAMPLE_SIZE` | `200` | preflight |
| `ALLOW_EXCEEDS_MAX_UID_COUNT` | `False` | preflight |
| `UID_PAGE_SIZE` | `10_000` | UID |
| `UID_MAX_RETRY_ATTEMPTS` | `5` | UID |
| `UID_REQUEST_DELAY_SECONDS` | `None` | UID |
| `UID_FETCH_TIMEOUT_SECONDS` | `30.0` | UID |
| `UID_REQUEST_DEADLINE_SECONDS` | `300.0` | UID |
| `XML_BATCH_SIZE` | `100` | XML |
| `XML_MAX_RETRY_ATTEMPTS` | `5` | XML |
| `XML_MAX_CONCURRENT_REQUESTS` | `4` | XML |
| `XML_REUSE_HTTP_CONNECTION` | `False` | XML |
| `XML_ENABLE_BATCH_RESUME` | `True` | XML |
| `XML_PURGE_BATCH_WORKSPACE_ON_SUCCESS` | `True` | XML |
| `SEQUENCE_LINE_WIDTH` | `60` | derived FASTA |
| `UPDATE_LATEST_DIRECTORY` | `True` | todos |
| todos os `*_SNAPSHOT_MODE` | `SnapshotMode.reuse_latest_or_create` | todos |

**Valores efetivamente registrados na execução real (manifestos):** UID `page_size=10000`, `request_delay_seconds=0.1` (o `None` do notebook virou o default do módulo), `fetch_timeout_seconds=30.0`, `request_deadline_seconds=300.0`; XML `batch_size=100`, `max_concurrent_requests=4`, `max_retry_attempts=5`, `reuse_http_connection=false`. `[EXEC-LOCAL]`

## 7.2 Como ESearch / History / UID / EFetch se conectam

`[CÓDIGO]` + `[EXEC-LOCAL]`:

1. **Preflight** — `esearch(db=protein, term=Q, retmax=0, usehistory="y")` retorna `Count`, `QueryTranslation`, `WebEnv`, `QueryKey`. Se `Count>0`: `efetch(rettype=uilist, retmax=200, WebEnv=…, query_key=…)` → 200 UIDs; `efetch(id=…, rettype=gb, retmode=xml)` → QC. O WebEnv do preflight é descartado.
2. **UID snapshot** — estratégia `esearch_history_efetch_uilist`: faz **o seu próprio** `esearch(usehistory=y)` (novo `WebEnv`, `MCID_6a94c3ecfbcf6bb364054dbf`) e depois 6 `efetch(rettype=uilist, retstart=k·10000, retmax=10000)` para paginar os 52.473 UIDs. Deduplica e ordena. `esearch_request_count=1`, `efetch_request_count=6`.
3. **XML snapshot** — lê `protein_uids.txt` do snapshot de UID; monta 525 lotes de 100; `efetch(id=<100 ids>, rettype=gp, retmode=xml)` por lote, até 4 em paralelo; consolida num XML só. Cada lote tem `xml_payload_sha256`; se `enable_batch_resume`, lotes já baixados em `.batch_workspace/` são reusados (0 reusos nesta execução — foi uma corrida única).
4. **Metadata snapshot** — lê o XML consolidado, achata em CSV de 148 colunas, roda QC.

Observação: a **cadeia de identidade** entre UID→XML→metadata é por `sha256` de arquivo e de manifesto (Seção 6.5), não por `WebEnv` (que é efêmero e difere entre preflight e UID).

## 7.3 Mecanismos de robustez — o que é PRÉ‑Fase‑A e o que é Fase A

| Mecanismo | Onde | Origem |
|---|---|---|
| Retry com backoff exponencial | `_run_with_retries` no preflight (`5·2^k`, 5 tentativas) | **Fase A** (`ncbi_esearch_preflight.py`) `[CÓDIGO]` |
| Retry/backoff no UID e XML (`retry_backoff_initial 0.1 → max 30s ×2`) | `request_policy` dos snapshots NCBI | **pré‑Fase‑A** `[EXEC-LOCAL]` |
| `fetch_timeout_seconds` (30s por requisição) | UID e XML | pré‑Fase‑A `[EXEC-LOCAL]` |
| `batch_deadline_seconds` / `request_deadline_seconds` (300s) | UID e XML | pré‑Fase‑A `[EXEC-LOCAL]` |
| Response validation (`response_validation`, `truncated_response` nos contadores) | XML | pré‑Fase‑A `[EXEC-LOCAL]` |
| Rate limiting (`max_request_starts_per_second: 8.0`) | XML | pré‑Fase‑A `[EXEC-LOCAL]` |
| Circuit breaker (`failure_threshold: 3`, `cooldown: 60s`) | XML | pré‑Fase‑A `[EXEC-LOCAL]` |
| Resumable XML batch workspace (`.batch_workspace/`, `enable_batch_resume`) | XML | pré‑Fase‑A (`24cba76 (feat) add resumable NCBI XML batch workspace`) `[CÓDIGO]` |
| Telemetria de requisições (latências, `failure_counts`) | UID e XML | pré‑Fase‑A (`50016f0`) `[EXEC-LOCAL]` |
| Guarda `max_uid_count` + relatório de preflight materializado para auditoria | preflight | **Fase A** `[CÓDIGO]` |
| Snapshot imutável + `latest/` + `reuse_latest_or_create` | todos | pré‑Fase‑A (padrão do projeto), **aplicado** aos 4 novos estágios pela Fase A `[CÓDIGO]` |

**Uso real do batch resume nesta execução:** `reused_batch_count = 0`, todos os 525 lotes com `reused_from_workspace: false` — a corrida XML foi única e completada sem retomada. O workspace foi então purgado (`purge_batch_workspace_on_success=True`). `[EXEC-LOCAL]`

## 7.4 Erros transitórios observados na execução real do XML

`[EXEC-LOCAL]` (telemetria do manifesto XML, estágio `xml_batches`):

```
request_count        = 525
retry_count          = 9
failure_counts:
  http_5xx           = 6      ← respostas 5xx do NCBI (ex. 502 Bad Gateway)
  truncated_response  = 3      ← leitura incompleta do corpo (IncompleteRead)
  http_429 / rate_limit / timeout / deadline / response_validation / other = 0
reused_batch_count   = 0
wall_seconds_total   ≈ 1170.5  (~19,5 min)
```

`[INFERÊNCIA]` As 9 tentativas extras (6 5xx + 3 respostas truncadas) foram **todas recuperadas pelo retry/backoff** — o resultado final tem `consolidated_record_count = 52473 == requested_protein_uid_count`, `batch_count = 525` completos, e o QC do metadata confirma `row_count_matches_source_xml = true`. Não há evidência de perda de registro.

---

# 8. Query reference recall — auditoria detalhada

## 8.1 Objetivo científico

Responder: **"o texto da query `(PIWI OR Argonaute) [All Fields] AND (Bacteria OR Archaea) [Organism]` recupera as pAgo/PIWI‑RE que já conhecemos?"** — estratificado por clado MID‑PIWI (`LONG_A`/`LONG_B`/`SHORT`) e por família (`PIWI_RE`). É um **teste de cobertura da query**, não uma afirmação de atividade bioquímica (por isso cada linha carrega `reference_label_evidence`). `[VERSIONADO]` (docstring do módulo e cabeçalho de `query_recall_reference_set_curation_notes.md`).

O estágio **não filtra** o dataset — só produz `summary.csv` e `detail.csv`. `[CÓDIGO]`

## 8.2 Painel de 21 referências

`[VERSIONADO]` (`tests/fixtures/query_recall_reference_set.csv`):

| # | accession | nome | ago_family | clade | seq_len | evidência | status |
|---|---|---|---|---|---|---|---|
| 1 | WP_011174533.1 | TtAgo | PAGO | LONG_A | 685 | EXPERIMENTAL | verified |
| 2 | WP_010870838.1 | MjAgo | PAGO | LONG_A | 713 | EXPERIMENTAL | verified |
| 3 | WP_011011654.1 | PfAgo | PAGO | LONG_A | 770 | EXPERIMENTAL | verified |
| 4 | WP_014295921.1 | MpAgo | PAGO | LONG_A | 639 | EXPERIMENTAL | verified |
| 5 | WP_072865986.1 | MhAgo | PAGO | LONG_A | 640 | EXPERIMENTAL | verified |
| 6 | WP_010880937.1 | AaAgo | PAGO | LONG_A | 706 | EXPERIMENTAL | verified |
| 7 | WP_011378069.1 | SeAgo | PAGO | LONG_A | 735 | EXPERIMENTAL | verified |
| 8 | WP_060384876.1 | TpsAgo | PAGO | LONG_A | 685 | LITERATURE_PHYLOGENETIC | provisional |
| 9 | **ABP72561.1** | RsAgo | PAGO | LONG_B | 777 | EXPERIMENTAL | verified |
| 10 | WP_005580376.1 | NgAgo | PAGO | LONG_B | 887 | LITERATURE_PHYLOGENETIC | provisional |
| 11 | WP_010878815.1 | AfAgo | PAGO | SHORT | 427 | EXPERIMENTAL | verified |
| 12 | WP_010942012.1 | GsAgo | PAGO | SHORT | 473 | EXPERIMENTAL | provisional |
| 13 | WP_109649955.1 | MapAgo | PAGO | SHORT | 507 | EXPERIMENTAL | provisional |
| 14 | WP_012735993.1 | SiAgo | PAGO | SHORT | 459 | LITERATURE_PHYLOGENETIC | provisional |
| 15 | WP_014597637.1 | PsPIWI‑RE | PIWI_RE | UNRESOLVED | 783 | EXPERIMENTAL | verified |
| 16 | WP_027844734.1 | Mtes‑pPIWI‑RE | PIWI_RE | UNRESOLVED | 1052 | CURATED_COMPUTATIONAL | verified |
| 17 | WP_017749591.1 | Shof‑pPIWI‑RE | PIWI_RE | UNRESOLVED | 925 | CURATED_COMPUTATIONAL | verified |
| 18 | WP_012408073.1 | Npun‑pPIWI‑RE | PIWI_RE | UNRESOLVED | 919 | CURATED_COMPUTATIONAL | verified |
| 19 | WP_012163488.1 | Amar‑pPIWI‑RE | PIWI_RE | UNRESOLVED | 1031 | CURATED_COMPUTATIONAL | verified |
| 20 | WP_015099491.1 | Sesp‑pPIWI‑RE | PIWI_RE | UNRESOLVED | 830 | CURATED_COMPUTATIONAL | verified |
| 21 | WP_012851864.1 | Tcur‑pPIWI‑RE | PIWI_RE | UNRESOLVED | 965 | CURATED_COMPUTATIONAL | verified |

**Contagem por estrato:** `LONG_A = 8`, `LONG_B = 2`, `SHORT = 4`, `PIWI_RE = 7`. Total 21. (`overall` = as 21.) `[VERSIONADO]`

**Níveis de evidência presentes:** `EXPERIMENTAL` (13), `LITERATURE_PHYLOGENETIC` (3), `CURATED_COMPUTATIONAL` (6 — todos PIWI‑RE), `DATABASE_ANNOTATION` (0). `[VERSIONADO]`

## 8.3 `NOT_EVALUABLE`

`_recall(reference_count, recovered_count)`: se `reference_count == 0` → `(None, "NOT_EVALUABLE")`; senão `(k/n, "EVALUABLE")`. No manifesto, `None` vira **JSON null** (não `0.0`). `[CÓDIGO]` No painel atual **todos os 5 estratos têm ≥1 referência**, então todos são `EVALUABLE` na execução real. `[EXEC-LOCAL]` O comportamento `NOT_EVALUABLE` está coberto por `test_empty_stratum_is_not_evaluable_not_zero`. `[VERSIONADO]`

## 8.4 Hierarquia de matching

`[CÓDIGO]` — testada uma por uma:

```
para cada referência:
  1. EXACT_ACCESSION_VERSION  ── se reference.accession ∈ retrieved_versioned_set
  2. SAME_BASE_ACCESSION      ── senão, se base(reference.accession) ∈ retrieved_by_bare
                                  (base = accession sem sufixo .\d+ ; representante = menor (acc,uid) lexicográfico)
  3. SEQUENCE_SHA256          ── senão, se reference.sequence_sha256 ∈ sequence_sha256_to_hits
                                  (representante = hits[0] após sort ; sequence_match_count = len(hits))
  4. NONE                     ── senão
recovered                 = método ≠ NONE
recovered_exact_accession = método ∈ {EXACT_ACCESSION_VERSION, SAME_BASE_ACCESSION}
```

Prioridade EXACT/SAME_BASE **sobre** SEQUENCE_SHA256 é testada em `test_exact_accession_takes_priority_over_sequence_hash` e `test_same_base_accession_takes_priority_over_sequence_hash`. `[VERSIONADO]`

## 8.5 Normalização de sequência e SHA‑256

- `SEQUENCE_NORMALIZATION = "strip_all_whitespace_then_uppercase"` (constante do módulo). `[CÓDIGO]`
- `normalize_protein_sequence(v) = "".join(str(v).split()).upper()`. `[CÓDIGO]`
- `protein_sequence_sha256(v) = sha256(normalize(v).encode("utf-8")).hexdigest()`. `[CÓDIGO]`
- A **mesma** normalização é aplicada às sequências recuperadas (`gbseq__sequence`) antes de hashear. `[CÓDIGO]`
- Reversível‑free: nenhuma substituição de resíduos, gaps ou ambiguidades. `[CÓDIGO]`

## 8.6 Múltiplas sequências idênticas & `sequence_match_count`

`sequence_sha256_to_hits[digest]` é uma **lista ordenada** (`hits.sort()`) de `(accession, protein_uid)`. Em caso de match por sequência: `matched_accession, matched_protein_uid = hits[0]` (menor accession lexicográfico → determinístico) e `sequence_match_count = len(hits)`. `[CÓDIGO]` Testado em `test_multiple_accessions_same_hash_resolved_deterministically` (3 accessions, roda 2×, incl. shuffle → sempre `AA_1.1`, count 3). `[VERSIONADO]`

## 8.7 `exact_accession_recall` vs `retrieval_equivalent_recall`

| Reading | Conta como recuperado | Interpretação |
|---|---|---|
| `exact_accession_recall` (`stratum_exact_recall`) | só EXACT + SAME_BASE | "a query trouxe **o mesmo accession**" — leitura estrita, comparável a benchmarks de accession |
| `retrieval_equivalent_recall` (`stratum_equivalent_recall`) | EXACT + SAME_BASE + SEQUENCE_SHA256 | "a query trouxe **a mesma proteína**, mesmo que sob outro accession" — leitura biológica |

Ambos são gravados no manifesto e impressos no notebook CELL 9. `[CÓDIGO]`

## 8.8 `matching_strategy` / `matching_strategy_sha256`

Ver Seção 6.7. Payload literal `strategy_version: "2.0"`; sha `3460b048fc6de363ddf9282c2943a44c284e51d2c48092ca14426600b2871a08`. `[CÓDIGO]` Pinado no teste `test_matching_strategy_sha256_is_pinned`. `[VERSIONADO]` Gravado e verificado no manifesto do recall. `[EXEC-LOCAL]`

## 8.9 Invalidação de snapshots antigos

`_validate_loaded_query_reference_recall_payload` rejeita, cada um independentemente: `artifact_type` errado; `snapshot_format_version != "1.1"`; `matching_strategy_sha256` divergente; qualquer `output_files[*].sha256` que não bate. `latest_*_is_available` também rejeita se `query_recall_reference_set_csv_sha256` ou `source_metadata_manifest_sha256` divergirem. `[CÓDIGO]` Evidência prática: o snapshot de recall `00-22-58Z` (v1.0, sem strategy sha, CSV `e239679e…`) **não** foi reusado; o notebook gerou `01-30-37Z`. `[EXEC-LOCAL]` (Seção 6.6).

## 8.10 Caso RsAgo `ABP72561.1` → `A4WYU7.1`

`[EXEC-LOCAL]` (`reference_recall_detail.csv`, linha 9) + `[VERSIONADO]` (notas de curadoria):

- `ABP72561.1` (RsAgo, GenBank, *Cereibacter sphaeroides* ATCC 17025, 777 aa; sem registro RefSeq `WP_`) **não** está entre os 52.473 accessions recuperados.
- `A4WYU7.1` (Swiss‑Prot, mesmo Identical Protein Group, sequência **byte‑idêntica**, 777 aa) **está** — `protein_uid = 2500461169`.
- Resultado no `detail.csv`: `match_method = SEQUENCE_SHA256`, `matched_accession = A4WYU7.1`, `matched_protein_uid = 2500461169`, `sequence_match_count = 1`, `recovered = True`, `recovered_exact_accession = False`.
- Efeito nos números: `stratum_exact_recall.long_b = 0.5` (1/2), `stratum_equivalent_recall.long_b = 1.0` (2/2); `overall` `0.952` (20/21) vs `1.0` (21/21).
- O hash `cbdb6bb64718c9e8ca78a34ac8445eff1556cb87b5ad687026373ed401c5fb36` da fixture foi derivado **offline** da própria `protein_metadata.csv` da execução real, usando a sequência de `A4WYU7.1` como fonte (mesma IPG, byte‑idêntica). `[VERSIONADO]` (notas de curadoria, seção "How the 21 hashes were derived").
- `ABP72561.1` permanece como accession da linha na fixture — é um accession correto para RsAgo.

## 8.11 Limitações científicas remanescentes do recall

`[INFERÊNCIA]` a partir de `[VERSIONADO]` (notas de curadoria) + `[CÓDIGO]`:

1. **O painel tem 21 proteínas fortemente enriquecidas** (quase todas experimentais ou de literatura). `recall = 21/21` mede recuperação **neste painel**, não a sensibilidade da query sobre o universo desconhecido de pAgos.
2. **Viés de anotação embutido**: proteínas conhecidas tendem a estar bem anotadas no RefSeq/CDD, exatamente o que `[All Fields]` indexa. O painel é o melhor caso para a query.
3. **PIWI‑RE `CURATED_COMPUTATIONAL` (6/7)** foi selecionado *usando* modelos de perfil PIWI‑RE (Pfam `pPIWI_RE_X`/`MID_pPIWI_RE`/`RNaseH_pPIWI_RE`). Isso **não** é circular para a pergunta "a query **textual** recupera?" (é para uma futura validação de HMM PIWI‑RE — explicitamente vedado nas notas).
4. **Sem rota sequence‑based**: descoberta por homologia (HMM/PSI‑BLAST sobre RefSeq) é trabalho de fase futura, documentado.
5. **`clade` das referências LONG/SHORT** vem de literatura filogenética, não de um placement feito neste projeto; 6 das 21 são `provisional`.

---

# 9. Technical prefilter

## 9.1 Política integral (`build_technical_prefilter_policy_payload`)

`[CÓDIGO]` / `[EXEC-LOCAL]` (`technical_prefilter_policy` no manifesto real):

```
policy_kind            = "pago_technical_prefilter"
policy_version         = "1.0"
allowed_residues       = A C D E F G H I K L M N P Q R S T V W Y      (20 canônicos)
tolerated_ambiguous    = B J O U X Z                                  (Asx, Xle, pyrrolisina, selenocisteína, qualquer, Glx)
tolerated_symbols      = * -                                          (terminador, gap)
length_warning_min     = 200
length_warning_max     = 2000
protein_uid_column     = "protein_uid"
sequence_column        = "gbseq__sequence"
length_column          = "gbseq__length"
exclusion_reasons      = [drop_unprocessable_record, drop_technical_duplicate,
                          drop_missing_sequence, drop_invalid_sequence_characters]
policy_sha256          = 032fbc727b68ceb97cdda00ca5764db3ddeebc0868058d0c73883c23771fa523
```

## 9.2 O que **pode** excluir (4 causas, primeira‑regra‑vence nesta prioridade)

`[CÓDIGO]`:

| Prioridade | Decisão | Condição |
|---|---|---|
| 1 (maior) | `drop_unprocessable_record` | `protein_uid` nulo/NaN ou string vazia após strip |
| 2 | `drop_technical_duplicate` | `protein_uid` repetido (`duplicated(keep="first")`), **e** não vazio |
| 3 | `drop_missing_sequence` | sequência normalizada (`"".join(str(v).split()).upper()`) == `""` |
| 4 | `drop_invalid_sequence_characters` | há caractere fora de (`allowed ∪ ambiguous ∪ symbols`) após normalização |
| — | `retain` | nenhuma das acima |

## 9.3 O que **propositalmente NÃO** exclui

`[CÓDIGO]` (constantes, docstring, `notes` da política) + `[VERSIONADO]` (`test_retains_records_regardless_of_annotation_text_or_length`):

- **Comprimento da sequência** — fora de `[200, 2000]` só liga `length_warning = True`; o registro é retido. "A varredura de domínio, não o prefiltro, decide se um comprimento é compatível com uma pAgo."
- **Texto de anotação** — `gbseq__definition`, `feature__*__note`, `product`, etc. **nunca** são lidos. Um registro anotado como "SAM‑dependent methyltransferase" ou "transposase" é retido (testado explicitamente).
- **Biologia / família / domínios** — nenhuma inferência.
- **Resíduos ambíguos/não‑canônicos** `BJOUXZ` e símbolos `*-` — retidos (testado: `"MKTXBZUO*-"` → retido).

## 9.4 Casos específicos

`[CÓDIGO]` / `[VERSIONADO]`:

| Caso | Tratamento |
|---|---|
| `protein_uid` ausente/vazio | `drop_unprocessable_record` (prioridade máxima — mesmo que a sequência também esteja vazia) |
| `protein_uid` duplicado | 1ª ocorrência retida, demais `drop_technical_duplicate` |
| sequência vazia / só espaços | `drop_missing_sequence` |
| sequência com dígitos, `@`, etc. | `drop_invalid_sequence_characters` |
| comprimento 3 aa, ou 5000 aa | **retido** com `length_warning = True` |
| coluna `gbseq__length` ausente | `length_warning` fica `True` para todas (via `numeric_length.isna()`) — ainda retém |

## 9.5 Por que 52.473 / 52.473 retidos é consistente com o design

`[INFERÊNCIA]` a partir de `[EXEC-LOCAL]` + `[CÓDIGO]`:

- A fonte é `protein_metadata.csv` já **pós‑QC** do metadata snapshot, cujo `qc_report.json` confirma `protein_uid_has_no_duplicates = true` e `protein_uid_has_no_empty_values = true`. Logo `drop_unprocessable_record` e `drop_technical_duplicate` são **necessariamente 0**.
- O preflight amostrou 200 registros e encontrou `sample_records_missing_sequence = 0` e `sample_records_with_extractable_uid = 200`. Registros do `db=protein` do NCBI com `rettype=gp` praticamente sempre trazem `GBSeq_sequence`. `drop_missing_sequence = 0` no dataset completo.
- Sequências de proteína do NCBI usam o alfabeto IUPAC — coberto por `allowed ∪ BJOUXZ ∪ *-`. `drop_invalid_sequence_characters = 0`.
- Portanto o design ("só técnico, sem comprimento, sem texto") **prevê** que um dataset NCBI limpo passe inteiro. Os 8.299 registros com `length_warning = True` (comprimento fora de [200,2000]; ver Seção 13) foram **todos retidos** — exatamente a intenção.

`counts_by_decision` real: `retain: 52473`, os 4 drops: `0`. `excluded_technical.csv` vazio. `[EXEC-LOCAL]`

---

# 10. Derived FASTA

## 10.1 Por que foi criado

`[CÓDIGO]` (docstrings, `test_resolve_produces_honest_provenance_manifest`):

O SWeeP a jusante consome um snapshot de FASTA. Antes da Fase A só existia `ncbi_protein_fasta_snapshot` (o proteoma NCBI inteiro). A Fase A precisa embedar **o subconjunto retido pelo prefiltro**. Reusar `ncbi_protein_fasta_snapshot` para isso seria desonesto quanto à proveniência (o artefato não teria como dizer "sou um recorte de X pela regra Y"). Daí um `artifact_type` **próprio** — `derived_protein_fasta_snapshot` — com metadados de transformação.

## 10.2 De qual estágio deriva

Deriva de **dois** pais `[CÓDIGO]`:
- **metadata snapshot** (fonte das sequências) — `source_metadata_csv_sha256`, `source_metadata_manifest_sha256`;
- **selection snapshot** = o `pago_technical_prefilter` (fonte da lista de UIDs) — `derived_from_artifact_type = "pago_technical_prefilter"`, `derived_from_manifest_sha256`, `source_selection_manifest_sha256`.

No notebook, `record_selection_rule = "pago_technical_prefilter.retained"` e `record_selection_config_sha256 = prefilter_manifest["technical_prefilter_policy_sha256"]`. `[VERSIONADO]` (CELL 11).

## 10.3 Como garante que só os UIDs retidos entram

`[CÓDIGO]`:
1. `_load_selection_snapshot` lê `latest/retained_protein_uids.txt` do prefilter e **exige** `manifest["artifact_type"] == "pago_technical_prefilter"` (senão `RuntimeError`).
2. `build_derived_fasta_selection` recebe **exatamente** essa lista como `selected_protein_uids`.
3. UID pedido ausente do metadata → `RuntimeError` (salvo `drop_missing_uids`, que o notebook **não** liga).
4. `resolved_uid_count` e `fasta_record_count` no manifesto; na execução real `requested = resolved = fasta_record_count = 52473`, `missing_uid_count = 0`. `[EXEC-LOCAL]`

## 10.4 Como preserva ordem

`record_order = "as_selected"` (default e usado no notebook): os registros saem **na ordem da lista de UIDs retidos**, que por sua vez é a ordem do metadata. `"sorted_by_uid"` existe como alternativa mas não é usado. `[CÓDIGO]` / `[VERSIONADO]`

## 10.5 Validação por hash

`[CÓDIGO]` (`_validate_loaded_derived_protein_fasta_payload`, executado a cada `load`):
- `fasta_file_sha256` vs bytes do arquivo;
- **reparse das deflines** → `len(uids)` vs `fasta_record_count`;
- **recomputa** `compute_record_ids_sha256_for_order(uids, record_order)` vs `source_record_ids_sha256` → pega **reordenação** mesmo se o hash de arquivo fosse remendado (`RuntimeError "record-id set … changed"`);
- `selection_report_file_sha256` vs bytes.

Testado em `test_validation_detects_reordered_or_tampered_fasta` (reverte os registros do FASTA; 1º pega por "mismatch" de hash; depois, com `fasta_file_sha256` remendado, pega por "record-id set"). `[VERSIONADO]`

Na execução real: `source_record_ids_sha256 = 7bea84d9…` **==** o `sha256` de `retained_protein_uids.txt` do prefilter. `[EXEC-LOCAL]`

## 10.6 `artifact_type`

`derived_protein_fasta_snapshot` (`snapshot_format_version = "1.0"`). Teste `test_resolve_produces_honest_provenance_manifest` assere explicitamente `manifest["artifact_type"] != "ncbi_protein_fasta_snapshot"`. `[VERSIONADO]`

## 10.7 Adaptação no caminho do SWeeP

`[CÓDIGO]` (ver Seção 4.9):
- `sweep_genes_snapshot.SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES = ("ncbi_protein_fasta_snapshot", "derived_protein_fasta_snapshot")`, passado como `allowed_artifact_types` ao `load_fasta_snapshot_by_directory`.
- `ncbi_fasta_snapshot.load_*` ganhou `allowed_artifact_types` (default só o tipo NCBI → o caminho `07_sweep_*` existente **não muda**).
- Manifesto do SWeeP passa a registrar `source_fasta_artifact_type` e `source_fasta_record_ids_sha256`; `_source_fasta_snapshot_identity_matches` compara esse fingerprint **quando ambos os lados o têm** (snapshots NCBI legados não têm → tolerado, sem quebra).
- Teste `test_sweep_supported_types_accept_derived_fasta_but_not_arbitrary`: aceita os 2 tipos, rejeita um 3º arbitrário. `[VERSIONADO]`
- **Nenhum notebook de SWeeP foi alterado na Fase A** — a habilitação está pronta no módulo, mas o consumo (notebooks 11–16 do plano) é fase futura. `[CÓDIGO]` / `[INFERÊNCIA]`

---

# 11. Notebook 10 — auditoria célula por célula

`notebooks/10_dataset_audit.ipynb` — 14 células (índice 0–13): 1 markdown + 13 código. `nbformat 4.5`, kernel `python3`, `language_info.version "3.14.0"`. `[VERSIONADO]`

| # | Tipo | Objetivo | Módulos/chamadas | Lê | Produz | Rede? | Falhas possíveis |
|---|---|---|---|---|---|---|---|
| 0 | markdown | Contexto: "annotation‑enriched candidate set", 2 readings de recall, "orquestração apenas" | — | — | — | não | — |
| 1 | código | Imports + `importlib.reload` dos 6 módulos de snapshot; liga os `resolve_*` a nomes locais | `importlib`, `src.pago_pipeline.*_snapshot` | — | símbolos no namespace | não | `ImportError` se um módulo tiver erro de sintaxe **no disco**; `reload` **não** resolve imports transitivos já cacheados (ver 11.1) |
| 2 | código | Carrega `.env` (`find_dotenv`, `load_dotenv override=True`); resolve `PROJECT_ROOT` subindo até achar `src/`; lê `NCBI_EMAIL`, `NCBI_API_KEY` | `dotenv` | `.env` | `PROJECT_ROOT`, `NCBI_EMAIL`, `NCBI_API_KEY` | não | `FileNotFoundError` se não há `.env`; `ValueError` se `NCBI_EMAIL` ausente |
| 3 | código | Constantes de configuração (Seção 7.1) | — | — | `DATASET_NAME`, `SEARCH_QUERY`, todos os `*_MODE`, tamanhos | não | — |
| 4 | código | Deriva os 7 `*_SNAPSHOT_ROOT_DIRECTORY` sob `data/…{01-raw,02-intermediate,03-features}` + `QUERY_RECALL_REFERENCE_SET_CSV_PATH` (`tests/fixtures/…`); imprime | `pathlib` | — | os roots | não | — |
| 5 | código | **Preflight** | `resolve_ncbi_esearch_preflight_snapshot(mode=reuse_latest_or_create, …, max_uid_count=250000, sample_size=200, allow_exceeds_max_uid_count=False)` | `esearch_preflight__…/latest` (se existe) | snapshot de preflight | **sim** se precisa criar; **não** se reusa | `RuntimeError` se `exceeds_max_uid_count`; `RuntimeError` do `_run_with_retries` |
| 6 | código | **UID snapshot** | `resolve_ncbi_protein_uid_snapshot(mode, page_size=10000, retry=5, timeout=30, deadline=300, …)` | `protein_uid_snapshots__…/latest` | `protein_uids.txt` | sim/não | erros de rede após retries; deadline |
| 7 | código | **XML snapshot** | `resolve_ncbi_protein_xml_snapshot(mode, batch_size=100, retry=5, concurrency=4, reuse_http=False, batch_resume=True, purge=True)` | UID snapshot | `protein_records.xml` | sim/não | 5xx/timeout após retries; circuit breaker |
| 8 | código | **Metadata snapshot** | `resolve_ncbi_protein_metadata_snapshot(mode, source_xml_snapshot_root_directory=…)` | XML snapshot | `protein_metadata.csv`, `qc_report.json` | não (só transforma) | falha de QC → `save_*` não materializa (comportamento do módulo pré‑existente) |
| 9 | código | **Query reference recall** | `resolve_query_reference_recall_snapshot(mode, source_metadata_snapshot_root_directory=…, query_recall_reference_set_csv_path=…)`; imprime `matching_strategy_sha256`, os 2 readings por estrato, e linhas `SEQUENCE_SHA256`; `display(summary)`, `display(detail)` | fixture CSV + metadata snapshot | snapshot de recall | não | `RuntimeError` (colunas); `FileNotFoundError` (sem metadata) |
| 10 | código | **Technical prefilter** | `resolve_pago_technical_prefilter_snapshot(mode, source_metadata_snapshot_root_directory=…)`; `assert policy_sha256 == build_technical_prefilter_policy_sha256()`; imprime exclusões **ignorando `retain`** | metadata snapshot | snapshot de prefilter | não | `AssertionError` se a política do disco divergir do código; `FileNotFoundError` |
| 11 | código | **Derived FASTA** | `resolve_derived_protein_fasta_snapshot(mode, source_metadata_…, source_selection_snapshot_root_directory=TECHNICAL_PREFILTER_…, selection_artifact_type=ARTIFACT_TYPE, selection_uid_list_file_name=DEFAULT_RETAINED_PROTEIN_UIDS_FILE_NAME, record_selection_rule="pago_technical_prefilter.retained", record_selection_config_sha256=prefilter policy sha, dataset_kind="annotation_enriched_proteome", record_order="as_selected")` | metadata + prefilter snapshots | `protein_sequences.fasta`, `selection_report.json` | não | `RuntimeError` (uid ausente / seq vazia / artifact_type de seleção errado) |
| 12 | código | **Audit summary** — dict com query, `translated_query`, `ncbi_count`, `uids_retrieved`, `xml_records`, `metadata_rows`, prefilter in/retained/excluded + `counts_by_decision`, `matching_strategy_sha256`, `exact_recovered_count`, `equivalent_recovered_count`, 10 recalls por estrato/reading, `derived_fasta_records`; imprime; distribuição de comprimento de `retained_records["gbseq__length"]` + contagem `length_warning` | payloads das células 5–11 | `audit_summary` (em memória) | não | `KeyError` se um payload não tiver a chave esperada |
| 13 | código | Lista as variáveis expostas aos notebooks a jusante | — | texto | não | — |

## 11.1 O problema `importlib.reload` e por que "Restart Kernel" resolve

`[CÓDIGO]` (CELL 1 faz `importlib.reload` de 6 módulos `*_snapshot`) + `[INFERÊNCIA]`:

- CELL 1 recarrega **apenas** os 6 módulos `*_snapshot`. Não recarrega os módulos de **lógica** que eles importam (`query_reference_recall`, `pago_technical_prefilter`, `derived_protein_fasta`, `ncbi_esearch_preflight`), nem as dependências transitivas (`ncbi_snapshot`, `storage`, `metadata_to_fasta`).
- Quando o `query_reference_recall.py` foi **reescrito** no commit `0013c6d` (nova API: `ReferenceMatchMethod`, `build_matching_strategy_*`, campos do dataclass trocados), um kernel Jupyter que já tinha `query_reference_recall` carregado em memória continuava com a **versão antiga** do módulo de lógica, mesmo após `reload` do `_snapshot`. Sintoma típico: `ImportError`/`AttributeError` (`cannot import name 'build_matching_strategy_payload'`, ou `QueryReferenceRecallResult` sem `stratum_exact_recall`), ou o `_snapshot` recarregado ligando‑se ao objeto antigo.
- **Restart Kernel** limpa `sys.modules` inteiro: o próximo `import` lê **todos** os `.py` do disco de novo, na versão atual. Isso resolve.

**Distinção crucial `[INFERÊNCIA]`:** era um problema de **estado do kernel** (módulos stale em `sys.modules`), **não** um problema do código no disco. O código no disco estava correto — prova: `python -m unittest discover -s tests -q` (processo novo, `sys.modules` limpo) passa 203/203, e a reexecução do notebook **após restart** produziu o snapshot de recall `01-30-37Z` consistente. `reload` seletivo é uma conveniência frágil quando a assinatura/símbolos de um módulo de lógica mudam; a mitigação robusta é reiniciar o kernel após mudanças estruturais.

---

# 12. Testes

## 12.1 Inventário (arquivos adicionados pela Fase A)

`[VERSIONADO]` (contagem por `unittest -v … | grep " ... ok"`):

| Arquivo | Nº testes | Comportamentos protegidos | Fixtures / mocks |
|---|---|---|---|
| `test_ncbi_esearch_preflight.py` | 4 | `parse_uilist_text`; `summarize_sample_xml` (contagens de seq/uid); flag `exceeds_max_uid_count`; fluxo completo com History + amostra | XML GBSet embutido; `mock.patch` de `Entrez` e `_configured_ncbi_entrez_urlopen` |
| `test_ncbi_esearch_preflight_snapshot.py` | 2 | materializa report + reusa (não rechama `run_*`); **levanta** se `result_count > max_uid_count`, mas materializa o report; `allow_exceeds_max_uid_count=True` retorna payload | `mock.patch` de `run_ncbi_esearch_preflight` com `_fake_result`; `tempfile` |
| `test_pago_technical_prefilter.py` | 7 | retém apesar de texto/comprimento; exclui só as 4 causas técnicas (com prioridade); tolera `BJOUXZ*-`; partição total e disjunta; erro por coluna faltando; `policy_sha256` estável e sensível a parâmetro; `counts` somam ao input | DataFrames sintéticos |
| `test_pago_technical_prefilter_snapshot.py` | 4 | persiste 4 saídas + manifesto; reusa snapshot congelado; **indisponível** quando o metadata muda; `_validate_*` rejeita `artifact_type` e `hash mismatch` | metadata snapshot fake em `tempfile`; `mock.patch` de `save_*` |
| `test_derived_protein_fasta.py` | 7 | preserva ordem pedida; `sorted_by_uid`; `record_ids_sha256` == helper; uid ausente levanta salvo `drop_missing_uids`; duplicata na lista levanta; metadata com uid duplicado levanta; sequência vazia levanta | DataFrames sintéticos |
| `test_derived_protein_fasta_snapshot.py` | 4 | manifesto de proveniência honesto (`artifact_type != ncbi_…`; `derived_from_*`; `dataset_kind`); **detecta FASTA reordenado/adulterado** (2 caminhos); indisponível quando a seleção muda; reusa sem re‑salvar | metadata + selection snapshots fake; `mock.patch` de `save_*` |
| `test_fasta_snapshot_artifact_types.py` | 2 | default de `ncbi_fasta_snapshot` ainda rejeita `derived_protein_fasta_snapshot`; `SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES` aceita os 2 tipos e rejeita um 3º | dirs de snapshot FASTA fake |
| `test_query_reference_recall.py` | 17 | `matching_strategy_sha256` pinado; normalização; EXACT; SAME_BASE; accession diferente + sequência idêntica (equivalent‑only); sequência diferente não recuperada; **prioridade EXACT/SAME_BASE sobre hash**; múltiplos hits determinísticos; recall estratificado + 2 readings; estrato vazio `NOT_EVALUABLE` (não `0.0`); erro por coluna; fixture bem‑formada; tamanhos de estrato (8/2/4/7); 21 hashes `^[0-9a-f]{64}$` e distintos; linha RsAgo; RsAgo recuperado via `A4WYU7.1` | fixture CSV real + DataFrames sintéticos; `protein_sequence_sha256` para gerar colisões controladas |
| `test_query_reference_recall_snapshot.py` | 4 | recall contra a fixture real (LONG_A 2/8 no cenário; PIWI_RE `EVALUABLE` 0.0); recuperação por identidade de sequência sob accession alias (end‑to‑end no snapshot); reuso e **invalidação ao mudar o CSV** (append de 1 linha de 12 campos); **snapshot de metodologia antiga não é reusado** (patch de `matching_strategy_sha256` → `"0"*64`) | metadata snapshot fake com `gbseq__sequence`; CSV de referência custom em `tempfile`; `mock.patch` de `save_*` |
| **Total Fase A** | **51** | | |

## 12.2 Edge cases cobertos

- estrato de referência com **zero** linhas → `None`/`NOT_EVALUABLE`, nunca `0.0`;
- estrato presente mas 0 recuperados → `EVALUABLE` 0.0 (distinto do anterior);
- sequência recuperada sob **múltiplos** accessions com o mesmo hash → representante determinístico (menor lexicográfico) + `sequence_match_count`;
- shuffle da ordem dos registros recuperados → mesmo resultado;
- FASTA reordenado com hash de arquivo remendado → ainda pego pelo `source_record_ids_sha256`;
- prefiltro: `protein_uid` vazio **e** sequência vazia ao mesmo tempo → vence `drop_unprocessable_record` (prioridade);
- append de 1 linha ao CSV de referência → snapshot de recall fica indisponível;
- `result_count > max_uid_count` → **report materializado** para auditoria **e** `RuntimeError`.

## 12.3 Regressões protegidas explicitamente

- `matching_strategy_sha256 == 3460b048…` pinado (qualquer mudança de metodologia quebra o teste deliberadamente);
- `technical_prefilter_policy_sha256` estável e sensível a parâmetro;
- `DEFAULT_FASTA_SNAPSHOT_ARTIFACT_TYPES == ("ncbi_protein_fasta_snapshot",)` (o caminho NCBI não regrediu);
- `snapshot_format_version` do recall deve ser `"1.1"` (v1.0 rejeitada);
- `clade` nunca aceita `"PIWI_RE"`; toda linha `ago_family==PIWI_RE` tem `clade==UNRESOLVED`.

## 12.4 O que "203 tests / OK" significa — e o que **não** prova

`python -m unittest discover -s tests -q` → `Ran 203 tests in ~4.4s / OK`. `[CÓDIGO]` (rodado nesta auditoria).

**Prova:**
- os 203 testes (51 da Fase A + 152 pré‑existentes) passam num processo Python novo (sem estado de kernel), no HEAD atual, no ambiente `.venv` (Python 3.14);
- a lógica pura da Fase A se comporta como especificado nos cenários testados;
- os `_snapshot` criam/validam/reusam corretamente com metadata/seleção **sintéticos** em `tempfile`;
- os contratos de invalidação (mudança de política, de CSV, de metodologia, de metadata) funcionam.

**NÃO prova:**
- que a **execução real contra o NCBI** produziu números corretos — os testes **mockam** toda a rede (`Entrez`, `_configured_ncbi_entrez_urlopen`, `run_ncbi_esearch_preflight`, `save_*`); os 52.473 vêm de manifestos **não versionados**;
- que os snapshots reais em `data/…` são íntegros — `unittest` não os toca (isso é `verify_raw_data.py`, e mesmo esse só cobre `data/01-raw`);
- que a curadoria biológica das 21 referências está **cientificamente correta** — os testes verificam **forma** (regex, contagens, unicidade de hash), não se `WP_011174533.1` é de fato um TtAgo long‑A;
- cobertura de código (não há gate de coverage);
- que o notebook 10 roda ponta a ponta num kernel Jupyter (nenhum teste executa `.ipynb`);
- reprodutibilidade do SWeeP/PCA a jusante (fora do escopo Fase A);
- ausência de regressão em cenários não testados.

---

# 13. Execução científica real

## 13.1 Natureza da evidência

**Nenhum** artefato numérico da execução completa está versionado. `[VERSIONADO]` (`git ls-files data/` só lista os snapshots da query **antiga** `q_891f443d754c`, de abril/2026; `.gitignore` ignora `data/02-intermediate/**` e `data/03-features/**` inteiros, e os dirs `__annotation_enriched_candidate_set` de `data/01-raw` estão **untracked**).

Os números abaixo vêm de manifestos e CSVs em `data/…` no working tree do usuário — **`[EXEC-LOCAL]`**. A auditoria recomputou hashes de arquivo na hora e confirmou consistência interna, mas isso **não é** verificação versionada.

## 13.2 Números da execução real (confirmados nos artefatos em disco)

| Métrica | Valor | Fonte (`[EXEC-LOCAL]`) |
|---|---|---|
| NCBI `Count` (preflight) | **52 473** | `esearch_preflight__…/latest/manifest.json:result_count` |
| `exceeds_max_uid_count` | `false` (teto 250 000) | idem |
| amostra do preflight | 200/200 com sequência, 200/200 uid extraível, `sample_fetch_error: null` | idem |
| `translated_query` | `(PIWI[All Fields] OR Argonaute[All Fields]) AND (("Bacteria"[Organism] OR "Bacteria Latreille et al. 1825"[Organism]) OR "Archaea"[Organism])` | idem |
| UIDs (raw = normalized = ncbi_reported) | **52 473** | `protein_uid_snapshots__…/latest/manifest.json` |
| `protein_uids_sha256` | `949bcaab…b253` | idem |
| XML `consolidated_record_count` | **52 473** | `protein_xml_snapshots__…/latest/manifest.json` |
| XML `batch_count` | 525 (× 100) | idem |
| `xml_file_sha256` | `77739043…6169` | idem |
| metadata `row_count` | **52 473** | `protein_metadata_csv__…/latest/manifest.json:` (`row_count` via prefilter `source_metadata_row_count`) + `qc_report.json` |
| metadata `column_count` | 148 | idem |
| `csv_file_sha256` | `b39096f3…1510` | idem |
| metadata QC | 5/5 checagens `true` | `qc_report.json` |
| prefilter `input_record_count` | **52 473** | `pago_technical_prefilter/latest/manifest.json` |
| prefilter `retained_record_count` | **52 473** | idem |
| prefilter `excluded_record_count` | **0** | idem (`counts_by_decision`: retain 52473, drops 0) |
| `technical_prefilter_policy_sha256` | `032fbc72…a523` | idem |
| `retained_protein_uids.txt` sha256 | `7bea84d9…8ef0` | idem + recomputado nesta auditoria ✔ |
| distribuição de comprimento (retidos) | min 26, p25 233, mediana 385, média 476,58, p75 703, max 5 433 | recomputado de `retained_records.csv` nesta auditoria |
| `length_warning == True` | **8 299** (de 52 473) | recomputado idem |
| recall `reference_count` | **21** | `query_reference_recall__…/latest/manifest.json` |
| recall `exact_recovered_count` | **20** | idem |
| recall `equivalent_recovered_count` | **21** | idem |
| `stratum_exact_recall` | overall 0.9524, LONG_A 1.0, LONG_B 0.5, SHORT 1.0, PIWI_RE 1.0 | idem |
| `stratum_equivalent_recall` | overall 1.0, LONG_A 1.0, LONG_B 1.0, SHORT 1.0, PIWI_RE 1.0 | idem |
| `stratum_recall_status` | 5× `EVALUABLE` | idem |
| `matching_strategy_sha256` | `3460b048…1a08` | idem |
| `query_recall_reference_set_csv_sha256` | `1da47d3f…a70b` (== fixture no HEAD ✔) | idem + recomputado ✔ |
| derived FASTA `fasta_record_count` | **52 473** | `derived_protein_fasta__…/latest/manifest.json` |
| `fasta_file_sha256` | `f7bac76b…ce43` | idem + recomputado ✔ |
| `source_record_ids_sha256` | `7bea84d9…8ef0` (== `retained_protein_uids.txt` ✔) | idem |

**Confirmação dos números do enunciado:** os 6 valores citados (`NCBI Count / UIDs / XML records / metadata rows / prefilter retained / derived FASTA` = **52 473** cada) **conferem** nos artefatos em disco. `[EXEC-LOCAL]` Não são versionados.

## 13.3 Erros transitórios no XML fetch (com evidência)

`[EXEC-LOCAL]` (telemetria `xml_batches` do manifesto XML):

| Categoria | Contagem | Interpretação |
|---|---|---|
| `http_5xx` | **6** | respostas 5xx do NCBI (compatível com HTTP 502 Bad Gateway) |
| `truncated_response` | **3** | corpo lido de forma incompleta (compatível com `IncompleteRead`) |
| `retry_count` (total do estágio) | **9** | = 6 + 3, todas seguidas de nova tentativa |
| `http_429` / `rate_limit` / `timeout` / `deadline` / `response_validation` / `other` | 0 | — |

Todas recuperadas: `request_count = 525` para `batch_count = 525`, `consolidated_record_count = 52 473`, `reused_batch_count = 0`. Não há terminologia literal "502" ou "IncompleteRead" gravada no manifesto — os contadores são `http_5xx` e `truncated_response`. `wall_seconds_total ≈ 1 170 s`. `[EXEC-LOCAL]`

## 13.4 O que é execução real vs comportamento testado/mocado

| Item | Real (`[EXEC-LOCAL]`) | Testado/mocado (`[VERSIONADO]`/`[CÓDIGO]`) |
|---|---|---|
| Contagem 52 473 e cadeia de hashes | sim (manifestos em disco) | não (testes usam ≤5 registros sintéticos) |
| Retry recuperando 6 5xx + 3 truncados | sim | testes só cobrem o retry do **preflight**, com mock |
| `recall exact 20 / equivalent 21`, RsAgo via `A4WYU7.1` | sim | a **lógica** do matcher é testada; o **valor** 21/21 vem da execução |
| Prefiltro 0 exclusões | sim | testes forçam exclusões com dados sintéticos |
| Notebook 10 rodando ponta a ponta | sim (kernel PyCharm do usuário) | nenhum teste executa o notebook |
| Integridade dos snapshots reais | recomputada nesta auditoria (não versionada) | `unittest` não toca `data/…` |

---

# 14. Verificações finais

## 14.1 `python -m unittest discover -s tests -q`

`[CÓDIGO]` (rodado nesta auditoria, `.venv` Python 3.14):

```
Ran 203 tests in 4.378s
OK
```

Processo novo, `sys.modules` limpo, sem rede (tudo mocado). Cobre os 51 testes da Fase A + 152 pré‑existentes.

## 14.2 `python scripts/verify_raw_data.py`

`[CÓDIGO]` (rodado nesta auditoria):

```
Verified 16 raw data file hashes from 13 manifests.
```

**O que o script faz:** `rglob("manifest.json")` sob `data/01-raw`; para cada manifesto, se tem `protein_uids_sha256` recomputa `sha256_of_lines` do `.txt` e compara; se tem `xml_file_sha256` recomputa `sha256_of_file` do `.xml` e compara. `[CÓDIGO]`

**Ponto de atenção:** `rglob` **não respeita o Git** — os "13 manifestos" incluem os dirs `__annotation_enriched_candidate_set` de `data/01-raw` que **não estão versionados**. Ou seja, `verify_raw_data.py` verifica o **working tree**, misturando snapshots versionados (query antiga) e não versionados (Fase A). Não cobre `data/02-intermediate` nem `data/03-features` (recall, prefilter, derived FASTA ficam de fora).

## 14.3 Funcionamento de snapshot reuse

`[CÓDIGO]` (observado no output do `unittest` e nos prints dos `resolve_*` durante a suíte):

```
Latest technical prefilter snapshot is available. Reusing frozen snapshot.
Found a compatible immutable PCA/KMeans snapshot. Reusing frozen snapshot.
Latest query reference recall snapshot is available. Reusing frozen snapshot.
Found a compatible immutable SWeeP Genes snapshot. Reusing frozen snapshot.
```

Os testes `test_reuse_latest_or_create_reuses_frozen_snapshot` / `…reuses_without_resaving` / `test_reuse_and_invalidation_on_reference_set_change` confirmam: com snapshot válido, `save_*` **não** é chamado; ao mudar um insumo (metadata, CSV de referência, política, metodologia), `latest_*_is_available` retorna `False` e `save_*` roda.

## 14.4 HEAD final

`[CÓDIGO]`:

```
branch : (feat)-siepe-ready-project
HEAD   : 0013c6d24d90070ee297ba6af7d12d125ebd1732
origin : 0013c6d…  (sincronizado; git status: working tree com apenas untracked data/ + output/ + tmp/)
```

`git status` no início da sessão mostrava só `?? output/` e `?? tmp/`; após o trabalho de `0013c6d`, os untracked são os 3 dirs `data/01-raw/*__annotation_enriched_candidate_set/` (produzidos pela execução) + `pdfs/output/`. Nenhum arquivo rastreado modificado. `[CÓDIGO]` / `[EXEC-LOCAL]`

## 14.5 Evidência versionada × execução local não versionada

| Categoria | Conteúdo |
|---|---|
| **Versionado** (`[VERSIONADO]`) | Todo o código (`src/pago_pipeline/*.py`), todos os 9 arquivos de teste, `tests/fixtures/query_recall_reference_set.csv` (+ notas), `notebooks/10_dataset_audit.ipynb`, `README.md`, `.gitignore`, `.gitattributes`. Estado congelado em `0013c6d`, empurrado para `origin`. |
| **Execução local NÃO versionada** (`[EXEC-LOCAL]`) | Tudo em `data/01-raw/*__annotation_enriched_candidate_set/`, `data/02-intermediate/{protein_metadata_csv,query_reference_recall}__…/`, `data/02-intermediate/derived_protein_fasta__annotation_enriched_proteome/`, `data/03-features/pago_technical_prefilter/`. Inclui **todos** os números 52 473, os hashes da cadeia de proveniência, os 2 snapshots de recall, a telemetria de retry. Existe apenas no disco do usuário; some se `data/` for limpo (é regenerável reexecutando o notebook 10). |

Consequência prática: **um clone limpo do repositório no HEAD atual não reproduz os números da Fase A sem reexecutar o notebook 10 contra o NCBI** (o que exige `.env` com `NCBI_EMAIL` e conectividade). Os testes passam num clone limpo; os artefatos científicos, não.

---

# 15. Decisões científicas e metodológicas

Cada decisão: **problema → alternativas → escolha → motivo → consequência → limitação residual.** Selos indicam onde a decisão está materializada.

### D1 — Chamar o dataset de `annotation_enriched_candidate_set`
- **Problema:** como nomear o produto da 2ª obtenção NCBI sem sugerir que é o universo de pAgos.
- **Alternativas:** "pago_universe", "bacterial_archaeal_argonautes", "candidate_set" sem qualificador.
- **Escolha:** `annotation_enriched_candidate_set` (dataset) / `annotation_enriched_proteome` (dataset_kind do FASTA). `[VERSIONADO]` (CELL 3, README, markdown do notebook).
- **Motivo:** o nome carrega o viés — a query `[All Fields]` só acha o que já foi anotado com "PIWI"/"Argonaute".
- **Consequência:** todo artefato a jusante herda o nome nos paths e manifestos; qualquer leitor vê o qualificador.
- **Limitação residual:** nome longo; `pago_technical_prefilter` e `sweep`/`pca` compartilham roots não‑sufixados (`data/03-features/pago_technical_prefilter`), então a distinção do dataset se perde nesses estágios.

### D2 — Não tratar o dataset como universo de pAgos
- **Problema:** a tentação de rodar estatística populacional sobre 52 473 "pAgos".
- **Alternativas:** assumir cobertura ~completa; assumir cobertura desconhecida sem medir.
- **Escolha:** medir explicitamente com `query_reference_recall` estratificado + documentar que é fase futura a rota sequence‑based. `[VERSIONADO]` / `[CÓDIGO]`
- **Motivo:** integridade epistemológica — não confundir "recuperável por texto" com "existente".
- **Consequência:** o recall é relatado como cobertura de painel, não como sensibilidade.
- **Limitação residual:** o painel de 21 é pequeno e enriquecido; o recall 21/21 **não** extrapola (ver Seção 17).

### D3 — Não filtrar por texto de anotação
- **Problema:** seria fácil dropar "SAM‑methyltransferase" etc. já no prefiltro.
- **Alternativas:** lista negra de termos; classificador de anotação.
- **Escolha:** o prefiltro **nunca** lê texto de anotação. `[CÓDIGO]` (`notes` da política; `test_retains_records_regardless_of_annotation_text_or_length`).
- **Motivo:** anotação textual é ruidosa e circular; a exclusão de "óbvios não‑pAgo" deve ser feita por **domínio/HMM** (fase futura), com evidência estrutural.
- **Consequência:** o proteoma retido contém falsos‑positivos textuais (transposases, metiltransferases co‑anotadas), a serem removidos pela varredura de domínio.
- **Limitação residual:** o FASTA derivado é maior e mais ruidoso do que o conjunto final de pAgos; o SWeeP/PCA sobre ele mistura pAgos e não‑pAgos.

### D4 — Não excluir por comprimento
- **Problema:** pAgos "long" têm ~700–900 aa, "short" ~400–500; 26 aa ou 5 433 aa "não parecem" pAgo.
- **Alternativas:** cortar fora de [200, 2000]; cortar por percentil.
- **Escolha:** comprimento fora da banda só liga `length_warning = True`; **nunca** exclui. `[CÓDIGO]`
- **Motivo:** "a varredura de domínio, não o prefiltro, decide se um comprimento é compatível com uma pAgo"; fragmentos e fusões podem ser biologicamente reais.
- **Consequência:** 8 299/52 473 registros retidos têm `length_warning`. `[EXEC-LOCAL]`
- **Limitação residual:** ruído adicional para o embedding; `length_warning` precisa ser propagado e considerado a jusante (está no `retained_records.csv`).

### D5 — Separar `ago_family` (`PIWI_RE`) de `pago_clade`
- **Problema:** PIWI‑RE às vezes é tratado como um "clado" de pAgo.
- **Alternativas:** `clade = "PIWI_RE"`; ignorar PIWI‑RE.
- **Escolha:** `PIWI_RE` é valor de `ago_family`; toda linha PIWI‑RE tem `clade = UNRESOLVED`; o estrato PIWI‑RE seleciona por `ago_family`. `[VERSIONADO]` (`RECALL_STRATA`, `98e66c6`, teste).
- **Motivo:** PIWI‑RE é uma **família** divergente (Burroughs et al. 2013), não um clado dentro da árvore MID‑PIWI de pAgos long/short.
- **Consequência:** o recall reporta um estrato `piwi_re_reference_recall` separado; o teste proíbe `"PIWI_RE"` na coluna `clade`.
- **Limitação residual:** a ontologia completa (`ago_family` ∈ {PAGO, PIWI_RE, UNRESOLVED}, `pago_clade` ∈ {LONG_A, LONG_B, SHORT, UNRESOLVED}) só existe no **plano**; a Fase A só materializa o suficiente para o recall.

### D6 — Painel de referência curado (21), com `reference_label_evidence`
- **Problema:** medir recall exige um "gabarito" de pAgos conhecidas; quão grande e quão confiável?
- **Alternativas:** usar só as ~10 experimentais clássicas; usar centenas de hits de BLAST.
- **Escolha:** 21 entradas, cada uma com `reference_label_source` (citação) e `reference_label_evidence` ∈ {EXPERIMENTAL, LITERATURE_PHYLOGENETIC, CURATED_COMPUTATIONAL, DATABASE_ANNOTATION}; PIWI‑RE expandido de 0 → 1 experimental (PsPIWI‑RE) → +6 computacionais. `[VERSIONADO]` (fixture, notas de curadoria, commits `d21f25a`, `0fd283a`, `21793bc`).
- **Motivo:** transparência sobre a força de cada rótulo; evitar tanto o painel minúsculo (estratos vazios → `NOT_EVALUABLE`) quanto o painel inflado com rótulos fracos.
- **Consequência:** todos os 5 estratos ficam `EVALUABLE`; os 6 PIWI‑RE `CURATED_COMPUTATIONAL` carregam a ressalva de circularidade para HMM futuro.
- **Limitação residual:** 21 é pequeno; `LONG_B` tem só 2; 6/21 são `provisional`.

### D7 — Equivalência por identidade de sequência (matcher `SEQUENCE_SHA256`)
- **Problema:** `ABP72561.1` (RsAgo) não está entre os accessions recuperados, embora sua sequência esteja (sob `A4WYU7.1`). Matching só por accession contava como miss.
- **Alternativas:** trocar o accession na fixture para `A4WYU7.1` (rejeitado pelo usuário — `ABP72561.1` é o accession correto do RsAgo experimental); resolver IPG online em tempo de execução; aceitar o miss.
- **Escolha:** 3ª camada de matching por `sha256(sequência normalizada)`, offline e determinística; hashes das 21 referências derivados **offline** da `protein_metadata.csv` local (RsAgo via alias byte‑idêntico `A4WYU7.1`). `[CÓDIGO]` / `[VERSIONADO]` (`0013c6d`, notas de curadoria).
- **Motivo:** medir recuperação **biológica** (a proteína) sem depender de rede em execução normal e sem falsificar a fixture.
- **Consequência:** `retrieval_equivalent_recall` = 21/21; `matching_strategy_sha256` no manifesto; snapshot `1.0` invalidado.
- **Limitação residual:** o hash é sensível a **qualquer** diferença de 1 resíduo (isoformas, variantes de comprimento não casam); depende da fixture ter o hash correto (derivado de uma execução específica).

### D8 — Preservar dois readings de recall
- **Problema:** um único número esconde se a recuperação foi "mesmo accession" ou "mesma proteína, outro accession".
- **Alternativas:** só `retrieval_equivalent_recall` (some o miss de accession); só `exact_accession_recall` (esconde a recuperação real).
- **Escolha:** manter os dois (`stratum_exact_recall` **e** `stratum_equivalent_recall`) no manifesto e no notebook. `[CÓDIGO]`
- **Motivo:** honestidade — o `exact` 20/21 (LONG_B 0.5) é reportado, não escondido; o `equivalent` 21/21 dá a leitura biológica.
- **Consequência:** o notebook CELL 9 imprime uma tabela de 2 colunas; `audit_summary` tem 10 chaves de recall.
- **Limitação residual:** exige explicar a diferença a cada leitura do relatório.

### D9 — Snapshots imutáveis + `reuse_latest_or_create`
- **Problema:** reexecutar o notebook não deve rebaixar reprodutibilidade nem custar re‑download.
- **Alternativas:** sempre recomputar; cache sem verificação de identidade.
- **Escolha:** diretório imutável (`mkdir exist_ok=False`) + `latest/` substituível + `latest_*_is_available` verificando hashes de todos os insumos + rollback por `rmtree` em falha parcial. `[CÓDIGO]`
- **Motivo:** cada execução deixa um registro congelado e auditável; o reuso só acontece se **tudo** que define o resultado bate.
- **Consequência:** a 2ª execução do notebook reusou preflight/uid/xml/metadata/prefilter/derived e só recomputou o recall (metodologia mudou). `[EXEC-LOCAL]`
- **Limitação residual:** os snapshots imutáveis são grandes e regeneráveis; a política de quais versionar é manual (`.gitignore` re‑inclui por nome). Nenhum snapshot da Fase A foi versionado ainda.

### D10 — Notebook = orquestração; lógica em `src/pago_pipeline/`
- **Problema:** manter os notebooks testáveis e o comportamento reproduzível.
- **Escolha:** o notebook 10 só chama `resolve_*` e imprime; toda a lógica (com `unittest`) vive nos módulos. `[VERSIONADO]` (markdown do notebook, README).
- **Consequência:** 51 testes cobrem a lógica; o notebook não tem teste próprio.
- **Limitação residual:** bugs de **fiação** do notebook (argumento errado, ordem de célula) não são pegos por `unittest` — só pela execução real (foi o caso do `importlib.reload`, Seção 11.1).

---

# 16. Problemas encontrados e lições

Formato: **sintoma → causa → investigação → solução → teste/proteção → lição.**

### P1 — PIWI‑RE inicialmente sem referências (estrato vazio)
- **Sintoma:** o estrato `PIWI_RE` do recall ficava vazio → `NOT_EVALUABLE`.
- **Causa:** a curadoria inicial só tinha pAgos long/short; PIWI‑RE não tinha nenhuma entrada com accession NCBI rastreável e caracterização.
- **Investigação:** busca por PIWI‑RE experimentalmente caracterizada com accession NCBI; o usuário encontrou PsPIWI‑RE `WP_014597637.1` (*Pseudomonas stutzeri* DSM 4166, Huang et al. 2022, mutantes D525A/D610A/R639A/E718A → implica ≥718 aa; EFetch confirmou 783 aa).
- **Solução:** commit `0fd283a` adiciona `WP_014597637.1` como `PIWI_RE`/`EXPERIMENTAL`/`verified`, `clade=UNRESOLVED`; o estrato passa a ser selecionado por `ago_family`.
- **Teste/proteção:** `test_committed_reference_set_is_well_formed` exige ≥1 linha PIWI‑RE; `test_stratum_sizes_match_curation` fixa `PIWI_RE=7`.
- **Lição:** um estrato vazio deve ser `NOT_EVALUABLE` (não `0.0`) — mas a meta é preenchê‑lo com a melhor evidência disponível, não deixá‑lo vazio "de propósito".

### P2 — Curadoria do conjunto PIWI‑RE (base accession vs accession.version)
- **Sintoma:** o usuário forneceu 5 accessions PIWI‑RE `CURATED_COMPUTATIONAL` como **base** (sem `.N`).
- **Causa:** `.1` não pode ser presumido; o RefSeq pode estar em `.2`, `.3`…
- **Investigação:** para cada base, EFetch `rettype=acc` → accession.version corrente; ESummary confirmou "pPIWI_RE module domain-containing protein" / "RNaseH domain-containing protein", 830–1052 aa.
- **Solução:** commit `21793bc` adiciona as 6 linhas com `.1` **verificadas** (não presumidas); notas de curadoria registram o método.
- **Teste/proteção:** `test_committed_reference_set_is_well_formed` exige `accession` casando `\.\d+$`.
- **Lição:** identificadores de referência sempre com versão explícita e verificada; documentar a resolução.

### P3 — Identificadores históricos (GI) da literatura PIWI‑RE
- **Sintoma:** Burroughs et al. 2013 cita números GI, que foram descontinuados como identificador primário.
- **Causa:** GI → accession não é 1:1 estável ao longo do tempo (registros suprimidos, dead records).
- **Investigação:** GI 269125748 → YP_003299118.1 → **WP_012851864.1** (cadeia de substituição NCBI viva) → **incluído**. GI 228927677 → ZP_04090728.1 (dead, sem substituição) → **excluído**. GI 119855142 → YP_935747.1 → WP_011767947.1 (suprimido, 477 aa) → **excluído**. GI 158336201 → é a REase associada, não a PIWI‑RE → **excluído**.
- **Solução:** só o caso com cadeia inequívoca e registro vivo entrou; notas de curadoria têm a tabela.
- **Teste/proteção:** documental (notas de curadoria); nenhum teste automatizado sobre resolução de GI.
- **Lição:** identificadores históricos precisam de uma cadeia de substituição rastreável e um registro **vivo**; na dúvida, excluir.

### P4 — RsAgo `ABP72561.1` como falso miss
- **Sintoma:** primeira execução do recall: `overall = 0.952`, `LONG_B = 0.5`; `ABP72561.1` reportado como único miss.
- **Causa:** o matcher só comparava accession.version; `ABP72561.1` (GenBank) não está entre os 52 473, mas sua sequência está sob `A4WYU7.1` (Swiss‑Prot, mesma IPG).
- **Investigação (mecânica, por identidade de sequência — não por nome/organismo):** IPG de `ABP72561.1` = {`A4WYU7.1`, `XLG71013.1`}, sem `WP_`. `A4WYU7.1` está nos 52 473 (`protein_uid 2500461169`), sequência byte‑idêntica (777 aa, `sha256 cbdb6bb6…`). RsAgo também presente como cadeias PDB 5AWH/6D8A/6D8F/6D8P/6D92. **Conclusão A** (existe accession alternativo recuperado com sequência idêntica).
- **Solução:** commit `0013c6d` — matcher `SEQUENCE_SHA256`; `sequence_sha256`/`sequence_length` nas 21 referências (derivados offline); dois readings; `ABP72561.1` mantido na fixture.
- **Teste/proteção:** `test_rsago_is_recovered_by_sequence_identity_under_alias_accession`, `test_rsago_row_carries_the_investigated_sequence_hash`, `test_reference_sequence_hashes_are_well_formed`.
- **Lição:** "miss" num benchmark de accession pode ser artefato do matcher; investigar por **identidade de sequência / IPG**, não por organismo/nome; e reportar **os dois** readings.

### P5 — Notação científica em `protein_uid` no `detail.csv`
- **Sintoma:** `protein_uid` (inteiro grande) sendo renderizado como `2.5e+09` ao carregar o `detail.csv`.
- **Causa:** `pd.read_csv` inferia `float64` para a coluna de identificador numérico.
- **Investigação:** direta.
- **Solução:** commit `34c0ca6` — o loader `load_query_reference_recall_snapshot_by_directory` passa `dtype={"accession":"string", "sequence_sha256":"string", "match_method":"string", "matched_accession":"string", "matched_protein_uid":"string"}`.
- **Teste/proteção:** implícita (o teste de snapshot lê `detail` e usa `.set_index("accession")` / compara strings).
- **Lição:** identificadores são strings, nunca números — fixar `dtype` na leitura.

### P6 — `retain` aparecendo sob "exclusions" no notebook
- **Sintoma:** CELL 10 listava `retain: 52473` na seção "technical exclusions by reason".
- **Causa:** o loop iterava `counts_by_decision.items()` sem filtrar a decisão `retain`.
- **Solução:** commit `34c0ca6` — `for decision, count in …: if decision == 'retain': continue`.
- **Teste/proteção:** nenhuma automatizada (é apresentação de notebook).
- **Lição:** `retain` não é exclusão; a apresentação deve distinguir a decisão‑padrão das causas de drop.

### P7 — `ImportError` por módulo/kernel stale no Jupyter
- **Sintoma:** ao reexecutar o notebook após `0013c6d`, `ImportError`/`AttributeError` em `query_reference_recall` (símbolos novos ausentes, dataclass com campos antigos).
- **Causa:** CELL 1 só faz `importlib.reload` dos 6 `*_snapshot`, não dos módulos de lógica nem das dependências transitivas; o kernel mantinha a versão pré‑`0013c6d` de `query_reference_recall` em `sys.modules`.
- **Investigação:** `python -m unittest` (processo novo) passava 203/203 → o código no disco estava correto; logo era estado do kernel.
- **Solução:** **Restart Kernel** (limpa `sys.modules`); a reexecução produziu o snapshot de recall `01-30-37Z` consistente.
- **Teste/proteção:** a suíte roda em processo novo por construção; nenhuma proteção contra staleness de kernel (é limitação intrínseca do `reload` seletivo).
- **Lição:** problema de **estado do kernel** ≠ problema de **código no disco**. Após mudança estrutural num módulo de lógica (assinatura, símbolos, dataclass), reiniciar o kernel; `importlib.reload` seletivo só é seguro para mudanças internas que não alteram a interface.

### P8 — Falhas transitórias do NCBI no XML fetch
- **Sintoma:** durante a corrida XML real, 6 respostas 5xx + 3 respostas truncadas.
- **Causa:** instabilidade transitória do serviço E‑utilities do NCBI sob carga (525 lotes, concorrência 4).
- **Investigação:** telemetria do manifesto XML: `failure_counts.http_5xx = 6`, `truncated_response = 3`, `retry_count = 9`.
- **Solução:** o mecanismo pré‑existente de retry/backoff (`0.1 → 30s ×2`, 5 tentativas) + circuit breaker recuperou **todas**; `consolidated_record_count = 52 473`, QC do metadata `row_count_matches_source_xml = true`.
- **Teste/proteção:** o retry do **preflight** é testado com mock; o retry de UID/XML é pré‑Fase‑A (fora do escopo de teste desta fase).
- **Lição:** o E‑utilities falha de forma transitória em escala; o design de retry/telemetria pré‑existente é adequado e deixa rastro auditável dos incidentes.

### P9 — `git push` bloqueado pelo classificador de permissões (Fase 0)
- **Sintoma:** o `git push` inicial da branch foi recusado pelo ambiente do agente.
- **Causa:** política de permissões do projeto no Claude Code.
- **Solução:** o usuário ajustou as permissões do projeto; o `push` seguinte funcionou. Nenhuma tentativa de contornar o bloqueio.
- **Lição:** operações de escrita remota dependem de autorização explícita do usuário; não contornar.

### P10 — Erro de sintaxe em `derived_protein_fasta_snapshot.py` (`_source_identity_matches`)
- **Sintoma:** `SyntaxError` — `if not _source_identity_matches(...)` fechava com `)` em vez de `):`, seguido de `return True` na lógica errada.
- **Causa:** edição manual malfeita.
- **Solução:** corrigido para `):` + `return False` + `return True`.
- **Teste/proteção:** `py_compile` + a suíte (`test_derived_protein_fasta_snapshot.py`).
- **Lição:** compilar cada módulo após edição; a suíte pega o resto.

### P11 — Fixture CSV ignorada pelo `.gitignore`
- **Sintoma:** `git add tests/fixtures/query_recall_reference_set.csv` não adicionava (regra `*.csv`).
- **Causa:** `.gitignore` tem `*.csv` global.
- **Solução:** commit `599a111` adiciona `!tests/fixtures/` + `!tests/fixtures/**` (e `!src/pago_pipeline/resources/**`).
- **Teste/proteção:** a existência dos testes que leem a fixture falha em CI se ela não estiver versionada.
- **Lição:** dados curados versionados precisam de negação explícita quando há regra global de tipo.

---

# 17. O que a Fase A permite concluir

## 17.1 Podemos afirmar

`[EXEC-LOCAL]` + `[CÓDIGO]`:

1. A query `(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])` retornou **52 473** UIDs de proteína no NCBI em 2026‑08‑30/31; o NCBI a traduziu incluindo `"Bacteria Latreille et al. 1825"` e `"Archaea"` como `[Organism]`.
2. Os 52 473 registros foram baixados (525 lotes), achatados num CSV de 148 colunas, e passaram 5/5 checagens de QC (sem `protein_uid` duplicado/vazio; contagem casa com o XML).
3. O prefiltro **técnico** reteve **todos** os 52 473 (0 exclusões técnicas) — consistente com o design e com a limpeza da fonte NCBI.
4. 8 299 dos retidos têm comprimento fora de [200, 2000] (`length_warning`), e mesmo assim foram retidos.
5. O FASTA derivado tem exatamente **52 473** registros, na ordem da seleção, com cadeia de proveniência (hash) completa e verificável até o CSV de metadata.
6. **Neste painel de 21 pAgos/PIWI‑RE conhecidas**, a query recuperou:
   - por **sequência** (`retrieval_equivalent_recall`): **21/21** (100%) — todos os estratos 1.0;
   - por **accession exato** (`exact_accession_recall`): **20/21** (95,2%); o único não‑exato é RsAgo `ABP72561.1`, recuperado sob o alias byte‑idêntico `A4WYU7.1`.
7. A infraestrutura de retry/telemetria pré‑existente absorveu 6 respostas 5xx + 3 truncadas do NCBI sem perda de registro.
8. Todo o código e os testes da Fase A (51 testes novos, 203 no total) passam num processo limpo no HEAD `0013c6d`, empurrado para `origin`.

## 17.2 NÃO podemos afirmar

1. **Que os 52 473 registros são pAgos.** São proteínas de bactéria/arqueia *co‑anotadas* com "PIWI" ou "Argonaute" em qualquer campo indexado — incluem transposases, metiltransferases, REases e fragmentos co‑mencionados. A separação pAgo/não‑pAgo é fase futura (HMM de domínio).
2. **Que `recall = 21/21` é a sensibilidade da query sobre o universo de pAgos.** O painel tem 21 proteínas, fortemente enriquecidas em entradas experimentais e bem anotadas — o melhor caso para uma busca textual. O número **não** estima quantas pAgos divergentes (mal anotadas, ou anotadas com outra terminologia) a query perde.
3. **Que a query recupera pAgos que não usam a terminologia "PIWI"/"Argonaute" em nenhum campo.** Por construção, `[All Fields]` não as acha; uma rota por homologia de sequência (HMM/PSI‑BLAST) é necessária e está fora da Fase A.
4. **Que os clados `LONG_A`/`LONG_B`/`SHORT` das referências estão corretos por análise deste projeto.** Vêm de literatura filogenética; 6/21 são `provisional`; nenhum placement foi feito aqui.
5. **Que os 6 PIWI‑RE `CURATED_COMPUTATIONAL` servem para validar um detector PIWI‑RE.** Foram selecionados *usando* modelos de perfil PIWI‑RE — usá‑los para pontuar um HMM PIWI‑RE seria circular (vedado nas notas de curadoria).
6. **Que a execução real é reprodutível a partir do repositório sozinho.** Os artefatos científicos não são versionados; um clone limpo precisa reexecutar o notebook 10 contra o NCBI.
7. **Nada sobre atividade catalítica, sistema (SPARTA/SPARSA/SPARDA), arquitetura de domínio ou estrutura** — a Fase A não toca nessas dimensões.

## 17.3 Continua para fases futuras (do plano, **não** implementado)

- Camada de referência de HMM (Pfam bundle, APAZ, clade seeds, árvore MID‑PIWI) com validação BUILD/CALIBRATION/FINAL_HOLDOUT — **Fase B**.
- Varredura de domínio (`pago_pfam_domain_scan`, `pago_apaz_scan`, `pago_clade_hmm_triage`) — **Fase C**.
- Detecção/família/anotação de proteína, placement filogenético de clado (EPA‑ng) — **Fase C**.
- Reconstrução de locus, vizinhança genômica (IPG), `system_class` — **Fase D**.
- Integração das 4 tabelas, redundância (MMseqs2), SWeeP→PCA→Clustering→plot 3D sobre (a) proteoma pré‑filtrado e (b) pAgos classificadas — **Fase E**.
- Pré‑filtro do "óbvio não‑pAgo" (ex. SAM‑metiltransferase) — parte da Fase C (por domínio), **não** do prefiltro técnico da Fase A.

---

# 18. Limitações e dívida técnica

> Estas são observações do estado atual — **não** requisitos retroativos da Fase A. Vários itens são deliberadamente adiados para fases futuras.

## 18.1 Limitações científicas

- **Painel de recall pequeno (21) e enriquecido** — `LONG_B` com só 2 entradas; 6/21 `provisional`; nenhum `DATABASE_ANNOTATION`. `[VERSIONADO]`
- **`clade` das referências não é derivado neste projeto** — vem de literatura; sem placement filogenético próprio. `[VERSIONADO]`
- **O recall mede cobertura textual, não sensibilidade real** — ver Seção 17.2. `[INFERÊNCIA]`
- **Nenhuma remoção de não‑pAgos** — o dataset "annotation‑enriched" contém falsos‑positivos textuais que só a Fase C removerá. `[CÓDIGO]` (design do prefiltro)
- **`sequence_sha256` é exato** — variantes/isoformas de comprimento diferente não casam; só pega identidade byte‑a‑byte. `[CÓDIGO]`

## 18.2 Limitações computacionais

- **Artefatos da execução real não versionados** — `data/02-intermediate/**` e `data/03-features/**` inteiramente ignorados; os dirs de `01-raw` untracked. Um clone limpo não tem os números. `[VERSIONADO]` (`.gitignore`)
- **`verify_raw_data.py` cobre só `data/01-raw`** e usa `rglob` (não respeita Git) — não verifica recall/prefilter/derived FASTA, e mistura versionado com não‑versionado. `[CÓDIGO]`
- **`XML fetch` levou ~19,5 min** para 52 473 (525 lotes, concorrência 4). `[EXEC-LOCAL]`
- **`importlib.reload` seletivo no notebook 10** é frágil a mudanças de interface nos módulos de lógica (exige Restart Kernel). `[CÓDIGO]` / `[INFERÊNCIA]`
- **`retained_records.csv` reescreve as 148 colunas + 3** para 52 473 linhas (~arquivo grande de metadados repetidos). `[EXEC-LOCAL]`

## 18.3 Limitações do painel

- Contagens desbalanceadas (`LONG_A 8` vs `LONG_B 2`).
- PIWI‑RE dominado por `CURATED_COMPUTATIONAL` (6/7).
- `NgAgo` (`WP_005580376.1`) é entrada de literatura com tétrade catalítica degradada e alegações originais disputadas — marcada `provisional`. `[VERSIONADO]` (notas)
- `SiAgo` da fixture é a cepa M.16.4, não a REY15A caracterizada — marcada `provisional`. `[VERSIONADO]` (notas)

## 18.4 Possíveis circularidades

- **PIWI‑RE `CURATED_COMPUTATIONAL` × futuro HMM PIWI‑RE** — explicitamente vedado nas notas; não afeta a pergunta textual da Fase A. `[VERSIONADO]`
- **Viés de anotação × query textual** — proteínas conhecidas estão bem anotadas exatamente nos campos que `[All Fields]` indexa; o recall alto é em parte tautológico. `[INFERÊNCIA]`

## 18.5 Decisões provisórias

- `record_selection_config_sha256` no derived FASTA recebe o `technical_prefilter_policy_sha256` — é uma escolha de "a config que produziu a seleção", razoável mas não é o sha da lista de UIDs em si (esse é `source_record_ids_sha256`). `[CÓDIGO]`
- `pago_technical_prefilter` grava em root **não sufixado** por dataset (`data/03-features/pago_technical_prefilter`) — se um 2º dataset rodar o prefiltro, o `latest/` colide (a invalidação por `source_metadata_manifest_sha256` protege a correção, mas o histórico de snapshots mistura datasets). `[CÓDIGO]` / `[INFERÊNCIA]`

## 18.6 Problemas de nomenclatura

- `annotation_enriched_candidate_set` (dataset) vs `annotation_enriched_proteome` (dataset_kind) — dois nomes para o mesmo conceito em níveis diferentes. `[VERSIONADO]`
- `query_reference_recall` sem sufixo `_final` (o plano previa `query_reference_recall_final` "após o snapshot completo"); a implementação usa o nome curto. `[CÓDIGO]` vs plano.
- Manifesto do derived FASTA tem **três** campos redundantes apontando para o mesmo pai (`derived_from_artifact_type`, `source_selection_artifact_type`; `derived_from_manifest_sha256`, `source_selection_manifest_sha256`). `[CÓDIGO]`

## 18.7 Melhorias futuras deliberadamente NÃO feitas na Fase A

- Rota de descoberta por homologia de sequência (HMM/PSI‑BLAST) — Fase futura documentada.
- Consumo do `derived_protein_fasta_snapshot` pelo SWeeP nos notebooks — o módulo está pronto, os notebooks não (Fase E).
- Versionar (curar) os snapshots reais da Fase A por nome no `.gitignore`.
- Um teste de fumaça que execute o notebook 10 (nbval/papermill).
- Resolução automatizada de IPG/GI (mantida offline e manual por decisão D7).
- `catalytic_site_status`, `system_class`, `architecture_*`, `pago_clade` por placement — Fases C/D.

---

# 19. Glossário técnico (inicial)

> Técnico agora; será reescrito pedagogicamente depois.

| Termo | Definição no contexto da Fase A |
|---|---|
| **NCBI** | National Center for Biotechnology Information; hospeda as bases de sequências (`protein`, `nuccore`, …) e o serviço E‑utilities. |
| **Entrez / E‑utilities** | API HTTP do NCBI para busca e recuperação. No código, acessada via `Bio.Entrez` (Biopython). |
| **ESearch** | endpoint que resolve uma query textual num conjunto de resultados: retorna `Count` e, com `usehistory=y`, handles de History. |
| **EFetch** | endpoint que baixa registros. `rettype` controla o formato: `uilist` (lista de UIDs), `gp`/`gb` + `retmode=xml` (GenPept/GenBank XML), `acc` (accession.version), `ipg` (Identical Protein Group), `fasta`. |
| **History API / WebEnv / QueryKey** | mecanismo do NCBI para paginar um resultado grande sem reenviar a query: `ESearch` guarda o conjunto no servidor e devolve `WebEnv` (sessão) + `QueryKey` (índice); `EFetch` os referencia. Efêmeros. |
| **UID** | *Unique Identifier* — no `db=protein`, um GI numérico (string de dígitos), o identificador interno estável do registro. |
| **accession** | identificador textual público de um registro (ex. `WP_011174533`). |
| **accession.version** | accession + sufixo de versão (`.1`, `.2`…); muda quando a sequência é revisada. A "base accession" é a parte sem `.\d+`. |
| **RefSeq** | coleção curada do NCBI; accessions de proteína começam com `WP_` (multiespécie), `NP_`, `YP_`. |
| **GenBank** | coleção primária (submissões); accessions de proteína como `ABP72561`, `KXK13845`. |
| **IPG (Identical Protein Group)** | agrupamento NCBI de todos os accessions cuja sequência de aminoácidos é idêntica. Usado para reconhecer que `ABP72561.1` e `A4WYU7.1` são a mesma proteína. |
| **CDD** | *Conserved Domain Database* do NCBI; anota regiões de domínio conservado nos registros — e essas anotações **são indexadas** por `[All Fields]`, o que torna a query textual mais completa do que só o nome do produto. |
| **Pfam** | base de famílias de domínio proteico (modelos HMM). Ex.: `PIWI PF02171`, `pPIWI_RE_X PF13111`. Relevante para as fases futuras. |
| **pAgo** | *prokaryotic Argonaute* — proteína Argonaute de bactéria/arqueia. Clados MID‑PIWI: `LONG_A`, `LONG_B` (com domínios N/PAZ/MID/PIWI), `SHORT` (só MID‑PIWI). |
| **PIWI‑RE** | família divergente (*PIWI‑RE module*), definida por Burroughs, Iyer & Aravind 2013; não é um clado de pAgo — na ontologia do projeto é `ago_family=PIWI_RE`, `clade=UNRESOLVED`. |
| **[All Fields] / [Organism]** | qualificadores de campo do Entrez. `PIWI[All Fields]` casa o termo em qualquer campo indexado (nome, definição, domínio CDD, …); `Bacteria[Organism]` restringe pela taxonomia. |
| **API** | interface programática; aqui, a HTTP do NCBI. |
| **HTTP / HTTP 5xx / HTTP 429** | protocolo web; `5xx` = erro do servidor (ex. 502 Bad Gateway); `429` = *Too Many Requests* (rate limit). |
| **XML** | formato de marcação hierárquico; o GenPept vem como `<GBSet><GBSeq>…`. |
| **CSV** | *Comma‑Separated Values*; o metadata achatado (`protein_metadata.csv`, 148 colunas). |
| **FASTA** | formato de sequência: linha `>` (defline/cabeçalho) + linhas de sequência. O derived FASTA usa deflines `>protein_uid=…|accession=…|length=…|organism=… <def>`. |
| **defline** | a linha de cabeçalho `>` de um registro FASTA. |
| **DataFrame** | tabela em memória do pandas; a unidade que a lógica pura da Fase A manipula. |
| **snapshot** | diretório imutável (`snapshots/<timestamp>__q_<hash>/`) com as saídas de um estágio + `manifest.json`, mais uma cópia `latest/` substituível. |
| **manifest** (`manifest.json`) | metadados de um snapshot: `artifact_type`, versão, hashes das saídas, hashes dos snapshots‑pai (proveniência), parâmetros. |
| **SHA‑256** | função de hash criptográfico de 256 bits (64 hex). Usada para: identidade de arquivo (`sha256_of_file`), identidade de conjunto de linhas (`sha256_of_lines`), fingerprint de sequência (`protein_sequence_sha256`), fingerprint de política/estratégia. |
| **proveniência** | a cadeia de "de onde veio": cada manifesto grava o hash do arquivo e do manifesto de cada pai, permitindo reconstruir query → preflight → UID → XML → metadata → {recall, prefilter → derived FASTA}. |
| **reprodutibilidade** | capacidade de obter o mesmo resultado; aqui garantida por snapshots imutáveis + verificação de identidade por hash + `reuse_latest_or_create`. Condicional para a rede (o NCBI pode mudar). |
| **retry** | repetir uma requisição que falhou. |
| **backoff** | esperar entre tentativas, com intervalo crescente. Preflight: `5·2^k` s. UID/XML: `0.1 → 30 s`, multiplicador 2. |
| **timeout** | limite de tempo de uma requisição individual (`fetch_timeout_seconds = 30`). |
| **deadline** | limite de tempo de uma operação inteira / lote (`request_deadline_seconds = 300`, `batch_deadline_seconds = 300`). |
| **circuit breaker** | após N falhas consecutivas (`failure_threshold = 3`), para de tentar por um período (`cooldown = 60 s`) para não martelar um serviço caído. |
| **rate limiting** | limitar a taxa de novas requisições (`max_request_starts_per_second = 8.0`). |
| **batch / lote** | grupo de itens processados juntos; XML: 100 UIDs por `EFetch`. |
| **resumable batch workspace** (`.batch_workspace/`) | diretório temporário onde cada lote XML baixado é guardado; permite retomar uma corrida interrompida sem rebaixar tudo. Purgado no sucesso. |
| **recall** | fração de itens relevantes conhecidos que foram recuperados. Aqui: referências do painel recuperadas ÷ total do painel, por estrato. |
| **benchmark** | conjunto de referência para medir desempenho; aqui `query_recall_reference_set.csv`. |
| **false negative / falso miss** | item que **deveria** casar e não casou; aqui, RsAgo contado como miss por o matcher só olhar accession (resolvido com `SEQUENCE_SHA256`). |
| **`EVALUABLE` / `NOT_EVALUABLE`** | status de um estrato de recall: `NOT_EVALUABLE` quando há 0 referências (recall = `None`, não `0.0`). |
| **`match_method`** | como uma referência casou: `EXACT_ACCESSION_VERSION`, `SAME_BASE_ACCESSION`, `SEQUENCE_SHA256`, `NONE`. |
| **`matching_strategy_sha256`** | hash do payload que descreve a metodologia de matching (`strategy_version 2.0`); gravado no manifesto e verificado no reuso. |
| **`length_warning`** | flag booleana ligada quando o comprimento está fora de [200, 2000]; **informativa**, nunca exclui. |
| **`artifact_type`** | rótulo do tipo de snapshot (`ncbi_esearch_preflight`, `query_reference_recall`, `pago_technical_prefilter`, `derived_protein_fasta_snapshot`, …); validado no `load`. |
| **`snapshot_format_version`** | versão do esquema do manifesto; um `load` rejeita versão diferente da esperada (recall exige `"1.1"`). |
| **fixture** | dado de teste versionado (`tests/fixtures/query_recall_reference_set.csv`). |
| **mock** (`unittest.mock`) | substituto de um objeto real num teste; aqui usado para não tocar a rede (`Entrez`, `run_ncbi_esearch_preflight`) e para checar que `save_*` não foi chamado no reuso. |
| **unit test** | teste de uma unidade de lógica isolada; framework `unittest`. |
| **regression test** | teste que fixa um valor/comportamento para detectar mudança não‑intencional (ex. `matching_strategy_sha256` pinado). |
| **kernel** (Jupyter) | o processo Python que executa as células do notebook; mantém estado (variáveis, `sys.modules`) entre células. |
| **module** (Python) | um arquivo `.py` importável; `sys.modules` cacheia os já carregados. |
| **import / `importlib.reload`** | `import` carrega um módulo (do cache se já visto); `reload` reexecuta o `.py` do disco **para aquele módulo**, sem recarregar seus dependentes. |
| **`SnapshotMode`** | enum: `create_new` (sempre recomputa), `reuse_latest` (erro se não há snapshot válido), `reuse_latest_or_create` (reusa se válido, senão computa) — o modo usado em todo o notebook 10. |
| **atomic write** | escrever num arquivo temporário e depois `replace` para o destino, de modo que um leitor nunca vê um arquivo meia‑escrito. |
| **namespace package** | pacote Python sem `__init__.py`; o projeto importa como `from src.pago_pipeline.X import …`. |
| **Git LFS** | *Large File Storage*; no projeto, só `data/01-raw/protein_xml_snapshots/**/*.xml` (query antiga) está sob LFS (`.gitattributes`). |
| **working tree** | os arquivos no disco do repositório, versionados ou não. Os artefatos da Fase A vivem no working tree mas fora do controle de versão. |

---

# 20. Notas de entrega

- **Escopo respeitado:** somente leitura. Nenhuma alteração no repositório; nenhum `commit`; nenhum `push`; Fase B não iniciada.
- **Comandos executados nesta auditoria** (todos read‑only ou de verificação): `git log/show/diff/rev-parse/ls-files/status/check-ignore`; leitura de arquivos; `python -c` (inspeção de módulos, releitura de manifestos, recomputo de hashes, `describe()` de CSV); `python -m unittest discover -s tests -q`; `python scripts/verify_raw_data.py`. As duas últimas são idempotentes e não escrevem no repositório.
- **HEAD ao final:** `0013c6d24d90070ee297ba6af7d12d125ebd1732` (inalterado).
- **Arquivo desta auditoria:** gerado fora do repositório (scratchpad da sessão).

## Sumário dos SHA‑256 citados

| Rótulo | Valor | Selo |
|---|---|---|
| fixture `query_recall_reference_set.csv` (HEAD) | `1da47d3f36db8f5a328f446804262e4d93ab017a49d7122602703ad3dea7a70b` | `[VERSIONADO]` |
| `matching_strategy_sha256` | `3460b048fc6de363ddf9282c2943a44c284e51d2c48092ca14426600b2871a08` | `[CÓDIGO]` |
| `technical_prefilter_policy_sha256` | `032fbc727b68ceb97cdda00ca5764db3ddeebc0868058d0c73883c23771fa523` | `[EXEC-LOCAL]` |
| metadata `csv_file_sha256` | `b39096f339a4079750b0880d49ae7a9221d07a16943bf476b1e27786a6541510` | `[EXEC-LOCAL]` |
| metadata `manifest_sha256` | `d6ab54d551432082db5b25502b1aacd1495f6a4eeb2de7892cf6735f56c80994` | `[EXEC-LOCAL]` |
| `protein_uids_sha256` | `949bcaabc6cae0a60238bcbd4ae88608dbaf6db1fa3557db5cd0ec96df5cb253` | `[EXEC-LOCAL]` |
| `xml_file_sha256` | `777390433da45c16c10bf0958dca429f03132e517dd10b13311a192268c96169` | `[EXEC-LOCAL]` |
| `retained_protein_uids.txt` / `derived.source_record_ids_sha256` | `7bea84d983c58bff09b9dc6c539377982a939c2ce9545b3a35e4bb96acde8ef0` | `[EXEC-LOCAL]` |
| derived `fasta_file_sha256` | `f7bac76b982bfe71d621416df493dc77ecda6fd3ef29b9b314cf21838954ce43` | `[EXEC-LOCAL]` |
| RsAgo `sequence_sha256` | `cbdb6bb64718c9e8ca78a34ac8445eff1556cb87b5ad687026373ed401c5fb36` | `[VERSIONADO]` |
| snapshot de recall antigo — CSV sha | `e239679ed2f0c949f2589a3b91f74d9da09952faf24b92e6734e367b714d46bf` | `[EXEC-LOCAL]` |

_Fim da auditoria._





