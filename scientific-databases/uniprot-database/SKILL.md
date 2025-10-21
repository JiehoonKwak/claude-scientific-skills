---
name: uniprot-database
description: "Access and query the UniProt protein database for comprehensive protein information retrieval. Use this skill when you need to search for proteins by name, gene symbol, accession number, or functional terms; retrieve protein sequences in FASTA format; access detailed protein annotations including function, structure, interactions, and pathways; map protein identifiers between different databases (Ensembl, RefSeq, PDB, KEGG, GO terms); query Swiss-Prot (reviewed) or TrEMBL (unreviewed) protein entries; perform batch operations on multiple proteins; download protein datasets; analyze protein families and domains; investigate protein-protein interactions; explore evolutionary relationships through UniRef clusters; access protein structure predictions from AlphaFoldDB; retrieve Gene Ontology annotations; analyze protein modifications and post-translational sites; or work with protein sequences for bioinformatics analysis. This skill provides programmatic access to UniProtKB, UniRef, UniParc, and related databases through REST API endpoints with support for various output formats (JSON, TSV, FASTA, XML)."
---

# UniProt Database

## Overview

Interact with UniProt (Universal Protein Resource), the world's leading comprehensive and freely accessible resource for protein sequence and functional information. This skill enables programmatic access to UniProtKB (Swiss-Prot and TrEMBL), UniRef, UniParc, and other UniProt databases through the REST API.

**Key databases:**
- **UniProtKB/Swiss-Prot**: Manually annotated, reviewed, high-quality protein entries
- **UniProtKB/TrEMBL**: Automatically annotated, unreviewed protein entries
- **UniRef**: Clustered sets of sequences for similarity searches
- **UniParc**: Comprehensive archive of all protein sequences

## Core Capabilities

### 1. Searching for Proteins

Search UniProt using natural language queries or structured search syntax.

**Common search patterns:**
```python
# Search by protein name
query = "insulin AND organism_name:\"Homo sapiens\""

# Search by gene name
query = "gene:BRCA1 AND reviewed:true"

# Search by accession
query = "accession:P12345"

# Search by sequence length
query = "length:[100 TO 500]"

# Search by taxonomy
query = "taxonomy_id:9606"  # Human proteins

# Search by GO term
query = "go:0005515"  # Protein binding
```

Use the API search endpoint: `https://rest.uniprot.org/uniprotkb/search?query={query}&format={format}`

**Supported formats:** JSON, TSV, Excel, XML, FASTA, RDF, TXT

### 2. Retrieving Individual Protein Entries

Retrieve specific protein entries by accession number.

**Accession number formats:**
- Classic: P12345, Q1AAA9, O15530 (6 characters: letter + 5 alphanumeric)
- Extended: A0A022YWF9 (10 characters for newer entries)

**Retrieve endpoint:** `https://rest.uniprot.org/uniprotkb/{accession}.{format}`

Example: `https://rest.uniprot.org/uniprotkb/P12345.fasta`

### 3. Batch Retrieval and ID Mapping

Map protein identifiers between different database systems and retrieve multiple entries efficiently.

**ID Mapping workflow:**
1. Submit mapping job to: `https://rest.uniprot.org/idmapping/run`
2. Check job status: `https://rest.uniprot.org/idmapping/status/{jobId}`
3. Retrieve results: `https://rest.uniprot.org/idmapping/results/{jobId}`

**Supported databases for mapping:**
- UniProtKB AC/ID
- Gene names
- Ensembl, RefSeq, EMBL
- PDB, AlphaFoldDB
- KEGG, GO terms
- And many more (see `/references/id_mapping_databases.md`)

**Limitations:**
- Maximum 100,000 IDs per job
- Results stored for 7 days

### 4. Streaming Large Result Sets

For large queries that exceed pagination limits, use the stream endpoint:

`https://rest.uniprot.org/uniprotkb/stream?query={query}&format={format}`

The stream endpoint returns all results without pagination, suitable for downloading complete datasets.

### 5. Customizing Retrieved Fields

Specify exactly which fields to retrieve for efficient data transfer.

**Common fields:**
- `accession` - UniProt accession number
- `id` - Entry name
- `gene_names` - Gene name(s)
- `organism_name` - Organism
- `protein_name` - Protein names
- `sequence` - Amino acid sequence
- `length` - Sequence length
- `go_*` - Gene Ontology annotations
- `cc_*` - Comment fields (function, interaction, etc.)
- `ft_*` - Feature annotations (domains, sites, etc.)

**Example:** `https://rest.uniprot.org/uniprotkb/search?query=insulin&fields=accession,gene_names,organism_name,length,sequence&format=tsv`

See `/references/api_fields.md` for complete field list.

## Python Implementation

For programmatic access, use the provided helper script `scripts/uniprot_client.py` which implements:

- `search_proteins(query, format)` - Search UniProt with any query
- `get_protein(accession, format)` - Retrieve single protein entry
- `map_ids(ids, from_db, to_db)` - Map between identifier types
- `batch_retrieve(accessions, format)` - Retrieve multiple entries
- `stream_results(query, format)` - Stream large result sets

**Alternative Python packages:**
- **Unipressed**: Modern, typed Python client for UniProt REST API
- **bioservices**: Comprehensive bioinformatics web services client

## Query Syntax Examples

**Boolean operators:**
```
kinase AND organism_name:human
(diabetes OR insulin) AND reviewed:true
cancer NOT lung
```

**Field-specific searches:**
```
gene:BRCA1
accession:P12345
organism_id:9606
taxonomy_name:"Homo sapiens"
annotation:(type:signal)
```

**Range queries:**
```
length:[100 TO 500]
mass:[50000 TO 100000]
```

**Wildcards:**
```
gene:BRCA*
protein_name:kinase*
```

See `/references/query_syntax.md` for comprehensive syntax documentation.

## Best Practices

1. **Use reviewed entries when possible**: Filter with `reviewed:true` for Swiss-Prot (manually curated) entries
2. **Specify format explicitly**: Choose the most appropriate format (FASTA for sequences, TSV for tabular data, JSON for programmatic parsing)
3. **Use field selection**: Only request fields you need to reduce bandwidth and processing time
4. **Handle pagination**: For large result sets, implement proper pagination or use the stream endpoint
5. **Cache results**: Store frequently accessed data locally to minimize API calls
6. **Rate limiting**: Be respectful of API resources; implement delays for large batch operations
7. **Check data quality**: TrEMBL entries are computational predictions; Swiss-Prot entries are manually reviewed

## Resources

### scripts/
`uniprot_client.py` - Python client with helper functions for common UniProt operations including search, retrieval, ID mapping, and streaming.

### references/
- `api_fields.md` - Complete list of available fields for customizing queries
- `id_mapping_databases.md` - Supported databases for ID mapping operations
- `query_syntax.md` - Comprehensive query syntax with advanced examples
- `api_examples.md` - Code examples in multiple languages (Python, curl, R)

## Additional Resources

- **API Documentation**: https://www.uniprot.org/help/api
- **Interactive API Explorer**: https://www.uniprot.org/api-documentation
- **REST Tutorial**: https://www.uniprot.org/help/uniprot_rest_tutorial
- **Query Syntax Help**: https://www.uniprot.org/help/query-fields
- **SPARQL Endpoint**: https://sparql.uniprot.org/ (for advanced graph queries)
