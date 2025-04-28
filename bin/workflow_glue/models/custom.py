"""Extended models for the workflow."""
from dataclasses import dataclass, field

from workflow_glue.models.common import Sample, WorkflowResult


@dataclass
class FastqStats:
    """Read statistics."""

    n_seqs: int | None = field(
        metadata={
            "title": "Number of reads",
            "description": "The number of sequencing reads"})
    n_bases: int | None = field(
        metadata={
            "title": "Number of bases",
            "description": "The number of bases"})
    min_length: int | None = field(
        metadata={
            "title": "Minimum read length",
            "description": "The minimum read length"})
    max_length: int | None = field(
        metadata={
            "title": "Maximum read length",
            "description": "The maximum read length"})
    mean_quality: float | None = field(
        metadata={
            "title": "Mean read quality",
            "description": "The mean read quality"})


@dataclass
class AMRVariants:
    """AMR associated variant information."""

    start: int | None = field(
        metadata={
            "title": "Start position",
            "description": """Start position of the detected gene or
            variant on the assembly contig"""})
    end: int | None = field(
        metadata={
            "title": "End position",
            "description": """End position of the detected gene or
            variant on the assembly contig"""})
    gene: str | None = field(
        metadata={
            "title": "Gene",
            "description": """Gene associated with the predicted antimicrobial
            resistance"""})
    database: str | None = field(
        metadata={
            "title": "Database",
            "description": """Database with which the associated predicted antimicrobial
            resistance was identified"""})
    drugs: list[str] | None = field(
        metadata={
            "title": "Antimicrobial drugs",
            "description": """Antimicrobials with predicted resistance associated with
            the variant or gene"""})
    aa: str | None = field(
        metadata={
            "title": "Amino acid mutation",
            "description": "Amino acid mutation"})
    nuc: str | None = field(
        metadata={
            "title": "Nucleotide mutation",
            "description": "Nucleotide mutation"})
    coverage: str | None = field(
        default=None,
        metadata={
            "title": "Coverage",
            "description": "How much of the gene is covered by the assembly"})
    identity: str | None = field(
        default=None,
        metadata={
            "title": "Identity",
            "description": """How much of the gene in the assembly is identical
            to the database gene entry"""})
    contig: str | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Assembly contig on which the variant was detected"})
    pmids: str | None = field(
        default=None,
        metadata={
            "title": "Pubmed identifiers",
            "description": "PMID or accession number for reference paper"})


@dataclass
class SequenceTypeSchema:
    """MLST schema and allele variant identified for sample."""

    schema_identifier: str | None = field(
        default=None,
        metadata={
            "title": "Schema identifier",
            "description": "Identifier for the schema used to classify the sequence"})
    allele_variant: str | None = field(
        default=None,
        metadata={
            "title": "Allele variant",
            "description": "Identifier for the specific allele variant detected"})


@dataclass
class Serotype:
    """Salmonella serotyping results."""

    predicted_serotype: str | None = field(
        default=None,
        metadata={
            "title": "Predicted serotype",
            "description": """Predicited serological typing idenitifer for the sample
            (Salmonella only)"""})
    predicted_antigenic_profile: str | None = field(
        default=None,
        metadata={
            "title": "Predicted antigenic profile",
            "description": """Predicted antigenic profile for the sample
            using the O, H1 and H2 antigens identified (Salmonella only)"""})
    o_antigen_prediction: str | None = field(
        default=None,
        metadata={
            "title": "O antigen prediction",
            "description": """Predicted O antigen found in the sample
            (Salmonella only)"""})
    h1_antigen_prediction: str | None = field(
        default=None,
        metadata={
            "title": "H1 antigen prediction",
            "description": """Predicted H1 antigen found in the sample
            (Salmonella only)"""})
    h2_antigen_prediction: str | None = field(
        default=None,
        metadata={
            "title": "H2 antigen prediction",
            "description": """Predicted H2 antigen found in the sample
            (Salmonella only)"""})


@dataclass
class Annotation:
    """Region of interest identified within assembly."""

    contig: str | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Assembly contig on which region was detected"})
    ID: str | None = field(
        default=None,
        metadata={
            "title": "Identifier",
            "description": "Unique identifier for the annotation from prokka"})
    start: int | None = field(
        default=None,
        metadata={
            "title": "Start position",
            "description": "Start position of the region on the assembly contig"})
    end: int | None = field(
        default=None,
        metadata={
            "title": "End position",
            "description": "End position of the region on the assembly contig"})
    strand: str | None = field(
        default=None,
        metadata={
            "title": "Strand",
            "description": "Which strand the region was identified on"})
    gene: str | None = field(
        default=None,
        metadata={
            "title": "Gene",
            "description": "Gene name of the region of interest, if applicable"})
    product: str | None = field(
        default=None,
        metadata={
            "title": "Product",
            "description": "Product name of the gene or region, if applicable"})
    ec_number: str | None = field(
        default=None,
        metadata={
            "title": "EC number",
            "description": "Identifier from the enzyme consortium catalogue"})


@dataclass
class Variant:
    """Variants identified in assembly compared to reference."""

    contig: str | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Contig on which the variant was detected"})
    pos: int | None = field(
        default=None,
        metadata={
            "title": "Position",
            "description": ""})
    ref: str | None = field(
        default=None,
        metadata={
            "title": "Reference",
            "description": "The reference nucleotide at the variant site"})
    alt: str | None = field(
        default=None,
        metadata={
            "title": "Alternate allele",
            "description": "The alternate nucleotide at the variant site"})
    depth: int | None = field(
        default=None,
        metadata={
            "title": "Depth",
            "description": "The frequency of the variant in the sample"})


@dataclass
class Coverage:
    """Coverage summary information for each contig in assembly."""

    counts: int | None = field(
        default=None,
        metadata={
            "title": "Counts",
            "description": "Number of reads that map to the contig"})
    median: float | None = field(
        default=None,
        metadata={
            "title": "Median coverage",
            "description": "Median coverage"})
    mean: float | None = field(
        default=None,
        metadata={
            "title": "Mean coverage",
            "description": "Mean coverage"})
    minimum: int | None = field(
        default=None,
        metadata={
            "title": "Minimum coverage",
            "description": "Minimum coverage"})
    maximum: int | None = field(
        default=None,
        metadata={
            "title": "Maximum coverage",
            "description": "Maximum coverage"})


@dataclass
class AntimicrobialResistance:
    """The antimicrobial resistance results for the sample."""

    detected_variants: list[AMRVariants] | None = field(
        metadata={
            "title": "Detected variants",
            "description": """A list of all predicted antimicrobial resistance
            genes or variants identified in the sample"""})


@dataclass
class MLST:
    """Multi-locus sequence typing results for the sample."""

    detected_species: str | None = field(
        default=None,
        metadata={
            "title": "Detected species",
            "description": """The detected species and MLST scheme
            used to classify the sample"""})
    sequence_type: str | None = field(
        default=None,
        metadata={
            "title": "Sequence type",
            "description": "The sequence type assigned to the sample"})
    typing_schema: list[SequenceTypeSchema] | None = field(
        default=None,
        metadata={
            "title": "Typing schema",
            "description": """The MLST schema alleles and variants
            identified in the sample"""})


@dataclass
class Contig:
    """Summary statistics for contig in assembly."""

    name: str | None = field(
        default=None,
        metadata={
            "title": "Name",
            "description": "The name of the contig"})
    length: int | None = field(
        default=None,
        metadata={
            "title": "Contig length",
            "description": "The length of the contig in base pairs"})
    coverage: Coverage | None = field(
        default=None,
        metadata={
            "title": "Coverage",
            "description": "Mean coverage of the contig in the sample"})


@dataclass
class Assembly:
    """Draft genome assembly statistics of the sample."""

    reference: str | None = field(
        default=None,
        metadata={
            "title": "Reference name",
            "description": """Name of the reference used in the assembly
            process. Null for de-novo"""})

    annotations: list[Annotation] | None = field(
        default=None,
        metadata={
            "title": "Annotations",
            "description": """Array of regions of interest identified within
            the assembly"""})
    variants: list[Variant] | None = field(
        default=None,
        metadata={
            "title": "Variants",
            "description": "A list of variants identified in the assembly"})
    contig: list[Contig] | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Summary statistics for each contig in the assembly"})


@dataclass
class ResultsContents:
    """Results for a sample."""

    antimicrobial_resistance: AntimicrobialResistance
    assembly: Assembly
    sequence_typing: MLST
    serotyping: Serotype
    fastq: FastqStats


@dataclass
class Sample(Sample):
    """A sample sheet entry and its corresponding checks and related results."""

    sample_pass: bool | None = field(
        default=None,
        metadata={
            "title": "Sample pass",
            "description": "If true the sample has passed workflow checks"})


@dataclass
class WorkflowResult(WorkflowResult):
    """Definition for results that will be returned by this workflow."""

    samples: list[Sample] = field(
        metadata={
            "title": "Samples",
            "description": "Samples in this workflow instance"})

    def __post_init__(self):
        """Determine overall status for the workflow given the sample results."""
        self.workflow_pass = all(
            sample.sample_pass for sample in self.samples)
