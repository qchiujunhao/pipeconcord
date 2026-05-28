from __future__ import annotations

from biocompare.comparators.bam_stats import BAMStatsComparator
from biocompare.comparators.bed import BEDComparator
from biocompare.comparators.counts import CountsComparator
from biocompare.comparators.deg import DEGComparator
from biocompare.comparators.fasta import FASTAComparator
from biocompare.comparators.table import TableComparator
from biocompare.comparators.vcf import VCFComparator
from biocompare.core.registry import ComparatorRegistry


def register_builtin_comparators(registry: type[ComparatorRegistry] = ComparatorRegistry) -> None:
    registry.register(DEGComparator)
    registry.register(CountsComparator)
    registry.register(BEDComparator)
    registry.register(FASTAComparator)
    registry.register(VCFComparator)
    registry.register(BAMStatsComparator)
    registry.register(TableComparator)


__all__ = ["BAMStatsComparator", "BEDComparator", "CountsComparator", "DEGComparator", "FASTAComparator", "TableComparator", "VCFComparator", "register_builtin_comparators"]
