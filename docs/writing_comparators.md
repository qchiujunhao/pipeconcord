# Writing Comparators

Custom comparators subclass `Comparator` and return a `ConcordanceReport`.

```python
from biocompare.comparators.base import Comparator
from biocompare.core.report import ConcordanceReport


class MyComparator(Comparator):
    name = "my-format"
    supported_types = ("my-format",)

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        return kwargs.get("file_type") == "my-format"

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        return ConcordanceReport(
            comparator=self.__class__.__name__,
            file_a=file_a,
            file_b=file_b,
            overall_concordance=1.0,
            metrics={"example": 1.0},
        )
```

Register it through `pyproject.toml`:

```toml
[project.entry-points."biocompare.comparators"]
my_format = "mypackage.compare:MyComparator"
```

