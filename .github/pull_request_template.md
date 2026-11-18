## Problem and resulting capability

Catalogue ID, related Issue, supported formulation/version and scientific scope.

## Scientific and numerical evidence

Reference/DOI, equations, adaptations, parameters and validation provenance. Distinguish numerical verification, source reproduction and independent validation.

## Checks

Commands/results, preserved existing behavior and remaining limitations. Update catalogue JSON and regenerate documentation when status or evidence changes.

## Contribution checklist

- [ ] Exact primary formulation/version, scope and relationship to existing criteria identified.
- [ ] Required inputs, units, assumptions and departures documented.
- [ ] Current plain/smooth CA scope or explicit future-extension status stated.
- [ ] External source files/revisions, variant correspondence, license permissions and retained notices reviewed before any reuse/translation; availability has not been mistaken for WG4.3 integration.
- [ ] Meaningful numerical checks and independently computed reference result provided, or scientific gaps stated.
- [ ] Implementation completeness and validation evidence classified separately.
- [ ] If the PR can alter numerical outputs, bumped `IMPLEMENTATION_REVISION` in `src/multifat/models/references.py` and refreshed affected `examples/` data manifests.
- [ ] Catalogue disposition and actual Issue links updated; generated docs pass `--check`.
