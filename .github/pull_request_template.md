# Summary

Describe the change and why it is needed.

## Type of Change

- [ ] Bug fix (wrong answer or error)
- [ ] Stan language or library feature
- [ ] Inference method or diagnostic
- [ ] Performance
- [ ] Documentation
- [ ] Refactor or maintenance
- [ ] Release preparation

## Stan Parity

Name the part of the Stan reference manual, Stan Math or CmdStan the change
follows, and how the expected values in the tests were obtained. If behavior
intentionally differs from Stan, say why.

## Verification

List the checks you ran.

```text
# commands and results
```

## Documentation

- [ ] README updated if needed
- [ ] Help pages or roxygen updated if needed
- [ ] Vignettes updated if needed
- [ ] NEWS.md updated for user-facing changes

## Checklist

- [ ] The change is focused and reviewable
- [ ] Relevant tests pass, including gradient tests against finite differences
- [ ] No local build artifacts, `.Rcheck` output, tarballs, or caches are committed
- [ ] The pull request follows the Code of Conduct
