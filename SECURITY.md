# Security Policy

## Supported Versions

Security fixes are considered for the active development branch and the most
recent release.

| Version | Supported |
|---------|-----------|
| `0.1.x` | Yes |
| `< 0.1.0` | No |

## Reporting a Vulnerability

Please do not report security vulnerabilities in public GitHub issues.

Report suspected vulnerabilities by email to `adamrsophie@gmail.com`, or use
GitHub private vulnerability reporting if it is enabled for the repository.

Include as much detail as you can:

- affected version or commit,
- operating system and R version,
- minimal reproduction steps,
- any proof-of-concept code.

## Response Process

The maintainer will acknowledge reports as soon as practical, investigate the
issue, and coordinate a fix if the report is confirmed. Public disclosure should
wait until a fix or mitigation is available, unless there is an overriding user
safety reason to disclose earlier.

## Scope

Security-sensitive areas include:

- compilation of Stan programs: a program is translated into R code and
  evaluated, so a Stan program, a data or initial value file, or a CSV file
  that makes bayescompR run R code it did not ask for would be a
  vulnerability (user functions supplied through `functions =` are R code by
  design),
- arbitrary file writes or path traversal, for example through `#include`,
  output file names in `bayescomp_write_csv()` or data files,
- unbounded memory or time on small inputs, for example in the parser,
- exposure of private data in examples, tests, or logs.

A wrong posterior, density or gradient that does not create a
software-security vulnerability should be opened as a normal GitHub issue
instead.
