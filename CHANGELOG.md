# Changelog

All notable changes to this project will be documented in this file.

The format is based on a mixture of [Keep a Changelog] and [Common Changelog].
This project adheres to [Semantic Versioning], with the exception that minor releases may include breaking changes.

## [Unreleased]

## [3.0.0] - 2025-05-08

_If you are upgrading: please see [`UPGRADING.md`](UPGRADING.md#300)._

### Added

- ✨ Support Qiskit 2.0+ ([#610]) ([**@burgholzer**])

### Changed

- **Breaking**: 🚚 Move MQT QMAP to the [munich-quantum-toolkit] GitHub organization
- **Breaking**: ♻️ Use the `mqt-core` Python package for handling circuits ([#418]) ([**@burgholzer**])
- **Breaking**: ⬆️ Bump minimum required CMake version to `3.24.0` ([#621]) ([**@burgholzer**])
- **Breaking**: ♻️ Adopt new `NAComputation` in NASP tool ([#608]) ([**@ystade**])
- **Breaking**: ♻️ Isolate NALAC from the main library ([#608], [#609]) ([**@ystade**])
- 📝 Rework existing project documentation ([#614]) ([**@burgholzer**])

### Removed

- **Breaking**: 🔥 Remove teleportation support for the heuristic mapping ([#621]) ([**@burgholzer**])
- **Breaking**: 🔥 Remove support for `BackendV1` Qiskit backends ([#610]) ([**@burgholzer**])
- **Breaking**: 🔥 Remove support for `.real`, `.qc`, `.tfc`, and `GRCS` files ([#621]) ([**@burgholzer**])
- **Breaking**: 🔥 Remove `yaml-cpp` dependency ([#608]) ([**@ystade**])

## [2.8.0] - 2024-11-18

_📚 Refer to the [GitHub Release Notes] for previous changelogs._

<!-- Version links -->

[unreleased]: https://github.com/munich-quantum-toolkit/qcec/compare/v3.0.0...HEAD
[3.0.0]: https://github.com/munich-quantum-toolkit//compare/v2.8.0...v3.0.0
[2.8.0]: https://github.com/munich-quantum-toolkit//releases/tag/v2.8.0

<!-- PR links -->

[#621]: https://github.com/munich-quantum-toolkit//pulls/621
[#614]: https://github.com/munich-quantum-toolkit//pulls/614
[#610]: https://github.com/munich-quantum-toolkit//pulls/610
[#609]: https://github.com/munich-quantum-toolkit//pulls/609
[#608]: https://github.com/munich-quantum-toolkit//pulls/608
[#418]: https://github.com/munich-quantum-toolkit//pulls/418

<!-- Contributor -->

[**@burgholzer**]: https://github.com/burgholzer
[**@ystade**]: https://github.com/ystade

<!-- General links -->

[Keep a Changelog]: https://keepachangelog.com/en/1.1.0/
[Common Changelog]: https://common-changelog.org
[Semantic Versioning]: https://semver.org/spec/v2.0.0.html
[GitHub Release Notes]: https://github.com/munich-quantum-toolkit/qmap/releases
[munich-quantum-toolkit]: https://github.com/munich-quantum-toolkit
