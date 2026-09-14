# M8.10 author package: landing manifest

Sixteen files were pinned by SHA-256 in [`../../tasks/m8_10_task_details.md`](../../tasks/m8_10_task_details.md). Fourteen land byte-identical. Two land as deterministic privacy redactions of the pinned originals, produced by [`redact.py`](redact.py), whose docstring states the rules and the gates. With `redact.py` and this manifest, the directory holds eighteen files.

| file | status | pinned SHA-256 | landed SHA-256 | substitutions |
| --- | --- | --- | --- | --- |
| `m810_core.py` | BYTE-IDENTICAL | `9339d87117a4fe524b913f5fe63972f9a19dfaa0d708f059323bde7553f7ecda` | `9339d87117a4fe524b913f5fe63972f9a19dfaa0d708f059323bde7553f7ecda` | - |
| `m810_selftest.py` | BYTE-IDENTICAL | `d0570144b051b1435c1194773ed10c3158e3b294e38107db541f2e8421ff9251` | `d0570144b051b1435c1194773ed10c3158e3b294e38107db541f2e8421ff9251` | - |
| `m810_parent.py` | BYTE-IDENTICAL | `bf4aeb556be3d73ad8adbadf8109c251ec32f8cce93fd90d1f249a6c5a38d965` | `bf4aeb556be3d73ad8adbadf8109c251ec32f8cce93fd90d1f249a6c5a38d965` | - |
| `m810_levels.py` | BYTE-IDENTICAL | `9d00d7ab200247ab806b45b93a80847b011fccb13467addfdaaf1d4623c41de4` | `9d00d7ab200247ab806b45b93a80847b011fccb13467addfdaaf1d4623c41de4` | - |
| `m810_exact.py` | BYTE-IDENTICAL | `8436ef9ffc465e6a91878b4b994c1385766167cc8ebc4047d482cf65c72463a9` | `8436ef9ffc465e6a91878b4b994c1385766167cc8ebc4047d482cf65c72463a9` | - |
| `m810_gates.py` | BYTE-IDENTICAL | `70c6490296272f9be60aec1ce8e45772bd9e93afe746eb13f5a51b57bbf58d1a` | `70c6490296272f9be60aec1ce8e45772bd9e93afe746eb13f5a51b57bbf58d1a` | - |
| `m810_zeros.py` | BYTE-IDENTICAL | `531a5f3e87c063e11e6d770d84d259f691e4697a741ebea02201a345acadc21b` | `531a5f3e87c063e11e6d770d84d259f691e4697a741ebea02201a345acadc21b` | - |
| `m810_results.py` | BYTE-IDENTICAL | `8fe498ff0cf59917bb346ba35a028fcef1557f6cdb6f1edf8d55c4e19d3f4d62` | `8fe498ff0cf59917bb346ba35a028fcef1557f6cdb6f1edf8d55c4e19d3f4d62` | - |
| `RESULTS.md` | BYTE-IDENTICAL | `9986195a4199d76db970c052bf5ab93b9ad373f3b190bc6dcb31226821a742ae` | `9986195a4199d76db970c052bf5ab93b9ad373f3b190bc6dcb31226821a742ae` | - |
| `POST_DERIVATION_NOTE.md` | BYTE-IDENTICAL | `4f6b96a0330ea982bfcfd29e48356897a5a3e8fb114d540b8eb169ade4007620` | `4f6b96a0330ea982bfcfd29e48356897a5a3e8fb114d540b8eb169ade4007620` | - |
| `audit_dryrun.py` | PRIVACY-REDACTED | `e0cbae2dea1750bf0692ea8f6698e7451b983cf0261ff7489b3f2a439ec4d371` | `f2414f480ac0e44b9ee5bf2eebc09e928b6ceb67af1f9b34aa0f9445ff356e96` | P1 0, P2 2 |
| `compare_dryrun.py` | BYTE-IDENTICAL | `dc45e596664a6fdb3d7bf8c955a97af54ea2c2b37c01bbc10b1e8034acaff4cd` | `dc45e596664a6fdb3d7bf8c955a97af54ea2c2b37c01bbc10b1e8034acaff4cd` | - |
| `dryrun/BRIEF.md` | BYTE-IDENTICAL | `4918be4276d2a5e3849bfdbe4e997df5ab4c05f081bb00aea1702c8299032541` | `4918be4276d2a5e3849bfdbe4e997df5ab4c05f081bb00aea1702c8299032541` | - |
| `dryrun/worklist.md` | BYTE-IDENTICAL | `b08d23c6ceba1ba18cc9ba53b3801d3654e11f495c4b9945e97e353692348248` | `b08d23c6ceba1ba18cc9ba53b3801d3654e11f495c4b9945e97e353692348248` | - |
| `dryrun/RETURN.md` | BYTE-IDENTICAL | `c048a1268a76ea9efcdafa8f1d6614351d8d8b92160312639ad24490b897a2c7` | `c048a1268a76ea9efcdafa8f1d6614351d8d8b92160312639ad24490b897a2c7` | - |
| `dryrun/transcript.jsonl` | PRIVACY-REDACTED | `a96cb211671baaaf72e8237e52e6244d02ae9c8a53c29b6e734eaa2296d6ddd7` | `b24ead8ba42a200d5fcb67d8ba85564e613308d0e522f92397c2966ba1ad07c9` | P1 3, P2 199, T1 19 |

T1 replaces a whole line with a marker carrying the SHA-256 of the removed line; P1 and P2 are the path rules. G3 and G4 checked the files as written against the pinned originals, before any landed code ran and again on the finished tree: every startup and rate-limit record is replaced by the marker carrying its hash, and every other line, with the path rules reversed, equals its original exactly. On the landed files, `audit_dryrun.py` reproduces the containment audit (55 tool calls, CLEAN) and `compare_dryrun.py` reproduces 103 matches with no mismatch, both exiting 0.

G2 found zero residuals in the landed UTF-8 text of all eighteen files, this manifest and `redact.py` included, so a home path the rules could not reach, such as a JSON-escaped one, would have refused the write. The sensitive-pattern list G2 applied was supplied at run time, so this repository names none of its entries; each entry matched the pinned transcript, so none was dead. The list's SHA-256 is `dd9ad7785a42ea2a6b00f4202c54cfd6ea584cd28c61e8152f56dcc8a8d1d330`. `check_no_local.py` does not read `.jsonl` files, so G2 is what covers the transcript.
