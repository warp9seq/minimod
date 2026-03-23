# Major changes from v0.4.0 to v0.5.0

- For <=v0.4.0, when non CG contexts were requested, both minimod view and freq were reporting modifications from erroneous contexts. We have fixed and tested this issue from >= v0.5.0.

- For <=v0.4.0, we allowed primary, secondary, supplementary alignments when viewing and calculating frequencies. From >=v0.5.0, we only consider primary and supplementary alignments by default. We introduced --allow-secondary option to enable considering secondary alignments. Minimod still errors out if hard-clipping is found.

- For <=v0.4.0, when --insertions option is used, erroneous reference positions were included in both freq and view outputs. We have fixed it in >=v0.5.0.

- For <=v0.4.0, status of skipped bases (low probability C+m. or unknown C+m? in MM tag) were not handled. We have fixed it in >=v0.5.0. When skipped as low probability, minimod consider them to have 0 modification probability. Otherwise, when skipped as unknown, they are ignored.

- For >=v0.5.0, mod_prob in the output of view subtool will slightly change due to a change on how probability conversion from unsigned 8-bit N to float probability p. Previously, it was p = N/255.0 and now we use p = (N+0.5)/256.0 as the conversion function.

- For >=v0.5.0, when a context is provided (ex: -c h[CG]) option, we output modifications where the modified read base matches aligned reference base. However, when the context is * (ex: -c a[*]), this comparison is ignored and modifications at both matched and mismatched positions are output. Note that only the modified base is compared with reference base, not the whole context.

- For >=v0.5.0, we print a warning if a certain modification and a context is untested.  For versions before this, only m[CG] for DNA (genome) and a[A] for transcriptome were tested.